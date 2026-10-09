#include "wfmat/kernel.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

#include <boost/multiprecision/cpp_bin_float.hpp>

namespace wfmat {

namespace {

using quad = boost::multiprecision::cpp_bin_float_quad;

struct Vec3 {
    double x = 0.0, y = 0.0, t = 0.0;

    friend Vec3 operator+(Vec3 a, Vec3 b) { return {a.x + b.x, a.y + b.y, a.t + b.t}; }
    friend Vec3 operator-(Vec3 a, Vec3 b) { return {a.x - b.x, a.y - b.y, a.t - b.t}; }
    friend Vec3 operator*(double s, Vec3 a) { return {s * a.x, s * a.y, s * a.t}; }
};

double dot3(Vec3 a, Vec3 b) { return a.x * b.x + a.y * b.y + a.t * b.t; }
Vec3 cross3(Vec3 a, Vec3 b) { return {a.y * b.t - a.t * b.y, a.t * b.x - a.x * b.t, a.x * b.y - a.y * b.x}; }
double norm3(Vec3 a) { return std::sqrt(dot3(a, a)); }

// The Minkowski form of Laguerre geometry, Q(u, v) = ux vx + uy vy - ut vt.
double mink(Vec3 u, Vec3 v) { return u.x * v.x + u.y * v.y - u.t * v.t; }

// A site equation q Q(z) + a.z + b = 0 with q in {0, 1} [section 7]:
//   line:   n.x - t - c = 0
//   circle: |x|^2 - t^2 - 2 c.x + 2 sigma R t + |c|^2 - R^2 = 0
struct Equation {
    bool quadratic = false;
    Vec3 a;
    double b = 0.0;
};

Equation equation(const Site& s)
{
    if (s.is_line()) return {false, {s.n.x, s.n.y, -1.0}, -s.c};
    return {true,
            {-2.0 * s.centre.x, -2.0 * s.centre.y, 2.0 * s.sigma * s.R},
            dot(s.centre, s.centre) - s.R * s.R};
}

int sign(double v) { return (v > 0.0) - (v < 0.0); }

// Allowed |d(p) - t| for an accepted point: eps_geom, plus the rounding error of evaluating the
// distance to a circle, which grows with its radius (a nearly flat fillet can have R = 10^6).
double residual_slack(const Site& s, Vec2 p, const Tolerances& tol)
{
    if (s.is_line()) return tol.geom;
    return tol.geom + 8.0 * std::numeric_limits<double>::epsilon() * (s.R + dist(p, s.centre));
}

// Newton on F_i(x, y, t) = d_i(x, y) - t, i = 1..3, in arithmetic T.
template <class T>
struct Newton {
    const std::array<Site, 3>& sites;

    // Value and gradient of d_i at (x, y) in T.
    void eval(const Site& s, const T& x, const T& y, T& d, T& gx, T& gy) const
    {
        using std::sqrt;
        if (s.is_line()) {
            d = T(s.n.x) * x + T(s.n.y) * y - T(s.c);
            gx = T(s.n.x);
            gy = T(s.n.y);
            return;
        }
        const T vx = x - T(s.centre.x), vy = y - T(s.centre.y);
        const T r = sqrt(vx * vx + vy * vy);
        d = T(s.sigma) * (T(s.R) - r);
        if (r == 0) {
            gx = gy = T(0);
        } else {
            gx = -T(s.sigma) * vx / r;
            gy = -T(s.sigma) * vy / r;
        }
    }

    // Runs at most max_iter steps; returns false when the Jacobian is singular.
    bool run(T& x, T& y, T& t, int max_iter, const T& stop) const
    {
        using std::abs;
        for (int it = 0; it < max_iter; ++it) {
            T f[3], g[3][3];
            for (int i = 0; i < 3; ++i) {
                T d, gx, gy;
                eval(sites[i], x, y, d, gx, gy);
                f[i] = d - t;
                g[i][0] = gx;
                g[i][1] = gy;
                g[i][2] = T(-1);
            }
            const T det = g[0][0] * (g[1][1] * g[2][2] - g[1][2] * g[2][1]) -
                          g[0][1] * (g[1][0] * g[2][2] - g[1][2] * g[2][0]) +
                          g[0][2] * (g[1][0] * g[2][1] - g[1][1] * g[2][0]);
            if (det == 0) return false;
            // Cramer's rule for J delta = -f.
            T delta[3];
            for (int k = 0; k < 3; ++k) {
                T m[3][3];
                for (int i = 0; i < 3; ++i)
                    for (int j = 0; j < 3; ++j) m[i][j] = j == k ? -f[i] : g[i][j];
                delta[k] = (m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1]) -
                            m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0]) +
                            m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])) /
                           det;
            }
            x += delta[0];
            y += delta[1];
            t += delta[2];
            if (abs(delta[0]) + abs(delta[1]) + abs(delta[2]) <= stop * (T(1) + abs(x) + abs(y) + abs(t)))
                return true;
        }
        return true;
    }
};

// Infinity-norm condition number of the Jacobian of the unsquared system at (p, t).
double jacobian_condition(const std::array<Site, 3>& sites, Vec2 p)
{
    double g[3][3];
    for (int i = 0; i < 3; ++i) {
        const Vec2 grad = sites[i].gradient(p);
        g[i][0] = grad.x;
        g[i][1] = grad.y;
        g[i][2] = -1.0;
    }
    double inv[3][3];
    const double det = g[0][0] * (g[1][1] * g[2][2] - g[1][2] * g[2][1]) -
                       g[0][1] * (g[1][0] * g[2][2] - g[1][2] * g[2][0]) +
                       g[0][2] * (g[1][0] * g[2][1] - g[1][1] * g[2][0]);
    if (det == 0.0) return std::numeric_limits<double>::infinity();
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j) {
            const int r0 = (j + 1) % 3, r1 = (j + 2) % 3, c0 = (i + 1) % 3, c1 = (i + 2) % 3;
            inv[i][j] = (g[r0][c0] * g[r1][c1] - g[r0][c1] * g[r1][c0]) / det;
        }
    double na = 0.0, ni = 0.0;
    for (int i = 0; i < 3; ++i) {
        na = std::max(na, std::abs(g[i][0]) + std::abs(g[i][1]) + std::abs(g[i][2]));
        ni = std::max(ni, std::abs(inv[i][0]) + std::abs(inv[i][1]) + std::abs(inv[i][2]));
    }
    return na * ni;
}

Root3 make_root(const std::array<Site, 3>& sites, Vec2 p, double t, double condition)
{
    Root3 r;
    r.p = p;
    r.t = t;
    r.condition = condition;
    for (int i = 0; i < 3; ++i) r.feet[i] = sites[i].foot(p, t);
    return r;
}

// Condition beyond which a double-precision root is refined in binary128.
constexpr double quad_threshold = 1e4;

} // namespace

double Site::distance(Vec2 x) const
{
    return is_line() ? dot(n, x) - c : sigma * (R - dist(x, centre));
}

Vec2 Site::gradient(Vec2 x) const
{
    if (is_line()) return n;
    const Vec2 v = x - centre;
    const double r = norm(v);
    return r == 0.0 ? Vec2{} : (-sigma / r) * v;
}

Vec2 Site::foot(Vec2 x, double t) const
{
    return x - t * gradient(x);
}

double Site::offset_radius(double t) const
{
    return R - sigma * t;
}

Site site_of(const Segment& s)
{
    if (!s.is_arc()) {
        const Vec2 n = perp(unit(s.b - s.a));
        return Site::line(n, dot(n, s.a));
    }
    return Site::circle(s.c, s.R, s.sweep > 0.0 ? 1 : -1);
}

ConicKind bisector_kind(const Site& a, const Site& b, double eps)
{
    if (a.is_line() && b.is_line())
        return std::abs(cross(a.n, b.n)) <= eps && dot(a.n, b.n) < 0.0 ? ConicKind::plateau_line : ConicKind::line;
    if (a.is_line() || b.is_line()) return ConicKind::parabola;
    if (dist(a.centre, b.centre) <= eps) return ConicKind::plateau_arc;
    if (a.sigma != b.sigma) return ConicKind::ellipse;
    return std::abs(a.R - b.R) <= eps ? ConicKind::line : ConicKind::hyperbola;
}

Root3 refine(const std::array<Site, 3>& sites, Root3 root, bool quad_precision)
{
    if (quad_precision) {
        quad x = root.p.x, y = root.p.y, t = root.t;
        Newton<quad>{sites}.run(x, y, t, 60, quad(std::numeric_limits<quad>::epsilon()) * 4);
        root = make_root(sites, {static_cast<double>(x), static_cast<double>(y)}, static_cast<double>(t), 1.0);
        root.refined_quad = true;
    } else {
        double x = root.p.x, y = root.p.y, t = root.t;
        Newton<double>{sites}.run(x, y, t, 8, 4 * std::numeric_limits<double>::epsilon());
        const bool was_quad = root.refined_quad;
        root = make_root(sites, {x, y}, t, 1.0);
        root.refined_quad = was_quad;
    }
    root.condition = jacobian_condition(sites, root.p);
    return root;
}

Solve3 solve_three(const Site& a, const Site& b, const Site& c, double t_now, const Tolerances& tol, double t_max)
{
    const std::array<Site, 3> sites{a, b, c};
    const std::array<Equation, 3> eq{equation(a), equation(b), equation(c)};
    Solve3 out;

    // Candidate (x, y, t) with the condition of its algebraic solve.
    std::vector<std::pair<Vec3, double>> candidates;

    // The pivot is the circle with the smallest coefficients: subtracting it from a nearly flat arc
    // (huge centre and radius) leaves a well-scaled plane, while pivoting on the flat arc would put
    // its cancellation into the quadratic.
    auto pivot = eq.end();
    double pivot_size = std::numeric_limits<double>::infinity();
    for (std::size_t i = 0; i < 3; ++i) {
        if (!eq[i].quadratic) continue;
        const double size = dot(sites[i].centre, sites[i].centre) + sites[i].R * sites[i].R;
        if (size < pivot_size) {
            pivot = eq.begin() + static_cast<std::ptrdiff_t>(i);
            pivot_size = size;
        }
    }
    if (pivot == eq.end()) {
        // Three lines: a unique solution unless two are parallel.
        const Vec3 &a0 = eq[0].a, &a1 = eq[1].a, &a2 = eq[2].a;
        const double det = dot3(a0, cross3(a1, a2));
        const double scale = norm3(a0) * norm3(a1) * norm3(a2);
        if (std::abs(det) <= 1e-12 * scale) {
            out.degenerate = true;
            return out;
        }
        const Vec3 z = (1.0 / det) * ((-eq[0].b) * cross3(a1, a2) + (-eq[1].b) * cross3(a2, a0) +
                                      (-eq[2].b) * cross3(a0, a1));
        candidates.emplace_back(z, scale / std::abs(det));
    } else {
        // Subtract the pivot circle from the other circles: two linear rows whose solutions form
        // the line z0 + s w; substituting into the pivot gives Q(w) s^2 + B s + C = 0 [KN-01].
        Vec3 rows[2];
        double rhs[2];
        int k = 0;
        for (const Equation& e : eq) {
            if (&e == &*pivot) continue;
            rows[k] = e.quadratic ? e.a - pivot->a : e.a;
            rhs[k] = -(e.quadratic ? e.b - pivot->b : e.b);
            ++k;
        }
        const Vec3 w = cross3(rows[0], rows[1]);
        const double ww = dot3(w, w);
        if (ww <= 1e-24 * dot3(rows[0], rows[0]) * dot3(rows[1], rows[1])) {
            out.degenerate = true;
            return out;
        }
        const Vec3 z0 = (1.0 / ww) * (rhs[0] * cross3(rows[1], w) + rhs[1] * cross3(w, rows[0]));
        const double qa = mink(w, w);
        const double qb = 2.0 * mink(z0, w) + dot3(pivot->a, w);
        const double qc = mink(z0, z0) + dot3(pivot->a, z0) + pivot->b;
        if (std::abs(qa) <= 1e-14 * ww) {
            if (qb != 0.0) candidates.emplace_back(z0 + (-qc / qb) * w, 1.0);
        } else {
            const double size = qb * qb + 4.0 * std::abs(qa * qc);
            double disc = qb * qb - 4.0 * qa * qc;
            if (disc < 0.0 && disc >= -1e-10 * size) disc = 0.0;  // a tangency blurred by rounding
            if (disc >= 0.0) {
                const double sq = std::sqrt(disc);
                const double condition = std::sqrt(size) / std::max(sq, 1e-16 * std::sqrt(size));
                const double q = -0.5 * (qb + std::copysign(sq, qb));
                if (q == 0.0) {
                    candidates.emplace_back(z0, condition);
                } else {
                    candidates.emplace_back(z0 + (q / qa) * w, condition);
                    if (sq > 0.0) candidates.emplace_back(z0 + (qc / q) * w, condition);
                }
            }
        }
    }

    for (const auto& [z, condition] : candidates) {
        // Squaring admits negative radii; reject those before refining.
        bool radii_ok = true;
        for (const Site& s : sites)
            if (!s.is_line() && s.offset_radius(z.t) < -std::max(tol.geom, 1e-8 * s.R)) radii_ok = false;
        if (!radii_ok) continue;

        Root3 r = refine(sites, make_root(sites, {z.x, z.y}, z.t, condition), false);
        r.condition = std::max(r.condition, condition);
        // A root well in the past cannot move into the window by refinement: its double error is
        // about condition * eps. Skip the binary128 pass for it.
        // Likewise a root far beyond t_max.
        if (std::isfinite(r.t) && r.t + 1e-15 * r.condition * (1.0 + std::abs(r.t)) < t_now - tol.time) continue;
        if (std::isfinite(r.t) && r.t - 1e-15 * r.condition * (1.0 + std::abs(r.t)) > t_max) continue;
        if (r.condition > quad_threshold) {
            const double cond = r.condition;
            r = refine(sites, r, true);
            r.condition = std::max(r.condition, cond);
        }

        if (!is_finite(r.p) || !std::isfinite(r.t)) continue;
        if (r.t < t_now - tol.time || r.t > t_max) continue;
        bool ok = true;
        for (const Site& s : sites) {
            if (!s.is_line() && s.offset_radius(r.t) < -tol.geom) ok = false;
            if (std::abs(s.distance(r.p) - r.t) > residual_slack(s, r.p, tol)) ok = false;
        }
        if (!ok) continue;
        const bool duplicate = std::any_of(out.roots.begin(), out.roots.end(), [&](const Root3& q) {
            return dist(q.p, r.p) <= tol.geom && std::abs(q.t - r.t) <= tol.time;
        });
        if (!duplicate) out.roots.push_back(r);
    }
    std::sort(out.roots.begin(), out.roots.end(), [](const Root3& x, const Root3& y) {
        return x.t != y.t ? x.t < y.t : (x.p.x != y.p.x ? x.p.x < y.p.x : x.p.y < y.p.y);
    });
    return out;
}

std::optional<Root3> earliest(const Solve3& s)
{
    if (s.roots.empty()) return std::nullopt;
    return s.roots.front();
}

std::vector<Vec2> offset_intersections(const Site& a, const Site& b, double t)
{
    std::vector<Vec2> out;
    auto radius = [t](const Site& s) { return s.offset_radius(t); };
    if (!a.is_line() && radius(a) < 0.0) return out;
    if (!b.is_line() && radius(b) < 0.0) return out;

    if (a.is_line() && b.is_line()) {
        const double det = cross(a.n, b.n);
        if (det == 0.0) return out;
        const double ka = a.c + t, kb = b.c + t;
        out.push_back({(ka * b.n.y - kb * a.n.y) / det, (a.n.x * kb - b.n.x * ka) / det});
        return out;
    }
    if (a.is_line() != b.is_line()) {
        const Site& line = a.is_line() ? a : b;
        const Site& circ = a.is_line() ? b : a;
        const double rho = radius(circ);
        const double h = dot(line.n, circ.centre) - (line.c + t);
        if (std::abs(h) > rho) {
            if (std::abs(h) - rho > 1e-14 * std::max(1.0, rho)) return out;
        }
        const double w = std::sqrt(std::max(0.0, rho * rho - h * h));
        const Vec2 f = circ.centre - h * line.n;
        out.push_back(f - w * perp(line.n));
        if (w > 0.0) out.push_back(f + w * perp(line.n));
    } else {
        const double ra = radius(a), rb = radius(b);
        const Vec2 e = b.centre - a.centre;
        const double d = norm(e);
        if (d == 0.0) return out;
        const double slack = 1e-14 * std::max(1.0, ra + rb);
        if (d > ra + rb + slack || d < std::abs(ra - rb) - slack) return out;
        // ra - along = (rb - (d - ra)) (rb + (d - ra)) / 2d avoids the cancellation in
        // ra^2 - along^2 when ra is large (a nearly flat arc).
        const double gap = d - ra;
        const double short_by = (rb - gap) * (rb + gap) / (2.0 * d);
        const double along = ra - short_by;
        const double h = std::sqrt(std::max(0.0, short_by * (ra + along)));
        const Vec2 base = a.centre + (along / d) * e;
        out.push_back(base - (h / d) * perp(e));
        if (h > 0.0) out.push_back(base + (h / d) * perp(e));
    }
    std::sort(out.begin(), out.end(), [&](Vec2 p, Vec2 q) { return branch_of(a, b, p) < branch_of(a, b, q); });
    return out;
}

int branch_of(const Site& a, const Site& b, Vec2 p)
{
    if (a.is_line() && b.is_line()) return 0;
    if (a.is_line() != b.is_line()) {
        const Site& line = a.is_line() ? a : b;
        const Site& circ = a.is_line() ? b : a;
        return sign(cross(line.n, p - circ.centre));
    }
    const Vec2 e = b.centre - a.centre;
    if (e == Vec2{}) return 0;
    return sign(cross(e, p - a.centre));
}

std::optional<Vec2> shock_position(const Site& a, const Site& b, double t, int branch)
{
    const std::vector<Vec2> pts = offset_intersections(a, b, t);
    if (pts.size() == 1) return pts.front();
    for (Vec2 p : pts)
        if (branch_of(a, b, p) == branch) return p;
    return std::nullopt;
}

std::optional<Vec2> shock_velocity(const Site& a, const Site& b, Vec2 p)
{
    const Vec2 ga = a.gradient(p), gb = b.gradient(p);
    const double det = cross(ga, gb);
    if (std::abs(det) <= 1e-14) return std::nullopt;
    return Vec2{(gb.y - ga.y) / det, (ga.x - gb.x) / det};
}

std::optional<double> regular_hit_time(Vec2 q0, Vec2 m, const Site& s, double t_now, const Tolerances& tol)
{
    double t;
    if (s.is_line()) {
        const double den = dot(s.n, m) - 1.0;
        if (std::abs(den) <= 1e-15) return std::nullopt;
        t = (s.c - dot(s.n, q0)) / den;
    } else {
        // |q0 - c + t m|^2 = (R - sigma t)^2 with |m| = 1: the t^2 terms cancel.
        const Vec2 v = q0 - s.centre;
        const double den = 2.0 * (dot(v, m) + s.sigma * s.R);
        if (std::abs(den) <= 1e-15) return std::nullopt;
        const double nv = norm(v);
        t = (s.R - nv) * (s.R + nv) / den;
        if (s.offset_radius(t) < -tol.geom) return std::nullopt;
    }
    if (t < t_now - tol.time) return std::nullopt;
    if (std::abs(s.distance(q0 + t * m) - t) > residual_slack(s, q0 + t * m, tol)) return std::nullopt;
    return t;
}

std::vector<Contact> contacts(const Site& a, const Site& b, double t_now, const Tolerances& tol)
{
    std::vector<Contact> out;
    auto accept = [&](double t) { return t >= t_now - tol.time; };

    if (a.is_line() && b.is_line()) {
        if (std::abs(cross(a.n, b.n)) > tol.ang || dot(a.n, b.n) >= 0.0) return out;
        // n_a.x = c_a + t and n_b.x = c_b + t with n_b = -n_a.
        const double t = -0.5 * (a.c + b.c);
        if (!accept(t)) return out;
        const Vec2 p = (a.c + t) * a.n;
        out.push_back({t, p, {p - t * a.n, p - t * b.n}, true});
        return out;
    }

    if (!a.is_line() && !b.is_line() && dist(a.centre, b.centre) <= tol.geom) {
        // Concentric circles meet along a whole circle, and only when one is convex and one concave.
        if (a.sigma == b.sigma) return out;
        const double r = (b.sigma * b.R - a.sigma * a.R) / (b.sigma - a.sigma);
        const double t = a.sigma * (a.R - r);
        if (r < 0.0 || !accept(t)) return out;
        const Vec2 foot_dir{1.0, 0.0};
        out.push_back({t, a.centre, {a.centre + a.R * foot_dir, a.centre + b.R * foot_dir}, true});
        return out;
    }

    // Common normal: the axis x(s) = o + s u. Along it each distance function is
    // d(s) = alpha + beta s + gamma |s - k| (lines: gamma = 0; circles: alpha = sigma R, gamma = -sigma).
    Vec2 o, u;
    if (a.is_line() || b.is_line()) {
        const Site& line = a.is_line() ? a : b;
        const Site& circ = a.is_line() ? b : a;
        o = circ.centre;
        u = line.n;
    } else {
        o = a.centre;
        u = unit(b.centre - a.centre);
    }
    struct Profile {
        double alpha, beta, gamma, k;
    };
    auto profile = [&](const Site& s) -> Profile {
        if (s.is_line()) return {dot(s.n, o) - s.c, dot(s.n, u), 0.0, 0.0};
        return {s.sigma * s.R, 0.0, -static_cast<double>(s.sigma), dot(s.centre - o, u)};
    };
    const Profile pa = profile(a), pb = profile(b);

    std::vector<double> kinks;
    if (pa.gamma != 0.0) kinks.push_back(pa.k);
    if (pb.gamma != 0.0) kinks.push_back(pb.k);
    std::sort(kinks.begin(), kinks.end());
    std::vector<std::pair<double, double>> intervals;
    double lo = -std::numeric_limits<double>::infinity();
    for (double k : kinks) {
        intervals.emplace_back(lo, k);
        lo = k;
    }
    intervals.emplace_back(lo, std::numeric_limits<double>::infinity());

    for (const auto& [s0, s1] : intervals) {
        const double mid = std::isinf(s0) ? s1 - 1.0 : std::isinf(s1) ? s0 + 1.0 : 0.5 * (s0 + s1);
        auto linear = [mid](const Profile& p) {
            const double sg = mid >= p.k ? 1.0 : -1.0;
            return std::pair{p.alpha - p.gamma * sg * p.k, p.beta + p.gamma * sg};
        };
        const auto [a0, a1] = linear(pa);
        const auto [b0, b1] = linear(pb);
        if (std::abs(a1 - b1) <= 1e-15) continue;
        const double s = (b0 - a0) / (a1 - b1);
        const double slack = 1e-12 * (1.0 + std::abs(s));
        if (s < s0 - slack || s > s1 + slack) continue;
        const double t = a0 + a1 * s;
        const Vec2 p = o + s * u;
        if (!accept(t)) continue;
        if ((!a.is_line() && a.offset_radius(t) < -tol.geom) || (!b.is_line() && b.offset_radius(t) < -tol.geom))
            continue;
        // The fronts must approach each other: opposed normals at the meeting point.
        const Vec2 ga = a.gradient(p), gb = b.gradient(p);
        if (ga == Vec2{} || gb == Vec2{} || dot(ga, gb) > -1.0 + 1e-9) continue;
        if (std::abs(a.distance(p) - t) > residual_slack(a, p, tol) ||
            std::abs(b.distance(p) - t) > residual_slack(b, p, tol))
            continue;
        const bool duplicate = std::any_of(out.begin(), out.end(), [&](const Contact& q) {
            return dist(q.p, p) <= tol.geom && std::abs(q.t - t) <= tol.time;
        });
        if (!duplicate) out.push_back({t, p, {a.foot(p, t), b.foot(p, t)}, false});
    }
    std::sort(out.begin(), out.end(), [](const Contact& x, const Contact& y) { return x.t < y.t; });
    return out;
}

} // namespace wfmat
