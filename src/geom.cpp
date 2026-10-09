#include "wfmat/geom.hpp"

#include <algorithm>
#include <array>

namespace wfmat {

namespace {

constexpr double two_pi = 2.0 * std::numbers::pi;

// Start angle and width of an arc as a counter-clockwise angular interval.
std::pair<double, double> ccw_interval(const Segment& s)
{
    return s.sweep > 0.0 ? std::pair{s.phi0, s.sweep} : std::pair{s.phi0 + s.sweep, -s.sweep};
}

// Measure of the intersection of the counter-clockwise intervals [a0, a0 + wa] and [b0, b0 + wb].
double circular_overlap(double a0, double wa, double b0, double wb)
{
    const double d = wrap_2pi(b0 - a0);
    const double first = std::max(0.0, std::min(wa, d + wb) - d);
    const double second = std::max(0.0, std::min(wa, d - two_pi + wb));
    return first + second;
}

} // namespace

double wrap_2pi(double angle)
{
    double r = std::fmod(angle, two_pi);
    if (r < 0.0) r += two_pi;
    return r >= two_pi ? 0.0 : r;
}

void Box::add(Vec2 p)
{
    lo = {std::min(lo.x, p.x), std::min(lo.y, p.y)};
    hi = {std::max(hi.x, p.x), std::max(hi.y, p.y)};
}

void Box::add(const Box& b)
{
    if (b.empty()) return;
    add(b.lo);
    add(b.hi);
}

Segment make_segment(Vec2 a, Vec2 b, double bulge)
{
    Segment s;
    s.a = a;
    s.b = b;
    s.bulge = bulge;
    if (bulge == 0.0) return s;
    const Vec2 d = b - a;
    const double len = norm(d);
    s.sweep = 4.0 * std::atan(bulge);
    // The centre lies on the chord's perpendicular bisector at signed distance
    // (len / 2) cot(sweep / 2) = (len / 2) (1 - bulge^2) / (2 bulge) to the left of the chord.
    s.c = 0.5 * (a + b) + (0.5 * (1.0 - bulge * bulge) / (2.0 * bulge)) * perp(d);
    s.R = len * (1.0 + bulge * bulge) / (4.0 * std::abs(bulge));
    s.phi0 = std::atan2(a.y - s.c.y, a.x - s.c.x);
    return s;
}

double Segment::length() const
{
    return is_arc() ? R * std::abs(sweep) : dist(a, b);
}

Vec2 Segment::point(double s) const
{
    if (s <= 0.0) return a;
    if (s >= 1.0) return b;
    if (!is_arc()) return a + s * (b - a);
    return c + R * polar(phi0 + s * sweep);
}

Vec2 Segment::start_tangent() const
{
    if (!is_arc()) return unit(b - a);
    return (sweep > 0.0 ? 1.0 : -1.0) * unit(perp(a - c));
}

Vec2 Segment::end_tangent() const
{
    if (!is_arc()) return unit(b - a);
    return (sweep > 0.0 ? 1.0 : -1.0) * unit(perp(b - c));
}

Box Segment::bbox() const
{
    Box box;
    box.add(a);
    box.add(b);
    if (is_arc()) {
        static constexpr std::array<Vec2, 4> axes{{{1, 0}, {0, 1}, {-1, 0}, {0, -1}}};
        for (Vec2 u : axes)
            if (angle_in_range(c + u)) box.add(c + R * u);
    }
    return box;
}

bool Segment::angle_in_range(Vec2 p, double slack) const
{
    const double ang = std::atan2(p.y - c.y, p.x - c.x);
    const double d = sweep > 0.0 ? wrap_2pi(ang - phi0) : wrap_2pi(phi0 - ang);
    return d <= std::abs(sweep) + slack || d >= two_pi - slack;
}

Vec2 Segment::closest_point(Vec2 p) const
{
    if (!is_arc()) {
        const Vec2 d = b - a;
        const double dd = dot(d, d);
        const double u = dd > 0.0 ? std::clamp(dot(p - a, d) / dd, 0.0, 1.0) : 0.0;
        return u <= 0.0 ? a : u >= 1.0 ? b : a + u * d;
    }
    const Vec2 v = p - c;
    const double r = norm(v);
    if (r == 0.0) return midpoint();
    if (angle_in_range(p)) return c + (R / r) * v;
    return dist(p, a) <= dist(p, b) ? a : b;
}

double Segment::cap_area() const
{
    return is_arc() ? 0.5 * R * R * (sweep - std::sin(sweep)) : 0.0;
}

Segment Segment::reversed() const
{
    return make_segment(b, a, -bulge);
}

Intersections intersect(const Segment& s, const Segment& t, double eps)
{
    Intersections out;
    auto add = [&](Vec2 p) {
        for (Vec2 q : out.points)
            if (dist(p, q) <= eps) return;
        out.points.push_back(p);
    };
    auto add_if_on_both = [&](Vec2 p) {
        if (s.distance(p) <= eps && t.distance(p) <= eps) add(p);
    };

    // End points touching the other segment cover T-junctions, shared vertices and near misses.
    for (Vec2 e : {s.a, s.b})
        if (t.distance(e) <= eps) add(e);
    for (Vec2 e : {t.a, t.b})
        if (s.distance(e) <= eps) add(e);

    if (!s.is_arc() && !t.is_arc()) {
        const Vec2 d1 = s.b - s.a, d2 = t.b - t.a;
        const double l1 = norm(d1), l2 = norm(d2);
        const double den = cross(d1, d2);
        if (std::abs(den) > 1e-14 * l1 * l2) {
            add_if_on_both(s.a + (cross(t.a - s.a, d2) / den) * d1);
        } else if (std::abs(cross(d1, t.a - s.a)) / l1 <= eps) {
            const Vec2 u = d1 / l1;
            const double u0 = dot(t.a - s.a, u), u1 = dot(t.b - s.a, u);
            if (std::min(l1, std::max(u0, u1)) - std::max(0.0, std::min(u0, u1)) > eps) out.overlap = true;
        }
        return out;
    }

    if (s.is_arc() != t.is_arc()) {
        const Segment& line = s.is_arc() ? t : s;
        const Segment& arc = s.is_arc() ? s : t;
        const Vec2 u = unit(line.b - line.a);
        const Vec2 f = line.a + dot(arc.c - line.a, u) * u;
        const double h = dist(arc.c, f);
        if (h <= arc.R + eps) {
            const double w = std::sqrt(std::max(0.0, arc.R * arc.R - h * h));
            add_if_on_both(f + w * u);
            add_if_on_both(f - w * u);
        }
        return out;
    }

    const Vec2 e = t.c - s.c;
    const double dcc = norm(e);
    if (dcc <= eps) {
        if (std::abs(s.R - t.R) <= eps) {
            const auto [a0, wa] = ccw_interval(s);
            const auto [b0, wb] = ccw_interval(t);
            if (circular_overlap(a0, wa, b0, wb) * s.R > eps) out.overlap = true;
        }
        return out;
    }
    if (dcc > s.R + t.R + eps || dcc < std::abs(s.R - t.R) - eps) return out;
    const double along = (dcc * dcc + s.R * s.R - t.R * t.R) / (2.0 * dcc);
    const double h = std::sqrt(std::max(0.0, s.R * s.R - along * along));
    const Vec2 base = s.c + (along / dcc) * e;
    const Vec2 off = (h / dcc) * perp(e);
    add_if_on_both(base + off);
    add_if_on_both(base - off);
    return out;
}

} // namespace wfmat
