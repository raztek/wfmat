#include "wfmat/mat.hpp"

#include <algorithm>
#include <cmath>

#include "engine.hpp"

namespace wfmat {

namespace {

Site to_caller(const Site& s, const Transform& xf)
{
    if (s.is_line()) return Site::line(s.n, xf.length_from_unit(s.c) + dot(s.n, xf.centre));
    return Site::circle(xf.from_unit(s.centre), xf.length_from_unit(s.R), s.sigma);
}

} // namespace

Vec2 MatEdge::unit_point(double t) const
{
    if (t <= d_.t0) return d_.p0;
    if (t >= d_.t1) return d_.p1;
    if (auto p = shock_position(d_.a, d_.b, t, d_.branch)) return *p;
    // Rounding just past a tangency (the apex of an E3 edge): the nearer end point.
    return t - d_.t0 < d_.t1 - t ? d_.p0 : d_.p1;
}

// [OUT-02] A plateau edge is a segment, or an arc about the common centre of its two sites.
Vec2 MatEdge::plateau_point(double u) const
{
    if (d_.kind == ConicKind::plateau_line) return d_.p0 + u * (d_.p1 - d_.p0);
    const Vec2 c = d_.a.centre;
    const double a0 = std::atan2(d_.p0.y - c.y, d_.p0.x - c.x);
    const double a1 = std::atan2(d_.p1.y - c.y, d_.p1.x - c.x);
    const double sweep = d_.branch * wrap_2pi(d_.branch * (a1 - a0));
    return c + dist(d_.p0, c) * polar(a0 + u * sweep);
}

Vec2 MatEdge::point(double r) const
{
    if (is_plateau()) return d_.xf.from_unit(d_.p0);
    return d_.xf.from_unit(unit_point(d_.xf.length_to_unit(r)));
}

Vec2 MatEdge::point_at(double u) const
{
    return d_.xf.from_unit(unit_point_at(u));
}

double MatEdge::radius_at(double u) const
{
    return d_.xf.length_from_unit(t_at(u));
}

Vec2 MatEdge::tangent_at(double u) const
{
    return unit_tangent_at(u);  // a similarity keeps directions
}

// Where the shock velocity is undefined (the apex of an E3 edge, where the two fronts are
// tangent), the bisector's tangent is still normal to the difference of the two gradients.
Vec2 MatEdge::unit_tangent_at(double u) const
{
    if (d_.kind == ConicKind::plateau_line) return unit(d_.p1 - d_.p0);
    if (d_.kind == ConicKind::plateau_arc) return d_.branch * unit(perp(plateau_point(u) - d_.a.centre));
    const Vec2 p = unit_point(t_at(u));
    if (auto v = shock_velocity(d_.a, d_.b, p)) return unit(*v);
    const double h = 1e-6;
    const Vec2 chord = unit_point(t_at(std::min(u + h, 1.0))) - unit_point(t_at(std::max(u - h, 0.0)));
    const Vec2 n = d_.a.gradient(p) - d_.b.gradient(p);
    if (norm(n) == 0.0) return unit(chord);
    const Vec2 tan = unit(perp(n));
    return dot(tan, chord) >= 0.0 ? tan : -1.0 * tan;
}

// [OUT-05] Each piece is fixed by its end points, their tangents and one interior point: the
// control point is where the end tangents meet, and the weight follows from the interior point's
// barycentric coordinates (a, b, c) in the control triangle, w = b / (2 sqrt(a c)).
std::vector<RationalQuadBezier> MatEdge::bezier() const
{
    std::vector<RationalQuadBezier> out;
    const Transform& xf = d_.xf;
    auto straight = [&](double ua, double ub) {
        const Vec2 a = unit_point_at(ua), b = unit_point_at(ub);
        out.push_back({xf.from_unit(a), xf.from_unit(0.5 * (a + b)), xf.from_unit(b), 1.0, ua, ub});
    };
    if (d_.kind == ConicKind::line || d_.kind == ConicKind::plateau_line) {
        straight(0.0, 1.0);
        return out;
    }
    auto piece = [&](auto&& self, double ua, double ub, int depth) -> void {
        const Vec2 a = unit_point_at(ua), b = unit_point_at(ub);
        const Vec2 ta = unit_tangent_at(ua), tb = unit_tangent_at(ub);
        const double sin_turn = cross(ta, tb), cos_turn = dot(ta, tb);
        const double um = 0.5 * (ua + ub);
        if (std::abs(sin_turn) <= 1e-12 * (1.0 + std::abs(cos_turn)) && cos_turn > 0.0) {
            straight(ua, ub);
            return;
        }
        auto split = [&] {
            self(self, ua, um, depth + 1);
            self(self, um, ub, depth + 1);
        };
        if (cos_turn < 0.0 && depth < 30) return split();  // turning by a right angle or more
        const double lambda = cross(b - a, tb) / sin_turn, mu = cross(ta, b - a) / sin_turn;
        const Vec2 c = a + lambda * ta;
        const Vec2 m = unit_point_at(um);
        // m = a + beta (c - a) + gamma (b - a)
        const double det = cross(c - a, b - a);
        const double beta = cross(m - a, b - a) / det, gamma = cross(c - a, m - a) / det;
        const double alpha = 1.0 - beta - gamma;
        if ((lambda <= 0.0 || mu <= 0.0 || alpha <= 0.0 || beta <= 0.0 || gamma <= 0.0) && depth < 30) return split();
        const double w = beta > 0.0 && alpha * gamma > 0.0 ? beta / (2.0 * std::sqrt(alpha * gamma)) : 1.0;
        out.push_back({xf.from_unit(a), xf.from_unit(c), xf.from_unit(b), w, ua, ub});
    };
    piece(piece, 0.0, 1.0, 0);
    return out;
}

void MatEdge::flatten(double chord_tol, std::vector<Vec2>& out) const
{
    const Vec2 start = point_at(0.0), end = point_at(1.0);
    out.push_back(start);
    if (d_.kind == ConicKind::line || d_.kind == ConicKind::plateau_line) {
        out.push_back(end);
        return;
    }
    auto deviation = [](Vec2 p, Vec2 a, Vec2 b) {
        const Vec2 d = b - a;
        const double len = norm(d);
        return len > 0.0 ? std::abs(cross(d, p - a)) / len : dist(p, a);
    };
    // At least four pieces, so that a symmetric arc cannot hide behind a chord through its middle;
    // u is not arc length, so the quarter points are tested as well as the middle, against a
    // margin for the deviation between them.
    auto step = [&](auto&& self, double ua, double ub, Vec2 a, Vec2 b, int depth) -> void {
        const double um = 0.5 * (ua + ub);
        const Vec2 m = point_at(um);
        const double dev = std::max({deviation(m, a, b), deviation(point_at(0.5 * (ua + um)), a, b),
                                     deviation(point_at(0.5 * (um + ub)), a, b)});
        if (depth < 2 || (dev > 0.7 * chord_tol && depth < 40)) {
            self(self, ua, um, a, m, depth + 1);
            self(self, um, ub, m, b, depth + 1);
        } else {
            out.push_back(b);
        }
    };
    step(step, 0.0, 1.0, start, end, 0);
}

std::pair<Vec2, Vec2> MatEdge::feet_at(double u) const
{
    const double t = t_at(u);
    const Vec2 p = unit_point_at(u);
    return {d_.xf.from_unit(d_.a.foot(p, t)), d_.xf.from_unit(d_.b.foot(p, t))};
}

Result<MedialAxis> compute_mat(const Region& region, const Options& options)
{
    auto prepared = prepare(region, options);
    if (!prepared) return tl::unexpected(prepared.error());
    return compute_mat(*prepared, options);
}

Result<MedialAxis> compute_mat(const PreparedRegion& region, const Options& options)
{
    detail::Engine engine(region, options);
    if (auto r = engine.run(); !r) return tl::unexpected(r.error());
    const Transform& xf = region.xf;

    std::vector<MatSite> sites = std::move(engine.site_info);
    for (std::size_t i = 0; i < sites.size(); ++i) sites[i].geometry = to_caller(engine.sites[i], xf);

    std::vector<MatVertex> vertices;
    vertices.reserve(engine.mat_vertices.size());
    for (const auto& v : engine.mat_vertices) {
        MatVertex out{xf.from_unit(v.p), xf.length_from_unit(v.t), v.kind, v.contacts, {}};
        for (MatContact& c : out.contacts) c.foot = xf.from_unit(c.foot);
        vertices.push_back(std::move(out));
    }

    std::vector<MatEdge> edges;
    edges.reserve(engine.mat_edges.size());
    for (const MatEdge::Data& d : engine.mat_edges) {
        const auto id = static_cast<EdgeId>(edges.size());
        edges.emplace_back(d);
        vertices[d.v0].edges.push_back(id);
        vertices[d.v1].edges.push_back(id);
    }

    // [OUT-01] Incident edges counter-clockwise, by the direction in which each leaves the vertex.
    for (VertexId vi = 0; vi < vertices.size(); ++vi) {
        MatVertex& v = vertices[vi];
        auto angle = [&](EdgeId e) {
            const MatEdge& edge = edges[e];
            const Vec2 q = edge.v0() == vi ? edge.point_at(1e-3) : edge.point_at(1.0 - 1e-3);
            return std::atan2(q.y - v.p.y, q.x - v.p.x);
        };
        std::sort(v.edges.begin(), v.edges.end(), [&](EdgeId x, EdgeId y) {
            const double ax = angle(x), ay = angle(y);
            return ax != ay ? ax < ay : x < y;
        });
    }

    RunStats stats = engine.stats;
    stats.max_residual = xf.length_from_unit(stats.max_residual);
    return MedialAxis(std::move(vertices), std::move(edges), std::move(sites), stats, xf);
}

namespace {

// [ALG-06] One piece of the front at t as bulge-form vertices in caller units: its start point
// and the bulge to the next piece's start; a whole circle is two half arcs.
void append_piece(const Site& s, const detail::Engine::FrontPiece& piece, double t, const Transform& xf, double eps,
                  std::vector<BulgeVertex>& out)
{
    if (s.is_line()) {
        if (dist(piece.a, piece.b) > eps) out.push_back({xf.from_unit(piece.a), 0.0});
        return;
    }
    const double r = s.offset_radius(t);
    if (r <= eps) return;
    const Vec2 c = s.centre;
    if (piece.full) {
        const Vec2 q = c + (c - piece.a);
        out.push_back({xf.from_unit(piece.a), static_cast<double>(s.sigma)});
        out.push_back({xf.from_unit(q), static_cast<double>(s.sigma)});
        return;
    }
    // Convex arcs run counter-clockwise about the centre, concave arcs and points clockwise.
    auto angle = [&](Vec2 q) { return std::atan2(q.y - c.y, q.x - c.x); };
    const double sweep = s.sigma * wrap_2pi(s.sigma * (angle(piece.b) - angle(piece.a)));
    if (std::abs(sweep) * r <= eps) return;
    out.push_back({xf.from_unit(piece.a), std::tan(sweep / 4.0)});
}

} // namespace

Result<std::vector<Loop>> inward_offset(const Region& region, double distance, const Options& options)
{
    auto prepared = prepare(region, options);
    if (!prepared) return tl::unexpected(prepared.error());
    return inward_offset(*prepared, distance, options);
}

Result<std::vector<Loop>> inward_offset(const PreparedRegion& region, double distance, const Options& options)
{
    if (!(distance >= 0.0) || !std::isfinite(distance))
        return make_error(ErrorCode::invalid_input, "ALG-06", "the offset distance must be finite and not negative");
    const double t = region.xf.length_to_unit(distance);
    detail::Engine engine(region, options);
    if (auto r = engine.run(t); !r) return tl::unexpected(r.error());
    const double eps = options.tol.len;
    std::vector<Loop> loops;
    for (const auto& front : engine.front_at(t)) {
        Loop loop;
        for (const auto& piece : front) append_piece(engine.sites[piece.site], piece, t, region.xf, eps, loop.vertices);
        if (!loop.vertices.empty()) loops.push_back(std::move(loop));
    }
    return loops;
}

} // namespace wfmat
