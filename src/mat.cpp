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
    if (d_.kind == ConicKind::plateau_line) return unit(d_.p1 - d_.p0);
    if (d_.kind == ConicKind::plateau_arc) return d_.branch * unit(perp(plateau_point(u) - d_.a.centre));
    const Vec2 p = unit_point(t_at(u));
    if (auto v = shock_velocity(d_.a, d_.b, p)) return unit(*v);
    const double h = 1e-6;
    return unit(unit_point(t_at(std::min(u + h, 1.0))) - unit_point(t_at(std::max(u - h, 0.0))));
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

} // namespace wfmat
