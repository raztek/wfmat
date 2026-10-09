// MAT checks shared by the engine tests: the [VER-04] property checks and the [VER-03] comparison
// with Boost.Polygon's Voronoi diagram of the boundary segments.
#pragma once

#include <algorithm>
#include <climits>
#include <cmath>
#include <cstdint>
#include <iterator>
#include <numbers>
#include <numeric>
#include <random>
#include <string>
#include <utility>
#include <vector>

#include <boost/polygon/point_data.hpp>
#include <boost/polygon/segment_data.hpp>
#include <boost/polygon/voronoi.hpp>

#include "wfmat/mat.hpp"
#include "wfmat/validate.hpp"

namespace wfmat::test {

// [VER-04] Residual of the radius function, at least two feet inside every edge, a tree whose
// vertex degrees match their kinds. Tolerances are relative to the half-diagonal L. Returns an
// empty string when everything holds, otherwise the first failure.
inline std::string check_mat(const PreparedRegion& region, const MedialAxis& mat, double rel_tol = 1e-9,
                             int samples = 8)
{
    const Transform& xf = region.xf;
    const double tol = rel_tol * xf.scale;
    const auto vertices = mat.vertices();
    const auto edges = mat.edges();
    auto where = [](const char* what, std::size_t i) { return std::string(what) + " " + std::to_string(i) + ": "; };

    if (vertices.size() != edges.size() + 1)
        return "not a tree: " + std::to_string(vertices.size()) + " vertices, " + std::to_string(edges.size()) + " edges";
    std::vector<std::size_t> parent(vertices.size());
    std::iota(parent.begin(), parent.end(), std::size_t{0});
    auto find = [&](std::size_t x) {
        while (parent[x] != x) x = parent[x] = parent[parent[x]];
        return x;
    };
    for (const MatEdge& e : edges) {
        const std::size_t a = find(e.v0()), b = find(e.v1());
        if (a == b) return "the MAT has a cycle";
        parent[a] = b;
    }

    for (std::size_t i = 0; i < vertices.size(); ++i) {
        const MatVertex& v = vertices[i];
        const std::size_t deg = v.edges.size();
        const bool deg_ok = v.kind == VertexKind::corner || v.kind == VertexKind::curvature_end ? deg == 1
                            : v.kind == VertexKind::junction                                     ? deg >= 3
                                                                                                 : deg == 2;
        if (!deg_ok) return where("vertex", i) + "degree " + std::to_string(deg) + " does not match its kind";
        const DiskCheck c = check_disk(region, xf.to_unit(v.p), xf.length_to_unit(v.r), 1e-7);
        if (c.residual * xf.scale > tol) return where("vertex", i) + "radius residual " + std::to_string(c.residual);
    }

    for (std::size_t i = 0; i < edges.size(); ++i) {
        const MatEdge& e = edges[i];
        if (e.r1() < e.r0() - tol) return where("edge", i) + "radius decreases";
        if (dist(e.point_at(0.0), vertices[e.v0()].p) > tol || dist(e.point_at(1.0), vertices[e.v1()].p) > tol)
            return where("edge", i) + "end points do not match its vertices";
        for (int k = 0; k < samples; ++k) {
            const double u = (k + 0.5) / samples;
            const Vec2 p = e.point_at(u);
            const DiskCheck c = check_disk(region, xf.to_unit(p), xf.length_to_unit(e.radius_at(u)), 1e-7);
            if (!c.inside) return where("edge", i) + "sample outside the region";
            if (c.residual * xf.scale > tol)
                return where("edge", i) + "radius residual " + std::to_string(c.residual) + " at u = " + std::to_string(u);
            // The disk touches the boundary at the edge's two feet, one on each site.
            const auto [fl, fr] = e.feet_at(u);
            for (Vec2 f : {fl, fr}) {
                const double on_boundary = boundary_distance(region.segments, xf.to_unit(f), 1e-7).distance;
                if (on_boundary * xf.scale > tol || std::abs(dist(p, f) - e.radius_at(u)) > tol)
                    return where("edge", i) + "a foot is not a contact point at u = " + std::to_string(u);
            }
            // Feet that are far apart must also be found by brute force (near a flat corner they
            // merge within the foot tolerance).
            if (c.feet < 2 && dist(fl, fr) > 1e-6 * xf.scale)
                return where("edge", i) + "fewer than two feet at u = " + std::to_string(u);
        }
    }
    return {};
}

using IPoint = std::pair<std::int32_t, std::int32_t>;

// [VER-03] Boost.Polygon's Voronoi diagram of the polygon's segments, restricted to the interior,
// without secondary edges (the spokes from reflex corners), against our MAT with its E3 vertices
// (minima of r inside a Voronoi edge) contracted. Vertices are matched by position within tol (in
// caller units) and the edge sets must then agree. Returns an empty string on success.
inline std::string compare_with_voronoi(const std::vector<IPoint>& polygon, const PreparedRegion& region,
                                        const MedialAxis& mat, double tol)
{
    namespace bp = boost::polygon;
    std::vector<bp::segment_data<std::int32_t>> segs;
    const std::size_t n = polygon.size();
    for (std::size_t i = 0; i < n; ++i) {
        const IPoint a = polygon[i], b = polygon[(i + 1) % n];
        segs.emplace_back(bp::point_data<std::int32_t>(a.first, a.second), bp::point_data<std::int32_t>(b.first, b.second));
    }
    bp::voronoi_diagram<double> vd;
    bp::construct_voronoi(segs.begin(), segs.end(), &vd);

    enum class Where { outside, boundary, inside };
    const auto& bverts = vd.vertices();
    std::vector<Where> where(bverts.size());
    for (std::size_t i = 0; i < bverts.size(); ++i) {
        const Vec2 p{bverts[i].x(), bverts[i].y()};
        bool on_corner = false;
        for (const IPoint& q : polygon)
            if (dist(p, {double(q.first), double(q.second)}) <= tol) on_corner = true;
        where[i] = on_corner ? Where::boundary
                   : contains(region, region.xf.to_unit(p), region.xf.length_to_unit(tol)) ? Where::inside
                                                                                           : Where::outside;
    }
    auto index = [&](const bp::voronoi_vertex<double>* v) { return static_cast<std::size_t>(v - &bverts[0]); };

    std::vector<std::pair<std::size_t, std::size_t>> bedges;
    std::vector<int> bdegree(bverts.size(), 0);
    for (const auto& e : vd.edges()) {
        if (e.is_infinite() || e.is_secondary() || &e > e.twin()) continue;
        const std::size_t i0 = index(e.vertex0()), i1 = index(e.vertex1());
        const Where w0 = where[i0], w1 = where[i1];
        bool interior = false;
        if (w0 == Where::outside || w1 == Where::outside)
            interior = false;
        else if (w0 == Where::inside || w1 == Where::inside)
            interior = true;
        else if (!e.is_curved()) {
            const Vec2 mid{0.5 * (bverts[i0].x() + bverts[i1].x()), 0.5 * (bverts[i0].y() + bverts[i1].y())};
            interior = contains(region, region.xf.to_unit(mid), region.xf.length_to_unit(tol));
        }
        if (!interior) continue;
        bedges.emplace_back(std::min(i0, i1), std::max(i0, i1));
        ++bdegree[i0];
        ++bdegree[i1];
    }
    std::vector<std::size_t> bused;
    for (std::size_t i = 0; i < bverts.size(); ++i)
        if (bdegree[i] > 0) bused.push_back(i);

    // Our graph with the E3 vertices contracted.
    const auto verts = mat.vertices();
    const auto edges = mat.edges();
    std::vector<std::size_t> ours;
    for (std::size_t i = 0; i < verts.size(); ++i)
        if (verts[i].kind != VertexKind::extremum_min) ours.push_back(i);
    if (ours.size() != bused.size())
        return "vertex count " + std::to_string(ours.size()) + " vs Voronoi " + std::to_string(bused.size());

    std::vector<std::size_t> to_boost(verts.size(), SIZE_MAX);
    std::vector<bool> taken(bverts.size(), false);
    for (std::size_t i : ours) {
        double best = INFINITY;
        std::size_t arg = SIZE_MAX;
        for (std::size_t j : bused) {
            const double d = dist(verts[i].p, {bverts[j].x(), bverts[j].y()});
            if (d < best) {
                best = d;
                arg = j;
            }
        }
        if (best > tol)
            return "vertex " + std::to_string(i) + " at (" + std::to_string(verts[i].p.x) + ", " +
                   std::to_string(verts[i].p.y) + ") has no Voronoi vertex within tolerance (nearest " +
                   std::to_string(best) + ")";
        if (taken[arg]) return "two vertices match the same Voronoi vertex";
        taken[arg] = true;
        to_boost[i] = arg;
    }

    std::vector<std::pair<std::size_t, std::size_t>> mine;
    auto other_end = [&](std::size_t edge, std::size_t v) {
        return edges[edge].v0() == v ? edges[edge].v1() : edges[edge].v0();
    };
    for (std::size_t k = 0; k < edges.size(); ++k) {
        const std::size_t a = edges[k].v0(), b = edges[k].v1();
        if (verts[a].kind == VertexKind::extremum_min || verts[b].kind == VertexKind::extremum_min) continue;
        mine.emplace_back(std::min(to_boost[a], to_boost[b]), std::max(to_boost[a], to_boost[b]));
    }
    for (std::size_t i = 0; i < verts.size(); ++i) {
        if (verts[i].kind != VertexKind::extremum_min) continue;
        const auto& es = verts[i].edges;
        if (es.size() != 2) return "an E3 vertex without two edges";
        const std::size_t a = other_end(es[0], i), b = other_end(es[1], i);
        if (to_boost[a] == SIZE_MAX || to_boost[b] == SIZE_MAX) return "an edge joins two E3 vertices";
        mine.emplace_back(std::min(to_boost[a], to_boost[b]), std::max(to_boost[a], to_boost[b]));
    }
    std::sort(mine.begin(), mine.end());
    std::sort(bedges.begin(), bedges.end());
    if (mine != bedges) {
        std::vector<std::pair<std::size_t, std::size_t>> only_mine, only_boost;
        std::set_difference(mine.begin(), mine.end(), bedges.begin(), bedges.end(), std::back_inserter(only_mine));
        std::set_difference(bedges.begin(), bedges.end(), mine.begin(), mine.end(), std::back_inserter(only_boost));
        auto text = [&](const std::pair<std::size_t, std::size_t>& e) {
            // Appended piecewise: GCC 12.3 reports a false -Wrestrict for "(" + std::string&&.
            std::string s = "(";
            s += std::to_string(bverts[e.first].x());
            s += ", ";
            s += std::to_string(bverts[e.first].y());
            s += ")-(";
            s += std::to_string(bverts[e.second].x());
            s += ", ";
            s += std::to_string(bverts[e.second].y());
            s += ")";
            return s;
        };
        return "edge sets differ: " + std::to_string(mine.size()) + " edges vs Voronoi " + std::to_string(bedges.size()) +
               (only_mine.empty() ? "" : "; only ours " + text(only_mine.front())) +
               (only_boost.empty() ? "" : "; only Voronoi " + text(only_boost.front()));
    }
    return {};
}

// [VER-05] A random star-shaped polygon with integer coordinates in [-R, R]^2, counter-clockwise,
// without collinear neighbours or parallel edges (those make plateaus and clusters, M5).
inline std::vector<IPoint> random_star_polygon(std::mt19937_64& rng, int n, std::int32_t R = 1'000'000)
{
    std::uniform_real_distribution<double> angle(0.0, 2.0 * std::numbers::pi), radius(0.2, 1.0);
    for (;;) {
        std::vector<double> a(n);
        for (double& x : a) x = angle(rng);
        std::sort(a.begin(), a.end());
        bool spaced = true;
        for (int i = 0; i < n; ++i) {
            const double gap = i + 1 < n ? a[i + 1] - a[i] : a[0] + 2.0 * std::numbers::pi - a[i];
            if (gap < 0.2 / n || gap > std::numbers::pi * 0.95) spaced = false;
        }
        if (!spaced) continue;
        std::vector<IPoint> p(n);
        for (int i = 0; i < n; ++i) {
            const double r = R * radius(rng);
            p[i] = {static_cast<std::int32_t>(std::lround(r * std::cos(a[i]))),
                    static_cast<std::int32_t>(std::lround(r * std::sin(a[i])))};
        }
        bool ok = true;
        for (int i = 0; i < n && ok; ++i) {
            const auto d = [&](int k) {
                const IPoint u = p[k % n], v = p[(k + 1) % n];
                return std::pair<std::int64_t, std::int64_t>{std::int64_t(v.first) - u.first,
                                                             std::int64_t(v.second) - u.second};
            };
            const auto di = d(i);
            for (int j = i + 1; j < n; ++j) {
                const auto dj = d(j);
                if (di.first * dj.second - di.second * dj.first == 0) ok = false;
            }
        }
        if (ok) return p;
    }
}

inline Region to_region(const std::vector<IPoint>& p)
{
    Region region;
    for (const IPoint& q : p) region.outer.vertices.push_back({{double(q.first), double(q.second)}, 0.0});
    return region;
}

} // namespace wfmat::test
