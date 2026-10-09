// MAT checks shared by the engine tests: the [VER-04] property checks and the [VER-03] comparison
// with Boost.Polygon's Voronoi diagram of the boundary segments.
#pragma once

#include <algorithm>
#include <climits>
#include <cmath>
#include <cstdio>
#include <cstdint>
#include <iterator>
#include <numbers>
#include <numeric>
#include <random>
#include <string>
#include <utility>
#include <vector>

#include <boost/geometry.hpp>
#include <boost/geometry/index/rtree.hpp>
#include <boost/polygon/point_data.hpp>
#include <boost/polygon/segment_data.hpp>
#include <boost/polygon/voronoi.hpp>

#include "wfmat/mat.hpp"
#include "wfmat/validate.hpp"

namespace wfmat::test {

// [VER-04] Residual of the radius function, at least two feet inside every edge, a tree whose
// vertex degrees match their kinds. Tolerances are relative to the half-diagonal L. Returns an
// empty string when everything holds, otherwise the first failure.
inline std::string sci(double x)
{
    char buf[32];
    std::snprintf(buf, sizeof buf, "%.3g", x);
    return buf;
}

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
        // A disk's MAT is its centre alone, a curvature end of degree 0.
        const std::size_t ends = vertices.size() == 1 ? 0 : 1;
        const bool deg_ok = v.kind == VertexKind::corner || v.kind == VertexKind::curvature_end ? deg == ends
                            : v.kind == VertexKind::junction                                     ? deg >= 3
                                                                                                 : deg == 2;
        if (!deg_ok) return where("vertex", i) + "degree " + std::to_string(deg) + " does not match its kind";
        const DiskCheck c = check_disk(region, xf.to_unit(v.p), xf.length_to_unit(v.r), 1e-7);
        if (c.residual * xf.scale > tol) return where("vertex", i) + "radius residual " + sci(c.residual * xf.scale);
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
                return where("edge", i) + "radius residual " + sci(c.residual * xf.scale) + " at u = " + std::to_string(u);
            // The disk touches the boundary at the edge's two feet, one on each site.
            const auto [fl, fr] = e.feet_at(u);
            for (Vec2 f : {fl, fr}) {
                const double on_boundary = boundary_distance(region.segments, xf.to_unit(f), 1e-7).distance;
                if (on_boundary * xf.scale > tol || std::abs(dist(p, f) - e.radius_at(u)) > tol)
                    return where("edge", i) + "a foot is not a contact point at u = " + std::to_string(u) + " (off by " +
                           sci(std::max(on_boundary * xf.scale, std::abs(dist(p, f) - e.radius_at(u)))) + ")";
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
// (minima of r inside a Voronoi edge) contracted, and with coincident Voronoi vertices merged.
// Vertices are matched by position within tol (in caller units) and the edge sets must then agree.
// Returns an empty string on success.
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
    }
    // A degenerate input (four or more sites on one circle) gives coincident Voronoi vertices
    // joined by zero-length edges: merge them, as our cluster resolution does.
    std::vector<std::size_t> rep(bverts.size());
    std::iota(rep.begin(), rep.end(), std::size_t{0});
    auto root = [&](std::size_t x) {
        while (rep[x] != x) x = rep[x] = rep[rep[x]];
        return x;
    };
    for (const auto& [i0, i1] : bedges)
        if (dist({bverts[i0].x(), bverts[i0].y()}, {bverts[i1].x(), bverts[i1].y()}) <= tol) rep[root(i1)] = root(i0);
    std::vector<std::pair<std::size_t, std::size_t>> merged;
    for (const auto& [i0, i1] : bedges) {
        const std::size_t r0 = root(i0), r1 = root(i1);
        if (r0 == r1) continue;
        merged.emplace_back(std::min(r0, r1), std::max(r0, r1));
        ++bdegree[r0];
        ++bdegree[r1];
    }
    bedges = std::move(merged);
    std::vector<std::size_t> bused;
    for (std::size_t i = 0; i < bverts.size(); ++i)
        if (bdegree[i] > 0) bused.push_back(i);

    // Our graph with the E3 vertices contracted, except where a degenerate contact coincides with a
    // Voronoi vertex (the contact at the end of an arc site, where the diagram changes sites).
    const auto verts = mat.vertices();
    const auto edges = mat.edges();
    std::vector<char> contracted(verts.size(), 0);
    for (std::size_t i = 0; i < verts.size(); ++i) {
        if (verts[i].kind != VertexKind::extremum_min) continue;
        contracted[i] = std::none_of(bused.begin(), bused.end(), [&](std::size_t j) {
            return dist(verts[i].p, {bverts[j].x(), bverts[j].y()}) <= tol;
        });
    }
    std::vector<std::size_t> ours;
    for (std::size_t i = 0; i < verts.size(); ++i)
        if (!contracted[i]) ours.push_back(i);
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
        if (contracted[a] || contracted[b]) continue;
        mine.emplace_back(std::min(to_boost[a], to_boost[b]), std::max(to_boost[a], to_boost[b]));
    }
    for (std::size_t i = 0; i < verts.size(); ++i) {
        if (!contracted[i]) continue;
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

// [VER-03] Sampling oracle for regions with arcs: Boost.Polygon's Voronoi diagram of the boundary
// flattened to a polyline with `chords` chords per half turn of arc. Distances are in caller units.
struct SampledComparison {
    double ours_to_voronoi = 0.0;  // largest distance from a point of our MAT to the Voronoi edges
    double voronoi_to_ours = 0.0;  // largest distance from a genuine Voronoi point to our MAT
    Vec2 worst_ours, worst_voronoi;  // where they occur
};

namespace detail {

namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;

inline double segment_distance(Vec2 p, Vec2 a, Vec2 b)
{
    const Vec2 ab = b - a;
    const double len2 = dot(ab, ab);
    const double s = len2 > 0.0 ? std::clamp(dot(p - a, ab) / len2, 0.0, 1.0) : 0.0;
    return dist(p, a + s * ab);
}

// Nearest-piece queries over a set of line pieces, through an R-tree of their boxes.
class PieceIndex {
public:
    explicit PieceIndex(const std::vector<std::pair<Vec2, Vec2>>& pieces) : pieces_(pieces)
    {
        std::vector<Item> items;
        items.reserve(pieces.size());
        for (std::size_t i = 0; i < pieces.size(); ++i) {
            const auto& [a, b] = pieces[i];
            items.emplace_back(BBox(BPoint(std::min(a.x, b.x), std::min(a.y, b.y)), BPoint(std::max(a.x, b.x), std::max(a.y, b.y))), i);
        }
        tree_ = Tree(items.begin(), items.end());
    }

    double distance(Vec2 p) const
    {
        double best = INFINITY;
        const BPoint q(p.x, p.y);
        for (auto it = tree_.qbegin(bgi::nearest(q, static_cast<unsigned>(tree_.size()))); it != tree_.qend(); ++it) {
            if (bg::distance(q, it->first) > best) break;
            const auto& [a, b] = pieces_[it->second];
            best = std::min(best, segment_distance(p, a, b));
        }
        return best;
    }

private:
    using BPoint = boost::geometry::model::point<double, 2, boost::geometry::cs::cartesian>;
    using BBox = boost::geometry::model::box<BPoint>;
    using Item = std::pair<BBox, std::size_t>;
    using Tree = boost::geometry::index::rtree<Item, boost::geometry::index::quadratic<16>>;
    const std::vector<std::pair<Vec2, Vec2>>& pieces_;
    Tree tree_;
};

} // namespace detail

inline SampledComparison compare_with_sampled_voronoi(const PreparedRegion& region, const MedialAxis& mat, int chords)
{
    namespace bp = boost::polygon;
    const Transform& xf = region.xf;
    const double grid = 1e8;  // integer grid per unit-frame length

    // Flatten in the unit frame, then round to the integer grid. Chord angles of at most
    // pi / (chords sqrt(max(R, 1))) keep every sagitta R phi^2 / 8 within (pi / chords)^2 / 8,
    // however flat the arc.
    std::vector<Vec2> pts;
    std::vector<std::int32_t> owner;  // the arc a point lies inside of, or -1 at a region corner or on a line
    std::vector<std::int32_t> chord;  // the arc the chord from a point flattens, or -1 for a line
    for (std::size_t i = 0; i < region.segments.size(); ++i) {
        const Segment& s = region.segments[i];
        const int m = s.is_arc() ? std::max(1, static_cast<int>(std::ceil(std::abs(s.sweep) * std::sqrt(std::max(s.R, 1.0)) *
                                                                         chords / std::numbers::pi)))
                                 : 1;
        for (int k = 0; k < m; ++k) {
            pts.push_back(s.point(double(k) / m));
            owner.push_back(k == 0 ? -1 : static_cast<std::int32_t>(i));
            chord.push_back(s.is_arc() ? static_cast<std::int32_t>(i) : -1);
        }
    }
    std::vector<bp::segment_data<std::int32_t>> isegs;
    std::vector<std::pair<Vec2, Vec2>> segs;  // the same, in caller units
    std::vector<bp::point_data<std::int32_t>> ip;
    std::vector<std::int32_t> arc_of_point, arc_of_chord;  // the arc a point or chord belongs to, or -1
    for (std::size_t k = 0; k < pts.size(); ++k) {
        const bp::point_data<std::int32_t> q(static_cast<std::int32_t>(std::lround(pts[k].x * grid)),
                                             static_cast<std::int32_t>(std::lround(pts[k].y * grid)));
        if (!ip.empty() && q == ip.back()) continue;
        ip.push_back(q);
        arc_of_point.push_back(owner[k]);
        arc_of_chord.push_back(chord[k]);
    }
    if (ip.size() > 1 && ip.front() == ip.back()) {
        ip.pop_back();
        arc_of_point.pop_back();
        arc_of_chord.pop_back();
    }
    auto caller = [&](const bp::point_data<std::int32_t>& q) { return xf.from_unit(Vec2{double(q.x()), double(q.y())} / grid); };
    for (std::size_t i = 0; i < ip.size(); ++i) {
        isegs.emplace_back(ip[i], ip[(i + 1) % ip.size()]);
        segs.emplace_back(caller(ip[i]), caller(ip[(i + 1) % ip.size()]));
    }
    bp::voronoi_diagram<double> vd;
    bp::construct_voronoi(isegs.begin(), isegs.end(), &vd);

    // Closest point of a cell's site to p, in caller units.
    auto site_foot = [&](const bp::voronoi_cell<double>& c, Vec2 p) {
        const auto& [a, b] = segs[c.source_index()];
        if (c.source_category() == bp::SOURCE_CATEGORY_SEGMENT_START_POINT) return a;
        if (c.source_category() == bp::SOURCE_CATEGORY_SEGMENT_END_POINT) return b;
        const Vec2 ab = b - a;
        return a + std::clamp(dot(p - a, ab) / dot(ab, ab), 0.0, 1.0) * ab;
    };
    auto inside = [&](Vec2 p) { return contains(region, xf.to_unit(p), 1e-9); };

    // Sample every finite primary Voronoi edge inside the region. Parabolic edges (a point site and
    // a segment site) are sampled through their parametrisation over the segment's line.
    const int per_edge = 8;
    std::vector<std::pair<Vec2, Vec2>> vpieces;
    std::vector<Vec2> genuine;
    const double min_angle = 4.0 * std::numbers::pi / chords;
    for (const auto& e : vd.edges()) {
        if (e.is_infinite() || e.is_secondary() || &e > e.twin()) continue;
        const Vec2 v0 = xf.from_unit(Vec2{e.vertex0()->x(), e.vertex0()->y()} / grid);
        const Vec2 v1 = xf.from_unit(Vec2{e.vertex1()->x(), e.vertex1()->y()} / grid);
        std::vector<Vec2> line{v0};
        if (e.is_curved()) {
            const auto& pc = e.cell()->contains_point() ? *e.cell() : *e.twin()->cell();
            const auto& sc = e.cell()->contains_point() ? *e.twin()->cell() : *e.cell();
            const Vec2 focus = site_foot(pc, v0);
            const auto& [a, b] = segs[sc.source_index()];
            const Vec2 d = unit(b - a);
            Vec2 n = perp(d);
            if (dot(focus - a, n) < 0.0) n = -n;
            const double px = dot(focus - a, d), py = dot(focus - a, n);
            const double x0 = dot(v0 - a, d), x1 = dot(v1 - a, d);
            for (int k = 1; k < per_edge; ++k) {
                const double x = x0 + (x1 - x0) * k / per_edge;
                line.push_back(a + x * d + ((x - px) * (x - px) + py * py) / (2.0 * py) * n);
            }
        } else {
            for (int k = 1; k < per_edge; ++k) line.push_back(v0 + (double(k) / per_edge) * (v1 - v0));
        }
        line.push_back(v1);
        bool in = true;
        for (std::size_t k = 1; k + 1 < line.size() && in; ++k) in = inside(line[k]);
        if (!in) continue;
        for (std::size_t k = 0; k + 1 < line.size(); ++k) vpieces.emplace_back(line[k], line[k + 1]);
        // Genuine axis points: their two feet are seen at an angle well above that of one chord.
        // Edges between chords of one arc, or a line and the first chord of an arc tangent to it,
        // approximate no part of the true axis.
        auto arc = [&](const bp::voronoi_cell<double>& c) {
            const std::size_t i = c.source_index();
            if (c.source_category() == bp::SOURCE_CATEGORY_SEGMENT_START_POINT) return arc_of_point[i];
            if (c.source_category() == bp::SOURCE_CATEGORY_SEGMENT_END_POINT) return arc_of_point[(i + 1) % ip.size()];
            return arc_of_chord[i];
        };
        const std::int32_t a0 = arc(*e.cell()), a1 = arc(*e.twin()->cell());
        if (a0 >= 0 && a0 == a1) continue;
        for (std::size_t k = 1; k + 1 < line.size(); ++k) {
            const Vec2 p = line[k];
            const Vec2 f0 = site_foot(*e.cell(), p), f1 = site_foot(*e.twin()->cell(), p);
            if (std::abs(signed_angle(f0 - p, f1 - p)) > min_angle) genuine.push_back(p);
        }
    }

    // Our MAT, flattened adaptively: edges are parametrised by radius, which is nearly constant
    // along the stretch around an E3 minimum, so uniform steps in u would leave long chords there.
    std::vector<std::pair<Vec2, Vec2>> opieces;
    std::vector<Vec2> ours;
    const double flat = 1e-7 * xf.scale;
    for (const MatEdge& e : mat.edges()) {
        std::vector<std::pair<double, Vec2>> stack{{1.0, e.point_at(1.0)}};
        double u = 0.0;
        Vec2 p = e.point_at(0.0);
        ours.push_back(p);
        while (!stack.empty()) {
            const auto [u1, q] = stack.back();
            const double um = 0.5 * (u + u1);
            const Vec2 m = e.point_at(um);
            if (u1 - u > 1.0 / 64 || (detail::segment_distance(m, p, q) > flat && u1 - u > 1e-9)) {
                stack.push_back({um, m});
                continue;
            }
            stack.pop_back();
            opieces.emplace_back(p, q);
            ours.push_back(q);
            u = u1;
            p = q;
        }
    }
    for (const MatVertex& v : mat.vertices()) {
        ours.push_back(v.p);
        opieces.emplace_back(v.p, v.p);
    }

    SampledComparison out;
    const detail::PieceIndex vindex(vpieces), oindex(opieces);
    for (Vec2 p : ours)
        if (const double d = vindex.distance(p); d > out.ours_to_voronoi) {
            out.ours_to_voronoi = d;
            out.worst_ours = p;
        }
    for (Vec2 p : genuine)
        if (const double d = oindex.distance(p); d > out.voronoi_to_ours) {
            out.voronoi_to_ours = d;
            out.worst_voronoi = p;
        }
    return out;
}

// [VER-05] A random star-shaped polygon with integer coordinates in [-R, R]^2, counter-clockwise,
// without collinear neighbours or parallel edges (those make plateaus and clusters, tested on their
// own in test_degenerate.cpp).
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
                return std::pair<std::int64_t, std::int32_t>{std::int64_t(v.first) - u.first,
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

// A polygon (counter-clockwise) with corner k rounded by an arc tangent to both of its edges, the
// tangent points at fraction[k] of the shorter half of the two edges (0 leaves the corner sharp).
// Convex corners get convex arcs, reflex corners concave ones.
inline Region filleted(const std::vector<Vec2>& polygon, const std::vector<double>& fraction)
{
    Region region;
    const std::size_t n = polygon.size();
    for (std::size_t k = 0; k < n; ++k) {
        const Vec2 prev = polygon[(k + n - 1) % n], p = polygon[k], next = polygon[(k + 1) % n];
        if (fraction[k] <= 0.0) {
            region.outer.vertices.push_back({p, 0.0});
            continue;
        }
        const Vec2 u = unit(p - prev), w = unit(next - p);
        const double d = fraction[k] * 0.5 * std::min(dist(prev, p), dist(p, next));
        const double turn = signed_angle(u, w);
        region.outer.vertices.push_back({p - d * u, bulge_from_sweep(turn)});
        region.outer.vertices.push_back({p + d * w, 0.0});
    }
    return region;
}

inline Region to_region(const std::vector<IPoint>& p)
{
    Region region;
    for (const IPoint& q : p) region.outer.vertices.push_back({{double(q.first), double(q.second)}, 0.0});
    return region;
}

// [VER-05] A random star polygon with every corner filleted (fraction 0 to 0.9 of the room).
inline Region random_filleted(std::mt19937_64& rng, int n)
{
    const auto poly = random_star_polygon(rng, n);
    std::vector<Vec2> p;
    for (const IPoint& q : poly) p.push_back({double(q.first), double(q.second)});
    std::uniform_real_distribution<double> f(0.05, 0.9);
    std::vector<double> fraction(p.size());
    for (double& x : fraction) x = f(rng);
    return filleted(p, fraction);
}

// [VER-05] A random arc chain: a star polygon whose edges are arcs with random bulges.
inline Region random_arc_chain(std::mt19937_64& rng, int n)
{
    Region region = to_region(random_star_polygon(rng, n));
    std::uniform_real_distribution<double> b(-0.25, 0.25);
    for (BulgeVertex& v : region.outer.vertices) v.bulge = b(rng);
    return region;
}

} // namespace wfmat::test
