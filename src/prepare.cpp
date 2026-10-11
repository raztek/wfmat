#include "wfmat/prepare.hpp"

#include <algorithm>
#include <cmath>
#include <iterator>
#include <numbers>
#include <string>
#include <utility>

#include "wfmat/validate.hpp"

#include <boost/geometry.hpp>
#include <boost/geometry/index/rtree.hpp>

namespace wfmat {

namespace {

namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;
using BPoint = bg::model::point<double, 2, bg::cs::cartesian>;
using BBox = bg::model::box<BPoint>;
using BoxItem = std::pair<BBox, std::uint32_t>;

constexpr double pi = std::numbers::pi;

BBox to_bbox(const Box& b)
{
    return {BPoint(b.lo.x, b.lo.y), BPoint(b.hi.x, b.hi.y)};
}

auto input_error(const char* requirement, std::string message, std::size_t loop, std::optional<std::size_t> segment)
{
    return make_error(ErrorCode::invalid_input, requirement, std::move(message), loop, segment);
}

double turning_angle(const Segment& s, const Segment& t)
{
    return signed_angle(s.end_tangent(), t.start_tangent());
}

// [IN-04] The loop is simple: segments meet only consecutively, at their shared vertex.
Result<void> check_simple(const std::vector<Segment>& segs, double eps, std::size_t loop)
{
    const std::size_t n = segs.size();
    std::vector<BoxItem> items;
    items.reserve(n);
    for (std::uint32_t i = 0; i < n; ++i) items.emplace_back(to_bbox(segs[i].bbox().inflated(eps)), i);
    const bgi::rtree<BoxItem, bgi::quadratic<16>> tree(items.begin(), items.end());

    std::vector<BoxItem> hits;
    for (std::uint32_t i = 0; i < n; ++i) {
        hits.clear();
        tree.query(bgi::intersects(items[i].first), std::back_inserter(hits));
        std::sort(hits.begin(), hits.end(), [](const BoxItem& x, const BoxItem& y) { return x.second < y.second; });
        for (const BoxItem& hit : hits) {
            const std::uint32_t j = hit.second;
            if (j <= i) continue;
            std::vector<Vec2> shared;
            if (j == i + 1) shared.push_back(segs[i].b);
            if (i == 0 && j == n - 1) shared.push_back(segs[0].a);
            const Intersections x = intersect(segs[i], segs[j], eps);
            bool bad = x.overlap;
            for (Vec2 p : x.points) {
                const bool at_shared = std::any_of(shared.begin(), shared.end(),
                                                   [&](Vec2 q) { return dist(p, q) <= 2.0 * eps; });
                if (!at_shared) bad = true;
            }
            if (bad)
                return input_error("IN-04",
                                   "the loop is not simple: segments " + std::to_string(i) + " and " +
                                       std::to_string(j) + (x.overlap ? " overlap" : " intersect"),
                                   loop, i);
        }
    }
    return {};
}

// [IN-12] Segments of different loops neither cross nor touch. loop_of[i] is the loop of segs[i];
// local[i] its index in that loop.
Result<void> check_disjoint(const std::vector<Segment>& segs, const std::vector<std::uint32_t>& loop_of,
                            const std::vector<std::uint32_t>& local, double eps)
{
    const std::size_t n = segs.size();
    std::vector<BoxItem> items;
    items.reserve(n);
    for (std::uint32_t i = 0; i < n; ++i) items.emplace_back(to_bbox(segs[i].bbox().inflated(eps)), i);
    const bgi::rtree<BoxItem, bgi::quadratic<16>> tree(items.begin(), items.end());

    std::vector<BoxItem> hits;
    for (std::uint32_t i = 0; i < n; ++i) {
        hits.clear();
        tree.query(bgi::intersects(items[i].first), std::back_inserter(hits));
        std::sort(hits.begin(), hits.end(), [](const BoxItem& x, const BoxItem& y) { return x.second < y.second; });
        for (const BoxItem& hit : hits) {
            const std::uint32_t j = hit.second;
            if (j <= i || loop_of[j] == loop_of[i]) continue;
            const Intersections x = intersect(segs[i], segs[j], eps);
            if (x.overlap || !x.points.empty())
                return input_error("IN-12",
                                   "loops " + std::to_string(loop_of[i]) + " and " + std::to_string(loop_of[j]) +
                                       " meet: segment " + std::to_string(local[i]) + " of loop " +
                                       std::to_string(loop_of[i]) + " and segment " + std::to_string(local[j]) +
                                       " of loop " + std::to_string(loop_of[j]),
                                   loop_of[i], local[i]);
        }
    }
    return {};
}

bool mergeable(const Segment& s, const Segment& t, const Tolerances& tol)
{
    if (!s.is_arc() && !t.is_arc()) return std::abs(turning_angle(s, t)) <= tol.ang;
    if (s.is_arc() && t.is_arc())
        return dist(s.c, t.c) <= tol.geom && std::abs(s.R - t.R) <= tol.geom && (s.sweep > 0.0) == (t.sweep > 0.0) &&
               std::abs(s.sweep + t.sweep) < 2.0 * pi - tol.ang;
    return false;
}

// [IN-09] Merges consecutive collinear lines and consecutive co-circular arcs.
void merge_neighbours(std::vector<Segment>& segs, std::vector<std::vector<std::uint32_t>>& sources,
                      const Tolerances& tol)
{
    bool changed = true;
    while (changed && segs.size() > 2) {
        changed = false;
        for (std::size_t i = 0; i < segs.size() && segs.size() > 2;) {
            const std::size_t j = (i + 1) % segs.size();
            if (!mergeable(segs[i], segs[j], tol)) {
                ++i;
                continue;
            }
            const Segment& s = segs[i];
            const Segment& t = segs[j];
            segs[i] = make_segment(s.a, t.b, s.is_arc() ? bulge_from_sweep(s.sweep + t.sweep) : 0.0);
            sources[i].insert(sources[i].end(), sources[j].begin(), sources[j].end());
            segs.erase(segs.begin() + static_cast<std::ptrdiff_t>(j));
            sources.erase(sources.begin() + static_cast<std::ptrdiff_t>(j));
            if (j < i) --i;
            changed = true;
        }
    }
}


// One input loop in the unit frame, after the per-loop checks.
struct LoopIn {
    std::vector<Segment> segs;
    double area = 0.0;
    Box bbox;
};

// [IN-03], [IN-06], [IN-04] for one loop, already in the unit frame.
Result<LoopIn> check_loop(const std::vector<BulgeVertex>& v, std::size_t loop, const Transform& xf, const Tolerances& tol)
{
    const std::size_t n = v.size();
    LoopIn out;
    out.segs.reserve(n);
    for (std::size_t i = 0; i < n; ++i)
        out.segs.push_back(make_segment(xf.to_unit(v[i].p), xf.to_unit(v[(i + 1) % n].p), v[i].bulge));
    const auto& segs = out.segs;

    // [IN-03] Minimum length and arc sweep.
    for (std::size_t i = 0; i < n; ++i) {
        const Segment& s = segs[i];
        if (s.length() <= tol.len)
            return input_error("IN-03", "segment " + std::to_string(i) + " is shorter than the minimum length", loop, i);
        if (s.is_arc() && !(std::abs(s.sweep) > tol.ang && std::abs(s.sweep) < 2.0 * pi))
            return input_error("IN-03", "arc " + std::to_string(i) + " has a sweep outside (eps_ang, 2 pi)", loop, i);
    }

    // [IN-06] No cusps.
    for (std::size_t i = 0; i < n; ++i) {
        const std::size_t j = (i + 1) % n;
        if (std::abs(turning_angle(segs[i], segs[j])) >= pi - tol.ang)
            return input_error("IN-06", "cusp at vertex " + std::to_string(j), loop, j);
    }

    if (auto simple = check_simple(segs, tol.geom, loop); !simple) return tl::unexpected(simple.error());

    for (const Segment& s : segs) {
        out.area += 0.5 * cross(s.a, s.b) + s.cap_area();
        out.bbox.add(s.bbox());
    }
    return out;
}

} // namespace

Result<PreparedRegion> prepare(const Region& region, const Options& options)
{
    const Tolerances& tol = options.tol;
    const std::size_t nloops = region.loops.size();
    if (nloops == 0) return make_error(ErrorCode::invalid_input, "IN-12", "a region needs at least one loop");

    Box vertex_box;
    for (std::size_t k = 0; k < nloops; ++k) {
        const std::vector<BulgeVertex>& v = region.loops[k].vertices;
        if (v.size() < 2)
            return make_error(ErrorCode::invalid_input, "IN-03", "a loop needs at least two vertices", k);
        for (std::size_t i = 0; i < v.size(); ++i) {
            if (!is_finite(v[i].p) || !std::isfinite(v[i].bulge))
                return input_error("IN-01", "non-finite coordinate or bulge at vertex " + std::to_string(i), k, i);
            vertex_box.add(v[i].p);
        }
    }
    const double rough_scale = 0.5 * norm(vertex_box.extent());
    if (!(rough_scale > 0.0)) return input_error("IN-03", "all vertices coincide", 0, 0);

    // Zero-length chords leave an arc undefined, so reject them before building segments.
    for (std::size_t k = 0; k < nloops; ++k) {
        const std::vector<BulgeVertex>& v = region.loops[k].vertices;
        const std::size_t n = v.size();
        for (std::size_t i = 0; i < n; ++i)
            if (dist(v[i].p, v[(i + 1) % n].p) <= tol.len * rough_scale)
                return input_error("IN-03", "segment " + std::to_string(i) + " is shorter than the minimum length", k, i);
    }

    // [IN-07] Bounding box of all loops (with arc extremes) to the unit frame.
    Box box;
    for (const Loop& loop : region.loops) {
        const std::vector<BulgeVertex>& v = loop.vertices;
        for (std::size_t i = 0; i < v.size(); ++i) box.add(make_segment(v[i].p, v[(i + 1) % v.size()].p, v[i].bulge).bbox());
    }
    PreparedRegion out;
    out.tol = tol;
    out.xf = Transform{box.center(), 0.5 * norm(box.extent())};

    std::vector<LoopIn> loops;
    loops.reserve(nloops);
    for (std::size_t k = 0; k < nloops; ++k) {
        auto loop = check_loop(region.loops[k].vertices, k, out.xf, tol);
        if (!loop) return tl::unexpected(loop.error());
        loops.push_back(std::move(*loop));
    }

    // [IN-12] Disjoint loops.
    if (nloops > 1) {
        std::vector<Segment> all;
        std::vector<std::uint32_t> loop_of, local;
        for (std::uint32_t k = 0; k < nloops; ++k)
            for (std::uint32_t i = 0; i < loops[k].segs.size(); ++i) {
                all.push_back(loops[k].segs[i]);
                loop_of.push_back(k);
                local.push_back(i);
            }
        if (auto disjoint = check_disjoint(all, loop_of, local, tol.len); !disjoint)
            return tl::unexpected(disjoint.error());
    }

    // [IN-13] Depth: the loops containing a vertex of the loop; disjointness makes one vertex enough.
    // The parent is the containing loop one level up.
    out.loops.resize(nloops);
    std::vector<std::vector<std::uint32_t>> containers(nloops);
    for (std::uint32_t k = 0; k < nloops; ++k) {
        const Vec2 p = loops[k].segs.front().a;
        for (std::uint32_t j = 0; j < nloops; ++j) {
            const Box& b = loops[j].bbox;
            if (j == k || p.x < b.lo.x || p.x > b.hi.x || p.y < b.lo.y || p.y > b.hi.y) continue;
            if (winding_number(loops[j].segs, p) != 0) containers[k].push_back(j);
        }
        out.loops[k].depth = static_cast<std::uint32_t>(containers[k].size());
        if (out.loops[k].depth > max_nesting_depth)
            return input_error("IN-13",
                               "loop " + std::to_string(k) + " is nested " + std::to_string(out.loops[k].depth) +
                                   " levels deep; at most 8 levels (depths 0 to 7) are allowed",
                               k, std::nullopt);
    }
    for (std::uint32_t k = 0; k < nloops; ++k)
        for (std::uint32_t j : containers[k])
            if (out.loops[j].depth + 1 == out.loops[k].depth) out.loops[k].parent = j;

    // [IN-08] Orientation is an input rule: counter-clockwise on even depths, clockwise on odd.
    for (std::uint32_t k = 0; k < nloops; ++k) {
        const bool even = out.loops[k].depth % 2 == 0;
        if ((loops[k].area > 0.0) != even)
            return input_error("IN-08",
                               "loop " + std::to_string(k) + " at depth " + std::to_string(out.loops[k].depth) +
                                   " must run " + (even ? "counter-clockwise" : "clockwise"),
                               k, std::nullopt);
    }

    // [IN-13] Components, in the input order of their outer loops.
    for (std::uint32_t k = 0; k < nloops; ++k)
        if (out.loops[k].depth % 2 == 0) {
            out.loops[k].component = static_cast<std::uint32_t>(out.components.size());
            out.components.push_back({k, {}});
        }
    for (std::uint32_t k = 0; k < nloops; ++k)
        if (out.loops[k].depth % 2 == 1) {
            const std::uint32_t c = out.loops[out.loops[k].parent].component;
            out.loops[k].component = c;
            out.components[c].holes.push_back(k);
        }

    // [IN-09], [IN-10], [ALG-01] Per loop, never across loops.
    for (std::uint32_t k = 0; k < nloops; ++k) {
        std::vector<Segment>& segs = loops[k].segs;
        std::vector<std::vector<std::uint32_t>> sources(segs.size());
        for (std::size_t i = 0; i < segs.size(); ++i) sources[i].assign(1, static_cast<std::uint32_t>(i));
        merge_neighbours(segs, sources, tol);

        const std::size_t m = segs.size();
        PreparedLoop& pl = out.loops[k];
        pl.first = static_cast<std::uint32_t>(out.segments.size());
        pl.count = static_cast<std::uint32_t>(m);
        pl.area = loops[k].area;
        out.area += loops[k].area;
        for (std::size_t i = 0; i < m; ++i) {
            const Segment& s = segs[i];
            const Segment& t = segs[(i + 1) % m];
            Join join;
            join.p = s.b;
            join.t_in = s.end_tangent();
            join.t_out = t.start_tangent();
            join.turn = signed_angle(join.t_in, join.t_out);
            if (join.turn > tol.ang) {
                join.kind = JoinKind::convex;
            } else if (join.turn < -tol.ang) {
                join.kind = JoinKind::reflex;
            } else {
                join.kind = JoinKind::tangent;
                join.turn = 0.0;
                join.t_in = join.t_out = unit(join.t_in + join.t_out);
            }
            out.joins.push_back(join);
            pl.bbox.add(s.bbox());
            out.segments.push_back(s);
            out.sources.push_back(std::move(sources[i]));
        }
        out.bbox.add(pl.bbox);
    }
    return out;
}

} // namespace wfmat
