#include "wfmat/prepare.hpp"

#include <algorithm>
#include <cmath>
#include <iterator>
#include <numbers>
#include <string>
#include <utility>

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

auto input_error(const char* requirement, std::string message, std::size_t segment)
{
    return make_error(ErrorCode::invalid_input, requirement, std::move(message), 0, segment);
}

double turning_angle(const Segment& s, const Segment& t)
{
    return signed_angle(s.end_tangent(), t.start_tangent());
}

// [IN-04] The loop is simple: segments meet only consecutively, at their shared vertex.
Result<void> check_simple(const std::vector<Segment>& segs, double eps)
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
                                   i);
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

} // namespace

Result<PreparedRegion> prepare(const Region& region, const Options& options)
{
    const Tolerances& tol = options.tol;
    const std::vector<BulgeVertex>& v = region.outer.vertices;
    const std::size_t n = v.size();
    if (n < 2) return make_error(ErrorCode::invalid_input, "IN-03", "a loop needs at least two vertices", 0);

    Box vertex_box;
    for (std::size_t i = 0; i < n; ++i) {
        if (!is_finite(v[i].p) || !std::isfinite(v[i].bulge))
            return input_error("IN-01", "non-finite coordinate or bulge at vertex " + std::to_string(i), i);
        vertex_box.add(v[i].p);
    }
    const double rough_scale = 0.5 * norm(vertex_box.extent());
    if (!(rough_scale > 0.0)) return input_error("IN-03", "all vertices coincide", 0);

    // Zero-length chords leave an arc undefined, so reject them before building segments.
    for (std::size_t i = 0; i < n; ++i)
        if (dist(v[i].p, v[(i + 1) % n].p) <= tol.len * rough_scale)
            return input_error("IN-03", "segment " + std::to_string(i) + " is shorter than the minimum length", i);

    // [IN-07] Bounding box (with arc extremes) to the unit frame.
    Box box;
    for (std::size_t i = 0; i < n; ++i) box.add(make_segment(v[i].p, v[(i + 1) % n].p, v[i].bulge).bbox());
    PreparedRegion out;
    out.tol = tol;
    out.xf = Transform{box.center(), 0.5 * norm(box.extent())};

    std::vector<Segment> segs;
    segs.reserve(n);
    for (std::size_t i = 0; i < n; ++i)
        segs.push_back(make_segment(out.xf.to_unit(v[i].p), out.xf.to_unit(v[(i + 1) % n].p), v[i].bulge));

    // [IN-03] Minimum length and arc sweep.
    for (std::size_t i = 0; i < n; ++i) {
        const Segment& s = segs[i];
        if (s.length() <= tol.len)
            return input_error("IN-03", "segment " + std::to_string(i) + " is shorter than the minimum length", i);
        if (s.is_arc() && !(std::abs(s.sweep) > tol.ang && std::abs(s.sweep) < 2.0 * pi))
            return input_error("IN-03", "arc " + std::to_string(i) + " has a sweep outside (eps_ang, 2 pi)", i);
    }

    // [IN-06] No cusps.
    for (std::size_t i = 0; i < n; ++i) {
        const std::size_t j = (i + 1) % n;
        if (std::abs(turning_angle(segs[i], segs[j])) >= pi - tol.ang)
            return input_error("IN-06", "cusp at vertex " + std::to_string(j), j);
    }

    if (auto simple = check_simple(segs, tol.geom); !simple) return tl::unexpected(simple.error());

    // [IN-08] Counter-clockwise orientation from the signed area.
    double area = 0.0;
    for (const Segment& s : segs) area += 0.5 * cross(s.a, s.b) + s.cap_area();
    std::vector<std::vector<std::uint32_t>> sources(n);
    if (area < 0.0) {
        std::vector<Segment> rev;
        rev.reserve(n);
        for (std::size_t k = 0; k < n; ++k) {
            rev.push_back(segs[n - 1 - k].reversed());
            sources[k] = {static_cast<std::uint32_t>(n - 1 - k)};
        }
        segs = std::move(rev);
        out.reversed = true;
        area = -area;
    } else {
        for (std::size_t k = 0; k < n; ++k) sources[k] = {static_cast<std::uint32_t>(k)};
    }
    out.area = area;

    merge_neighbours(segs, sources, tol);

    // [ALG-01] Corner classification; [IN-10] near-tangent joins become G1 with one shared tangent.
    const std::size_t m = segs.size();
    out.joins.reserve(m);
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
    }

    for (const Segment& s : segs) out.bbox.add(s.bbox());
    out.segments = std::move(segs);
    out.sources = std::move(sources);
    return out;
}

} // namespace wfmat
