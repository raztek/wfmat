#include "wfmat/validate.hpp"

#include <cmath>
#include <limits>
#include <numbers>

namespace wfmat {

std::size_t BoundaryDistance::distinct_feet() const
{
    std::size_t count = 0;
    for (const Foot& f : feet) count += f.whole_arc ? 2 : 1;
    return count;
}

BoundaryDistance boundary_distance(std::span<const Segment> segments, Vec2 p, double foot_tol)
{
    struct Candidate {
        double d;
        Foot foot;
    };
    std::vector<Candidate> candidates;
    candidates.reserve(segments.size());
    double best = std::numeric_limits<double>::infinity();
    for (std::uint32_t i = 0; i < segments.size(); ++i) {
        const Segment& s = segments[i];
        const Vec2 q = s.closest_point(p);
        const bool whole = s.is_arc() && dist(p, s.c) <= foot_tol;
        const double d = whole ? s.R : dist(p, q);
        candidates.push_back({d, Foot{i, q, whole}});
        best = std::min(best, d);
    }

    BoundaryDistance out;
    out.distance = best;
    for (const Candidate& c : candidates) {
        if (c.d > best + foot_tol) continue;
        bool duplicate = false;
        if (!c.foot.whole_arc)
            for (const Foot& f : out.feet)
                if (!f.whole_arc && dist(f.point, c.foot.point) <= foot_tol) duplicate = true;
        if (!duplicate) out.feet.push_back(c.foot);
    }
    return out;
}

namespace {

// Angle subtended at p by an arc: its chord's angle plus a full turn when p lies in the cap
// between chord and arc. A p on the chord line makes the chord angle ambiguous (+-pi), so the
// arc is split until p is clear of every sub-chord.
double arc_angle(const Segment& s, Vec2 p, int depth)
{
    const Vec2 chord = s.b - s.a;
    const double side = cross(chord, p - s.a);
    if (depth < 16 && std::abs(side) <= 1e-12 * dot(chord, chord)) {
        const Vec2 m = s.midpoint();
        const double half = bulge_from_sweep(0.5 * s.sweep);
        return arc_angle(make_segment(s.a, m, half), p, depth + 1) + arc_angle(make_segment(m, s.b, half), p, depth + 1);
    }
    double angle = signed_angle(s.a - p, s.b - p);
    const double arc_side = cross(chord, s.midpoint() - s.a);
    if (dist(p, s.c) < s.R && side * arc_side > 0.0) angle += (s.sweep > 0.0 ? 2.0 : -2.0) * std::numbers::pi;
    return angle;
}

} // namespace

int winding_number(std::span<const Segment> segments, Vec2 p)
{
    double total = 0.0;
    for (const Segment& s : segments) total += s.is_arc() ? arc_angle(s, p, 0) : signed_angle(s.a - p, s.b - p);
    return static_cast<int>(std::lround(total / (2.0 * std::numbers::pi)));
}

bool contains(const PreparedRegion& region, Vec2 p, double tol)
{
    if (boundary_distance(region.segments, p, tol).distance <= tol) return false;
    return winding_number(region.segments, p) != 0;
}

DiskCheck check_disk(const PreparedRegion& region, Vec2 p, double r, double foot_tol)
{
    const BoundaryDistance bd = boundary_distance(region.segments, p, foot_tol);
    return {std::abs(bd.distance - r), bd.distinct_feet(), contains(region, p, foot_tol)};
}

} // namespace wfmat
