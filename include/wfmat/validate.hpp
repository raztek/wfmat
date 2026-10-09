// Brute-force boundary distance and containment: the reference validator of [OUT-07] and [VER-04].
// Everything here is O(n) per query and deliberately independent of the propagation engine.
#pragma once

#include <cstddef>
#include <cstdint>
#include <span>
#include <vector>

#include "wfmat/geom.hpp"
#include "wfmat/prepare.hpp"

namespace wfmat {

// A closest boundary point (foot) of a query point.
struct Foot {
    std::uint32_t segment = 0;  // index into PreparedRegion::segments
    Vec2 point;
    bool whole_arc = false;     // the query point is the arc's centre: every arc point is a foot
};

struct BoundaryDistance {
    double distance = 0.0;
    std::vector<Foot> feet;  // all feet within foot_tol of the minimum, distinct to within foot_tol

    // Number of distinct feet; an arc whose every point is a foot counts as two.
    std::size_t distinct_feet() const;
};

// d(p) = dist(p, boundary), with every foot.
BoundaryDistance boundary_distance(std::span<const Segment> segments, Vec2 p, double foot_tol);

// Winding number of the closed chain of segments about p (p not on the boundary).
int winding_number(std::span<const Segment> segments, Vec2 p);

// p strictly inside the (counter-clockwise) region; points within tol of the boundary are outside.
bool contains(const PreparedRegion& region, Vec2 p, double tol);

// [OUT-07] check of a candidate inscribed disk (p, r) in the unit frame.
struct DiskCheck {
    double residual = 0.0;    // |d(p) - r|
    std::size_t feet = 0;     // distinct feet of p
    bool inside = false;      // p inside the region
};
DiskCheck check_disk(const PreparedRegion& region, Vec2 p, double r, double foot_tol);

} // namespace wfmat
