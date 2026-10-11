// Input regions in bulge form [IN-01], [API-02].
#pragma once

#include <span>
#include <variant>
#include <vector>

#include "wfmat/geom.hpp"
#include "wfmat/result.hpp"

namespace wfmat {

// Vertex i starts segment i, which ends at vertex i + 1 (cyclically): a line when bulge == 0,
// otherwise an arc with bulge = tan(sweep / 4).
struct BulgeVertex {
    Vec2 p;
    double bulge = 0.0;
};

struct Loop {
    std::vector<BulgeVertex> vertices;
};

// [SC-01], [IN-13] Any number of disjoint loops, in any order. A loop's nesting depth is the
// number of loops containing it, at most 7; loops of even depth run counter-clockwise and loops
// of odd depth clockwise [IN-08]. Region{{loop}} is a region of one loop.
struct Region {
    std::vector<Loop> loops;
};

// Explicit segment forms, converted to bulge form by loop_from_segments [IN-01].
struct LineSpec {
    Vec2 a, b;
};
struct ArcSpec {
    Vec2 centre;
    double radius = 0.0;
    double start_angle = 0.0;  // radians
    double sweep = 0.0;        // signed, positive counter-clockwise, |sweep| < 2 pi
};
using SegmentSpec = std::variant<LineSpec, ArcSpec>;

// Builds a loop from explicit segments in boundary order. Each segment must start where the
// previous one ends (and the last must end where the first starts) to within gap, measured
// in the caller's units.
Result<Loop> loop_from_segments(std::span<const SegmentSpec> segments, double gap);

} // namespace wfmat
