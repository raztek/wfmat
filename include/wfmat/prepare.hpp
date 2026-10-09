// Input validation and normalisation [IN-03..IN-10], corner classification [ALG-01].
#pragma once

#include <cstdint>
#include <vector>

#include "wfmat/geom.hpp"
#include "wfmat/options.hpp"
#include "wfmat/region.hpp"
#include "wfmat/result.hpp"

namespace wfmat {

// Similarity mapping caller coordinates to the unit frame [IN-07]: q = (p - centre) / scale.
struct Transform {
    Vec2 centre;
    double scale = 1.0;  // L, the half-diagonal of the caller-frame bounding box

    Vec2 to_unit(Vec2 p) const { return (p - centre) / scale; }
    Vec2 from_unit(Vec2 q) const { return centre + scale * q; }
    double length_to_unit(double d) const { return d / scale; }
    double length_from_unit(double d) const { return d * scale; }
};

enum class JoinKind { convex, tangent, reflex };

// The junction at the end of segment i and the start of segment i + 1 [ALG-01].
struct Join {
    Vec2 p;
    double turn = 0.0;   // signed turning angle, positive to the left (convex)
    JoinKind kind = JoinKind::convex;
    Vec2 t_in, t_out;    // unit tangents of the two segments at p; equal for tangent joins [IN-10]
};

// A validated, normalised region: one counter-clockwise loop in the unit frame.
struct PreparedRegion {
    Transform xf;
    std::vector<Segment> segments;                    // unit frame, counter-clockwise [IN-08]
    std::vector<Join> joins;                          // joins[i] follows segments[i]
    std::vector<std::vector<std::uint32_t>> sources;  // input segment indices merged into each [IN-09]
    bool reversed = false;                            // input loop was clockwise
    double area = 0.0;                                // unit frame
    Box bbox;                                         // unit frame
    Tolerances tol;
};

// Validates and normalises the region. On failure the error names the requirement violated
// and the input segment, numbered as in the caller's vertex list.
Result<PreparedRegion> prepare(const Region& region, const Options& options = {});

} // namespace wfmat
