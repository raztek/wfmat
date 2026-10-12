// Input validation and normalisation [IN-03..IN-10], corner classification [ALG-01].
#pragma once

#include <cstdint>
#include <span>
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

// [IN-13] Loops may be nested at most 8 levels: depths 0 to 7.
inline constexpr std::uint32_t max_nesting_depth = 7;

// One loop of a prepared region [IN-14]: its segments are segments[first, first + count), each
// followed by its join, the last joining back to the first.
struct PreparedLoop {
    static constexpr std::uint32_t no_loop = 0xffffffffu;
    std::uint32_t first = 0, count = 0;
    std::uint32_t depth = 0;       // number of loops containing it [IN-13]
    std::uint32_t component = 0;
    std::uint32_t parent = no_loop;  // the loop directly containing it
    double area = 0.0;             // signed, unit frame: positive counter-clockwise
    Box bbox;                      // unit frame
};

// A component of the region [IN-13]: an outer loop and the loops directly inside it.
struct PreparedComponent {
    std::uint32_t outer = 0;
    std::vector<std::uint32_t> holes;  // loop indices, ascending
};

// A validated, normalised region in the unit frame. Loops keep their input order and numbering;
// segments, joins and sources run loop by loop, so a brute-force query over all segments sees the
// whole boundary.
struct PreparedRegion {
    Transform xf;
    std::vector<Segment> segments;                    // unit frame, the material on the left [IN-08]
    std::vector<Join> joins;                          // joins[i] follows segments[i] within its loop
    std::vector<std::vector<std::uint32_t>> sources;  // input segment indices (within the loop) merged into each [IN-09]
    std::vector<PreparedLoop> loops;
    std::vector<PreparedComponent> components;        // in the input order of their outer loops
    double area = 0.0;                                // of the material, unit frame
    Box bbox;                                         // unit frame
    Tolerances tol;

    // The segments of one loop.
    std::span<const Segment> loop_segments(std::size_t loop) const
    {
        return std::span<const Segment>(segments).subspan(loops[loop].first, loops[loop].count);
    }
};

// Validates and normalises the region. On failure the error names the requirement violated,
// the input loop, and the segment, numbered as in the caller's vertex list.
Result<PreparedRegion> prepare(const Region& region, const Options& options = {});

} // namespace wfmat
