// JSON input in the v1 bulge schema [IN-02].
//
//   { "units": "mm", "outer": [[x, y, bulge], [x, y], ...] }
//
// A vertex without a bulge starts a line. "units" is optional and informational. Any other
// key is rejected, and so is "holes" with a pointer to v2 [IN-05].
#pragma once

#include <filesystem>
#include <string>
#include <string_view>

#include <span>
#include <vector>

#include "wfmat/mat.hpp"
#include "wfmat/region.hpp"
#include "wfmat/result.hpp"

namespace wfmat::io {

struct RegionDocument {
    Region region;
    std::string units;
};

Result<RegionDocument> parse_region_json(std::string_view text);
Result<RegionDocument> read_region_json(const std::filesystem::path& path);

// Writes the v1 schema; the output parses back to the same region bit for bit.
std::string write_region_json(const RegionDocument& doc);

// [OUT-06] The medial axis as JSON, in caller units:
//
//   { "units": "mm",
//     "sites":    [{"segment": 3, "corner": false, "kind": "line", "normal": [nx, ny], "offset": c}
//                  | {..., "kind": "circle", "centre": [x, y], "radius": R, "sigma": 1}],
//     "vertices": [{"p": [x, y], "r": r, "kind": "junction", "edges": [...],
//                   "contacts": [{"site": i, "foot": [x, y]}]}],
//     "edges":    [{"v0": i, "v1": j, "left": a, "right": b, "kind": "parabola", "r0": r0, "r1": r1,
//                   "bezier": [{"p0": [x, y], "p1": [x, y], "p2": [x, y], "w": w, "u": [u0, u1]}],
//                   "polyline": [[x, y], ...]}],
//     "offsets":  [{"distance": d, "loops": [[[x, y, bulge], ...], ...]}],
//     "stats":    {"collapses": n, ...} }
//
// A site's distance is normal . x - offset for a line, sigma (R - |x - centre|) for a circle. u
// runs from v0 (radius r0) to v1 (radius r1); "polyline" is present when chord_tol > 0.
struct OffsetLoops {
    double distance = 0.0;
    std::vector<Loop> loops;
};
struct MatJsonOptions {
    std::string units;
    double chord_tol = 0.0;  // flatten each edge to this tolerance; 0 leaves the polylines out
};
std::string write_mat_json(const MedialAxis& mat, std::span<const OffsetLoops> offsets = {},
                           const MatJsonOptions& options = {});

} // namespace wfmat::io
