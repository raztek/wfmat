// JSON input in the v2 bulge schema [IN-02].
//
//   { "units": "mm", "loops": [[[x, y, bulge], [x, y], ...], ...] }
//
// or the v1 form { "outer": [...], "holes": [[...], ...] }, where "holes" is optional and its
// loops simply follow "outer" in the loop numbering. A vertex without a bulge starts a line.
// "units" is optional and informational. Any other key is rejected.
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

// Writes the v2 schema ("loops"); the output parses back to the same region bit for bit.
std::string write_region_json(const RegionDocument& doc);

// [OUT-06] The medial axis as JSON, in caller units:
//
//   { "units": "mm",
//     "components": [{"outer": 0, "holes": [1, 2], "vertices": [first, count], "edges": [first, count]}],
//     "sites":    [{"loop": 0, "segment": 3, "corner": false, "kind": "line", "normal": [nx, ny], "offset": c}
//                  | {..., "kind": "circle", "centre": [x, y], "radius": R, "sigma": 1}],
//     "vertices": [{"p": [x, y], "r": r, "kind": "junction", "component": k, "edges": [...],
//                   "contacts": [{"site": i, "foot": [x, y]}]}],
//     "edges":    [{"v0": i, "v1": j, "left": a, "right": b, "component": k, "kind": "parabola",
//                   "r0": r0, "r1": r1,
//                   "bezier": [{"p0": [x, y], "p1": [x, y], "p2": [x, y], "w": w, "u": [u0, u1]}],
//                   "polyline": [[x, y], ...]}],
//     "offsets":  [{"distance": d, "components": [[[[x, y, bulge], ...], ...], ...]}],
//     "stats":    {"collapses": n, "merges": n, ...} }
//
// Input loops and segments are numbered as in the input [IN-02]; an offset lists, per component,
// its loops.
//
// A site's distance is normal . x - offset for a line, sigma (R - |x - centre|) for a circle. u
// runs from v0 (radius r0) to v1 (radius r1); "polyline" is present when chord_tol > 0.
struct OffsetLoops {
    double distance = 0.0;
    std::vector<std::vector<Loop>> components;  // as inward_offset returns them
};
struct MatJsonOptions {
    std::string units;
    double chord_tol = 0.0;  // flatten each edge to this tolerance; 0 leaves the polylines out
};
std::string write_mat_json(const MedialAxis& mat, std::span<const OffsetLoops> offsets = {},
                           const MatJsonOptions& options = {});

} // namespace wfmat::io
