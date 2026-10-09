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

} // namespace wfmat::io
