// Shared helpers for the test suite.
#pragma once

#include <filesystem>
#include <string>

#include <catch2/catch_test_macros.hpp>

#include "wfmat/io/json.hpp"
#include "wfmat/mat.hpp"
#include "wfmat/prepare.hpp"

namespace wfmat::test {

inline std::filesystem::path data_path(const std::string& name)
{
    return std::filesystem::path(WFMAT_TEST_DATA) / name;
}

inline std::filesystem::path output_path(const std::string& name)
{
    const std::filesystem::path dir(WFMAT_TEST_OUTPUT);
    std::filesystem::create_directories(dir);
    return dir / name;
}

inline Region load(const std::string& name)
{
    auto doc = io::read_region_json(data_path(name));
    REQUIRE(doc.has_value());
    return doc->region;
}

inline PreparedRegion prepared(const Region& region, const Options& options = {})
{
    auto r = prepare(region, options);
    if (!r) FAIL(to_string(r.error()));
    return *r;
}

// The inward offset with the loops of every component in one list.
inline Result<std::vector<Loop>> offset_loops(const PreparedRegion& region, double distance)
{
    auto comps = inward_offset(region, distance);
    if (!comps) return tl::unexpected(comps.error());
    std::vector<Loop> out;
    for (const auto& c : *comps) out.insert(out.end(), c.begin(), c.end());
    return out;
}

inline Region polygon(std::initializer_list<Vec2> points)
{
    Region region{{Loop{}}};
    for (Vec2 p : points) region.loops[0].vertices.push_back({p, 0.0});
    return region;
}

} // namespace wfmat::test
