#include <catch2/catch_test_macros.hpp>

#include "support.hpp"

using namespace wfmat;
using namespace wfmat::test;

TEST_CASE("the v1 schema parses", "[IN-02]")
{
    auto doc = io::parse_region_json(R"({"units": "mm", "outer": [[0,0,0],[100,0,0],[100,40,1],[60,40,0],[0,40,0]]})");
    REQUIRE(doc);
    CHECK(doc->units == "mm");
    REQUIRE(doc->region.outer.vertices.size() == 5);
    CHECK(doc->region.outer.vertices[2].p == Vec2{100, 40});
    CHECK(doc->region.outer.vertices[2].bulge == 1.0);
}

TEST_CASE("the bulge and the units are optional", "[IN-02]")
{
    auto doc = io::parse_region_json(R"({"outer": [[0,0],[1,0],[0,1,0.5]]})");
    REQUIRE(doc);
    CHECK(doc->units.empty());
    CHECK(doc->region.outer.vertices[0].bulge == 0.0);
    CHECK(doc->region.outer.vertices[2].bulge == 0.5);
}

TEST_CASE("malformed documents are rejected with a diagnostic", "[IN-02]")
{
    const char* bad[] = {
        "not json",
        "[1, 2]",
        R"({"units": "mm"})",
        R"({"outer": {"x": 1}})",
        R"({"outer": [[0,0],[1,0],[0,"1"]]})",
        R"({"outer": [[0,0],[1,0],[0]]})",
        R"({"outer": [[0,0],[1,0],[0,1,0,0]]})",
        R"({"outer": [[0,0],[1,0],[0,1]], "unit": "mm"})",
        R"({"outer": [[0,0],[1,0],[0,1]], "units": 3})",
    };
    for (const char* text : bad) {
        INFO(text);
        auto doc = io::parse_region_json(text);
        REQUIRE_FALSE(doc);
        CHECK(doc.error().requirement == "IN-02");
        CHECK(doc.error().code == ErrorCode::invalid_input);
    }
    auto doc = io::parse_region_json(R"({"outer": [[0,0],[1,0],[0,"1"]]})");
    CHECK(doc.error().segment == 2u);
}

TEST_CASE("regions with holes are rejected until v2", "[IN-05]")
{
    for (const char* text : {R"({"outer": [[0,0],[9,0],[0,9]], "holes": [[[1,1],[2,1],[1,2]]]})",
                             R"({"outer": [[[0,0],[9,0],[0,9]], [[1,1],[2,1],[1,2]]]})"}) {
        INFO(text);
        auto doc = io::parse_region_json(text);
        REQUIRE_FALSE(doc);
        CHECK(doc.error().requirement == "IN-05");
        CHECK(doc.error().message.find("v2") != std::string::npos);
    }
    CHECK(io::parse_region_json(R"({"outer": [[0,0],[9,0],[0,9]], "holes": []})").has_value());
}

TEST_CASE("written documents read back bit for bit", "[IN-02]")
{
    const io::RegionDocument doc{load("filleted-rectangle.json"), "mm"};
    auto back = io::parse_region_json(io::write_region_json(doc));
    REQUIRE(back);
    CHECK(back->units == "mm");
    REQUIRE(back->region.outer.vertices.size() == doc.region.outer.vertices.size());
    for (std::size_t i = 0; i < doc.region.outer.vertices.size(); ++i) {
        CHECK(back->region.outer.vertices[i].p == doc.region.outer.vertices[i].p);
        CHECK(back->region.outer.vertices[i].bulge == doc.region.outer.vertices[i].bulge);
    }
}

TEST_CASE("missing files are reported", "[IN-02]")
{
    auto doc = io::read_region_json(data_path("no-such-file.json"));
    REQUIRE_FALSE(doc);
    CHECK(doc.error().message.find("cannot open") != std::string::npos);
}
