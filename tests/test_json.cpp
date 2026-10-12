#include <catch2/catch_test_macros.hpp>

#include "support.hpp"

using namespace wfmat;
using namespace wfmat::test;

TEST_CASE("the v1 schema parses", "[IN-02]")
{
    auto doc = io::parse_region_json(R"({"units": "mm", "outer": [[0,0,0],[100,0,0],[100,40,1],[60,40,0],[0,40,0]]})");
    REQUIRE(doc);
    CHECK(doc->units == "mm");
    REQUIRE(doc->region.loops[0].vertices.size() == 5);
    CHECK(doc->region.loops[0].vertices[2].p == Vec2{100, 40});
    CHECK(doc->region.loops[0].vertices[2].bulge == 1.0);
}

TEST_CASE("the bulge and the units are optional", "[IN-02]")
{
    auto doc = io::parse_region_json(R"({"outer": [[0,0],[1,0],[0,1,0.5]]})");
    REQUIRE(doc);
    CHECK(doc->units.empty());
    CHECK(doc->region.loops[0].vertices[0].bulge == 0.0);
    CHECK(doc->region.loops[0].vertices[2].bulge == 0.5);
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

TEST_CASE("several loops read as \"loops\" or as \"outer\" with \"holes\"", "[IN-02]")
{
    const char* forms[] = {R"({"outer": [[0,0],[9,0],[0,9]], "holes": [[[1,1],[1,2],[2,1]]]})",
                           R"({"loops": [[[0,0],[9,0],[0,9]], [[1,1],[1,2],[2,1]]]})"};
    for (const char* text : forms) {
        INFO(text);
        auto doc = io::parse_region_json(text);
        REQUIRE(doc);
        REQUIRE(doc->region.loops.size() == 2);
        CHECK(doc->region.loops[1].vertices[1].p == Vec2{1, 2});
    }
    auto empty_holes = io::parse_region_json(R"({"outer": [[0,0],[9,0],[0,9]], "holes": []})");
    REQUIRE(empty_holes);
    CHECK(empty_holes->region.loops.size() == 1);

    const char* bad[] = {
        R"({"outer": [[[0,0],[9,0],[0,9]], [[1,1],[2,1],[1,2]]]})",
        R"({"loops": [[[0,0],[9,0],[0,9]]], "outer": [[0,0],[9,0],[0,9]]})",
        R"({"loops": []})",
        R"({"loops": [[[0,0],[9,0],[0,9]], [[1,1],[2,1],[1,"2"]]]})",
        R"({"holes": [[[1,1],[2,1],[1,2]]]})",
    };
    for (const char* text : bad) {
        INFO(text);
        auto doc = io::parse_region_json(text);
        REQUIRE_FALSE(doc);
        CHECK(doc.error().requirement == "IN-02");
    }
    auto doc = io::parse_region_json(bad[3]);
    CHECK(doc.error().loop == 1u);
    CHECK(doc.error().segment == 2u);
}

TEST_CASE("written documents read back bit for bit", "[IN-02]")
{
    const io::RegionDocument doc{load("filleted-rectangle.json"), "mm"};
    auto back = io::parse_region_json(io::write_region_json(doc));
    REQUIRE(back);
    CHECK(back->units == "mm");
    REQUIRE(back->region.loops[0].vertices.size() == doc.region.loops[0].vertices.size());
    for (std::size_t i = 0; i < doc.region.loops[0].vertices.size(); ++i) {
        CHECK(back->region.loops[0].vertices[i].p == doc.region.loops[0].vertices[i].p);
        CHECK(back->region.loops[0].vertices[i].bulge == doc.region.loops[0].vertices[i].bulge);
    }
}

TEST_CASE("missing files are reported", "[IN-02]")
{
    auto doc = io::read_region_json(data_path("no-such-file.json"));
    REQUIRE_FALSE(doc);
    CHECK(doc.error().message.find("cannot open") != std::string::npos);
}
