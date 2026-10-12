// M4.1: regions of several loops: disjointness, nesting, orientation and components.
#include <algorithm>
#include <cmath>
#include <fstream>
#include <numbers>
#include <string>

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include "support.hpp"
#include "wfmat/io/svg.hpp"
#include "wfmat/mat.hpp"
#include "wfmat/validate.hpp"

using namespace wfmat;
using namespace wfmat::test;
using Catch::Approx;

namespace {

constexpr double pi = std::numbers::pi;

// An axis-aligned square of half-width h about c, counter-clockwise or clockwise.
Loop square(Vec2 c, double h, bool ccw)
{
    Loop loop;
    const Vec2 corners[] = {{-h, -h}, {h, -h}, {h, h}, {-h, h}};
    for (int k = 0; k < 4; ++k) loop.vertices.push_back({c + corners[ccw ? k : 3 - k], 0.0});
    return loop;
}

// A circle of radius r about c as two half arcs.
Loop circle(Vec2 c, double r, bool ccw)
{
    const double b = ccw ? 1.0 : -1.0;
    return Loop{{{c + Vec2{r, 0}, b}, {c - Vec2{r, 0}, b}}};
}

// Concentric squares, alternating in orientation from counter-clockwise outside.
Region nested_squares(int levels)
{
    Region r;
    for (int k = 0; k < levels; ++k) r.loops.push_back(square({0, 0}, 100.0 - 10.0 * k, k % 2 == 0));
    return r;
}

std::size_t occurrences(const std::string& s, const std::string& what)
{
    std::size_t n = 0;
    for (auto at = s.find(what); at != std::string::npos; at = s.find(what, at + 1)) ++n;
    return n;
}

} // namespace

TEST_CASE("a square ring has one component with one hole", "[IN-02][IN-13][IN-14][IN-08]")
{
    const PreparedRegion r = prepared(load("square-ring.json"));
    REQUIRE(r.loops.size() == 2);
    CHECK(r.loops[0].depth == 0);
    CHECK(r.loops[1].depth == 1);
    CHECK(r.loops[1].parent == 0);
    CHECK(r.loops[0].parent == PreparedLoop::no_loop);
    REQUIRE(r.components.size() == 1);
    CHECK(r.components[0].outer == 0);
    CHECK(r.components[0].holes == std::vector<std::uint32_t>{1});
    CHECK(r.loops[1].first == 4);
    CHECK(r.loops[1].count == 4);
    CHECK(r.loop_segments(1).size() == 4);
    CHECK(r.loops[0].area > 0.0);
    CHECK(r.loops[1].area < 0.0);
    const double s2 = r.xf.scale * r.xf.scale;
    CHECK(r.area * s2 == Approx(100.0 * 100.0 - 40.0 * 40.0));

    // The hole's corners point into the material: reflex.
    for (std::uint32_t i = 0; i < 4; ++i) CHECK(r.joins[i].kind == JoinKind::convex);
    for (std::uint32_t i = 4; i < 8; ++i) CHECK(r.joins[i].kind == JoinKind::reflex);

    CHECK(contains(r, r.xf.to_unit({10, 50}), 1e-12));
    CHECK_FALSE(contains(r, r.xf.to_unit({50, 50}), 1e-12));
    CHECK_FALSE(contains(r, r.xf.to_unit({150, 50}), 1e-12));

    const std::string svg = io::render_svg(r);
    CHECK(occurrences(svg, "<path ") == r.segments.size());
    std::ofstream(output_path("square-ring.json.svg")) << svg;
}

TEST_CASE("a plate with an arc hole, a slot and an island", "[IN-02][IN-13][IN-14]")
{
    const PreparedRegion r = prepared(load("plate-with-holes.json"));
    REQUIRE(r.loops.size() == 4);
    CHECK(r.loops[1].depth == 1);
    CHECK(r.loops[2].depth == 1);
    CHECK(r.loops[3].depth == 2);
    CHECK(r.loops[3].parent == 2);
    REQUIRE(r.components.size() == 2);
    CHECK(r.components[0].outer == 0);
    CHECK(r.components[0].holes == std::vector<std::uint32_t>{1, 2});
    CHECK(r.components[1].outer == 3);
    CHECK(r.components[1].holes.empty());
    CHECK(r.loops[3].component == 1);
    // The circular hole keeps two half arcs joined tangentially, as the v1 disk does.
    CHECK(r.loops[1].count == 2);
    CHECK(r.joins[r.loops[1].first].kind == JoinKind::tangent);

    const double s2 = r.xf.scale * r.xf.scale;
    const double plate = 4000.0 + 200.0 * pi, hole = 25.0 * pi, slot = 15.0 * 20.0, island = pi;
    CHECK(r.area * s2 == Approx(plate - hole - slot + island));
    CHECK(contains(r, r.xf.to_unit({48, 15}), 1e-12));   // on the island
    CHECK_FALSE(contains(r, r.xf.to_unit({42, 15}), 1e-12));  // in the slot
    CHECK_FALSE(contains(r, r.xf.to_unit({20, 20}), 1e-12));  // in the circular hole
    std::ofstream(output_path("plate-with-holes.json.svg")) << io::render_svg(r);
}

TEST_CASE("loops may come in any order", "[IN-13][VER-04]")
{
    Region r = load("plate-with-holes.json");
    std::reverse(r.loops.begin(), r.loops.end());
    const PreparedRegion p = prepared(r);
    // Input loop 0 is now the island: components follow the input order of their outer loops.
    REQUIRE(p.components.size() == 2);
    CHECK(p.components[0].outer == 0);
    CHECK(p.loops[0].depth == 2);
    CHECK(p.components[1].outer == 3);
    CHECK(p.components[1].holes == std::vector<std::uint32_t>{1, 2});
    CHECK(p.area == Approx(prepared(load("plate-with-holes.json")).area).epsilon(1e-14));
}

TEST_CASE("separate outer loops are separate components", "[IN-13]")
{
    Region r;
    r.loops.push_back(square({0, 0}, 10, true));
    r.loops.push_back(square({30, 0}, 10, true));
    const PreparedRegion p = prepared(r);
    REQUIRE(p.components.size() == 2);
    CHECK(p.loops[0].depth == 0);
    CHECK(p.loops[1].depth == 0);
    CHECK(p.loops[1].component == 1);
}

TEST_CASE("a hole must run clockwise, an island counter-clockwise", "[IN-08]")
{
    Region r;
    r.loops.push_back(square({0, 0}, 100, true));
    r.loops.push_back(square({0, 0}, 50, true));
    auto bad = prepare(r);
    REQUIRE_FALSE(bad);
    CHECK(bad.error().requirement == "IN-08");
    CHECK(bad.error().loop == 1u);
    CHECK(bad.error().message.find("depth 1 must run clockwise") != std::string::npos);

    r.loops[1] = square({0, 0}, 50, false);
    r.loops.push_back(circle({0, 0}, 20, false));
    bad = prepare(r);
    REQUIRE_FALSE(bad);
    CHECK(bad.error().loop == 2u);
    CHECK(bad.error().message.find("depth 2 must run counter-clockwise") != std::string::npos);

    r.loops[2] = circle({0, 0}, 20, true);
    CHECK(prepare(r).has_value());
}

TEST_CASE("loops that cross or touch are rejected", "[IN-12]")
{
    SECTION("crossing")
    {
        Region r;
        r.loops.push_back(square({0, 0}, 10, true));
        r.loops.push_back(square({10, 0}, 3, false));
        auto bad = prepare(r);
        REQUIRE_FALSE(bad);
        CHECK(bad.error().requirement == "IN-12");
        CHECK(bad.error().loop == 0u);
    }
    SECTION("a hole touching the outer loop at one point")
    {
        Region r;
        r.loops.push_back(square({0, 0}, 10, true));
        r.loops.push_back(Loop{{{{10, 0}, 0.0}, {{5, -2}, 0.0}, {{5, 2}, 0.0}}});
        auto bad = prepare(r);
        REQUIRE_FALSE(bad);
        CHECK(bad.error().requirement == "IN-12");
    }
    SECTION("a circular hole tangent to the outer loop")
    {
        Region r;
        r.loops.push_back(square({0, 0}, 10, true));
        r.loops.push_back(circle({5, 0}, 5, false));
        auto bad = prepare(r);
        REQUIRE_FALSE(bad);
        CHECK(bad.error().requirement == "IN-12");
    }
    SECTION("two holes sharing an edge")
    {
        Region r;
        r.loops.push_back(square({0, 0}, 10, true));
        r.loops.push_back(square({-2, 0}, 2, false));
        r.loops.push_back(square({2, 0}, 2, false));
        auto bad = prepare(r);
        REQUIRE_FALSE(bad);
        CHECK(bad.error().requirement == "IN-12");
        CHECK(bad.error().loop == 1u);
    }
    SECTION("a thin wall is fine")
    {
        Region r;
        r.loops.push_back(square({0, 0}, 10, true));
        r.loops.push_back(square({0, 0}, 10 - 1e-6, false));
        CHECK(prepare(r).has_value());
    }
}

TEST_CASE("nesting is limited to 8 levels", "[IN-13]")
{
    const PreparedRegion eight = prepared(nested_squares(8));
    REQUIRE(eight.loops.size() == 8);
    for (std::uint32_t k = 0; k < 8; ++k) {
        CHECK(eight.loops[k].depth == k);
        CHECK(eight.loops[k].component == k / 2);
        if (k > 0) CHECK(eight.loops[k].parent == k - 1);
    }
    REQUIRE(eight.components.size() == 4);
    CHECK(eight.components[3].outer == 6);
    CHECK(eight.components[3].holes == std::vector<std::uint32_t>{7});

    auto nine = prepare(nested_squares(9));
    REQUIRE_FALSE(nine);
    CHECK(nine.error().requirement == "IN-13");
    CHECK(nine.error().loop == 8u);
}

TEST_CASE("loops are checked one by one", "[IN-03][IN-04][IN-01]")
{
    Region r;
    r.loops.push_back(square({0, 0}, 10, true));
    r.loops.push_back(Loop{{{{-5, -5}, 0.0}, {{-5, 5}, 0.0}, {{5, -5}, 0.0}, {{5, 5}, 0.0}}});  // a bow tie
    auto bad = prepare(r);
    REQUIRE_FALSE(bad);
    CHECK(bad.error().requirement == "IN-04");
    CHECK(bad.error().loop == 1u);

    r.loops[1] = Loop{{{{1, 1}, 0.0}}};
    bad = prepare(r);
    REQUIRE_FALSE(bad);
    CHECK(bad.error().requirement == "IN-03");
    CHECK(bad.error().loop == 1u);

    CHECK(prepare(Region{}).error().requirement == "IN-12");
}

TEST_CASE("the axis of a region with holes waits for M4.2", "[API-03]")
{
    const Region r = load("square-ring.json");
    auto mat = compute_mat(r);
    REQUIRE_FALSE(mat);
    CHECK(mat.error().code == ErrorCode::unsupported);
    CHECK(mat.error().requirement == "M4.2");
    auto offset = inward_offset(r, 5.0);
    REQUIRE_FALSE(offset);
    CHECK(offset.error().code == ErrorCode::unsupported);
}

TEST_CASE("written documents with several loops read back bit for bit", "[IN-02]")
{
    const io::RegionDocument doc{load("plate-with-holes.json"), "mm"};
    auto back = io::parse_region_json(io::write_region_json(doc));
    REQUIRE(back);
    REQUIRE(back->region.loops.size() == doc.region.loops.size());
    for (std::size_t k = 0; k < doc.region.loops.size(); ++k) {
        const auto& a = doc.region.loops[k].vertices;
        const auto& b = back->region.loops[k].vertices;
        REQUIRE(a.size() == b.size());
        for (std::size_t i = 0; i < a.size(); ++i) {
            CHECK(a[i].p == b[i].p);
            CHECK(a[i].bulge == b[i].bulge);
        }
    }
}
