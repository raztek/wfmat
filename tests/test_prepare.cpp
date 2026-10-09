#include <algorithm>
#include <cmath>
#include <numbers>

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include "support.hpp"

using namespace wfmat;
using namespace wfmat::test;
using Catch::Approx;

namespace {

constexpr double pi = std::numbers::pi;

Error rejected(const Region& region, const Options& options = {})
{
    auto r = prepare(region, options);
    REQUIRE_FALSE(r.has_value());
    CHECK(r.error().code == ErrorCode::invalid_input);
    CHECK(r.error().loop == 0u);
    return r.error();
}

std::size_t count(const PreparedRegion& r, JoinKind kind)
{
    std::size_t n = 0;
    for (const Join& j : r.joins) n += j.kind == kind;
    return n;
}

} // namespace

TEST_CASE("segments shorter than eps_len are rejected", "[IN-03]")
{
    const Error e = rejected(polygon({{0, 0}, {10, 0}, {10, 1e-12}, {10, 10}, {0, 10}}));
    CHECK(e.requirement == "IN-03");
    CHECK(e.segment == 1u);
}

TEST_CASE("arcs with a sweep below eps_ang are rejected", "[IN-03]")
{
    Region r = polygon({{0, 0}, {10, 0}, {10, 10}, {0, 10}});
    r.outer.vertices[1].bulge = 1e-12;
    const Error e = rejected(r);
    CHECK(e.requirement == "IN-03");
    CHECK(e.segment == 1u);
}

TEST_CASE("non-finite input is rejected", "[IN-01]")
{
    Region r = polygon({{0, 0}, {10, 0}, {10, 10}});
    r.outer.vertices[2].p.y = std::nan("");
    CHECK(rejected(r).requirement == "IN-01");
}

TEST_CASE("loops that are not simple are rejected", "[IN-04]")
{
    SECTION("bow tie")
    {
        const Error e = rejected(polygon({{0, 0}, {10, 10}, {10, 0}, {0, 10}}));
        CHECK(e.requirement == "IN-04");
        CHECK(e.segment == 0u);
    }
    SECTION("vertex touching a non-adjacent segment")
    {
        CHECK(rejected(polygon({{0, 0}, {10, 0}, {10, 10}, {5, 0.0}, {0, 10}})).requirement == "IN-04");
    }
    SECTION("collinear overlap of non-adjacent segments")
    {
        CHECK(rejected(polygon({{0, 0}, {10, 0}, {10, 5}, {5, 5}, {5, 0}, {2, 0}, {2, -5}, {0, -5}})).requirement ==
              "IN-04");
    }
    SECTION("arc crossing a line")
    {
        Region r = polygon({{0, 0}, {10, 0}, {10, 2}, {0, 2}});
        r.outer.vertices[2].bulge = -0.6;  // the top edge sags to y = -1, through the bottom edge
        const Error e = rejected(r);
        CHECK(e.requirement == "IN-04");
        CHECK(e.segment == 0u);
    }
    SECTION("a near miss closer than eps_geom counts as contact")
    {
        CHECK(rejected(polygon({{0, 0}, {10, 0}, {10, 10}, {5, 1e-11}, {0, 10}})).requirement == "IN-04");
    }
    SECTION("a gap wider than eps_geom is fine")
    {
        CHECK(prepare(polygon({{0, 0}, {10, 0}, {10, 10}, {5, 1e-6}, {0, 10}})).has_value());
    }
}

TEST_CASE("cusps are rejected", "[IN-06]")
{
    SECTION("two lines folding back")
    {
        CHECK(rejected(polygon({{0, 0}, {10, 0}})).requirement == "IN-06");
    }
    SECTION("a spike")
    {
        const Error e = rejected(polygon({{0, 0}, {10, 0}, {10, 10}, {20, 10}, {0, 10}}));
        CHECK(e.requirement == "IN-06");
    }
    SECTION("an arc leaving tangentially backwards")
    {
        // The line arrives at (10, 0) heading +x; the clockwise arc leaves it heading -x.
        Region r = polygon({{0, 0}, {10, 0}, {10, -4}, {0, -4}});
        r.outer.vertices[1].bulge = 1.0;
        const Error e = rejected(r);
        CHECK(e.requirement == "IN-06");
        CHECK(e.segment == 1u);
    }
}

TEST_CASE("normalisation maps the bounding box to the unit frame", "[IN-07]")
{
    const PreparedRegion r = prepared(load("spec-example.json"));
    // The semicircle on top reaches y = 60, so the box is [0, 100] x [0, 60].
    CHECK(r.xf.centre.x == Approx(50.0));
    CHECK(r.xf.centre.y == Approx(30.0));
    CHECK(r.xf.scale == Approx(0.5 * std::hypot(100.0, 60.0)));
    CHECK(0.5 * norm(r.bbox.extent()) == Approx(1.0));
    CHECK(norm(r.bbox.center()) == Approx(0.0).margin(1e-15));
    for (const Segment& s : r.segments) {
        CHECK(std::abs(s.a.x) <= 1.0);
        CHECK(std::abs(s.a.y) <= 1.0);
    }
    CHECK(r.xf.from_unit(r.xf.to_unit({3, 4})).x == Approx(3.0));
    CHECK(r.area * r.xf.scale * r.xf.scale == Approx(4000.0 + 200.0 * pi));
}

TEST_CASE("clockwise loops are reversed", "[IN-08]")
{
    const PreparedRegion ccw = prepared(polygon({{0, 0}, {4, 0}, {4, 2}, {0, 2}}));
    const PreparedRegion cw = prepared(polygon({{0, 0}, {0, 2}, {4, 2}, {4, 0}}));
    CHECK_FALSE(ccw.reversed);
    CHECK(cw.reversed);
    CHECK(cw.area == Approx(ccw.area));
    CHECK(cw.area > 0.0);
    REQUIRE(cw.sources.size() == 4);
    // Input segment 3 runs (4, 0) -> (0, 0); reversed it is the first counter-clockwise segment.
    CHECK(cw.sources[0] == std::vector<std::uint32_t>{3});
    CHECK(cw.xf.from_unit(cw.segments[0].a).x == Approx(0.0));
    CHECK(cw.xf.from_unit(cw.segments[0].b).x == Approx(4.0));
    CHECK(count(cw, JoinKind::convex) == 4);

    SECTION("arc bulges flip with the direction")
    {
        Region r = load("slot.json");
        std::reverse(r.outer.vertices.begin(), r.outer.vertices.end());
        // Reversing the vertex order moves each bulge to the previous vertex, negated.
        Region cw_slot;
        const auto& v = r.outer.vertices;
        for (std::size_t i = 0; i < v.size(); ++i) cw_slot.outer.vertices.push_back({v[i].p, -v[(i + 1) % v.size()].bulge});
        const PreparedRegion p = prepared(cw_slot);
        CHECK(p.reversed);
        CHECK(p.area * p.xf.scale * p.xf.scale == Approx(1200.0 + 100.0 * pi));
        CHECK(count(p, JoinKind::tangent) == 4);
        for (const Segment& s : p.segments)
            if (s.is_arc()) CHECK(s.sweep > 0.0);
    }
}

TEST_CASE("collinear lines and co-circular arcs are merged", "[IN-09]")
{
    SECTION("collinear lines")
    {
        const PreparedRegion r = prepared(polygon({{0, 0}, {5, 0}, {10, 0}, {10, 10}, {5, 10}, {0, 10}}));
        CHECK(r.segments.size() == 4);
        CHECK(r.sources[0] == std::vector<std::uint32_t>{0, 1});
        CHECK(count(r, JoinKind::convex) == 4);
    }
    SECTION("wrapping around the start of the loop")
    {
        const PreparedRegion r = prepared(polygon({{5, 0}, {10, 0}, {10, 10}, {0, 10}, {0, 0}}));
        CHECK(r.segments.size() == 4);
    }
    SECTION("co-circular arcs")
    {
        // A disk given as four quarter arcs merges into two arcs of 3 pi / 2 and pi / 2: one arc
        // cannot hold a full turn.
        const double b = bulge_from_sweep(pi / 2);
        Region disk;
        disk.outer.vertices = {{{1, 0}, b}, {{0, 1}, b}, {{-1, 0}, b}, {{0, -1}, b}};
        const PreparedRegion r = prepared(disk);
        REQUIRE(r.segments.size() == 2);
        CHECK(std::abs(r.segments[0].sweep) + std::abs(r.segments[1].sweep) == Approx(2 * pi));
        CHECK(count(r, JoinKind::tangent) == 2);
        CHECK(r.area * r.xf.scale * r.xf.scale == Approx(pi));
    }
    SECTION("a disk as two half arcs stays two arcs")
    {
        CHECK(prepared(load("disk.json")).segments.size() == 2);
    }
    SECTION("arcs on different circles are kept")
    {
        CHECK(prepared(load("slot.json")).segments.size() == 4);
    }
}

TEST_CASE("near-tangent joins snap to G1", "[IN-10]")
{
    // A line meeting a half circle with a tangent error of 1e-10 rad (below eps_ang = 1e-9).
    const double tilt = 1e-10;
    Region r;
    r.outer.vertices = {{{0, 0}, 0.0}, {{10, 0}, std::tan((pi + tilt) / 4)}, {{10, 4}, 0.0}, {{0, 4}, 0.0}};
    const PreparedRegion p = prepared(r);
    REQUIRE(p.joins[0].kind == JoinKind::tangent);
    CHECK(p.joins[0].turn == 0.0);
    CHECK(p.joins[0].t_in == p.joins[0].t_out);

    SECTION("a larger error is a real corner")
    {
        Options options;
        options.tol.ang = 1e-11;
        CHECK(prepared(r, options).joins[0].kind != JoinKind::tangent);
    }
}

TEST_CASE("tolerances are configurable", "[IN-11]")
{
    const Region r = polygon({{0, 0}, {10, 0}, {10, 1e-3}, {10, 10}, {0, 10}});
    CHECK(prepare(r).has_value());
    Options strict;
    strict.tol.len = 1e-3;
    CHECK(prepare(r, strict).error().requirement == "IN-03");
}

TEST_CASE("corners are classified by turning angle", "[ALG-01]")
{
    const PreparedRegion l = prepared(load("l-shape.json"));
    CHECK(count(l, JoinKind::convex) == 5);
    CHECK(count(l, JoinKind::reflex) == 1);
    for (const Join& j : l.joins) {
        if (j.kind == JoinKind::reflex) {
            CHECK(l.xf.from_unit(j.p).x == Approx(20.0));
            CHECK(l.xf.from_unit(j.p).y == Approx(20.0));
            CHECK(j.turn == Approx(-pi / 2));
        } else {
            CHECK(j.turn == Approx(pi / 2));
        }
    }

    const PreparedRegion spec = prepared(load("spec-example.json"));
    CHECK(count(spec, JoinKind::convex) == 3);
    CHECK(count(spec, JoinKind::tangent) == 1);
    CHECK(count(spec, JoinKind::reflex) == 1);
}
