#include <cmath>
#include <numbers>

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include "support.hpp"
#include "wfmat/validate.hpp"

using namespace wfmat;
using namespace wfmat::test;
using Catch::Approx;

namespace {

// Distance check of a caller-frame disk (p, r).
DiskCheck check(const PreparedRegion& region, Vec2 p, double r)
{
    return check_disk(region, region.xf.to_unit(p), region.xf.length_to_unit(r), 1e-9);
}

} // namespace

TEST_CASE("distance to the boundary and its feet", "[OUT-07][VER-04]")
{
    const PreparedRegion rect = prepared(load("rectangle.json"));
    const double L = rect.xf.scale;

    SECTION("points on the plateau of the known axis have two feet at r = 20")
    {
        for (double x : {20.0, 35.0, 50.0, 80.0}) {
            const DiskCheck c = check(rect, {x, 20}, 20);
            CHECK(c.residual * L < 1e-12);
            CHECK(c.feet == (x == 20.0 || x == 80.0 ? 3u : 2u));  // the plateau ends touch three sides
            CHECK(c.inside);
        }
    }
    SECTION("points on the corner bisectors have two feet")
    {
        const DiskCheck c = check(rect, {10, 10}, 10);
        CHECK(c.residual * L < 1e-12);
        CHECK(c.feet == 2u);
    }
    SECTION("points off the axis have one foot")
    {
        const DiskCheck c = check(rect, {50, 10}, 10);
        CHECK(c.residual * L < 1e-12);
        CHECK(c.feet == 1u);
    }
    SECTION("a wrong radius shows in the residual")
    {
        CHECK(check(rect, {50, 20}, 19).residual * L == Approx(1.0));
    }
}

TEST_CASE("arc centres see the whole arc", "[OUT-07][VER-04]")
{
    const PreparedRegion disk = prepared(load("disk.json"));
    const BoundaryDistance bd = boundary_distance(disk.segments, disk.xf.to_unit({0, 0}), 1e-9);
    CHECK(bd.distance * disk.xf.scale == Approx(10.0));
    CHECK(bd.distinct_feet() >= 2u);

    const PreparedRegion slot = prepared(load("slot.json"));
    const DiskCheck c = check(slot, {60, 10}, 10);  // a curvature end of the slot's axis
    CHECK(c.residual * slot.xf.scale < 1e-12);
    CHECK(c.feet >= 2u);
}

TEST_CASE("containment follows lines, convex arcs and concave arcs", "[VER-04]")
{
    const PreparedRegion slot = prepared(load("slot.json"));
    auto inside = [](const PreparedRegion& r, Vec2 p) { return contains(r, r.xf.to_unit(p), 1e-12); };
    CHECK(inside(slot, {30, 10}));
    CHECK(inside(slot, {69, 10}));    // in the right cap
    CHECK_FALSE(inside(slot, {71, 10}));
    CHECK(inside(slot, {-9, 10}));    // in the left cap
    CHECK_FALSE(inside(slot, {30, 21}));
    CHECK_FALSE(inside(slot, {30, 20}));  // on the boundary

    const PreparedRegion concave = prepared(load("concave-top.json"));
    CHECK(inside(concave, {50, 25}));
    CHECK_FALSE(inside(concave, {50, 30}));  // above the dip of the concave arc, which reaches y = 27.5
    CHECK(inside(concave, {1, 39}));
    CHECK(winding_number(concave.segments, concave.xf.to_unit({50, 10})) == 1);
    CHECK(winding_number(concave.segments, concave.xf.to_unit({50, 100})) == 0);
}
