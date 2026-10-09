#include <cmath>
#include <numbers>

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include "wfmat/geom.hpp"
#include "wfmat/region.hpp"

using namespace wfmat;
using Catch::Approx;

namespace {
constexpr double pi = std::numbers::pi;
bool near(Vec2 a, Vec2 b, double eps = 1e-12) { return dist(a, b) <= eps; }
} // namespace

TEST_CASE("bulge 1 is a counter-clockwise half circle", "[IN-01]")
{
    const Segment s = make_segment({1, 0}, {-1, 0}, 1.0);
    CHECK(s.is_arc());
    CHECK(s.sweep == Approx(pi));
    CHECK(s.R == Approx(1.0));
    CHECK(near(s.c, {0, 0}));
    CHECK(near(s.midpoint(), {0, 1}));
    CHECK(near(s.start_tangent(), {0, 1}));
    CHECK(near(s.end_tangent(), {0, -1}));
    CHECK(s.length() == Approx(pi));
    CHECK(s.cap_area() == Approx(pi / 2));
}

TEST_CASE("negative bulge is clockwise, centre on the right of the chord", "[IN-01]")
{
    const double sweep = -pi / 2;
    const Segment s = make_segment({0, 0}, {1, 1}, bulge_from_sweep(sweep));
    CHECK(s.sweep == Approx(sweep));
    CHECK(near(s.c, {1, 0}));
    CHECK(s.R == Approx(1.0));
    CHECK(near(s.point(0.5), Vec2{1, 0} + polar(3 * pi / 4)));
    CHECK(s.cap_area() == Approx(-(pi / 4 - 0.5)));
}

TEST_CASE("major arcs put the centre on the far side", "[IN-01]")
{
    const Segment s = make_segment({1, 0}, {0, 1}, bulge_from_sweep(-1.5 * pi));  // clockwise, 270 degrees
    CHECK(near(s.c, {0, 0}));
    CHECK(s.R == Approx(1.0));
    CHECK(near(s.midpoint(), polar(-0.75 * pi)));
    const Box b = s.bbox();
    CHECK(near(b.lo, {-1, -1}));
    CHECK(near(b.hi, {1, 1}));
}

TEST_CASE("explicit arcs convert to bulge form", "[IN-01]")
{
    const SegmentSpec specs[] = {ArcSpec{{0, 0}, 2.0, 0.0, pi}, LineSpec{{-2, 0}, {2, 0}}};
    auto loop = loop_from_segments(specs, 1e-12);
    REQUIRE(loop);
    REQUIRE(loop->vertices.size() == 2);
    CHECK(near(loop->vertices[0].p, {2, 0}));
    CHECK(loop->vertices[0].bulge == Approx(1.0));
    CHECK(near(loop->vertices[1].p, {-2, 0}));
    CHECK(loop->vertices[1].bulge == 0.0);

    const SegmentSpec gap[] = {LineSpec{{0, 0}, {1, 0}}, LineSpec{{1, 0.1}, {0, 0}}};
    auto bad = loop_from_segments(gap, 1e-12);
    REQUIRE_FALSE(bad);
    CHECK(bad.error().requirement == "IN-01");
    CHECK(bad.error().segment == 0u);
}

TEST_CASE("closest points on lines and arcs", "[geom]")
{
    const Segment line = make_segment({0, 0}, {2, 0}, 0.0);
    CHECK(near(line.closest_point({1, 5}), {1, 0}));
    CHECK(near(line.closest_point({-3, 1}), {0, 0}));

    const Segment arc = make_segment({1, 0}, {-1, 0}, 1.0);
    CHECK(near(arc.closest_point({0, 3}), {0, 1}));
    CHECK(near(arc.closest_point({0, -3}), {1, 0}));  // outside the angular range: nearer end point
    CHECK(near(arc.closest_point({0, 0}), {0, 1}));   // centre: the midpoint stands for every point
    CHECK(arc.distance({0, 0.25}) == Approx(0.75));
}

TEST_CASE("segment intersections", "[IN-04]")
{
    const double eps = 1e-10;
    const Segment a = make_segment({0, 0}, {2, 2}, 0.0);
    const Segment b = make_segment({0, 2}, {2, 0}, 0.0);
    auto x = intersect(a, b, eps);
    REQUIRE(x.points.size() == 1);
    CHECK(near(x.points[0], {1, 1}));

    SECTION("parallel lines only meet when collinear and overlapping")
    {
        CHECK(intersect(a, make_segment({0, 1}, {2, 3}, 0.0), eps).points.empty());
        CHECK(intersect(a, make_segment({1, 1}, {3, 3}, 0.0), eps).overlap);
        auto touch = intersect(a, make_segment({2, 2}, {3, 3}, 0.0), eps);
        CHECK_FALSE(touch.overlap);
        CHECK(touch.points.size() == 1);
    }
    SECTION("near misses within eps count as contact")
    {
        CHECK(intersect(make_segment({0, 0}, {2, 0}, 0.0), make_segment({1, 0.5e-10}, {1, 1}, 0.0), eps).points.size() == 1);
        CHECK(intersect(make_segment({0, 0}, {2, 0}, 0.0), make_segment({1, 1e-9}, {1, 1}, 0.0), eps).points.empty());
    }
    SECTION("line and arc")
    {
        const Segment arc = make_segment({1, 0}, {-1, 0}, 1.0);  // upper half of the unit circle
        CHECK(intersect(arc, make_segment({-2, 0.5}, {2, 0.5}, 0.0), eps).points.size() == 2);
        CHECK(intersect(arc, make_segment({-2, -0.5}, {2, -0.5}, 0.0), eps).points.empty());
        CHECK(intersect(arc, make_segment({-2, 1}, {2, 1}, 0.0), eps).points.size() == 1);  // tangent
    }
    SECTION("arc and arc")
    {
        const Segment upper = make_segment({1, 0}, {-1, 0}, 1.0);
        const Segment shifted = make_segment({2, 0}, {0, 0}, 1.0);  // upper half about (1, 0)
        auto xs = intersect(upper, shifted, eps);
        REQUIRE(xs.points.size() == 1);
        CHECK(near(xs.points[0], {0.5, std::sqrt(3.0) / 2}));
        CHECK_FALSE(xs.overlap);
        const Segment lower = make_segment({-1, 0}, {1, 0}, 1.0);
        auto halves = intersect(upper, lower, eps);
        CHECK_FALSE(halves.overlap);
        CHECK(halves.points.size() == 2);
        CHECK(intersect(upper, make_segment({0, 1}, {-1, 0}, bulge_from_sweep(pi / 2)), eps).overlap);
    }
}
