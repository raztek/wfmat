#include <string>

#include <catch2/catch_test_macros.hpp>

#include "support.hpp"
#include "wfmat/io/svg.hpp"

using namespace wfmat;
using namespace wfmat::test;

namespace {

std::size_t occurrences(const std::string& s, const std::string& what)
{
    std::size_t n = 0;
    for (auto pos = s.find(what); pos != std::string::npos; pos = s.find(what, pos + what.size())) ++n;
    return n;
}

} // namespace

TEST_CASE("the boundary renders in caller units", "[OUT-06]")
{
    const PreparedRegion r = prepared(load("spec-example.json"));
    const std::string svg = io::render_svg(r);
    CHECK(svg.rfind("<svg xmlns=\"http://www.w3.org/2000/svg\"", 0) == 0);
    CHECK(svg.find("</svg>") != std::string::npos);
    CHECK(occurrences(svg, "class=\"line\"") == 4);
    CHECK(occurrences(svg, "class=\"arc\"") == 1);
    // The semicircle from (100, 40) to (60, 40), radius 20, counter-clockwise, not a large arc.
    CHECK(svg.find("M100,40 A20,20 0 0 1 60,40") != std::string::npos);
    CHECK(occurrences(svg, "class=\"convex\"") == 3);
    CHECK(occurrences(svg, "class=\"tangent\"") == 1);
    CHECK(occurrences(svg, "class=\"reflex\"") == 1);
}

TEST_CASE("clockwise and large arcs set the SVG flags", "[OUT-06]")
{
    // Concave top edge: a clockwise arc from (100, 40) to (0, 40).
    const std::string concave = io::render_svg(prepared(load("concave-top.json")));
    CHECK(concave.find(" 0 0 0 0,40\"") != std::string::npos);

    const double b = bulge_from_sweep(1.5 * std::numbers::pi);
    Region r;
    r.outer.vertices = {{{0, -1}, b}, {{-1, 0}, 0.0}};  // three quarters of the unit circle plus a chord
    const std::string large = io::render_svg(prepared(r));
    CHECK(large.find(" 0 1 1 ") != std::string::npos);
}

TEST_CASE("disks, points and polylines are drawn", "[OUT-06]")
{
    io::SvgWriter w(Box{{0, 0}, {10, 10}});
    w.disk({5, 5}, 2);
    w.point({1, 1});
    const Vec2 pts[] = {{0, 0}, {1, 2}, {3, 4}};
    w.polyline(pts);
    const std::string svg = w.str();
    CHECK(svg.find("class=\"disk\" cx=\"5\" cy=\"5\" r=\"2\"") != std::string::npos);
    CHECK(svg.find("points=\"0,0 1,2 3,4\"") != std::string::npos);
    CHECK(svg.find("class=\"point\"") != std::string::npos);
    CHECK(svg.find("viewBox=\"-0.5 -10.5 11 11\"") != std::string::npos);
}
