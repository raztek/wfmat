// [VER-02] The canonical shapes load, validate and render: the M0 exit test.
#include <cmath>
#include <fstream>
#include <numbers>
#include <string>

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/generators/catch_generators_range.hpp>

#include "support.hpp"
#include "wfmat/io/svg.hpp"
#include "wfmat/validate.hpp"

using namespace wfmat;
using namespace wfmat::test;
using Catch::Approx;

namespace {

constexpr double pi = std::numbers::pi;

struct Shape {
    const char* file;
    std::size_t lines, arcs;
    std::size_t convex, tangent, reflex;
    double area;
    Vec2 interior;  // a point known to be inside
};

const Shape shapes[] = {
    {"disk.json", 0, 2, 0, 2, 0, 100 * pi, {0, 0}},
    {"rectangle.json", 4, 0, 4, 0, 0, 4000, {50, 20}},
    {"equilateral-triangle.json", 3, 0, 3, 0, 0, 2500 * std::sqrt(3.0), {50, 28}},
    {"regular-hexagon.json", 6, 0, 6, 0, 0, 3750 * std::sqrt(3.0), {0, 0}},
    {"slot.json", 2, 2, 0, 4, 0, 1200 + 100 * pi, {30, 10}},
    {"l-shape.json", 6, 0, 5, 0, 1, 2000, {10, 50}},
    {"filleted-rectangle.json", 4, 1, 3, 2, 0, 4000 - 100 * (1 - pi / 4), {95, 35}},
    {"spec-example.json", 4, 1, 3, 1, 1, 4000 + 200 * pi, {80, 50}},
    {"concave-top.json", 3, 1, 4, 0, 0, 0.0, {50, 20}},
};

std::size_t count(const PreparedRegion& r, JoinKind kind)
{
    std::size_t n = 0;
    for (const Join& j : r.joins) n += j.kind == kind;
    return n;
}

std::size_t occurrences(const std::string& s, const std::string& what)
{
    std::size_t n = 0;
    for (auto pos = s.find(what); pos != std::string::npos; pos = s.find(what, pos + what.size())) ++n;
    return n;
}

// The region transformed by p -> s R(theta) p + t, optionally with the loop reversed.
Region transformed(const Region& in, double s, double theta, Vec2 t, bool reverse)
{
    Region out;
    const auto& v = in.outer.vertices;
    const std::size_t n = v.size();
    for (std::size_t k = 0; k < n; ++k) {
        const std::size_t i = reverse ? n - 1 - k : k;
        const Vec2 p = v[i].p;
        const Vec2 q{s * (std::cos(theta) * p.x - std::sin(theta) * p.y) + t.x,
                     s * (std::sin(theta) * p.x + std::cos(theta) * p.y) + t.y};
        // Reversed, vertex i starts the reverse of input segment i - 1.
        const double bulge = reverse ? -v[(i + n - 1) % n].bulge : v[i].bulge;
        out.outer.vertices.push_back({q, bulge});
    }
    return out;
}

} // namespace

TEST_CASE("canonical shapes load, validate and render", "[VER-02][IN-02][IN-07][IN-08][ALG-01][OUT-06]")
{
    const Shape& shape = GENERATE(from_range(std::begin(shapes), std::end(shapes)));
    INFO(shape.file);
    const PreparedRegion r = prepared(load(shape.file));

    std::size_t arcs = 0;
    for (const Segment& s : r.segments) arcs += s.is_arc();
    CHECK(r.segments.size() - arcs == shape.lines);
    CHECK(arcs == shape.arcs);
    CHECK(count(r, JoinKind::convex) == shape.convex);
    CHECK(count(r, JoinKind::tangent) == shape.tangent);
    CHECK(count(r, JoinKind::reflex) == shape.reflex);
    if (shape.area > 0.0) CHECK(r.area * r.xf.scale * r.xf.scale == Approx(shape.area).epsilon(1e-13));
    CHECK(contains(r, r.xf.to_unit(shape.interior), 1e-12));
    CHECK(0.5 * norm(r.bbox.extent()) == Approx(1.0));

    const std::string svg = io::render_svg(r);
    CHECK(occurrences(svg, "<path ") == r.segments.size());
    CHECK(occurrences(svg, "<circle ") == r.joins.size());
    CHECK(svg.find("nan") == std::string::npos);
    std::ofstream(output_path(std::string(shape.file) + ".svg")) << svg;
}

TEST_CASE("validation is invariant under similarity and reversal", "[VER-04][IN-07][IN-08]")
{
    const Shape& shape = GENERATE(from_range(std::begin(shapes), std::end(shapes)));
    INFO(shape.file);
    const Region base = load(shape.file);
    const PreparedRegion r0 = prepared(base);
    const double area0 = r0.area * r0.xf.scale * r0.xf.scale;

    for (const bool reverse : {false, true}) {
        INFO("reverse " << reverse);
        const double s = 3.5, theta = 0.7;
        const PreparedRegion r = prepared(transformed(base, s, theta, {-1e3, 250}, reverse));
        CHECK(r.reversed == (r0.reversed != reverse));
        CHECK(r.segments.size() == r0.segments.size());
        CHECK(count(r, JoinKind::convex) == count(r0, JoinKind::convex));
        CHECK(count(r, JoinKind::tangent) == count(r0, JoinKind::tangent));
        CHECK(count(r, JoinKind::reflex) == count(r0, JoinKind::reflex));
        CHECK(r.area * r.xf.scale * r.xf.scale == Approx(area0 * s * s).epsilon(1e-12));
    }
}
