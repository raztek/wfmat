// [M7] Edge export and inward offsets: rational quadratic Bezier arcs and adaptive flattening
// [OUT-05], and the front snapshot as an inward offset [ALG-06].
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <limits>
#include <numbers>
#include <random>
#include <string>
#include <vector>

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <nlohmann/json.hpp>

#include "mat_checks.hpp"
#include "support.hpp"
#include "wfmat/io/json.hpp"
#include "wfmat/io/svg.hpp"
#include "wfmat/mat.hpp"
#include "wfmat/validate.hpp"

using namespace wfmat;
using namespace wfmat::test;
using Catch::Approx;

namespace {

constexpr double pi = std::numbers::pi;

int shape_count()
{
    if (const char* s = std::getenv("WFMAT_RANDOM_SHAPES")) return std::max(1, std::atoi(s));
    return 300;
}

double segment_distance(Vec2 p, Vec2 a, Vec2 b)
{
    const Vec2 d = b - a;
    const double len2 = dot(d, d);
    const double s = len2 > 0.0 ? std::clamp(dot(p - a, d) / len2, 0.0, 1.0) : 0.0;
    return dist(p, a + s * d);
}

// Every Bezier piece lies on the bisector of the edge's two sites, between the edge's radii, and
// the pieces run end to end from v0 to v1.
std::string check_bezier(const MedialAxis& mat, double scale)
{
    const double tol = 1e-8 * scale;
    for (const MatEdge& e : mat.edges()) {
        const auto pieces = e.bezier();
        if (pieces.empty()) return "an edge without Bezier pieces";
        if (dist(pieces.front().p0, e.point_at(0.0)) > tol || dist(pieces.back().p2, e.point_at(1.0)) > tol)
            return "Bezier ends off the edge ends";
        for (std::size_t i = 0; i < pieces.size(); ++i) {
            const RationalQuadBezier& b = pieces[i];
            if (!(b.w > 0.0)) return "a Bezier weight that is not positive: " + sci(b.w);
            if (i > 0 && dist(pieces[i - 1].p2, b.p0) > tol) return "Bezier pieces that do not join";
            const Site& l = mat.sites()[e.left()].geometry;
            const Site& r = mat.sites()[e.right()].geometry;
            const double lo = std::min(e.r0(), e.r1()) - tol, hi = std::max(e.r0(), e.r1()) + tol;
            for (double s : {0.125, 0.25, 0.5, 0.75, 0.875}) {
                const Vec2 q = b.point(s);
                const double dl = l.distance(q), dr = r.distance(q);
                if (std::abs(dl - dr) > tol) return "a Bezier point off the bisector by " + sci(std::abs(dl - dr));
                if (dl < lo || dl > hi) return "a Bezier point outside the edge's radii";
            }
        }
    }
    return {};
}

// Every point of the edge is within the chord tolerance of its polyline.
std::string check_flatten(const MedialAxis& mat, double chord_tol)
{
    for (const MatEdge& e : mat.edges()) {
        std::vector<Vec2> poly;
        e.flatten(chord_tol, poly);
        if (poly.size() < 2) return "a polyline with fewer than two points";
        if (poly.front() != e.point_at(0.0) || poly.back() != e.point_at(1.0)) return "polyline ends off the edge";
        for (int k = 0; k <= 256; ++k) {
            const Vec2 p = e.point_at(k / 256.0);
            double d = std::numeric_limits<double>::infinity();
            for (std::size_t i = 0; i + 1 < poly.size(); ++i) d = std::min(d, segment_distance(p, poly[i], poly[i + 1]));
            if (d > 1.01 * chord_tol + 1e-12)
                return "an edge point " + sci(d) + " from its polyline (tolerance " + sci(chord_tol) + ", u = " +
                       std::to_string(k / 256.0) + ", kind " + std::to_string(int(e.kind())) + ", " +
                       std::to_string(poly.size()) + " points)";
        }
    }
    return {};
}

// The start and the middle of every piece of every offset loop are at the offset distance from the
// boundary, inside the region.
std::string check_offset(const PreparedRegion& prep, const std::vector<Loop>& loops, double distance)
{
    const double tol = 1e-9;
    const double t = prep.xf.length_to_unit(distance);
    for (const Loop& loop : loops) {
        const auto& v = loop.vertices;
        for (std::size_t i = 0; i < v.size(); ++i) {
            const Vec2 a = v[i].p, b = v[(i + 1) % v.size()].p;
            const Vec2 d = b - a;
            const Vec2 right = norm(d) > 0.0 ? Vec2{d.y, -d.x} / norm(d) : Vec2{};
            const Vec2 mid = 0.5 * (a + b) + (v[i].bulge * 0.5 * norm(d)) * right;
            for (Vec2 p : {a, mid}) {
                const Vec2 q = prep.xf.to_unit(p);
                const double dq = boundary_distance(prep.segments, q, 1e-9).distance;
                if (std::abs(dq - t) > tol) return "an offset point at distance " + sci(dq) + " instead of " + sci(t);
                if (t > tol && !contains(prep, q, tol / 2)) return "an offset point outside the region";
            }
        }
    }
    return {};
}

double loop_area(const Loop& loop)
{
    Region r;
    r.outer = loop;
    auto prep = prepare(r);
    REQUIRE(prep);
    return prep->area * prep->xf.scale * prep->xf.scale;
}

} // namespace

TEST_CASE("Bezier export lies on the bisectors of the canonical shapes", "[M7][OUT-05]")
{
    for (const char* name : {"rectangle.json", "l-shape.json", "disk.json", "slot.json", "filleted-rectangle.json",
                             "concave-top.json", "regular-hexagon.json", "spec-example.json"}) {
        INFO(name);
        const PreparedRegion prep = prepared(load(name));
        auto mat = compute_mat(prep);
        REQUIRE(mat);
        CHECK(check_bezier(*mat, prep.xf.scale).empty());
        CHECK(check_flatten(*mat, 1e-3 * prep.xf.scale).empty());
    }
}

TEST_CASE("Bezier pieces of a line edge, a parabola and a circle plateau", "[M7][OUT-05]")
{
    // The L-shape's reflex corner gives parabolic edges; the slot a plateau line.
    const PreparedRegion prep = prepared(load("l-shape.json"));
    auto mat = compute_mat(prep);
    REQUIRE(mat);
    bool parabola = false;
    for (const MatEdge& e : mat->edges()) {
        const auto pieces = e.bezier();
        if (e.kind() == ConicKind::line) {
            REQUIRE(pieces.size() == 1);
            CHECK(pieces[0].w == 1.0);
        }
        if (e.kind() == ConicKind::parabola) {
            parabola = true;
            for (const auto& b : pieces) CHECK(b.w == Approx(1.0).epsilon(1e-9));
        }
    }
    CHECK(parabola);
}

TEST_CASE("Bezier export and flattening on random shapes with arcs", "[M7][OUT-05][VER-05]")
{
    std::mt19937_64 rng(20261010);
    std::uniform_int_distribution<int> size(3, 40);
    const int n = shape_count();
    for (int k = 0; k < n; ++k) {
        const Region region = k % 2 ? random_filleted(rng, size(rng)) : random_arc_chain(rng, size(rng));
        auto prep = prepare(region);
        if (!prep) continue;
        auto mat = compute_mat(*prep);
        REQUIRE(mat);
        const std::string bez = check_bezier(*mat, prep->xf.scale);
        const std::string flat = check_flatten(*mat, 1e-4 * prep->xf.scale);
        INFO("shape " << k << ": " << bez << flat);
        CHECK(bez.empty());
        CHECK(flat.empty());
        if (!bez.empty() || !flat.empty()) break;
    }
}

TEST_CASE("Inward offsets of the canonical shapes", "[M7][ALG-06]")
{
    SECTION("rectangle: a smaller rectangle, then nothing")
    {
        const PreparedRegion prep = prepared(load("rectangle.json"));  // 100 x 40
        auto loops = inward_offset(prep, 5.0);
        REQUIRE(loops);
        REQUIRE(loops->size() == 1);
        CHECK((*loops)[0].vertices.size() == 4);
        CHECK(loop_area((*loops)[0]) == Approx(90.0 * 30.0).epsilon(1e-12));
        CHECK(check_offset(prep, *loops, 5.0).empty());
        auto none = inward_offset(prep, 25.0);
        REQUIRE(none);
        CHECK(none->empty());
    }
    SECTION("disk: a concentric circle")
    {
        const PreparedRegion prep = prepared(load("disk.json"));  // radius 10
        auto loops = inward_offset(prep, 4.0);
        REQUIRE(loops);
        REQUIRE(loops->size() == 1);
        CHECK(loop_area((*loops)[0]) == Approx(pi * 36.0).epsilon(1e-12));
        CHECK(check_offset(prep, *loops, 4.0).empty());
    }
    SECTION("L-shape: an arc about the reflex corner")
    {
        const PreparedRegion prep = prepared(load("l-shape.json"));
        auto loops = inward_offset(prep, 5.0);
        REQUIRE(loops);
        REQUIRE(loops->size() == 1);
        int arcs = 0;
        for (const BulgeVertex& v : (*loops)[0].vertices) arcs += v.bulge != 0.0;
        CHECK(arcs == 1);
        CHECK(check_offset(prep, *loops, 5.0).empty());
    }
    SECTION("slot: a smaller slot")
    {
        const PreparedRegion prep = prepared(load("slot.json"));  // 60 between centres, radius 10
        auto loops = inward_offset(prep, 5.0);
        REQUIRE(loops);
        REQUIRE(loops->size() == 1);
        CHECK(loop_area((*loops)[0]) == Approx(60.0 * 10.0 + pi * 25.0).epsilon(1e-12));
        CHECK(check_offset(prep, *loops, 5.0).empty());
    }
    SECTION("distance zero is the boundary; a negative distance is refused")
    {
        const PreparedRegion prep = prepared(load("l-shape.json"));
        auto loops = inward_offset(prep, 0.0);
        REQUIRE(loops);
        REQUIRE(loops->size() == 1);
        CHECK(loop_area((*loops)[0]) == Approx(2000.0).epsilon(1e-12));
        auto bad = inward_offset(prep, -1.0);
        REQUIRE_FALSE(bad);
        CHECK(bad.error().requirement == "ALG-06");
    }
}

TEST_CASE("Inward offsets split a region into several loops", "[M7][ALG-06]")
{
    // A dumbbell: two 40 x 40 squares joined by a 10-wide neck; at 8 the neck is gone.
    const Region region = polygon(
        {{0, 0}, {40, 0}, {40, 15}, {60, 15}, {60, 0}, {100, 0}, {100, 40}, {60, 40}, {60, 25}, {40, 25}, {40, 40}, {0, 40}});
    const PreparedRegion prep = prepared(region);
    auto loops = inward_offset(prep, 8.0);
    REQUIRE(loops);
    CHECK(loops->size() == 2);
    CHECK(check_offset(prep, *loops, 8.0).empty());
}

TEST_CASE("Inward offsets of random shapes", "[M7][ALG-06][VER-05]")
{
    std::mt19937_64 rng(20261011);
    std::uniform_int_distribution<int> size(3, 40);
    std::uniform_real_distribution<double> fraction(0.05, 0.95);
    const int n = shape_count();
    for (int k = 0; k < n; ++k) {
        const Region region = k % 3 == 0   ? random_filleted(rng, size(rng))
                              : k % 3 == 1 ? random_arc_chain(rng, size(rng))
                                           : to_region(random_star_polygon(rng, size(rng)));
        auto prep = prepare(region);
        if (!prep) continue;
        auto mat = compute_mat(*prep);
        REQUIRE(mat);
        double rmax = 0.0;
        for (const MatVertex& v : mat->vertices()) rmax = std::max(rmax, v.r);
        const double d = fraction(rng) * rmax;
        auto loops = inward_offset(*prep, d);
        REQUIRE(loops);
        INFO("shape " << k << " at " << d << " of " << rmax);
        CHECK_FALSE(loops->empty());
        const std::string why = check_offset(*prep, *loops, d);
        CHECK(why.empty());
        if (!why.empty()) break;
    }
}

TEST_CASE("The MAT and its offsets as JSON and SVG", "[M7][OUT-06]")
{
    const PreparedRegion prep = prepared(load("spec-example.json"));
    auto mat = compute_mat(prep);
    REQUIRE(mat);
    auto loops = inward_offset(prep, 3.0);
    REQUIRE(loops);
    const std::vector<io::OffsetLoops> offsets{{3.0, *loops}};
    const auto doc = nlohmann::json::parse(io::write_mat_json(*mat, offsets, {"mm", 0.5}));
    CHECK(doc["units"] == "mm");
    REQUIRE(doc["vertices"].size() == mat->vertices().size());
    REQUIRE(doc["edges"].size() == mat->edges().size());
    CHECK(doc["sites"].size() == mat->sites().size());
    for (std::size_t i = 0; i < mat->edges().size(); ++i) {
        const auto& e = doc["edges"][i];
        CHECK(e["bezier"].size() == mat->edges()[i].bezier().size());
        CHECK(e["polyline"].size() >= 2);
        CHECK(e["r0"].get<double>() == mat->edges()[i].r0());
    }
    REQUIRE(doc["offsets"].size() == 1);
    CHECK(doc["offsets"][0]["loops"].size() == loops->size());
    CHECK(doc["stats"]["clusters"].get<std::size_t>() == mat->stats().clusters);
    CHECK(nlohmann::json::parse(io::write_mat_json(*mat)).contains("offsets") == false);

    io::SvgWriter w(io::caller_bbox(prep));
    w.loops(*loops);
    w.medial_axis(*mat);
    const std::string svg = w.str();
    CHECK(svg.find("class=\"offset\"") != std::string::npos);
    CHECK(svg.find(" A") != std::string::npos);  // the offset's arcs
    CHECK(svg.find("style=\"stroke:#") != std::string::npos);
}
