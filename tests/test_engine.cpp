// [M2] The polygon engine: line and reflex-point sites, events E1a, E1b, E2, E3, all-pairs
// scheduling, checked by the [VER-04] properties and against Boost.Polygon's Voronoi diagram.
#include <cmath>
#include <cstdlib>
#include <random>
#include <string>
#include <vector>

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include "mat_checks.hpp"
#include "support.hpp"
#include "wfmat/mat.hpp"

using namespace wfmat;
using namespace wfmat::test;
using Catch::Approx;

namespace {

MedialAxis mat_of(const Region& region, const Options& options = {})
{
    auto mat = compute_mat(region, options);
    if (!mat) FAIL(to_string(mat.error()));
    return *mat;
}

std::size_t count(const MedialAxis& mat, VertexKind kind)
{
    std::size_t n = 0;
    for (const MatVertex& v : mat.vertices()) n += v.kind == kind;
    return n;
}

int random_polygon_count()
{
    if (const char* s = std::getenv("WFMAT_RANDOM_POLYGONS")) return std::max(1, std::atoi(s));
    return 10'000;
}

} // namespace

TEST_CASE("Scalene triangle: three corner edges meeting at the incentre", "[M2][EV-01][EV-08][VER-02]")
{
    const Region region = polygon({{0, 0}, {7, 0}, {2, 5}});
    const MedialAxis mat = mat_of(region);
    REQUIRE(mat.vertices().size() == 4);
    REQUIRE(mat.edges().size() == 3);
    CHECK(count(mat, VertexKind::corner) == 3);
    CHECK(count(mat, VertexKind::junction) == 1);
    CHECK(mat.stats().annihilations == 1);

    // Incentre: the side-length-weighted mean of the corners; inradius = area / semiperimeter.
    const Vec2 A{0, 0}, B{7, 0}, C{2, 5};
    const double a = dist(B, C), b = dist(C, A), c = dist(A, B);
    const Vec2 incentre = (a * A + b * B + c * C) / (a + b + c);
    const double inradius = 17.5 / (0.5 * (a + b + c));
    for (const MatVertex& v : mat.vertices()) {
        if (v.kind != VertexKind::junction) continue;
        CHECK(dist(v.p, incentre) < 1e-12);
        CHECK(v.r == Approx(inradius).epsilon(1e-13));
        CHECK(v.contacts.size() == 3);
    }
    CHECK(check_mat(prepared(region), mat).empty());
}

TEST_CASE("Convex quadrilateral: two junctions joined by one edge", "[M2][EV-01][VER-04]")
{
    const Region region = polygon({{0, 0}, {10, 0}, {9, 4}, {1, 5}});
    const MedialAxis mat = mat_of(region);
    CHECK(mat.vertices().size() == 6);
    CHECK(count(mat, VertexKind::junction) == 2);
    CHECK(mat.stats().collapses == 1);
    CHECK(check_mat(prepared(region), mat).empty());
}

TEST_CASE("A reflex corner gives parabolic edges and transitions", "[M2][EV-02][VER-02]")
{
    // An L-like hexagon without parallel sides.
    const Region region = polygon({{0, 0}, {61, 2}, {58, 21}, {22, 19}, {19, 63}, {1, 60}});
    const MedialAxis mat = mat_of(region);
    std::size_t parabolas = 0;
    for (const MatEdge& e : mat.edges()) parabolas += e.kind() == ConicKind::parabola;
    CHECK(parabolas >= 1);
    CHECK(mat.stats().transitions >= 1);
    CHECK(check_mat(prepared(region), mat).empty());

    // Every parabolic edge has the reflex corner as one site: its foot there is the corner.
    for (const MatEdge& e : mat.edges()) {
        if (e.kind() != ConicKind::parabola) continue;
        const auto [fl, fr] = e.feet_at(0.5);
        const bool at_corner = dist(fl, {22, 19}) < 1e-9 || dist(fr, {22, 19}) < 1e-9;
        CHECK(at_corner);
    }
}

TEST_CASE("A waist: two reflex corners facing each other give an E3 minimum", "[M2][EV-05][KN-05]")
{
    // A hexagon pinched by two reflex corners.
    const Region region = polygon({{0, 0}, {20, 6}, {41, 1}, {43, 30}, {21, 23}, {-2, 29}});
    const MedialAxis mat = mat_of(region);
    CHECK(mat.stats().contacts == 1);
    CHECK(count(mat, VertexKind::extremum_min) == 1);
    for (const MatVertex& v : mat.vertices()) {
        if (v.kind != VertexKind::extremum_min) continue;
        // The minimum lies halfway between the two reflex corners.
        CHECK(dist(v.p, Vec2{20.5, 14.5}) < 1e-9);
        CHECK(v.r == Approx(0.5 * dist({20, 6}, {21, 23})).epsilon(1e-12));
    }
    CHECK(check_mat(prepared(region), mat).empty());
}

TEST_CASE("Split: a shock reaching a far element", "[M2][EV-04]")
{
    // A wide pentagon whose apex shock runs into the base before the side elements collapse.
    const Region region = polygon({{0, 0}, {30, 1}, {29, 9}, {15, 4}, {1, 10}});
    const MedialAxis mat = mat_of(region);
    CHECK(mat.stats().splits + mat.stats().contacts >= 1);
    CHECK(check_mat(prepared(region), mat).empty());
}

TEST_CASE("Output is deterministic and independent of orientation and placement", "[M2][OUT-10][VER-04]")
{
    std::mt19937_64 rng(7);
    for (int k = 0; k < 20; ++k) {
        const auto poly = random_star_polygon(rng, 5 + k % 13);
        const Region region = to_region(poly);
        const MedialAxis m1 = mat_of(region), m2 = mat_of(region);
        REQUIRE(m1.vertices().size() == m2.vertices().size());
        for (std::size_t i = 0; i < m1.vertices().size(); ++i) {
            CHECK(m1.vertices()[i].p == m2.vertices()[i].p);
            CHECK(m1.vertices()[i].r == m2.vertices()[i].r);
        }

        // Reversed, rotated by 90 degrees and shifted: the same vertex set up to the transform.
        Region moved;
        for (auto it = poly.rbegin(); it != poly.rend(); ++it)
            moved.outer.vertices.push_back({{-double(it->second) + 12345.0, double(it->first) - 678.0}, 0.0});
        const MedialAxis m3 = mat_of(moved);
        REQUIRE(m3.vertices().size() == m1.vertices().size());
        for (const MatVertex& v : m1.vertices()) {
            const Vec2 q{-v.p.y + 12345.0, v.p.x - 678.0};
            double best = INFINITY;
            for (const MatVertex& w : m3.vertices())
                if (w.kind == v.kind) best = std::min(best, dist(w.p, q));
            CHECK(best < 1e-9 * 1e6);
        }
    }
}

TEST_CASE("Arcs and clusters are refused until their milestones", "[M2][API-03]")
{
    SECTION("an arc") {
        const auto mat = compute_mat(load("slot.json"));
        REQUIRE_FALSE(mat.has_value());
        CHECK(mat.error().code == ErrorCode::unsupported);
        CHECK(mat.error().requirement == "M3");
    }
    SECTION("a square: four collapses at the centre") {
        const auto mat = compute_mat(polygon({{0, 0}, {1, 0}, {1, 1}, {0, 1}}));
        REQUIRE_FALSE(mat.has_value());
        CHECK(mat.error().code == ErrorCode::unsupported);
    }
    SECTION("a rectangle: a plateau") {
        const auto mat = compute_mat(load("rectangle.json"));
        REQUIRE_FALSE(mat.has_value());
        CHECK(mat.error().code == ErrorCode::unsupported);
    }
}

TEST_CASE("Random polygons match the Boost.Polygon Voronoi diagram", "[M2][VER-03][VER-05][EV-12]")
{
    const int shapes = random_polygon_count();
    std::mt19937_64 rng(20261009);
    std::uniform_int_distribution<int> size(3, 40);
    int failures = 0;
    std::string first;
    for (int k = 0; k < shapes; ++k) {
        const auto poly = random_star_polygon(rng, size(rng));
        const Region region = to_region(poly);
        auto prep = prepare(region);
        if (!prep) continue;  // rounding made the polygon non-simple
        Options options;
        options.debug_checks = k % 16 == 0;
        auto mat = compute_mat(*prep, options);
        std::string why = mat ? check_mat(*prep, *mat, 1e-9, 4) : to_string(mat.error());
        if (why.empty()) why = compare_with_voronoi(poly, *prep, *mat, 1e-9 * prep->xf.scale);
        if (!why.empty()) {
            if (failures++ == 0) first = "shape " + std::to_string(k) + " (" + std::to_string(poly.size()) + " vertices): " + why;
        }
    }
    INFO(first);
    CHECK(failures == 0);
}
