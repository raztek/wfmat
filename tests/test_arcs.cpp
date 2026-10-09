// [M3] Arcs: convex and concave arc sites, tangent joins, E1c curvature ends, and loops of two
// elements; checked by the [VER-04] properties and the [VER-03] sampling oracle.
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <numbers>
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

} // namespace

TEST_CASE("Disk: one vertex at the centre", "[M3][EV-08][VER-02]")
{
    const Region region = load("disk.json");
    const MedialAxis mat = mat_of(region);
    REQUIRE(mat.vertices().size() == 1);
    CHECK(mat.edges().empty());
    CHECK(mat.vertices()[0].kind == VertexKind::curvature_end);
    CHECK(dist(mat.vertices()[0].p, {0, 0}) < 1e-12);
    CHECK(mat.vertices()[0].r == Approx(10.0).epsilon(1e-14));
    CHECK(check_mat(prepared(region), mat).empty());
}

TEST_CASE("Half disk: two corner edges meeting at a maximum", "[M3][EV-08][VER-02]")
{
    Region region;
    region.outer.vertices = {{{-10, 0}, 0.0}, {{10, 0}, 1.0}};
    const MedialAxis mat = mat_of(region);
    CHECK(mat.vertices().size() == 3);
    CHECK(count(mat, VertexKind::extremum_max) == 1);
    for (const MatVertex& v : mat.vertices())
        if (v.kind == VertexKind::extremum_max) {
            CHECK(dist(v.p, {0, 5}) < 1e-9);
            CHECK(v.r == Approx(5.0).epsilon(1e-12));
        }
    CHECK(check_mat(prepared(region), mat).empty());
}

TEST_CASE("Lens: two convex arcs", "[M3][EV-08][VER-04]")
{
    Region region;
    region.outer.vertices = {{{0, 0}, 0.3}, {{10, 1}, 0.45}};
    const MedialAxis mat = mat_of(region);
    CHECK(mat.vertices().size() == 3);
    CHECK(count(mat, VertexKind::extremum_max) == 1);
    CHECK(check_mat(prepared(region), mat).empty());
}

TEST_CASE("A concave arc opposite a line", "[M3][EV-05][VER-04]")
{
    const Region region = load("concave-top.json");
    const MedialAxis mat = mat_of(region);
    CHECK(mat.stats().contacts == 1);
    CHECK(check_mat(prepared(region), mat).empty());
}

TEST_CASE("Ice-cream cone: the corner edge ends at the arc centre", "[M3][EV-03][VER-02]")
{
    // Two lines from the apex (0, -20), tangent to the circle of radius 5 about the origin.
    const double R = 5.0, h = 20.0;
    const double s = R / h, c = std::sqrt(1.0 - s * s);  // sine and cosine of the half-angle
    const Vec2 t1{R * c, -R * s}, t2{-R * c, -R * s};
    Region region;
    region.outer.vertices = {{{0, -h}, 0.0}, {t1, bulge_from_sweep(std::numbers::pi + 2.0 * std::asin(s))}, {t2, 0.0}};
    const MedialAxis mat = mat_of(region);
    REQUIRE(mat.vertices().size() == 2);
    CHECK(count(mat, VertexKind::curvature_end) == 1);
    for (const MatVertex& v : mat.vertices())
        if (v.kind == VertexKind::curvature_end) {
            CHECK(dist(v.p, {0, 0}) < 1e-9);
            CHECK(v.r == Approx(R).epsilon(1e-12));
        }
    CHECK(check_mat(prepared(region), mat).empty());
}

TEST_CASE("Filleted quadrilateral: corner branches end at the fillet centres", "[M3][EV-03][VER-02]")
{
    const std::vector<Vec2> quad{{0, 0}, {10, 0}, {9, 4}, {1, 5}};
    const Region region = filleted(quad, {0.3, 0.5, 0.2, 0.4});
    const MedialAxis mat = mat_of(region);
    CHECK(count(mat, VertexKind::curvature_end) == 4);
    CHECK(mat.stats().curvature_ends == 4);
    CHECK(check_mat(prepared(region), mat).empty());
}

TEST_CASE("Filleted reflex corner: a concave arc", "[M3][EV-02][VER-04]")
{
    const std::vector<Vec2> hexagon{{0, 0}, {61, 2}, {58, 21}, {22, 19}, {19, 63}, {1, 60}};
    const Region region = filleted(hexagon, {0.0, 0.0, 0.0, 0.6, 0.0, 0.0});
    const MedialAxis mat = mat_of(region);
    CHECK(mat.stats().transitions >= 2);
    CHECK(check_mat(prepared(region), mat).empty());
}

namespace {

int random_shape_count()
{
    if (const char* s = std::getenv("WFMAT_RANDOM_SHAPES")) return std::max(1, std::atoi(s));
    return 10'000;
}

} // namespace

TEST_CASE("Random filleted polygons and arc chains pass the property checks", "[M3][VER-03][VER-04][VER-05]")
{
    const int shapes = random_shape_count();
    std::mt19937_64 rng(20261009);
    std::uniform_int_distribution<int> size(3, 30);
    int failures = 0, refused = 0, skipped = 0;
    std::string first;
    for (int k = 0; k < shapes; ++k) {
        Region region = k % 2 == 0 ? random_filleted(rng, size(rng)) : random_arc_chain(rng, size(rng));
        auto prep = prepare(region);
        while (!prep) {  // a bulge made the loop self-intersect: draw again
            ++skipped;
            region = random_arc_chain(rng, size(rng));
            prep = prepare(region);
        }
        Options options;
        options.debug_checks = k % 16 == 0;
        auto mat = compute_mat(*prep, options);
        std::string why;
        if (!mat && mat.error().code == ErrorCode::unsupported) {
            ++refused;
            why = to_string(mat.error());
        } else {
            why = mat ? check_mat(*prep, *mat, 1e-9, 4) : to_string(mat.error());
        }
        // [VER-03] Every hundredth shape also against the sampled Voronoi oracle.
        if (why.empty() && k % 100 == 0) {
            const SampledComparison c = compare_with_sampled_voronoi(*prep, *mat, 128);
            if (std::max(c.ours_to_voronoi, c.voronoi_to_ours) > 3e-3 * prep->xf.scale)
                why = "sampled Voronoi oracle off by " + sci(std::max(c.ours_to_voronoi, c.voronoi_to_ours) / prep->xf.scale) + " L";
        }
        if (!why.empty() && failures++ == 0) first = "shape " + std::to_string(k) + ": " + why;
        if (!why.empty() && std::getenv("WFMAT_DUMP")) {
            std::fprintf(stderr, "%d: %s\n", k, why.c_str());
            std::ofstream(output_path("failure-" + std::to_string(k) + ".json")) << io::write_region_json({region, ""});
        }
    }
    INFO(first);
    INFO(std::to_string(skipped) + " self-intersecting, " + std::to_string(refused) + " refused");
    CHECK(failures == 0);
}

TEST_CASE("The sampled Voronoi oracle converges to the axis", "[M3][VER-03]")
{
    std::vector<std::pair<std::string, Region>> shapes;
    shapes.emplace_back("concave top", load("concave-top.json"));
    shapes.emplace_back("filleted quadrilateral", filleted({{0, 0}, {10, 0}, {9, 4}, {1, 5}}, {0.3, 0.5, 0.2, 0.4}));
    shapes.emplace_back("filleted reflex corner",
                        filleted({{0, 0}, {61, 2}, {58, 21}, {22, 19}, {19, 63}, {1, 60}}, {0.5, 0.0, 0.3, 0.6, 0.0, 0.4}));
    Region lens;
    lens.outer.vertices = {{{0, 0}, 0.3}, {{10, 1}, 0.45}};
    shapes.emplace_back("lens", lens);
    std::mt19937_64 rng(3);
    for (int k = 0; k < 6; ++k) {
        std::uniform_int_distribution<int> size(4, 12);
        Region r = k % 2 == 0 ? random_filleted(rng, size(rng)) : random_arc_chain(rng, size(rng));
        if (prepare(r)) shapes.emplace_back("random " + std::to_string(k), r);
    }

    for (const auto& [name, region] : shapes) {
        INFO(name);
        const PreparedRegion prep = prepared(region);
        const MedialAxis mat = mat_of(region);
        const double L = prep.xf.scale;
        if (std::getenv("WFMAT_DUMP")) std::ofstream(output_path("oracle-" + name + ".json")) << io::write_region_json({region, ""});
        // The polyline is within (pi / chords)^2 L / 8 of the boundary, so the error shrinks about
        // 16-fold for 4 times the chords, down to the 1e-8 L grid of the integer coordinates. Near
        // the centre of a convex arc the Voronoi diagram of its chords converges at first order
        // only, with an error of about R pi / (4 chords): ask for a 3-fold reduction, and 3e-3 L at
        // 128 chords.
        const SampledComparison coarse = compare_with_sampled_voronoi(prep, mat, 32);
        const SampledComparison fine = compare_with_sampled_voronoi(prep, mat, 128);
        INFO("coarse " << coarse.ours_to_voronoi / L << ", " << coarse.voronoi_to_ours / L << "; fine "
                       << fine.ours_to_voronoi / L << ", " << fine.voronoi_to_ours / L);
        CHECK(fine.ours_to_voronoi < std::max(coarse.ours_to_voronoi / 3, 1e-7 * L));
        CHECK(fine.voronoi_to_ours < std::max(coarse.voronoi_to_ours / 3, 1e-7 * L));
        CHECK(fine.ours_to_voronoi < 3e-3 * L);
        CHECK(fine.voronoi_to_ours < 3e-3 * L);
    }
}
