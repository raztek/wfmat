// M4.2: the axis of regions with holes and islands [ALG-07], [EV-09], [RB-07], [OUT-08], [OUT-11].
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <numbers>
#include <random>
#include <string>

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include "mat_checks.hpp"
#include "support.hpp"
#include "wfmat/mat.hpp"

using namespace wfmat;
using namespace wfmat::test;
using Catch::Approx;

namespace {

constexpr double pi = std::numbers::pi;

Loop square(Vec2 c, double h, bool ccw)
{
    Loop loop;
    const Vec2 corners[] = {{-h, -h}, {h, -h}, {h, h}, {-h, h}};
    for (int k = 0; k < 4; ++k) loop.vertices.push_back({c + corners[ccw ? k : 3 - k], 0.0});
    return loop;
}

Loop circle(Vec2 c, double r, bool ccw)
{
    const double b = ccw ? 1.0 : -1.0;
    return Loop{{{c + Vec2{r, 0}, b}, {c - Vec2{r, 0}, b}}};
}

// A rectangle [-a, a] x [-b, b] with its corners rounded to radius f (f = 0 for sharp corners).
Loop rounded_rectangle(double a, double b, double f, bool ccw)
{
    const double q = std::tan(pi / 8.0);  // a quarter turn
    std::vector<BulgeVertex> v;
    if (f == 0.0) {
        v = {{{-a, -b}, 0}, {{a, -b}, 0}, {{a, b}, 0}, {{-a, b}, 0}};
    } else {
        v = {{{-a + f, -b}, 0}, {{a - f, -b}, q}, {{a, -b + f}, 0}, {{a, b - f}, q},
             {{a - f, b}, 0},   {{-a + f, b}, q}, {{-a, b - f}, 0}, {{-a, -b + f}, q}};
    }
    Loop loop{v};
    if (!ccw) {
        // Reverse: vertex k + 1 keeps the bulge of the segment that ended at it, negated.
        Loop r;
        const std::size_t n = v.size();
        for (std::size_t k = 0; k < n; ++k) {
            const std::size_t i = (n - k) % n, before = (i + n - 1) % n;
            r.vertices.push_back({v[i].p, -v[before].bulge});
        }
        loop = r;
    }
    return loop;
}

// Every MAT point sampled along the edges, with its radius.
std::vector<std::pair<Vec2, double>> samples(const MedialAxis& mat, int per_edge = 16)
{
    std::vector<std::pair<Vec2, double>> out;
    for (const MatEdge& e : mat.edges())
        for (int k = 0; k <= per_edge; ++k) {
            const double u = double(k) / per_edge;
            out.emplace_back(e.point_at(u), e.radius_at(u));
        }
    return out;
}

// The checks every region with holes passes: check_mat, the sampled Voronoi oracle, the same axis
// from the all-pairs broad phase, and a valid axis when every event goes through cluster resolution.
void check_region(const Region& region, int chords = 128)
{
    auto prep = prepare(region);
    REQUIRE(prep);
    auto mat = compute_mat(*prep);
    REQUIRE(mat);
    CHECK(check_mat(*prep, *mat) == "");
    const double L = prep->xf.scale;
    const SampledComparison c = compare_with_sampled_voronoi(*prep, *mat, chords);
    INFO("oracle " << c.ours_to_voronoi / L << ", " << c.voronoi_to_ours / L);
    CHECK(c.ours_to_voronoi < 3e-3 * L);
    CHECK(c.voronoi_to_ours < 3e-3 * L);

    Options all_pairs;
    all_pairs.broad_phase = BroadPhase::all_pairs;
    auto ap = compute_mat(*prep, all_pairs);
    REQUIRE(ap);
    CHECK(identical_mat(*mat, *ap) == "");

    Options clusters;
    clusters.resolve_all_events = true;
    auto cl = compute_mat(*prep, clusters);
    REQUIRE(cl);
    CHECK(check_mat(*prep, *cl) == "");
}

} // namespace

TEST_CASE("concentric annulus: a closed plateau circle halfway between the rims", "[RB-07][OUT-08]")
{
    const double Ro = 10.0, Ri = 4.0;
    const Region region{{circle({0, 0}, Ro, true), circle({0, 0}, Ri, false)}};
    auto prep = prepare(region);
    REQUIRE(prep);
    auto mat = compute_mat(*prep);
    REQUIRE(mat);
    CHECK(check_mat(*prep, *mat) == "");
    REQUIRE(mat->components().size() == 1);
    for (const MatEdge& e : mat->edges()) CHECK(e.kind() == ConicKind::plateau_arc);
    for (const auto& [p, r] : samples(*mat)) {
        CHECK(norm(p) == Approx(0.5 * (Ro + Ri)).epsilon(1e-9));
        CHECK(r == Approx(0.5 * (Ro - Ri)).epsilon(1e-9));
    }
    CHECK(mat->stats().merges == 1);
}

TEST_CASE("eccentric annulus: an ellipse with foci at the two centres", "[EV-09][OUT-08]")
{
    const double Ro = 10.0, Ri = 4.0, e = 2.0;
    const Vec2 c{e, 0.0};
    const Region region{{circle({0, 0}, Ro, true), circle(c, Ri, false)}};
    auto prep = prepare(region);
    REQUIRE(prep);
    auto mat = compute_mat(*prep);
    REQUIRE(mat);
    CHECK(check_mat(*prep, *mat) == "");
    double rmin = INFINITY, rmax = 0.0;
    for (const auto& [p, r] : samples(*mat)) {
        CHECK(norm(p) + dist(p, c) == Approx(Ro + Ri).epsilon(1e-9));
        CHECK(r == Approx(Ro - norm(p)).epsilon(1e-9));
        rmin = std::min(rmin, r);
        rmax = std::max(rmax, r);
    }
    CHECK(rmin == Approx(0.5 * (Ro - Ri - e)).epsilon(1e-9));
    CHECK(rmax == Approx(0.5 * (Ro - Ri + e)).epsilon(1e-9));
    CHECK(mat->stats().contacts + mat->stats().clusters >= 1);
    CHECK(mat->stats().merges == 1);
}

TEST_CASE("square ring: one cycle around the hole", "[EV-09][OUT-08]")
{
    check_region(load("square-ring.json"));
    auto mat = compute_mat(load("square-ring.json"));
    REQUIRE(mat);
    CHECK(mat->components()[0].holes.size() == 1);
}

TEST_CASE("off-centre square holes", "[EV-09][OUT-08]")
{
    check_region(Region{{square({0, 0}, 100, true), square({20, 10}, 30, false)}});
    check_region(Region{{square({0, 0}, 100, true), square({-35, 40}, 20, false), square({40, -30}, 25, false)}});
}

TEST_CASE("disk with k circular holes: k independent cycles", "[EV-09][OUT-08]")
{
    for (int k = 1; k <= 6; ++k) {
        INFO(k << " holes");
        Region region{{circle({0, 0}, 100.0, true)}};
        for (int i = 0; i < k; ++i) {
            const double a = 2.0 * pi * i / k + 0.3, d = 55.0 + 7.0 * (i % 2);
            region.loops.push_back(circle({d * std::cos(a), d * std::sin(a)}, 12.0 + 2.0 * (i % 3), false));
        }
        check_region(region);
        auto mat = compute_mat(region);
        REQUIRE(mat);
        CHECK(mat->components()[0].holes.size() == static_cast<std::size_t>(k));
        CHECK(mat->stats().merges == static_cast<std::size_t>(k));
    }
}

TEST_CASE("nested islands: one axis per component", "[ALG-07][OUT-11]")
{
    // Squares inside circular holes, three levels of islands, each slightly off centre.
    const double size[] = {100.0, 90.0, 60.0, 55.0, 36.0, 32.0};
    Region region;
    for (int k = 0; k < 6; ++k) {
        const Vec2 c{1.0 * k, 0.5 * k};
        region.loops.push_back(k % 2 == 0 ? square(c, size[k], true) : circle(c, size[k], false));
    }
    check_region(region);
    auto mat = compute_mat(region);
    REQUIRE(mat);
    REQUIRE(mat->components().size() == 3);
    for (const MatComponent& comp : mat->components()) CHECK(comp.holes.size() == 1);
}

TEST_CASE("plate with two holes and an island", "[ALG-07][OUT-08][OUT-11]")
{
    check_region(load("plate-with-holes.json"));
}

TEST_CASE("constant-width rings: a closed plateau chain of lines and arcs", "[RB-07][OUT-08]")
{
    // The outer loop is the offset of the inner one by w, sharp or filleted, so every point of the
    // mid offset is equidistant from both: the axis is that offset, all plateau.
    const double a = 40.0, b = 25.0, w = 10.0;
    for (const double f : {0.0, 5.0}) {
        INFO("inner fillet " << f);
        const Region region{{rounded_rectangle(a + w, b + w, f + w, true), rounded_rectangle(a, b, f, false)}};
        auto prep = prepare(region);
        REQUIRE(prep);
        auto mat = compute_mat(*prep);
        REQUIRE(mat);
        CHECK(check_mat(*prep, *mat) == "");
        double rmin = INFINITY, rmax = 0.0;
        for (const MatEdge& e : mat->edges()) {
            rmin = std::min({rmin, e.r0(), e.r1()});
            rmax = std::max({rmax, e.r0(), e.r1()});
        }
        for (const MatEdge& e : mat->edges()) CHECK(e.is_plateau());
        CHECK(rmin == Approx(0.5 * w));
        CHECK(rmax == Approx(0.5 * w));
    }
}

TEST_CASE("the inward offset of a region with holes", "[ALG-06]")
{
    const Region ring = load("square-ring.json");
    auto two = inward_offset(ring, 5.0);
    REQUIRE(two);
    REQUIRE(two->size() == 1);
    CHECK((*two)[0].size() == 2);
    auto none = inward_offset(ring, 20.0);
    REQUIRE(none);
    REQUIRE(none->size() == 1);
    CHECK((*none)[0].empty());
}

TEST_CASE("random regions with holes pass the property checks and the oracle", "[VER-05][VER-03][OUT-08]")
{
    // Polygons with polygonal holes, then filleted polygons with circular and filleted holes. Each
    // passes check_mat; every tenth also matches the sampled Voronoi oracle and the all-pairs run.
    std::mt19937_64 rng(20261012);
    const int shapes = 400;
    int failures = 0;
    std::size_t holes = 0;
    std::string first;
    for (int k = 0; k < shapes; ++k) {
        const Region region = random_holed(rng, k >= shapes / 2);
        holes += region.loops.size() - 1;
        auto prep = prepare(region);
        REQUIRE(prep);
        Options options;
        options.debug_checks = k % 16 == 0;
        auto mat = compute_mat(*prep, options);
        std::string why = mat ? check_mat(*prep, *mat, 1e-9, 4) : to_string(mat.error());
        if (why.empty() && k % 10 == 0) {
            const SampledComparison c = compare_with_sampled_voronoi(*prep, *mat, 128);
            if (std::max(c.ours_to_voronoi, c.voronoi_to_ours) > 3e-3 * prep->xf.scale)
                why = "sampled Voronoi oracle off by " +
                      sci(std::max(c.ours_to_voronoi, c.voronoi_to_ours) / prep->xf.scale) + " L";
            Options all_pairs;
            all_pairs.broad_phase = BroadPhase::all_pairs;
            auto ap = compute_mat(*prep, all_pairs);
            if (why.empty()) why = ap ? identical_mat(*mat, *ap) : "all pairs: " + to_string(ap.error());
        }
        if (!why.empty() && failures++ == 0) first = "shape " + std::to_string(k) + ": " + why;
        if (!why.empty() && std::getenv("WFMAT_DUMP")) {
            std::fprintf(stderr, "%d: %s\n", k, why.c_str());
            std::ofstream(output_path("holes-failure-" + std::to_string(k) + ".json"))
                << io::write_region_json({region, ""});
        }
    }
    INFO(first);
    CHECK(failures == 0);
    CHECK(holes > shapes * 3 / 2);  // most shapes have several holes
}
