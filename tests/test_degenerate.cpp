// [M5] Degeneracies: cluster formation and resolution and plateau edges [RB-02..RB-04], checked on
// the [VER-06] degeneracy suite, against Boost.Polygon's Voronoi diagram where the input is an
// integer polygon, and by sending every event of the random suites through cluster resolution.
#include <algorithm>
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

constexpr double pi = std::numbers::pi;

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

std::vector<const MatEdge*> plateaus(const MedialAxis& mat)
{
    std::vector<const MatEdge*> out;
    for (const MatEdge& e : mat.edges())
        if (e.is_plateau()) out.push_back(&e);
    return out;
}

int degenerate_shape_count()
{
    if (const char* s = std::getenv("WFMAT_DEGENERATE_SHAPES")) return std::max(1, std::atoi(s));
    return 2'000;
}

// Drops repeated points, so that equal neighbouring heights make no zero-length segment, and the
// middle points of collinear runs, which Boost.Polygon would take as sites of their own.
std::vector<IPoint> clean(const std::vector<IPoint>& p)
{
    std::vector<IPoint> out;
    for (const IPoint& q : p)
        if (out.empty() || out.back() != q) out.push_back(q);
    while (out.size() > 1 && out.front() == out.back()) out.pop_back();
    for (bool removed = true; removed && out.size() > 3;) {
        removed = false;
        for (std::size_t i = 0; i < out.size(); ++i) {
            const IPoint a = out[(i + out.size() - 1) % out.size()], b = out[i], c = out[(i + 1) % out.size()];
            const std::int64_t turn = std::int64_t(b.first - a.first) * (c.second - b.second) -
                                      std::int64_t(b.second - a.second) * (c.first - b.first);
            if (turn == 0) {
                out.erase(out.begin() + static_cast<std::ptrdiff_t>(i));
                removed = true;
                break;
            }
        }
    }
    return out;
}

// [VER-06] An x-monotone orthogonal polygon: columns of random widths between a random floor and
// ceiling, every pair of neighbours overlapping. Equal widths and heights everywhere make plateaus
// and simultaneous events.
std::vector<IPoint> random_columns(std::mt19937_64& rng, int n, bool flat_floor, std::int32_t unit = 1000)
{
    std::uniform_int_distribution<int> width(1, 3), floor(0, 2), height(1, 4);
    std::vector<int> x{0}, lo(n), hi(n);
    for (int i = 0; i < n; ++i) {
        x.push_back(x.back() + width(rng));
        lo[i] = flat_floor ? 0 : floor(rng);
        hi[i] = lo[i] + height(rng);
        // Neighbouring columns share at least one unit of height.
        if (i > 0) {
            lo[i] = std::min(lo[i], hi[i - 1] - 1);
            hi[i] = std::max(hi[i], lo[i - 1] + 1);
        }
    }
    std::vector<IPoint> p;
    for (int i = 0; i < n; ++i) {
        p.push_back({x[i] * unit, lo[i] * unit});
        p.push_back({x[i + 1] * unit, lo[i] * unit});
    }
    for (int i = n - 1; i >= 0; --i) {
        p.push_back({x[i + 1] * unit, hi[i] * unit});
        p.push_back({x[i] * unit, hi[i] * unit});
    }
    return clean(p);
}

// [VER-06] A star polygon with exact symmetry about both axes: a random chain in the first quadrant
// from the x axis to the y axis, mirrored.
std::vector<IPoint> random_symmetric(std::mt19937_64& rng, int k, std::int32_t R = 1'000'000)
{
    std::uniform_real_distribution<double> radius(0.3, 1.0);
    std::vector<double> a(k);
    for (int i = 0; i < k; ++i) a[i] = (i + 1) * (pi / 2) / (k + 1);
    auto at = [&](double angle, double r) {
        return IPoint{static_cast<std::int32_t>(std::lround(R * r * std::cos(angle))),
                      static_cast<std::int32_t>(std::lround(R * r * std::sin(angle)))};
    };
    std::vector<IPoint> quadrant{at(0.0, radius(rng))};
    for (double angle : a) quadrant.push_back(at(angle, radius(rng)));
    quadrant.push_back({0, static_cast<std::int32_t>(std::lround(R * radius(rng)))});
    // The upper half, mirrored about the y axis, then the lower half, mirrored about the x axis.
    std::vector<IPoint> p = quadrant;
    for (auto it = quadrant.rbegin() + 1; it != quadrant.rend(); ++it) p.push_back({-it->first, it->second});
    for (std::size_t i = p.size() - 2; i >= 1; --i) p.push_back({p[i].first, -p[i].second});
    return clean(p);
}

// A regular n-gon of circumradius 1.
Region regular_polygon(int n, double phase = 0.0)
{
    Region region;
    for (int k = 0; k < n; ++k) region.outer.vertices.push_back({polar(phase + 2 * pi * k / n), 0.0});
    return region;
}

// [VER-06] A polygon whose sides all touch the unit circle (a tangential polygon): every side
// collapses at the centre at t = 1.
Region random_tangential(std::mt19937_64& rng, int n)
{
    // Gaps between the tangent directions: at least 0.05, at most 0.8 pi, summing to 2 pi.
    std::uniform_real_distribution<double> weight(0.0, 1.0), phase(0.0, 2 * pi);
    for (;;) {
        std::vector<double> gap(n);
        double sum = 0.0;
        for (double& g : gap) sum += g = weight(rng);
        const double spare = 2 * pi - 0.05 * n;
        for (double& g : gap) g = 0.05 + spare * g / sum;
        if (*std::max_element(gap.begin(), gap.end()) > 0.8 * pi) continue;
        Region region;
        double a = phase(rng);
        for (double g : gap) {
            region.outer.vertices.push_back({polar(a + g / 2) / std::cos(g / 2), 0.0});
            a += g;
        }
        return region;
    }
}

// [VER-06] A gear: n trapezoidal teeth between roots on the unit circle, so that the reflex roots
// are cocircular and the outline has n-fold symmetry.
Region gear(int n, double tip = 1.3, double land = 0.35, double top = 0.2)
{
    Region region;
    const double step = 2 * pi / n;
    for (int k = 0; k < n; ++k) {
        const double c = k * step;
        region.outer.vertices.push_back({polar(c - land * step), 0.0});
        region.outer.vertices.push_back({tip * polar(c - top * step), 0.0});
        region.outer.vertices.push_back({tip * polar(c + top * step), 0.0});
        region.outer.vertices.push_back({polar(c + land * step), 0.0});
    }
    return region;
}

// [RB-04] Two MATs agree: the same vertices (kind, position, radius) and the same edges.
std::string same_mat(const MedialAxis& a, const MedialAxis& b, double tol)
{
    if (a.vertices().size() != b.vertices().size() || a.edges().size() != b.edges().size())
        return std::to_string(a.vertices().size()) + " vertices and " + std::to_string(a.edges().size()) +
               " edges vs " + std::to_string(b.vertices().size()) + " and " + std::to_string(b.edges().size());
    std::vector<std::size_t> match(a.vertices().size());
    std::vector<char> taken(b.vertices().size(), 0);
    for (std::size_t i = 0; i < a.vertices().size(); ++i) {
        const MatVertex& v = a.vertices()[i];
        std::size_t best = SIZE_MAX;
        for (std::size_t j = 0; j < b.vertices().size(); ++j) {
            const MatVertex& w = b.vertices()[j];
            if (taken[j] || w.kind != v.kind || dist(v.p, w.p) > tol || std::abs(v.r - w.r) > tol) continue;
            if (best == SIZE_MAX || dist(v.p, w.p) < dist(v.p, b.vertices()[best].p)) best = j;
        }
        if (best == SIZE_MAX)
            return "vertex " + std::to_string(i) + " at (" + std::to_string(v.p.x) + ", " + std::to_string(v.p.y) +
                   ") has no counterpart";
        taken[best] = 1;
        match[i] = best;
    }
    std::vector<std::pair<std::size_t, std::size_t>> ea, eb;
    for (const MatEdge& e : a.edges()) ea.emplace_back(std::minmax(match[e.v0()], match[e.v1()]));
    for (const MatEdge& e : b.edges()) eb.emplace_back(std::minmax(std::size_t{e.v0()}, std::size_t{e.v1()}));
    std::sort(ea.begin(), ea.end());
    std::sort(eb.begin(), eb.end());
    return ea == eb ? std::string{} : std::string("the edges differ");
}

// Runs a family of shapes, with the property checks, the Boost.Polygon comparison for integer
// polygons, the [RB-04] cross-check against cluster resolution of every event, and the [EV-12]
// check that the all-pairs broad phase processes the same events.
struct FamilyRun {
    int failures = 0;
    std::string first;
    void fail(const std::string& name, const std::string& why, const Region& region)
    {
        if (failures++ == 0) first = name + ": " + why;
        if (std::getenv("WFMAT_DUMP")) {
            std::fprintf(stderr, "%s: %s\n", name.c_str(), why.c_str());
            std::ofstream(output_path("degenerate-" + name + ".json")) << io::write_region_json({region, ""});
        }
    }
    void run(const std::string& name, const Region& region, const std::vector<IPoint>* integer_polygon = nullptr,
             bool debug = true)
    {
        auto prep = prepare(region);
        if (!prep) return fail(name, to_string(prep.error()), region);
        Options options;
        options.debug_checks = debug;
        auto mat = compute_mat(*prep, options);
        if (!mat) return fail(name, to_string(mat.error()), region);
        std::string why = check_mat(*prep, *mat, 1e-9, 4);
        if (why.empty() && integer_polygon)
            why = compare_with_voronoi(*integer_polygon, *prep, *mat, 1e-9 * prep->xf.scale);
        if (why.empty()) {
            options.resolve_all_events = true;
            auto all = compute_mat(*prep, options);
            why = all ? same_mat(*mat, *all, 1e-9 * prep->xf.scale) : to_string(all.error());
            if (!why.empty()) why = "cluster resolution of every event: " + why;
        }
        if (why.empty()) {
            options.resolve_all_events = false;
            options.broad_phase = BroadPhase::all_pairs;
            auto all = compute_mat(*prep, options);
            why = all ? identical_mat(*mat, *all) : to_string(all.error());
            if (!why.empty()) why = "all-pairs broad phase: " + why;
        }
        if (!why.empty()) fail(name, why, region);
    }
};

} // namespace

TEST_CASE("Rectangle: four corner edges and a plateau", "[M5][RB-03][OUT-02][VER-02]")
{
    const Region region = load("rectangle.json");  // 100 x 40
    const MedialAxis mat = mat_of(region);
    REQUIRE(mat.vertices().size() == 6);
    REQUIRE(mat.edges().size() == 5);
    CHECK(count(mat, VertexKind::corner) == 4);
    CHECK(count(mat, VertexKind::junction) == 2);
    const auto flat = plateaus(mat);
    REQUIRE(flat.size() == 1);
    CHECK(flat[0]->kind() == ConicKind::plateau_line);
    CHECK(flat[0]->r0() == Approx(20.0).epsilon(1e-12));
    CHECK(flat[0]->r1() == Approx(20.0).epsilon(1e-12));
    const Vec2 a = flat[0]->point_at(0.0), b = flat[0]->point_at(1.0);
    CHECK(dist(a, b) == Approx(60.0).epsilon(1e-12));
    CHECK(std::abs(a.y - 20.0) < 1e-12);
    CHECK(std::min(a.x, b.x) == Approx(20.0).epsilon(1e-12));
    CHECK(mat.stats().clusters == 1);
    CHECK(check_mat(prepared(region), mat).empty());
    const std::vector<IPoint> poly{{0, 0}, {100, 0}, {100, 40}, {0, 40}};
    CHECK(compare_with_voronoi(poly, prepared(region), mat, 1e-9 * 100).empty());
}

TEST_CASE("Regular polygons: one vertex of degree n at the centre", "[M5][RB-02][RB-03][VER-02][VER-06]")
{
    for (int n : {3, 4, 5, 6, 8, 12, 31, 64, 256, 1024}) {
        INFO("n = " << n);
        const Region region = regular_polygon(n, 0.1 * n);
        const MedialAxis mat = mat_of(region);
        REQUIRE(mat.vertices().size() == static_cast<std::size_t>(n + 1));
        CHECK(count(mat, VertexKind::corner) == static_cast<std::size_t>(n));
        for (const MatVertex& v : mat.vertices()) {
            if (v.kind != VertexKind::junction) continue;
            CHECK(v.edges.size() == static_cast<std::size_t>(n));
            CHECK(norm(v.p) < 1e-9);
            CHECK(v.r == Approx(std::cos(pi / n)).epsilon(1e-9));
        }
        CHECK(check_mat(prepared(region), mat, 1e-9, 2).empty());
    }
}

TEST_CASE("Slot: a plateau between two curvature ends", "[M5][RB-03][OUT-02][VER-02]")
{
    const Region region = load("slot.json");  // 60 x 20, end radius 10
    const MedialAxis mat = mat_of(region);
    REQUIRE(mat.vertices().size() == 2);
    REQUIRE(mat.edges().size() == 1);
    CHECK(count(mat, VertexKind::curvature_end) == 2);
    const MatEdge& e = mat.edges()[0];
    CHECK(e.kind() == ConicKind::plateau_line);
    CHECK(e.r0() == Approx(10.0).epsilon(1e-12));
    for (const MatVertex& v : mat.vertices()) {
        CHECK(std::abs(v.p.y - 10.0) < 1e-12);
        CHECK((std::abs(v.p.x) < 1e-12 || std::abs(v.p.x - 60.0) < 1e-12));
    }
    CHECK(check_mat(prepared(region), mat).empty());
}

TEST_CASE("Curved slot: a plateau arc between two curvature ends", "[M5][RB-03][OUT-02][VER-02]")
{
    // Centre line of radius 8 from 0 to 90 degrees, half-width 2, round ends.
    const double b = std::tan(pi / 8);
    Region region;
    region.outer.vertices = {{{10, 0}, b}, {{0, 10}, 1.0}, {{0, 6}, -b}, {{6, 0}, 1.0}};
    const MedialAxis mat = mat_of(region);
    REQUIRE(mat.vertices().size() == 2);
    REQUIRE(mat.edges().size() == 1);
    CHECK(count(mat, VertexKind::curvature_end) == 2);
    const MatEdge& e = mat.edges()[0];
    CHECK(e.kind() == ConicKind::plateau_arc);
    CHECK(e.r0() == Approx(2.0).epsilon(1e-12));
    for (double u : {0.0, 0.25, 0.5, 0.75, 1.0}) CHECK(norm(e.point_at(u)) == Approx(8.0).epsilon(1e-12));
    CHECK(norm(e.point_at(0.5) - Vec2{8 / std::sqrt(2.0), 8 / std::sqrt(2.0)}) < 1e-9);
    CHECK(check_mat(prepared(region), mat).empty());
}

TEST_CASE("Filleted rectangles: curvature ends, corner edges and a plateau", "[M5][RB-03][VER-02]")
{
    SECTION("one filleted corner") {
        const Region region = load("filleted-rectangle.json");
        const MedialAxis mat = mat_of(region);
        CHECK(count(mat, VertexKind::curvature_end) == 1);
        CHECK(plateaus(mat).size() == 1);
        CHECK(check_mat(prepared(region), mat).empty());
    }
    SECTION("all four corners") {
        for (double f : {0.2, 0.5, 0.9}) {
            INFO("fraction " << f);
            const Region region = filleted({{0, 0}, {100, 0}, {100, 40}, {0, 40}}, {f, f, f, f});
            const MedialAxis mat = mat_of(region);
            CHECK(check_mat(prepared(region), mat).empty());
            CHECK(plateaus(mat).size() == 1);
            CHECK(count(mat, VertexKind::curvature_end) == 4);
        }
    }
}

TEST_CASE("Corridors: plateaus ending at reflex corners", "[M5][RB-03][VER-03][VER-06]")
{
    const std::vector<std::pair<std::string, std::vector<IPoint>>> shapes = {
        {"T", {{0, 0}, {30, 0}, {30, 10}, {20, 10}, {20, 40}, {10, 40}, {10, 10}, {0, 10}}},
        {"plus", {{10, 0}, {20, 0}, {20, 10}, {30, 10}, {30, 20}, {20, 20}, {20, 30}, {10, 30}, {10, 20}, {0, 20}, {0, 10}, {10, 10}}},
        {"H", {{0, 0}, {10, 0}, {10, 15}, {20, 15}, {20, 0}, {30, 0}, {30, 40}, {20, 40}, {20, 25}, {10, 25}, {10, 40}, {0, 40}}},
        {"staircase", {{0, 0}, {40, 0}, {40, 10}, {30, 10}, {30, 20}, {20, 20}, {20, 30}, {10, 30}, {10, 40}, {0, 40}}},
        {"offset corridor", {{0, 0}, {50, 0}, {50, 30}, {35, 30}, {35, 10}, {15, 10}, {15, 30}, {0, 30}}},
        {"L, equal arms", {{0, 0}, {40, 0}, {40, 10}, {10, 10}, {10, 40}, {0, 40}}},
    };
    FamilyRun family;
    for (const auto& [name, poly] : shapes) {
        std::vector<IPoint> scaled;
        for (const IPoint& q : poly) scaled.push_back({q.first * 1000, q.second * 1000});
        family.run(name, to_region(scaled), &scaled);
    }
    INFO(family.first);
    CHECK(family.failures == 0);
}

TEST_CASE("Degeneracy suite", "[M5][RB-02][RB-03][RB-04][EV-12][VER-03][VER-06]")
{
    const int shapes = degenerate_shape_count();
    std::mt19937_64 rng(20261010);
    FamilyRun family;
    for (int k = 0; k < shapes; ++k) {
        const std::string id = std::to_string(k);
        const bool debug = k % 8 < 2;  // full invariants on a quarter of each family
        switch (k % 6) {
        case 0:
        case 1: {
            std::uniform_int_distribution<int> n(1, 24);
            const auto poly = random_columns(rng, n(rng), k % 6 == 0);
            family.run("columns-" + id, to_region(poly), &poly, debug);
            break;
        }
        case 2: {
            std::uniform_int_distribution<int> n(1, 8);
            const auto poly = random_symmetric(rng, n(rng));
            family.run("symmetric-" + id, to_region(poly), &poly, debug);
            break;
        }
        case 3: {
            std::uniform_int_distribution<int> n(3, 40);
            family.run("tangential-" + id, random_tangential(rng, n(rng)), nullptr, debug);
            break;
        }
        case 4: {
            std::uniform_int_distribution<int> n(3, 32);
            family.run("gear-" + id, gear(n(rng)), nullptr, debug);
            break;
        }
        default: {
            // A rectangle or a slot with rounded corners, at a random size and position.
            std::uniform_int_distribution<int> side(1, 50);
            const double w = side(rng), h = side(rng);
            std::uniform_real_distribution<double> f(0.0, 1.0);
            const double fr = k % 12 == 5 ? 0.0 : f(rng);
            family.run("rounded-" + id, filleted({{0, 0}, {w, 0}, {w, h}, {0, h}}, {fr, fr, fr, fr}), nullptr, debug);
            break;
        }
        }
    }
    INFO(family.first);
    CHECK(family.failures == 0);
}

TEST_CASE("Simple handlers agree with cluster resolution", "[M5][RB-04][VER-05]")
{
    const int shapes = degenerate_shape_count();
    std::mt19937_64 rng(4);
    std::uniform_int_distribution<int> size(3, 30);
    int failures = 0;
    std::string first;
    for (int k = 0; k < shapes; ++k) {
        Region region = k % 3 == 0   ? to_region(random_star_polygon(rng, size(rng)))
                        : k % 3 == 1 ? random_filleted(rng, size(rng))
                                     : random_arc_chain(rng, size(rng));
        auto prep = prepare(region);
        if (!prep) continue;
        Options options;
        auto simple = compute_mat(*prep, options);
        options.resolve_all_events = true;
        auto all = compute_mat(*prep, options);
        std::string why;
        if (!simple || !all)
            why = to_string(!simple ? simple.error() : all.error());
        else
            why = same_mat(*simple, *all, 1e-9 * prep->xf.scale);
        if (!why.empty() && failures++ == 0) first = "shape " + std::to_string(k) + ": " + why;
    }
    INFO(first);
    CHECK(failures == 0);
}
