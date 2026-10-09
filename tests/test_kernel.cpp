#include <algorithm>
#include <cmath>
#include <fstream>
#include <numbers>
#include <string>

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <nlohmann/json.hpp>

#include "support.hpp"
#include "wfmat/kernel.hpp"

using namespace wfmat;
using namespace wfmat::test;
using Catch::Approx;

namespace {

struct Fixture {
    std::string name;
    std::array<Site, 3> sites;
    double t_now = 0.0;
    bool degenerate = false;
    std::vector<std::array<double, 3>> roots;  // x, y, t
};

Site site_from_json(const nlohmann::json& j)
{
    if (j.at("kind") == "line") return Site::line({j["n"][0], j["n"][1]}, j["c"]);
    return Site::circle({j["centre"][0], j["centre"][1]}, j["R"], j["sigma"]);
}

std::vector<Fixture> load_fixtures()
{
    std::ifstream in(data_path("kernel/three_site.json"));
    REQUIRE(in);
    const auto doc = nlohmann::json::parse(in);
    std::vector<Fixture> out;
    for (const auto& c : doc.at("cases")) {
        Fixture f;
        f.name = c.at("name");
        for (int i = 0; i < 3; ++i) f.sites[i] = site_from_json(c.at("sites")[i]);
        f.t_now = c.at("t_now");
        f.degenerate = c.at("degenerate");
        for (const auto& r : c.at("roots")) f.roots.push_back({r["x"], r["y"], r["t"]});
        out.push_back(std::move(f));
    }
    return out;
}

bool close(double a, double b, double rel = 1e-11)
{
    return std::abs(a - b) <= rel * std::max(1.0, std::abs(b));
}

bool on_site(const Site& s, Vec2 foot)
{
    if (s.is_line()) return std::abs(dot(s.n, foot) - s.c) <= 1e-10;
    return std::abs(dist(foot, s.centre) - s.R) <= 1e-10 * std::max(1.0, s.R);
}

} // namespace

TEST_CASE("three-site solve matches the reference roots", "[VER-01][KN-01][KN-02]")
{
    const auto fixtures = load_fixtures();
    REQUIRE(fixtures.size() >= 60);
    std::size_t quad_refined = 0;
    for (const Fixture& f : fixtures) {
        INFO(f.name);
        const Solve3 s = solve_three(f.sites[0], f.sites[1], f.sites[2], f.t_now, Tolerances{});
        CHECK(s.degenerate == f.degenerate);
        REQUIRE(s.roots.size() == f.roots.size());
        for (std::size_t i = 0; i < f.roots.size(); ++i) {
            INFO("root " << i);
            const Root3& r = s.roots[i];
            CHECK(close(r.p.x, f.roots[i][0]));
            CHECK(close(r.p.y, f.roots[i][1]));
            CHECK(close(r.t, f.roots[i][2]));
            for (int k = 0; k < 3; ++k) {
                CHECK(on_site(f.sites[k], r.feet[k]));
                CHECK(std::abs(f.sites[k].distance(r.p) - r.t) <= 1e-12 * std::max(1.0, r.t));
            }
            quad_refined += r.refined_quad;
        }
        if (!f.roots.empty()) {
            const auto first = earliest(s);
            REQUIRE(first);
            CHECK(first->t == s.roots.front().t);
        } else {
            CHECK_FALSE(earliest(s));
        }
    }
    CHECK(quad_refined > 0);
}

TEST_CASE("the solve does not depend on the site order", "[KN-01]")
{
    for (const Fixture& f : load_fixtures()) {
        if (f.degenerate) continue;
        INFO(f.name);
        const Solve3 s = solve_three(f.sites[2], f.sites[0], f.sites[1], f.t_now, Tolerances{});
        REQUIRE(s.roots.size() == f.roots.size());
        for (std::size_t i = 0; i < f.roots.size(); ++i) CHECK(close(s.roots[i].t, f.roots[i][2]));
    }
}

TEST_CASE("roots before t_now are rejected", "[KN-01]")
{
    const auto fixtures = load_fixtures();
    const auto it = std::find_if(fixtures.begin(), fixtures.end(), [](const Fixture& f) { return f.roots.size() == 2; });
    REQUIRE(it != fixtures.end());
    const double between = 0.5 * (it->roots[0][2] + it->roots[1][2]);
    const Solve3 s = solve_three(it->sites[0], it->sites[1], it->sites[2], between, Tolerances{});
    REQUIRE(s.roots.size() == 1);
    CHECK(close(s.roots[0].t, it->roots[1][2]));
}

TEST_CASE("ill-conditioned roots are refined in binary128", "[KN-06][VER-01]")
{
    for (const Fixture& f : load_fixtures()) {
        if (f.name.find("near-tangent") == std::string::npos || f.name.find("e-10") == std::string::npos) continue;
        INFO(f.name);
        const Solve3 s = solve_three(f.sites[0], f.sites[1], f.sites[2], 0.0, Tolerances{});
        REQUIRE(s.roots.size() == 2);
        for (const Root3& r : s.roots) {
            CHECK(r.refined_quad);
            CHECK(r.condition > 1e4);
        }
    }

    SECTION("refinement restores a perturbed root")
    {
        const auto fixtures = load_fixtures();
        const Fixture& f = fixtures.front();
        Root3 r;
        r.p = {f.roots[0][0] + 1e-7, f.roots[0][1] - 1e-7};
        r.t = f.roots[0][2] + 1e-7;
        const Root3 d = refine(f.sites, r, false);
        const Root3 q = refine(f.sites, r, true);
        CHECK(close(d.t, f.roots[0][2], 1e-14));
        CHECK(close(q.t, f.roots[0][2], 1e-15));
        CHECK(q.refined_quad);
    }
}

TEST_CASE("bisector conics by site pair", "[KN-03]")
{
    const double eps = 1e-12;
    const Site lx = Site::line({1, 0}, 0), ly = Site::line({0, 1}, 0), ly2 = Site::line({0, -1}, -2);
    const Site convex = Site::circle({0, 0}, 3, 1), concave = Site::circle({5, 0}, 1, -1);
    CHECK(bisector_kind(lx, ly, eps) == ConicKind::line);
    CHECK(bisector_kind(ly, ly2, eps) == ConicKind::plateau_line);
    CHECK(bisector_kind(lx, convex, eps) == ConicKind::parabola);
    CHECK(bisector_kind(lx, Site::point({1, 1}), eps) == ConicKind::parabola);
    CHECK(bisector_kind(convex, concave, eps) == ConicKind::ellipse);
    CHECK(bisector_kind(concave, Site::circle({-5, 0}, 2, -1), eps) == ConicKind::hyperbola);
    CHECK(bisector_kind(Site::point({0, 0}), Site::point({1, 0}), eps) == ConicKind::line);
    CHECK(bisector_kind(convex, Site::circle({0, 0}, 1, -1), eps) == ConicKind::plateau_arc);
}

TEST_CASE("shock positions and the branch selector", "[KN-03]")
{
    SECTION("two lines: affine in t, speed 1 / sin(alpha / 2)")
    {
        const Site a = Site::line({0, 1}, 0), b = Site::line({1, 0}, 0);  // the corner at the origin
        for (double t : {0.0, 0.5, 2.0}) {
            const auto p = shock_position(a, b, t, 0);
            REQUIRE(p);
            CHECK(p->x == Approx(t));
            CHECK(p->y == Approx(t));
        }
        const auto v = shock_velocity(a, b, {1, 1});
        REQUIRE(v);
        CHECK(norm(*v) == Approx(1.0 / std::sin(std::numbers::pi / 4)));
    }
    SECTION("line and point: a parabola with two branches")
    {
        const Site line = Site::line({0, 1}, 0);
        const Site point = Site::point({0, 2});
        // At t = 2: y = 2 and x^2 + (y - 2)^2 = 4, so x = +-2.
        const auto pts = offset_intersections(line, point, 2.0);
        REQUIRE(pts.size() == 2);
        const auto right = shock_position(line, point, 2.0, branch_of(line, point, {2, 2}));
        REQUIRE(right);
        CHECK(right->x == Approx(2.0));
        CHECK(right->y == Approx(2.0));
        CHECK(branch_of(line, point, {2, 2}) == -branch_of(line, point, {-2, 2}));
        // Before the fronts touch (t < 1) there is no shock.
        CHECK(offset_intersections(line, point, 0.5).empty());
        // The branch is constant along the edge.
        for (double t : {1.5, 3.0, 10.0}) {
            const auto p = shock_position(line, point, t, -1);
            REQUIRE(p);
            CHECK(p->x > 0.0);
            CHECK(line.distance(*p) == Approx(t));
            CHECK(point.distance(*p) == Approx(t));
        }
    }
    SECTION("two circles")
    {
        const Site convex = Site::circle({0, 0}, 3, 1), concave = Site::circle({2, 0}, 0.5, -1);
        for (double t : {0.5, 1.0}) {
            for (Vec2 p : offset_intersections(convex, concave, t)) {
                CHECK(convex.distance(p) == Approx(t));
                CHECK(concave.distance(p) == Approx(t));
            }
        }
        // A convex circle's front is gone once t > R.
        CHECK(offset_intersections(convex, concave, 3.5).empty());
    }
}

TEST_CASE("regular vertices meet fronts in closed form", "[KN-04]")
{
    const Tolerances tol;
    const Vec2 q0{0, 0}, up{0, 1};
    CHECK(*regular_hit_time(q0, up, Site::line({0, -1}, -3), 0.0, tol) == Approx(1.5));       // y = 3 facing down
    CHECK(*regular_hit_time(q0, up, Site::point({0, 4}), 0.0, tol) == Approx(2.0));          // reflex corner
    CHECK(*regular_hit_time(q0, up, Site::circle({0, 1}, 5, 1), 0.0, tol) == Approx(3.0));   // convex arc above
    CHECK_FALSE(regular_hit_time(q0, up, Site::line({0, 1}, -1), 0.0, tol));                 // moving with the front
    CHECK_FALSE(regular_hit_time(q0, up, Site::point({0, 4}), 2.5, tol));                    // already past
}

TEST_CASE("contacts along common normals", "[KN-05]")
{
    const Tolerances tol;
    SECTION("antiparallel lines: a plateau at half the gap")
    {
        const auto c = contacts(Site::line({0, 1}, 0), Site::line({0, -1}, -2), 0.0, tol);
        REQUIRE(c.size() == 1);
        CHECK(c[0].plateau);
        CHECK(c[0].t == Approx(1.0));
        CHECK(c[0].p.y == Approx(1.0));
    }
    SECTION("non-parallel lines never make contact")
    {
        CHECK(contacts(Site::line({0, 1}, 0), Site::line({1, 0}, 0), 0.0, tol).empty());
    }
    SECTION("line and convex arc facing each other")
    {
        const auto c = contacts(Site::line({0, -1}, -4), Site::circle({0, 3}, 2, 1), 0.0, tol);
        REQUIRE(c.size() == 1);
        CHECK_FALSE(c[0].plateau);
        CHECK(c[0].t == Approx(1.5));
        CHECK(c[0].p.x == Approx(0.0).margin(1e-15));
        CHECK(c[0].p.y == Approx(2.5));
        CHECK(c[0].feet[0].y == Approx(4.0));
        CHECK(c[0].feet[1].y == Approx(1.0));
    }
    SECTION("line and reflex corner")
    {
        const auto c = contacts(Site::point({0, 3}), Site::line({0, 1}, 0), 0.0, tol);
        REQUIRE(c.size() == 1);
        CHECK(c[0].t == Approx(1.5));
        CHECK(c[0].p.y == Approx(1.5));
    }
    SECTION("two reflex corners")
    {
        const auto c = contacts(Site::point({0, 0}), Site::point({4, 0}), 0.0, tol);
        REQUIRE(c.size() == 1);
        CHECK(c[0].t == Approx(2.0));
        CHECK(c[0].p.x == Approx(2.0));
        CHECK(contacts(Site::point({0, 0}), Site::point({4, 0}), 2.5, tol).empty());
    }
    SECTION("convex and concave arcs: contacts on both sides of the centre line")
    {
        const Site outer = Site::circle({0, 0}, 4, 1), inner = Site::circle({1, 0}, 1, -1);
        const auto c = contacts(outer, inner, 0.0, tol);
        REQUIRE(c.size() == 2);
        // The gap from x = 2 to x = 4 closes at t = 1; the gap from x = -4 to x = 0 at t = 2.
        CHECK(c[0].t == Approx(1.0));
        CHECK(c[0].p.x == Approx(3.0));
        CHECK(c[1].t == Approx(2.0));
        CHECK(c[1].p.x == Approx(-2.0));
        for (const Contact& k : c) {
            CHECK(outer.distance(k.p) == Approx(k.t));
            CHECK(inner.distance(k.p) == Approx(k.t));
        }
    }
    SECTION("concentric convex and concave arcs: a plateau arc")
    {
        const auto c = contacts(Site::circle({0, 0}, 3, 1), Site::circle({0, 0}, 1, -1), 0.0, tol);
        REQUIRE(c.size() == 1);
        CHECK(c[0].plateau);
        CHECK(c[0].t == Approx(1.0));
    }
}

TEST_CASE("sites of prepared segments", "[KN-01]")
{
    const PreparedRegion r = prepared(load("spec-example.json"));
    for (std::size_t i = 0; i < r.segments.size(); ++i) {
        const Segment& s = r.segments[i];
        const Site site = site_of(s);
        // Every segment point is at distance 0, and the region side is positive.
        CHECK(site.distance(s.midpoint()) == Approx(0.0).margin(1e-14));
        const Vec2 inward = perp(s.is_arc() ? (s.sweep > 0 ? 1.0 : -1.0) * unit(perp(s.midpoint() - s.c)) : unit(s.b - s.a));
        CHECK(site.distance(s.midpoint() + 1e-3 * inward) > 0.0);
    }
}
