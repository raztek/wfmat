// [M6] The windowed broad phase [EV-11] against the all-pairs mode [EV-12], and heap compaction
// [PF-01]: both modes must process the same events in the same order.
#include <cstdlib>
#include <random>
#include <string>

#include <catch2/catch_test_macros.hpp>

#include "mat_checks.hpp"
#include "support.hpp"
#include "wfmat/mat.hpp"

using namespace wfmat;
using namespace wfmat::test;

namespace {

int shape_count()
{
    if (const char* s = std::getenv("WFMAT_RANDOM_SHAPES")) return std::max(1, std::atoi(s));
    return 2'000;
}

// Both broad phases on one region; an empty string when they agree bit for bit.
std::string compare_modes(const Region& region, bool debug)
{
    auto prep = prepare(region);
    if (!prep) return {};
    Options options;
    options.debug_checks = debug;
    auto windowed = compute_mat(*prep, options);
    options.broad_phase = BroadPhase::all_pairs;
    auto all = compute_mat(*prep, options);
    if (!windowed || !all) return to_string(!windowed ? windowed.error() : all.error());
    if (auto why = check_mat(*prep, *windowed, 1e-9, 4); !why.empty()) return why;
    return identical_mat(*windowed, *all);
}

} // namespace

TEST_CASE("The windowed broad phase finds the events of the all-pairs mode", "[M6][EV-11][EV-12][VER-05]")
{
    const int shapes = shape_count();
    std::mt19937_64 rng(6);
    std::uniform_int_distribution<int> size(3, 120);
    int failures = 0;
    std::string first;
    for (int k = 0; k < shapes; ++k) {
        const int n = size(rng);
        const Region region = k % 3 == 0   ? to_region(random_star_polygon(rng, n))
                              : k % 3 == 1 ? random_filleted(rng, n)
                                           : random_arc_chain(rng, n);
        if (auto why = compare_modes(region, k % 16 == 0); !why.empty() && failures++ == 0)
            first = "shape " + std::to_string(k) + ": " + why;
    }
    INFO(first);
    CHECK(failures == 0);
}

TEST_CASE("Canonical shapes: identical events in both broad phases", "[M6][EV-12][VER-02]")
{
    for (const char* name : {"concave-top.json", "disk.json", "equilateral-triangle.json", "filleted-rectangle.json", "l-shape.json",
                             "rectangle.json", "regular-hexagon.json", "slot.json", "spec-example.json"}) {
        INFO(name);
        CHECK(compare_modes(load(name), true).empty());
    }
}

TEST_CASE("Windows: several open on a large shape and none is needed in all-pairs mode", "[M6][EV-11]")
{
    std::mt19937_64 rng(11);
    auto prep = prepare(to_region(random_star_polygon(rng, 400)));
    REQUIRE(prep);
    auto windowed = compute_mat(*prep);
    REQUIRE(windowed);
    CHECK(windowed->stats().windows > 1);
    Options options;
    options.broad_phase = BroadPhase::all_pairs;
    auto all = compute_mat(*prep, options);
    REQUIRE(all);
    CHECK(all->stats().windows == 0);
    CHECK(identical_mat(*windowed, *all).empty());
}
