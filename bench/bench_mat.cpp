// [VER-07] Benchmarks: scaled shape families at 10^3 to 10^5 segments, tracked against the [PF-02]
// targets. Each family and size has a prepare benchmark and a compute_mat benchmark on the
// prepared region; the segment count is reported as a counter.
#include <cstdio>
#include <functional>
#include <string>

#include <benchmark/benchmark.h>

#include "shapes.hpp"
#include "wfmat/mat.hpp"
#include "wfmat/prepare.hpp"

namespace {

using wfmat::Region;

struct Family {
    const char* name;
    int per_unit;  // segments per unit of the generator's size argument
    std::function<Region(int)> make;
};

const Family families[] = {
    {"gear", 4, [](int n) { return wfmat::bench::gear(n); }},
    {"filleted_star", 2, [](int n) { return wfmat::bench::filleted_star(n); }},
    {"wavy", 1, [](int n) { return wfmat::bench::wavy_polygon(n); }},
    {"star", 1, [](int n) { return wfmat::bench::star_polygon(n); }},
};

void register_family(const Family& f, int segments, wfmat::BroadPhase broad_phase)
{
    const int n = segments / f.per_unit;
    const std::string size = std::to_string(segments);
    const bool all_pairs = broad_phase == wfmat::BroadPhase::all_pairs;
    if (!all_pairs) {
        benchmark::RegisterBenchmark(("prepare/" + std::string(f.name) + "/" + size).c_str(), [f, n](benchmark::State& state) {
            const Region region = f.make(n);
            for (auto _ : state) benchmark::DoNotOptimize(wfmat::prepare(region));
        })->Unit(benchmark::kMillisecond);
    }
    const std::string name = std::string(all_pairs ? "mat_all_pairs/" : "mat/") + f.name + "/" + size;
    benchmark::RegisterBenchmark(name.c_str(), [f, n, broad_phase](benchmark::State& state) {
        const auto prepared = wfmat::prepare(f.make(n));
        if (!prepared) {
            state.SkipWithError(wfmat::to_string(prepared.error()).c_str());
            return;
        }
        wfmat::Options options;
        options.broad_phase = broad_phase;
        for (auto _ : state) {
            auto mat = wfmat::compute_mat(*prepared, options);
            if (!mat) {
                state.SkipWithError(wfmat::to_string(mat.error()).c_str());
                return;
            }
            benchmark::DoNotOptimize(mat);
        }
        state.counters["segments"] = static_cast<double>(prepared->segments.size());
    })->Unit(benchmark::kMillisecond);
}

} // namespace

int main(int argc, char** argv)
{
    for (const Family& f : families) {
        for (const int segments : {1'000, 10'000, 100'000}) register_family(f, segments, wfmat::BroadPhase::windowed_rtree);
        register_family(f, 1'000, wfmat::BroadPhase::all_pairs);  // the [EV-12] reference, for scale
    }
    benchmark::Initialize(&argc, argv);
    if (benchmark::ReportUnrecognizedArguments(argc, argv)) return 1;
    benchmark::RunSpecifiedBenchmarks();
    benchmark::Shutdown();
    return 0;
}
