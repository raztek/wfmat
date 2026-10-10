# wfmat

Wavefront (grassfire) computation of the interior medial axis transform of a planar region
bounded by one closed loop of line segments and circular arcs.

The design and specification live in [`docs/spec.md`](docs/spec.md) (Markdown with LaTeX maths);
this repository is the master copy. `bash docs/build-pdf.sh` rebuilds
`docs/Wavefront-MAT-Spec-v<version>.pdf` (needs pandoc, Node and Playwright with Chromium).
Each approved spec version is tagged `spec-vX.Y`. The public API is documented in
[`docs/api.md`](docs/api.md).

## Status

- M0 (scaffold): input in bulge form from C++ or JSON, validation and normalisation, corner
  classification, SVG rendering and a brute-force distance validator.
- M1 (kernel): sites and their distance functions, the three-site solve, shock positions and
  velocities, regular-vertex hits, contact times, and binary128 refinement (`wfmat/kernel.hpp`).

- M2 (polygons): the propagation engine for line sites and reflex-corner point sites, with events
  E1a (collapse), E1b (transition), E2 (split) and E3 (contact), all-pairs scheduling, and the MAT
  as an embedded graph (`wfmat/mat.hpp`, `compute_mat`). It matches Boost.Polygon's Voronoi
  diagram on 10^4 random polygons in the test suite.
- M3 (arcs): convex and concave arc sites, tangent joins, E1c (a convex arc shrinking to its
  centre, a curvature end), and loops of two elements. 10^4 random filleted polygons and arc
  chains pass the property checks, and a sampled Voronoi oracle (Boost.Polygon on the boundary
  flattened to fine polylines) converges to the computed axis.

- M5 (degeneracies): simultaneous events are gathered into clusters and resolved by one generic
  procedure, so a square's four corners meet in one vertex of degree 4, and plateaus (the
  constant-radius edges of rectangles, slots and curved slots) are plateau edges. A degeneracy
  suite (orthogonal column polygons, shapes symmetric about both axes, tangential polygons,
  gears, rounded rectangles, regular polygons up to 1024 sides) passes the property checks and,
  for integer polygons, matches Boost.Polygon's Voronoi diagram. Sending every event through the
  cluster resolution gives the same axis as the simple handlers (`Options::resolve_all_events`).

- M6 (performance): non-local events come from a windowed R-tree broad phase
  (`BroadPhase::windowed_rtree`, the default), which processes exactly the events of the
  all-pairs mode, bit for bit, and is 10 to 200 times faster on 1,000 segments. The event heap
  drops stale events in bulk. The spec's speed targets are not met yet; section 9 of the spec lists
  the measured times, and milestone M8 is to close the gap.
- M7 (release, version 1.0.0): edges exported exactly as rational quadratic Bézier arcs and
  flattened to a chord tolerance, inward offsets from the front at any distance, MAT JSON output,
  an SVG with the axis coloured by radius and the offsets, CLI options for all of these, and the
  API documentation. Merging a project version of 1.0.0 or later tags `vX.Y.Z`.

`ErrorCode::unsupported` is now only the code for valid input that a later milestone handles;
no input currently gets it. See section 13 of the spec.

## Building

CMake 3.25 or newer and a C++20 compiler (GCC 12+, Clang 15+, MSVC 19.36+).

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build
ctest --test-dir build --output-on-failure
./build/wfmat-cli tests/data/l-shape.json --svg l-shape.svg --json l-shape.json --flatten 0.1 --offset 5
```

Dependencies are found as installed packages (a vcpkg manifest is provided: pass
`-DCMAKE_TOOLCHAIN_FILE=<vcpkg>/scripts/buildsystems/vcpkg.cmake`, or use system packages such as
Ubuntu's `libboost-dev nlohmann-json3-dev catch2`). Anything missing is fetched with
FetchContent; `-DWFMAT_FETCH_DEPS=ON` fetches everything. `tl::expected` falls back to the copy
in `third_party/`.

With Visual Studio 2026, set `VCPKG_ROOT` (the Visual Studio vcpkg component or a vcpkg clone)
and use the `vs2026` preset, from a terminal or by opening the folder in Visual Studio:

```bat
cmake --preset vs2026
cmake --build --preset vs2026-release
ctest --preset vs2026-release
```

| Dependency | Licence | Used by |
| --- | --- | --- |
| Boost.Geometry, Multiprecision, Math (header-only) | BSL-1.0 | core |
| tl::expected | CC0 | core |
| nlohmann/json | MIT | `wfmat_io` |
| Catch2 v3 | BSL-1.0 | tests |
| Boost.Polygon (header-only) | BSL-1.0 | tests (Voronoi oracle) |
| Google Benchmark | Apache-2.0 | benchmarks (optional, `-DWFMAT_BUILD_BENCHMARKS=ON`) |

The library is compiled with `-ffp-contract=off` (`/fp:precise` on MSVC) and never with
fast-math [LIB-03].

## Layout

| Path | Contents |
| --- | --- |
| `include/wfmat/` | public headers: `geom`, `region`, `options`, `result`, `prepare`, `kernel`, `mat`, `validate`, `io/json`, `io/svg` |
| `src/` | library sources (`wfmat` core with the engine in `engine.cpp` and cluster resolution in `cluster.cpp`, `wfmat_io` for JSON and SVG) |
| `tools/wfmat-cli.cpp` | command-line tool |
| `examples/basic.cpp` | the API example (`wfmat-example`) |
| `tests/` | Catch2 tests; `tests/data/` holds the canonical shapes and kernel fixtures |
| `bench/` | Google Benchmark suite (`wfmat-bench`) and its scaled shape families |
| `tools/fixtures/` | fixture generator for the kernel tests (SymPy, `gen_kernel_fixtures.py`) |
| `docs/` | specification sources and PDF, API documentation (`api.md`) |

## Tests and requirement IDs

Every test names the spec requirements it verifies as Catch2 tags, for example `[IN-04]`, so
`./build/tests/wfmat_tests "[IN-04]"` runs the tests for one requirement.

The Voronoi comparison runs on 10^4 random polygons, and the property checks on 10^4 random
shapes with arcs (about 30 s each in a release build); `WFMAT_RANDOM_POLYGONS=500` and
`WFMAT_RANDOM_SHAPES=500` set smaller counts for quick local runs. With `WFMAT_DUMP=1` the arc
test writes every failing shape to `build/tests/output/` as JSON, for `wfmat-cli`. The degeneracy
suite and the cross-check run 2000 shapes each; `WFMAT_DEGENERATE_SHAPES` sets the count, and
`WFMAT_DUMP=1` writes their failures there too. The broad-phase test compares both broad phases
on 2000 random shapes (`WFMAT_RANDOM_SHAPES`).

## Benchmarks

With Google Benchmark installed (`libbenchmark-dev` on Ubuntu), configure with
`-DWFMAT_BUILD_BENCHMARKS=ON` and run `./build/wfmat-bench`. It times `prepare` and
`compute_mat` on gears, filleted stars, wavy outlines and spiky stars at 1,000, 10,000 and 100,000
segments; `--benchmark_filter='/1000$'` keeps the small sizes.

## Licence

MIT; see [`LICENSE`](LICENSE). Third-party code keeps its own licence (`third_party/README.md`).
