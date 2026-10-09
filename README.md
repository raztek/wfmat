# wfmat

Wavefront (grassfire) computation of the interior medial axis transform of a planar region
bounded by one closed loop of line segments and circular arcs.

The design and specification live in [`docs/spec.md`](docs/spec.md) (Markdown with LaTeX maths);
this repository is the master copy. `bash docs/build-pdf.sh` rebuilds
`docs/Wavefront-MAT-Spec-v<version>.pdf` (needs pandoc, Node and Playwright with Chromium).
Each approved spec version is tagged `spec-vX.Y`.

## Status

- M0 (scaffold): input in bulge form from C++ or JSON, validation and normalisation, corner
  classification, SVG rendering and a brute-force distance validator.
- M1 (kernel): sites and their distance functions, the three-site solve, shock positions and
  velocities, regular-vertex hits, contact times, and binary128 refinement (`wfmat/kernel.hpp`).

- M2 (polygons): the propagation engine for line sites and reflex-corner point sites, with events
  E1a (collapse), E1b (transition), E2 (split) and E3 (contact), all-pairs scheduling, and the MAT
  as an embedded graph (`wfmat/mat.hpp`, `compute_mat`). It matches Boost.Polygon's Voronoi
  diagram on 10^4 random polygons in the test suite.

Not yet handled, and refused with `ErrorCode::unsupported`: circular arcs and tangent joins (M3),
and simultaneous events such as the four corners of a square meeting at its centre, or the plateau
of a rectangle (M5). See section 13 of the spec.

## Building

CMake 3.25 or newer and a C++20 compiler (GCC 12+, Clang 15+, MSVC 19.36+).

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build
ctest --test-dir build --output-on-failure
./build/wfmat-cli tests/data/equilateral-triangle.json --svg triangle.svg
```

Dependencies are found as installed packages (a vcpkg manifest is provided: pass
`-DCMAKE_TOOLCHAIN_FILE=<vcpkg>/scripts/buildsystems/vcpkg.cmake`, or use system packages such as
Ubuntu's `libboost-dev nlohmann-json3-dev catch2`). Anything missing is fetched with
FetchContent; `-DWFMAT_FETCH_DEPS=ON` fetches everything. `tl::expected` falls back to the copy
in `third_party/`.

| Dependency | Licence | Used by |
| --- | --- | --- |
| Boost.Geometry, Multiprecision, Math (header-only) | BSL-1.0 | core |
| tl::expected | CC0 | core |
| nlohmann/json | MIT | `wfmat_io` |
| Catch2 v3 | BSL-1.0 | tests |
| Boost.Polygon (header-only) | BSL-1.0 | tests (Voronoi oracle) |

The library is compiled with `-ffp-contract=off` (`/fp:precise` on MSVC) and never with
fast-math [LIB-03].

## Layout

| Path | Contents |
| --- | --- |
| `include/wfmat/` | public headers: `geom`, `region`, `options`, `result`, `prepare`, `kernel`, `mat`, `validate`, `io/json`, `io/svg` |
| `src/` | library sources (`wfmat` core with the engine in `engine.cpp`, `wfmat_io` for JSON and SVG) |
| `tools/wfmat-cli.cpp` | command-line tool |
| `tests/` | Catch2 tests; `tests/data/` holds the canonical shapes and kernel fixtures |
| `tools/fixtures/` | fixture generator for the kernel tests (SymPy, `gen_kernel_fixtures.py`) |
| `docs/` | specification sources and PDF |

## Tests and requirement IDs

Every test names the spec requirements it verifies as Catch2 tags, for example `[IN-04]`, so
`./build/tests/wfmat_tests "[IN-04]"` runs the tests for one requirement.

The Voronoi comparison runs on 10^4 random polygons (about 30 s in a release build);
`WFMAT_RANDOM_POLYGONS=500` sets a smaller count for quick local runs.

## Licence

MIT; see [`LICENSE`](LICENSE). Third-party code keeps its own licence (`third_party/README.md`).
