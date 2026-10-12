# Building and running wfmat with Visual Studio 2026

Step-by-step instructions for Windows. Verified for wfmat 2.0.0 with Visual Studio Community
2026 (18.10), its bundled CMake 4.3.1 and its bundled vcpkg: configure, Release build and all 95
tests pass. For other platforms and compilers see [Building](../README.md#building) in the README.
The input file format is described in [input-format.md](input-format.md).

## One-time setup

1. **Visual Studio components.** Install the "Desktop development with C++" workload, which
   includes "C++ CMake tools for Windows" and "vcpkg package manager".
2. **vcpkg baseline.** Visual Studio's bundled vcpkg only installs packages when the project pins
   a registry version. The repository's `vcpkg-configuration.json` does this (the vcpkg 2026.07.29
   release), so nothing needs to be added by hand.
3. **Dependencies.** The first configure has vcpkg build Boost (geometry, math, multiprecision,
   polygon), nlohmann-json, tl-expected and Catch2 from source. That takes several minutes and is
   cached afterwards (in `build\vs2026\vcpkg_installed` and vcpkg's binary cache).

## Getting or updating the code

```powershell
git clone https://github.com/raztek/wfmat.git
cd wfmat
```

To update an existing clone later:

```powershell
git fetch --tags origin
git pull --ff-only origin main
git describe --tags                   # the version tag, e.g. v2.0.0
```

Then rebuild and run the tests as in option A or B. The configure step only rebuilds dependencies
when `vcpkg.json` changed.

## Option A: inside the Visual Studio IDE

1. Start Visual Studio 2026, choose **Open a local folder** (or File > Open > Folder) and pick the
   `wfmat` folder. Visual Studio sees `CMakePresets.json` and configures.
2. In the toolbar, set the configure preset to **Visual Studio 2026 (x64, vcpkg)** and the build
   preset to **vs2026-release** (or **vs2026-debug** for debugging).
   - If configure fails with a message about `VCPKG_ROOT` or a missing `vcpkg.cmake`, set the
     variable once and restart Visual Studio (adjust the edition in the path if needed):
     `setx VCPKG_ROOT "C:\Program Files\Microsoft Visual Studio\18\Community\VC\vcpkg"`
   - After pulling a new version, use **Project > Delete Cache and Reconfigure** if the build
     behaves oddly.
3. **Build > Build All** (Ctrl+Shift+B).
4. **Tests:** open **Test > Test Explorer** and click **Run All**. All 95 tests should pass.
5. **Run the CLI:** in the startup item dropdown (green arrow) pick **wfmat-cli.exe**. To give it
   arguments, right-click `CMakeLists.txt` in Solution Explorer > **Add Debug Configuration** (or
   **Debug > Debug and Launch Settings for wfmat-cli**) and add, in `launch.vs.json`:

   ```json
   "args": [ "tests/data/plate-with-holes.json", "--svg", "out/plate-with-holes.svg",
             "--json", "out/plate-with-holes.mat.json", "--flatten", "0.1", "--offset", "3" ],
   "currentDir": "${workspaceRoot}"
   ```

   Create the `out` folder first. Press F5 or Ctrl+F5.

## Option B: from a terminal

Open **Developer PowerShell for VS 2026** (Start menu). It puts Visual Studio's `cmake` and
`ctest` on the PATH and sets `VCPKG_ROOT` to the bundled vcpkg. Then, in the `wfmat` folder:

```powershell
cmake --preset vs2026                 # configure (generates build\vs2026\wfmat.slnx)
cmake --build --preset vs2026-release # compile Release
ctest --preset vs2026-release         # 95 tests, about 80 s
```

The last line of `ctest` should read:

```text
100% tests passed, 0 tests failed out of 95
```

Run one requirement's tests only: `build\vs2026\tests\Release\wfmat_tests.exe "[IN-04]"`.

After the configure step you can also open `build\vs2026\wfmat.slnx` in Visual Studio and work
with it as a normal solution.

## Running wfmat-cli on the sample shapes

The executable is `build\vs2026\Release\wfmat-cli.exe` (`wfmat-cli --version` prints
`wfmat-cli 2.0.0`). Usage:

```text
wfmat-cli <input.json> [--svg out.svg] [--json out.json] [--flatten tol] [--offset d]...
wfmat-cli --version
```

- `--svg`: picture of the shape, the medial axis coloured by radius, and the offsets.
- `--json`: the medial axis transform, grouped by component (sites, vertices, edges as rational
  quadratic arcs).
- `--flatten tol`: also flatten edges to polylines within chord tolerance `tol`.
- `--offset d`: inward offset curve at distance `d`, drawn in the SVG and written to the JSON. May
  be repeated, e.g. `--offset 2 --offset 5`.

Distances and tolerances are in the input's units.

```powershell
mkdir out -Force
$cli = ".\build\vs2026\Release\wfmat-cli.exe"
& $cli tests\data\rectangle.json          --svg out\rectangle.svg          --json out\rectangle.mat.json          --flatten 0.1 --offset 5
& $cli tests\data\l-shape.json            --svg out\l-shape.svg            --json out\l-shape.mat.json            --flatten 0.1 --offset 5
& $cli tests\data\disk.json               --svg out\disk.svg               --json out\disk.mat.json               --flatten 0.1 --offset 5
& $cli tests\data\filleted-rectangle.json --svg out\filleted-rectangle.svg --json out\filleted-rectangle.mat.json --flatten 0.1 --offset 5
& $cli tests\data\square-ring.json        --svg out\square-ring.svg        --json out\square-ring.mat.json        --flatten 0.1 --offset 5
& $cli tests\data\plate-with-holes.json   --svg out\plate-with-holes.svg   --json out\plate-with-holes.mat.json   --flatten 0.1 --offset 3
```

Console output for the plate with holes (an outline, two holes, and an island inside one hole):

```text
tests\data\plate-with-holes.json: valid; 4 loops in 2 components; 13 segments (8 lines, 5 arcs) after merging 13 input segments
joins: 3 convex, 5 tangent, 5 reflex; area 4252.920307 mm^2
MAT: 24 vertices, 24 edges; events: 0 collapses, 5 transitions, 0 curvature ends, 0 splits, 2 contacts, 2 merges, 6 clusters, 6 loops vanished
offset 3: 3 loops
```

The other shapes:

| Shape | Loops / components | MAT | Notes |
| --- | --- | --- | --- |
| rectangle | 1 / 1 | 6 vertices, 5 edges | |
| l-shape | 1 / 1 | 10 vertices, 9 edges | 1 reflex corner |
| disk | 1 / 1 | 1 vertex, 0 edges | the axis of a circle is its centre |
| filleted-rectangle | 1 / 1 | 6 vertices, 5 edges | 1 curvature end |
| square-ring | 2 / 1 | 16 vertices, 16 edges | one cycle around the hole; offset 5 gives 2 loops |

Open the `.svg` files in a browser (drag them onto Edge or Chrome) to see the result.

More shapes in `tests\data`: `equilateral-triangle`, `regular-hexagon`, `concave-top`, `slot`,
`near-tangent-slit`, `fast-shock-triangle`, `spec-example`. To describe your own shapes, see
[input-format.md](input-format.md).
