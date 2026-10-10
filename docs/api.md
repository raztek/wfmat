# wfmat API

wfmat computes the interior medial axis transform (MAT) of a planar region bounded by one closed loop
of line segments and circular arcs, by simulating the inward grassfire front. This page documents the
public C++ API of version 1.0. The design and the requirement IDs cited below are in
[`spec.md`](spec.md).

- [Using the library](#using-the-library)
- [Input: regions in bulge form](#input-regions-in-bulge-form)
- [Computing the MAT](#computing-the-mat)
- [The medial axis](#the-medial-axis)
- [Edge geometry](#edge-geometry)
- [Inward offsets](#inward-offsets)
- [Errors](#errors)
- [Options and tolerances](#options-and-tolerances)
- [JSON and SVG](#json-and-svg)
- [Command-line tool](#command-line-tool)

## Using the library

wfmat is a C++20 CMake project. Add it as a subdirectory, or fetch it, and link the targets you need:

```cmake
include(FetchContent)
FetchContent_Declare(wfmat GIT_REPOSITORY https://github.com/raztek/wfmat.git GIT_TAG v1.0.0)
FetchContent_MakeAvailable(wfmat)
target_link_libraries(app PRIVATE wfmat::wfmat)       # the algorithm
target_link_libraries(app PRIVATE wfmat::io)          # JSON and SVG (needs nlohmann/json)
```

`#include "wfmat/wfmat.hpp"` brings in the whole core API; `wfmat/io/json.hpp` and
`wfmat/io/svg.hpp` hold the file formats. Everything is in namespace `wfmat`, the formats in
`wfmat::io`. `wfmat::version` is the library version as a string.

The library keeps no global state: runs on different regions may proceed on different threads. A
`MedialAxis` is an immutable value once returned.

A complete program is in [`examples/basic.cpp`](../examples/basic.cpp), built as `wfmat-example`.

## Input: regions in bulge form

```cpp
struct BulgeVertex { Vec2 p; double bulge = 0.0; };
struct Loop        { std::vector<BulgeVertex> vertices; };
struct Region      { Loop outer; };                  // holes arrive in v2
```

Vertex `i` starts segment `i`, which ends at vertex `i + 1` (the last segment closes the loop). The
segment is a line when `bulge == 0`, otherwise a circular arc with `bulge = tan(sweep / 4)`: 1 is a
counter-clockwise half circle, -1 a clockwise one [IN-01]. A full circle is two half arcs. The loop
may run either way round; it must be simple, with no segment shorter than `tol.len` and no cusp
corners [IN-03..IN-06].

`loop_from_segments` builds a loop from explicit `LineSpec{a, b}` and
`ArcSpec{centre, radius, start_angle, sweep}` pieces, checking that consecutive pieces meet to
within a gap you give.

Coordinates are in your units. Internally the region is mapped to a unit frame (bounding box
centred, half-diagonal 1) where every tolerance applies [IN-07]; every output is mapped back, so
you never see the unit frame.

## Computing the MAT

```cpp
Result<MedialAxis> compute_mat(const Region& region, const Options& options = {});

Result<PreparedRegion> prepare(const Region& region, const Options& options = {});
Result<MedialAxis>     compute_mat(const PreparedRegion& region, const Options& options = {});
```

`compute_mat(region)` validates and normalises the region, then runs the propagation. Call
`prepare` yourself to validate once and reuse the result for several runs or offsets, or to
inspect the normalised loop: its segments in the unit frame, the joins classified as convex, tangent
or reflex, and the transform `xf` to and from your units.

## The medial axis

```cpp
struct MatVertex { Vec2 p; double r; VertexKind kind;
                   std::vector<MatContact> contacts;  // site and foot point of each contact
                   std::vector<EdgeId> edges; };      // counter-clockwise around p

class MedialAxis {
    std::span<const MatVertex> vertices() const;
    std::span<const MatEdge>   edges() const;
    std::span<const MatSite>   sites() const;
    const RunStats&            stats() const;
};
```

The MAT is a tree embedded in the region [OUT-08]. Each vertex has its position, the radius of
its inscribed disk and its kind:

| Kind | Degree | Where |
| --- | --- | --- |
| `corner` | 1 | a convex corner of the boundary, `r = 0` |
| `curvature_end` | 1 | the centre of a convex arc, `r` = the arc's radius |
| `junction` | 3 or more | the disk touches three or more sites |
| `transition` | 2 | the contact passes from a site to its neighbour, so the edge's conic type changes |
| `extremum_min` | 2 | two fronts first touch: a local minimum of `r` along the axis |
| `extremum_max` | 2 | two shocks meet: a local maximum of `r` |

A site is a boundary piece seen as a source of distance: a line, an arc, or the point of a reflex
corner. `MatSite::geometry` holds it in your units, with the signed distance
`Site::distance(x)`; `segment` is the input segment it comes from, and `corner` marks a reflex
corner point. Each edge names its left and right sites, seen along its direction of travel.

`RunStats` counts the events by type, the clusters of simultaneous events resolved, and the largest
distance residual at an event point.

## Edge geometry

```cpp
class MatEdge {
    VertexId v0() const, v1() const;        // birth and death vertex
    SiteId left() const, right() const;
    ConicKind kind() const;                 // line, parabola, ellipse, hyperbola,
                                            // plateau_line, plateau_arc
    double r0() const, r1() const;          // r0 <= r1

    Vec2 point_at(double u) const;          // u in [0, 1], from v0 to v1
    double radius_at(double u) const;
    Vec2 tangent_at(double u) const;        // unit, along the direction of travel
    std::pair<Vec2, Vec2> feet_at(double u) const;  // contact points on the left and right sites
    Vec2 point(double r) const;             // by radius; not for plateau edges

    std::vector<RationalQuadBezier> bezier() const;
    void flatten(double chord_tol, std::vector<Vec2>& out) const;
};
```

Every edge is the path of one shock of the grassfire, so it is a conic arc parametrised by the
radius: `u` maps linearly onto `[r0, r1]` [OUT-02]. A plateau edge, where two parallel lines or two
concentric arcs meet all along at once, has `r0 == r1`; there `u` runs along its length instead.
Evaluation is in closed form from the two sites [OUT-04], so any `u` is exact to rounding.

`bezier()` exports the edge exactly as rational quadratic arcs [OUT-05]:

```cpp
struct RationalQuadBezier {
    Vec2 p0, p1, p2;  double w;             // w < 1 ellipse, 1 parabola or line, > 1 hyperbola
    double u0, u1;                          // the edge parameters it covers
    Vec2 point(double s) const;             // s in [0, 1]
};
```

The arcs run end to end from `v0` to `v1`, each turning by less than a right angle, so every
weight is positive and the arcs map directly onto CAD and PostScript-style rational curves. The
radius is not polynomial in the Bézier parameter; evaluate it with `radius_at` over `[u0, u1]`.

`flatten(chord_tol, out)` appends a polyline from `v0` to `v1` whose chords stay within `chord_tol`
(your units) of the edge. Lines and plateau lines give their two end points.

## Inward offsets

```cpp
Result<std::vector<Loop>> inward_offset(const Region& region, double distance, const Options& options = {});
Result<std::vector<Loop>> inward_offset(const PreparedRegion& region, double distance, const Options& options = {});
```

The inward offset at `distance` is the grassfire front at that time [ALG-06]: lines stay lines and
arcs stay arcs, so the result is exact, in the same bulge form as the input, counter-clockwise. A
narrow neck splits the offset into several loops; at `distance = 0` it is the boundary, and past
the largest inscribed radius the list is empty. A negative or non-finite distance is an
`invalid_input` error.

## Errors

Every entry point returns `Result<T>`, which is `tl::expected<T, Error>` [API-03]:

```cpp
struct Error {
    ErrorCode code;              // invalid_input, numerical_failure, invariant_violation, unsupported
    std::string requirement;     // the spec requirement that failed, e.g. "IN-04"
    std::string message;
    std::optional<std::size_t> loop, segment;  // segment numbered as in your vertex list
    std::string dump_path;
};
std::string to_string(const Error&);           // "IN-04: loop 0, segment 3: ..."
```

`invalid_input` means the region breaks the input contract (not simple, a segment too short, a
cusp, holes). `numerical_failure` and `invariant_violation` mean the run could not be completed
consistently; they are bugs or inputs at the edge of double precision, and worth reporting with the
input. `unsupported` is reserved for valid input that a later version handles; no input gets it
in 1.0.

## Options and tolerances

```cpp
struct Tolerances { double geom = 1e-10, len = 1e-9, ang = 1e-9, time = 1e-11; };
struct Options {
    Tolerances tol;                                       // in the unit frame
    BroadPhase broad_phase = BroadPhase::windowed_rtree;  // or all_pairs, a slow cross-check
    bool debug_checks = false;                            // full invariants after every event
    bool resolve_all_events = false;                      // every event through cluster resolution
};
```

The tolerances are relative to the region's size, because they apply in the unit frame [IN-11]:
`geom` for coincident points (events within ten times this distance form one cluster), `len` for
the shortest segment, `ang` for tangent joins and directions, `time` for simultaneous events. The
defaults suit double-precision input; there is rarely a reason to change them. The other options
are for testing: `all_pairs` must give the same axis as the default, bit for bit [EV-12].

## JSON and SVG

Input is the bulge JSON of [IN-02]:

```json
{ "units": "mm", "outer": [[0, 0], [100, 0, 1], [100, 40], [0, 40]] }
```

`io::read_region_json(path)` and `io::parse_region_json(text)` read it; `io::write_region_json`
writes it back bit for bit.

`io::write_mat_json(mat, offsets, {units, chord_tol})` writes the MAT [OUT-06]: the sites, the
vertices with their contacts, the edges with their conic kind, radii and Bézier arcs (and a
polyline each when `chord_tol > 0`), any inward offsets as bulge loops, and the run statistics. The
layout is spelled out in `wfmat/io/json.hpp`.

`io::SvgWriter` draws for inspection, in your units with y up: `boundary`, `joins`, `medial_axis`
(edges coloured from blue to red by radius), `loops` (offsets, dashed), `disk`, `polyline` and
`point`. `io::render_svg(region, mat)` draws the usual picture in one call.

## Command-line tool

```text
wfmat-cli <input.json> [--svg <output.svg>] [--json <output.json>]
          [--flatten <chord tolerance>] [--offset <distance>]...
wfmat-cli --version
```

It validates the region, reports its joins and the MAT's size and event counts, and writes the JSON
and SVG outputs above. `--offset` may be repeated; distances and tolerances are in the input's units.
