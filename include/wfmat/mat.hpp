// The medial axis transform as an embedded graph [OUT-01..OUT-05] and its entry point [API-02].
//
// Every edge is the trajectory of one shock of the grassfire over a radius interval, so it is
// parametrised by the radius r, which increases from the birth vertex v0 to the death vertex v1.
// All coordinates and radii are in the caller's units; the unit frame of [IN-07] is internal.
#pragma once

#include <cstddef>
#include <cstdint>
#include <span>
#include <utility>
#include <vector>

#include "wfmat/geom.hpp"
#include "wfmat/kernel.hpp"
#include "wfmat/options.hpp"
#include "wfmat/prepare.hpp"
#include "wfmat/region.hpp"
#include "wfmat/result.hpp"

namespace wfmat {

using SiteId = std::uint32_t;
using VertexId = std::uint32_t;
using EdgeId = std::uint32_t;
inline constexpr std::uint32_t no_id = 0xffffffffu;

// [section 2] MAT vertex kinds.
enum class VertexKind { corner, curvature_end, junction, transition, extremum_min, extremum_max };

// [OUT-03] A site in caller units and the boundary piece it comes from.
struct MatSite {
    Site geometry;
    std::uint32_t segment = 0;  // input segment index; a corner point names the segment ending
                                // at the corner in counter-clockwise order
    bool corner = false;        // a point site at a reflex corner
};

// A contact of the inscribed disk at a vertex: the site and the foot point on it.
struct MatContact {
    SiteId site = 0;
    Vec2 foot;
};

// [OUT-01]
struct MatVertex {
    Vec2 p;
    double r = 0.0;
    VertexKind kind = VertexKind::junction;
    std::vector<MatContact> contacts;
    std::vector<EdgeId> edges;  // counter-clockwise around p
};

// [OUT-02], [OUT-04], [OUT-05] One MAT edge. u in [0, 1] maps linearly onto [r0, r1]; on a plateau
// edge (r0 = r1) it maps linearly onto arc length instead.
class MatEdge {
public:
    // Construction data, filled in by the engine. Sites, points and radii are in the unit frame.
    struct Data {
        VertexId v0 = no_id, v1 = no_id;
        SiteId left = 0, right = 0;  // seen along the direction of travel
        ConicKind kind = ConicKind::line;
        int branch = 0;              // [KN-03] branch selector fixed at birth; plateau_arc: +1 when
                                     // the edge runs counter-clockwise about the common centre
        Site a, b;                   // left and right sites
        double t0 = 0.0, t1 = 0.0;   // radius interval
        Vec2 p0, p1;                 // end points
        Transform xf;
    };

    MatEdge() = default;
    explicit MatEdge(const Data& data) : d_(data) {}

    VertexId v0() const { return d_.v0; }
    VertexId v1() const { return d_.v1; }
    SiteId left() const { return d_.left; }
    SiteId right() const { return d_.right; }
    ConicKind kind() const { return d_.kind; }
    double r0() const { return d_.xf.length_from_unit(d_.t0); }
    double r1() const { return d_.xf.length_from_unit(d_.t1); }

    bool is_plateau() const { return d_.kind == ConicKind::plateau_line || d_.kind == ConicKind::plateau_arc; }
    Vec2 point(double r) const;                     // shock position at radius r in [r0, r1]; not
                                                    // for plateau edges (gives the start point)
    Vec2 point_at(double u) const;
    double radius_at(double u) const;
    Vec2 tangent_at(double u) const;                // unit, along the direction of travel
    std::pair<Vec2, Vec2> feet_at(double u) const;  // contact points on the left and right sites

private:
    double t_at(double u) const { return d_.t0 + u * (d_.t1 - d_.t0); }
    Vec2 unit_point(double t) const;
    Vec2 plateau_point(double u) const;
    Vec2 unit_point_at(double u) const { return is_plateau() ? plateau_point(u) : unit_point(t_at(u)); }

    Data d_;
};

// [OUT-03] Per-run statistics.
struct RunStats {
    std::size_t collapses = 0;          // E1a, shock + shock [EV-01]
    std::size_t transitions = 0;        // E1b, shock + regular [EV-02]
    std::size_t curvature_ends = 0;     // E1c, a convex arc shrinks to its centre [EV-03]
    std::size_t splits = 0;             // E2 [EV-04]
    std::size_t contacts = 0;           // E3 [EV-05]
    std::size_t annihilations = 0;      // loops that vanished [EV-08]
    std::size_t clusters = 0;           // resolved clusters [RB-03]
    std::size_t stale_events = 0;       // discarded at pop time [EV-10]
    double max_residual = 0.0;          // largest |d_i(p) - t| over event points, caller units
};

class MedialAxis {
public:
    MedialAxis() = default;
    MedialAxis(std::vector<MatVertex> vertices, std::vector<MatEdge> edges, std::vector<MatSite> sites,
               RunStats stats, Transform xf)
        : vertices_(std::move(vertices)), edges_(std::move(edges)), sites_(std::move(sites)), stats_(stats),
          xf_(xf)
    {
    }

    std::span<const MatVertex> vertices() const { return vertices_; }
    std::span<const MatEdge> edges() const { return edges_; }
    std::span<const MatSite> sites() const { return sites_; }
    const RunStats& stats() const { return stats_; }
    const Transform& transform() const { return xf_; }

private:
    std::vector<MatVertex> vertices_;
    std::vector<MatEdge> edges_;
    std::vector<MatSite> sites_;
    RunStats stats_;
    Transform xf_;
};

// [API-02] The interior MAT of a region bounded by lines and arcs (polygons since M2, arcs since
// M3, simultaneous events and plateaus since M5).
Result<MedialAxis> compute_mat(const Region& region, const Options& options = {});
Result<MedialAxis> compute_mat(const PreparedRegion& region, const Options& options = {});

} // namespace wfmat
