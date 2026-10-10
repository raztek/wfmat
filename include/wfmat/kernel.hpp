// Geometric kernel [KN-01..KN-06]: sites, their distance functions, and the small solvers that
// give every event time and every shock trajectory. Pure functions of their arguments, no state.
//
// Unknowns are z = (x, y, t): a point and the time (= offset distance = radius) at which the
// front reaches it. Every site i has a distance function d_i, and the front of site i at time t
// is the level set d_i(x) = t.
#pragma once

#include <array>
#include <limits>
#include <optional>
#include <vector>

#include "wfmat/geom.hpp"
#include "wfmat/options.hpp"

namespace wfmat {

enum class SiteKind { line, circle };

// A site of the grassfire [section 2]. Points are concave circles of radius zero.
struct Site {
    SiteKind kind = SiteKind::line;
    Vec2 n;            // line: unit normal into the region
    double c = 0.0;    // line: offset, d(x) = n.x - c
    Vec2 centre;       // circle
    double R = 0.0;    // circle radius
    int sigma = 1;     // circle: +1 convex (centre on the region side), -1 concave

    static Site line(Vec2 n, double c) { return {SiteKind::line, n, c, {}, 0.0, 1}; }
    static Site circle(Vec2 centre, double R, int sigma) { return {SiteKind::circle, {}, 0.0, centre, R, sigma}; }
    static Site point(Vec2 p) { return circle(p, 0.0, -1); }

    bool is_line() const { return kind == SiteKind::line; }
    bool is_point() const { return kind == SiteKind::circle && R == 0.0; }

    double distance(Vec2 x) const;          // d(x), signed
    Vec2 gradient(Vec2 x) const;            // unit gradient of d, away from the site's front
    Vec2 foot(Vec2 x, double t) const;      // x - t grad d(x): the contact point on the site
    double offset_radius(double t) const;   // circle: radius of the front at time t, R - sigma t
};

// Site of a boundary segment of a counter-clockwise loop (the region on the left).
Site site_of(const Segment& s);

// [section 2] Conic traced by a shock between two sites.
enum class ConicKind { line, parabola, ellipse, hyperbola, plateau_line, plateau_arc };
ConicKind bisector_kind(const Site& a, const Site& b, double eps);

// A point equidistant (distance t) from three sites, with its contact points.
struct Root3 {
    Vec2 p;
    double t = 0.0;
    std::array<Vec2, 3> feet;
    double condition = 1.0;  // condition estimate of the solve; large near tangency [KN-06]
    bool refined_quad = false;
};

struct Solve3 {
    std::vector<Root3> roots;  // accepted roots, earliest first
    bool degenerate = false;   // rank-deficient system: parallel lines, concentric circles, ... [KN-02]
};

// [KN-01] Three-site solve. Accepts roots with t_now - eps_t <= t <= t_max, non-negative circle
// radii and a unit-frame residual |d_i(p) - t| <= eps_geom for all three sites. Whether each foot
// lies on the live part of its front element is the engine's check. The engine passes the largest
// possible radius as t_max for collapses and the end of its event window [EV-11] for splits, which
// spares the refinement of irrelevant far roots.
Solve3 solve_three(const Site& a, const Site& b, const Site& c, double t_now, const Tolerances& tol,
                   double t_max = std::numeric_limits<double>::infinity());

// [KN-02] The earliest accepted root, if any.
std::optional<Root3> earliest(const Solve3& s);

// [KN-03] Points of the two offset fronts at time t, i.e. the candidate shock positions.
std::vector<Vec2> offset_intersections(const Site& a, const Site& b, double t);

// [KN-03] Branch selector: the side of the pair's axis that p lies on (-1, 0 or +1). The axis runs
// through both centres for two circles and along the normal through the centre for a line and a
// circle; line pairs have no axis (0).
int branch_of(const Site& a, const Site& b, Vec2 p);

// [KN-03] Shock position at time t on the given branch.
std::optional<Vec2> shock_position(const Site& a, const Site& b, double t, int branch);

// Shock velocity dp/dt at p: grad d_a . v = grad d_b . v = 1. None when the gradients are parallel.
std::optional<Vec2> shock_velocity(const Site& a, const Site& b, Vec2 p);

// [KN-04] Time at which a regular vertex q(t) = q0 + t m (|m| = 1) meets the front of a site:
// d(q0 + t m) = t, after t_now - eps_t. Linear in t for lines and circles alike.
std::optional<double> regular_hit_time(Vec2 q0, Vec2 m, const Site& s, double t_now, const Tolerances& tol);

// [KN-05] First contact of two fronts along a common normal (an E3 candidate).
struct Contact {
    double t = 0.0;
    Vec2 p;                  // meeting point; plateau: the midline point nearest the origin (lines)
                             // or the common centre (arcs)
    std::array<Vec2, 2> feet;
    bool plateau = false;    // antiparallel lines or concentric circles: a whole curve meets at once
};
std::vector<Contact> contacts(const Site& a, const Site& b, double t_now, const Tolerances& tol);

// [KN-06] Newton refinement of a three-site root on the unsquared equations d_i(p) = t, in double
// or in binary128. Returns the refined root and updates its condition estimate.
Root3 refine(const std::array<Site, 3>& sites, Root3 root, bool quad);

} // namespace wfmat
