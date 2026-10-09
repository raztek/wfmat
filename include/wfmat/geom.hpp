// Planar primitives: vectors, boxes and boundary segments (lines and circular arcs).
#pragma once

#include <cmath>
#include <cstdint>
#include <limits>
#include <numbers>
#include <optional>
#include <vector>

namespace wfmat {

struct Vec2 {
    double x = 0.0, y = 0.0;

    friend constexpr Vec2 operator+(Vec2 a, Vec2 b) { return {a.x + b.x, a.y + b.y}; }
    friend constexpr Vec2 operator-(Vec2 a, Vec2 b) { return {a.x - b.x, a.y - b.y}; }
    friend constexpr Vec2 operator-(Vec2 a) { return {-a.x, -a.y}; }
    friend constexpr Vec2 operator*(double s, Vec2 a) { return {s * a.x, s * a.y}; }
    friend constexpr Vec2 operator*(Vec2 a, double s) { return {s * a.x, s * a.y}; }
    friend constexpr Vec2 operator/(Vec2 a, double s) { return {a.x / s, a.y / s}; }
    friend constexpr bool operator==(Vec2, Vec2) = default;
};

constexpr double dot(Vec2 a, Vec2 b) { return a.x * b.x + a.y * b.y; }
constexpr double cross(Vec2 a, Vec2 b) { return a.x * b.y - a.y * b.x; }
constexpr Vec2 perp(Vec2 a) { return {-a.y, a.x}; }  // rotation by +90 degrees
inline double norm(Vec2 a) { return std::hypot(a.x, a.y); }
inline double dist(Vec2 a, Vec2 b) { return norm(a - b); }
inline Vec2 unit(Vec2 a) { return a / norm(a); }
inline Vec2 polar(double angle) { return {std::cos(angle), std::sin(angle)}; }
inline bool is_finite(Vec2 a) { return std::isfinite(a.x) && std::isfinite(a.y); }

// Signed angle from direction a to direction b, in (-pi, pi].
inline double signed_angle(Vec2 a, Vec2 b) { return std::atan2(cross(a, b), dot(a, b)); }

// Angle reduced to [0, 2 pi).
double wrap_2pi(double angle);

struct Box {
    Vec2 lo{std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity()};
    Vec2 hi{-std::numeric_limits<double>::infinity(), -std::numeric_limits<double>::infinity()};

    bool empty() const { return lo.x > hi.x || lo.y > hi.y; }
    void add(Vec2 p);
    void add(const Box& b);
    Box inflated(double d) const { return {lo - Vec2{d, d}, hi + Vec2{d, d}}; }
    Vec2 center() const { return 0.5 * (lo + hi); }
    Vec2 extent() const { return hi - lo; }
};

// A boundary segment from a to b: a line when bulge == 0, otherwise a circular arc with
// bulge = tan(sweep / 4) (the DXF LWPOLYLINE convention), sweep signed, positive counter-clockwise.
// Arc data (centre, radius, start angle, sweep) is derived once by make_segment.
struct Segment {
    Vec2 a, b;
    double bulge = 0.0;

    Vec2 c;             // arc centre (arcs only)
    double R = 0.0;     // arc radius (arcs only)
    double phi0 = 0.0;  // angle of a about c (arcs only)
    double sweep = 0.0; // signed sweep in (-2 pi, 2 pi), 0 for lines

    bool is_arc() const { return bulge != 0.0; }
    double length() const;
    Vec2 point(double s) const;          // s in [0, 1], arc-length proportional
    Vec2 start_tangent() const;          // unit direction of travel at a
    Vec2 end_tangent() const;            // unit direction of travel at b
    Vec2 midpoint() const { return point(0.5); }
    Box bbox() const;                    // tight, including arc extremes

    // True when the direction of p - c lies within the arc's angular range, widened by slack
    // radians at both ends. Arcs only.
    bool angle_in_range(Vec2 p, double slack = 0.0) const;

    // Closest point on the segment to p. For p at an arc centre every point is closest;
    // the arc midpoint is returned then.
    Vec2 closest_point(Vec2 p) const;
    double distance(Vec2 p) const { return dist(p, closest_point(p)); }

    // Signed area between the chord a->b and the arc (positive for counter-clockwise arcs).
    double cap_area() const;

    Segment reversed() const;
};

Segment make_segment(Vec2 a, Vec2 b, double bulge);

// Bulge of an arc with the given signed sweep.
inline double bulge_from_sweep(double sweep) { return std::tan(sweep / 4.0); }

// Points where two segments meet or come within eps of each other. Overlapping collinear or
// co-circular pieces (a common part longer than eps) set the overlap flag.
struct Intersections {
    std::vector<Vec2> points;
    bool overlap = false;
};
Intersections intersect(const Segment& s, const Segment& t, double eps);

} // namespace wfmat
