// Scaled shape families for the [VER-07] benchmarks: valid simple outlines with n segments.
#pragma once

#include <cmath>
#include <numbers>
#include <random>
#include <vector>

#include "wfmat/geom.hpp"
#include "wfmat/region.hpp"

namespace wfmat::bench {

inline constexpr double pi = std::numbers::pi;

// A star-shaped polygon with n vertices: evenly spaced directions with jitter, random radii.
inline Region star_polygon(int n, std::uint64_t seed = 1)
{
    std::mt19937_64 rng(seed);
    std::uniform_real_distribution<double> jitter(-0.3, 0.3), radius(0.5, 1.0);
    Region region;
    for (int k = 0; k < n; ++k) {
        const double a = 2 * pi * (k + jitter(rng)) / n;
        region.outer.vertices.push_back({radius(rng) * polar(a), 0.0});
    }
    return region;
}

// A wavy outline with n vertices: a smooth star r(a) = 1 + 0.25 sin 5a + 0.1 sin(13a + 1) with
// a little radial noise, so that most corners are nearly flat and a few are reflex.
inline Region wavy_polygon(int n, std::uint64_t seed = 1)
{
    std::mt19937_64 rng(seed);
    std::uniform_real_distribution<double> noise(-0.002, 0.002);
    Region region;
    for (int k = 0; k < n; ++k) {
        const double a = 2 * pi * k / n;
        const double r = 1.0 + 0.25 * std::sin(5 * a) + 0.1 * std::sin(13 * a + 1) + noise(rng);
        region.outer.vertices.push_back({r * polar(a), 0.0});
    }
    return region;
}

// The star polygon with every corner rounded by a tangent arc: 2n segments.
inline Region filleted_star(int n, std::uint64_t seed = 1)
{
    const Region star = star_polygon(n, seed);
    std::mt19937_64 rng(seed + 1);
    std::uniform_real_distribution<double> f(0.05, 0.9);
    Region region;
    const auto& v = star.outer.vertices;
    for (int k = 0; k < n; ++k) {
        const Vec2 prev = v[(k + n - 1) % n].p, p = v[k].p, next = v[(k + 1) % n].p;
        const Vec2 u = unit(p - prev), w = unit(next - p);
        const double d = f(rng) * 0.5 * std::min(dist(prev, p), dist(p, next));
        region.outer.vertices.push_back({p - d * u, bulge_from_sweep(signed_angle(u, w))});
        region.outer.vertices.push_back({p + d * w, 0.0});
    }
    return region;
}

// A gear with n trapezoidal teeth: 4n segments.
inline Region gear(int n)
{
    Region region;
    const double step = 2 * pi / n;
    for (int k = 0; k < n; ++k) {
        const double c = k * step;
        region.outer.vertices.push_back({polar(c - 0.35 * step), 0.0});
        region.outer.vertices.push_back({1.3 * polar(c - 0.2 * step), 0.0});
        region.outer.vertices.push_back({1.3 * polar(c + 0.2 * step), 0.0});
        region.outer.vertices.push_back({polar(c + 0.35 * step), 0.0});
    }
    return region;
}

} // namespace wfmat::bench
