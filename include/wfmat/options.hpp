// Run options and tolerances [IN-11], [API-02].
#pragma once

namespace wfmat {

// Tolerances in the normalised unit frame [IN-07].
struct Tolerances {
    double geom = 1e-10;   // coincidence of points and of site parameters
    double len = 1e-9;     // minimum segment length
    double ang = 1e-9;     // tangent and cusp detection, radians
    double time = 1e-11;   // event-time clustering
};

enum class BroadPhase { windowed_rtree, all_pairs };

struct Options {
    Tolerances tol;
    BroadPhase broad_phase = BroadPhase::windowed_rtree;
    bool debug_checks = false;  // full invariants after every event
    bool resolve_all_events = false;  // debug [RB-04]: send every event through cluster resolution
};

} // namespace wfmat
