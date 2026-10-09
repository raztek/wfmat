#include "wfmat/region.hpp"

#include <cmath>
#include <numbers>
#include <string>

namespace wfmat {

namespace {

struct Endpoints {
    Vec2 start, end;
    double bulge;
};

Endpoints endpoints(const SegmentSpec& spec)
{
    if (const auto* line = std::get_if<LineSpec>(&spec)) return {line->a, line->b, 0.0};
    const auto& arc = std::get<ArcSpec>(spec);
    return {arc.centre + arc.radius * polar(arc.start_angle),
            arc.centre + arc.radius * polar(arc.start_angle + arc.sweep), bulge_from_sweep(arc.sweep)};
}

} // namespace

Result<Loop> loop_from_segments(std::span<const SegmentSpec> segments, double gap)
{
    Loop loop;
    const std::size_t n = segments.size();
    for (std::size_t i = 0; i < n; ++i) {
        if (const auto* arc = std::get_if<ArcSpec>(&segments[i])) {
            if (!(arc->radius > 0.0) || !(std::abs(arc->sweep) > 0.0 && std::abs(arc->sweep) < 2.0 * std::numbers::pi))
                return make_error(ErrorCode::invalid_input, "IN-01",
                                  "arc " + std::to_string(i) + " needs a positive radius and 0 < |sweep| < 2 pi", 0, i);
        }
        const Endpoints e = endpoints(segments[i]);
        const Endpoints next = endpoints(segments[(i + 1) % n]);
        if (dist(e.end, next.start) > gap)
            return make_error(ErrorCode::invalid_input, "IN-01",
                              "segment " + std::to_string(i) + " does not end where segment " +
                                  std::to_string((i + 1) % n) + " starts",
                              0, i);
        loop.vertices.push_back({e.start, e.bulge});
    }
    return loop;
}

std::string to_string(const Error& e)
{
    std::string s = e.requirement.empty() ? std::string("error") : e.requirement;
    if (e.loop) s += ": loop " + std::to_string(*e.loop);
    if (e.segment) s += (e.loop ? ", segment " : ": segment ") + std::to_string(*e.segment);
    s += ": " + e.message;
    if (!e.dump_path.empty()) s += " (diagnostic dump: " + e.dump_path + ")";
    return s;
}

} // namespace wfmat
