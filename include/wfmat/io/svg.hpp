// SVG output for inspection [OUT-06]. Drawing coordinates are in the caller's units, y up.
#pragma once

#include <span>
#include <string>
#include <string_view>

#include "wfmat/geom.hpp"
#include "wfmat/mat.hpp"
#include "wfmat/prepare.hpp"

namespace wfmat::io {

struct SvgOptions {
    double width_px = 800.0;   // rendered width; height follows the aspect ratio
    double margin = 0.05;      // fraction of the larger extent added on every side
    bool show_joins = true;    // mark joins by kind: convex, tangent, reflex
};

class SvgWriter {
public:
    explicit SvgWriter(const Box& view, SvgOptions options = {});

    void boundary(const PreparedRegion& region);   // lines and arcs, mapped back from the unit frame
    void joins(const PreparedRegion& region);
    void medial_axis(const MedialAxis& mat, int samples_per_edge = 32);  // edges as polylines
    void disk(Vec2 centre, double r, std::string_view css_class = "disk");
    void polyline(std::span<const Vec2> points, std::string_view css_class = "curve");
    void point(Vec2 p, std::string_view css_class = "point");

    std::string str() const;

private:
    Box view_;
    SvgOptions options_;
    double marker_ = 0.0;
    std::string body_;
};

// Boundary and joins of a prepared region.
std::string render_svg(const PreparedRegion& region, const SvgOptions& options = {});

// Boundary, joins and the medial axis.
std::string render_svg(const PreparedRegion& region, const MedialAxis& mat, const SvgOptions& options = {});

// Caller-frame bounding box of a prepared region.
Box caller_bbox(const PreparedRegion& region);

} // namespace wfmat::io
