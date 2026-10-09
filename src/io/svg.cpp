#include "wfmat/io/svg.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <numbers>

namespace wfmat::io {

namespace {

std::string num(double v)
{
    char buf[32];
    std::snprintf(buf, sizeof buf, "%.10g", v);
    return buf;
}

std::string xy(Vec2 p)
{
    return num(p.x) + "," + num(p.y);
}

const char* join_class(JoinKind k)
{
    switch (k) {
    case JoinKind::convex: return "convex";
    case JoinKind::tangent: return "tangent";
    case JoinKind::reflex: return "reflex";
    }
    return "";
}

constexpr const char* style = R"(<style>
path, polyline { fill: none; vector-effect: non-scaling-stroke; stroke-linejoin: round; stroke-linecap: round; }
.line { stroke: #1f2937; stroke-width: 2; }
.arc { stroke: #2563eb; stroke-width: 2; }
.curve { stroke: #dc2626; stroke-width: 1.5; }
.disk { fill: #f59e0b22; stroke: #f59e0b; stroke-width: 1; vector-effect: non-scaling-stroke; }
.point { fill: #dc2626; }
.convex { fill: #16a34a; }
.tangent { fill: #ffffff; stroke: #2563eb; stroke-width: 1.5; vector-effect: non-scaling-stroke; }
.reflex { fill: #9333ea; }
</style>
)";

} // namespace

SvgWriter::SvgWriter(const Box& view, SvgOptions options) : view_(view), options_(options)
{
    const Vec2 e = view.extent();
    const double pad = options.margin * std::max(e.x, e.y);
    view_ = view.inflated(pad);
    marker_ = 0.006 * std::max(view_.extent().x, view_.extent().y);
}

void SvgWriter::boundary(const PreparedRegion& region)
{
    const Transform& xf = region.xf;
    for (const Segment& s : region.segments) {
        const Vec2 a = xf.from_unit(s.a), b = xf.from_unit(s.b);
        if (!s.is_arc()) {
            body_ += "<path class=\"line\" d=\"M" + xy(a) + " L" + xy(b) + "\"/>\n";
            continue;
        }
        // Inside the y-flipped group, sweep-flag 1 is the counter-clockwise direction.
        const double r = xf.length_from_unit(s.R);
        const int large = std::abs(s.sweep) > std::numbers::pi ? 1 : 0;
        const int ccw = s.sweep > 0.0 ? 1 : 0;
        body_ += "<path class=\"arc\" d=\"M" + xy(a) + " A" + num(r) + "," + num(r) + " 0 " + std::to_string(large) +
                 " " + std::to_string(ccw) + " " + xy(b) + "\"/>\n";
    }
}

void SvgWriter::joins(const PreparedRegion& region)
{
    for (const Join& j : region.joins)
        body_ += "<circle class=\"" + std::string(join_class(j.kind)) + "\" cx=\"" + num(region.xf.from_unit(j.p).x) +
                 "\" cy=\"" + num(region.xf.from_unit(j.p).y) + "\" r=\"" + num(marker_) + "\"/>\n";
}

void SvgWriter::disk(Vec2 centre, double r, std::string_view css_class)
{
    body_ += "<circle class=\"" + std::string(css_class) + "\" cx=\"" + num(centre.x) + "\" cy=\"" + num(centre.y) +
             "\" r=\"" + num(r) + "\"/>\n";
}

void SvgWriter::polyline(std::span<const Vec2> points, std::string_view css_class)
{
    body_ += "<polyline class=\"" + std::string(css_class) + "\" points=\"";
    for (std::size_t i = 0; i < points.size(); ++i) body_ += (i ? " " : "") + xy(points[i]);
    body_ += "\"/>\n";
}

void SvgWriter::point(Vec2 p, std::string_view css_class)
{
    disk(p, marker_, css_class);
}

std::string SvgWriter::str() const
{
    const Vec2 e = view_.extent();
    const double height_px = options_.width_px * e.y / e.x;
    std::string out = "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"" + num(options_.width_px) + "\" height=\"" +
                      num(height_px) + "\" viewBox=\"" + num(view_.lo.x) + " " + num(-view_.hi.y) + " " + num(e.x) +
                      " " + num(e.y) + "\">\n";
    out += style;
    out += "<rect x=\"" + num(view_.lo.x) + "\" y=\"" + num(-view_.hi.y) + "\" width=\"" + num(e.x) + "\" height=\"" +
           num(e.y) + "\" fill=\"#ffffff\"/>\n";
    out += "<g transform=\"scale(1,-1)\">\n" + body_ + "</g>\n</svg>\n";
    return out;
}

Box caller_bbox(const PreparedRegion& region)
{
    Box b;
    b.add(region.xf.from_unit(region.bbox.lo));
    b.add(region.xf.from_unit(region.bbox.hi));
    return b;
}

std::string render_svg(const PreparedRegion& region, const SvgOptions& options)
{
    SvgWriter w(caller_bbox(region), options);
    w.boundary(region);
    if (options.show_joins) w.joins(region);
    return w.str();
}

} // namespace wfmat::io
