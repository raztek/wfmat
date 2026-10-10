// wfmat-cli: validate a region given as bulge JSON, compute its MAT and inward offsets, and write
// them as JSON and SVG [API-03], [OUT-06]. Input the engine does not handle yet
// (ErrorCode::unsupported) is reported and rendered without its MAT.

#include <cstdio>
#include <cstdlib>
#include <vector>
#include <fstream>
#include <iostream>
#include <string>
#include <string_view>

#include "wfmat/io/json.hpp"
#include "wfmat/io/svg.hpp"
#include "wfmat/mat.hpp"
#include "wfmat/prepare.hpp"
#include "wfmat/version.hpp"

namespace {

int usage()
{
    std::cerr << "usage: wfmat-cli <input.json> [--svg <output.svg>] [--json <output.json>]\n"
                 "                 [--flatten <chord tolerance>] [--offset <distance>]...\n"
                 "       wfmat-cli --version\n"
                 "  --svg      the boundary, the MAT coloured by radius and the offsets\n"
                 "  --json     sites, vertices, edges as rational quadratic Bezier arcs, offsets, statistics\n"
                 "  --flatten  add each edge as a polyline within this tolerance to the JSON\n"
                 "  --offset   an inward offset at this distance; repeat for several\n"
                 "Distances and tolerances are in the input's units.\n";
    return 2;
}

bool number(const char* text, double& out)
{
    char* end = nullptr;
    out = std::strtod(text, &end);
    return end != text && *end == '\0';
}

bool write_file(const std::string& path, const std::string& text)
{
    std::ofstream out(path, std::ios::binary);
    out << text;
    if (!out) std::cerr << "cannot write " << path << "\n";
    return static_cast<bool>(out);
}

} // namespace

int main(int argc, char** argv)
{
    std::string input, svg_path, json_path;
    double chord_tol = 0.0;
    std::vector<double> distances;
    for (int i = 1; i < argc; ++i) {
        const std::string_view arg = argv[i];
        double value = 0.0;
        if (arg == "--version") {
            std::printf("wfmat-cli %s\n", wfmat::version);
            return 0;
        }
        if (arg == "--svg" && i + 1 < argc)
            svg_path = argv[++i];
        else if (arg == "--json" && i + 1 < argc)
            json_path = argv[++i];
        else if (arg == "--flatten" && i + 1 < argc && number(argv[i + 1], value) && value > 0.0) {
            chord_tol = value;
            ++i;
        } else if (arg == "--offset" && i + 1 < argc && number(argv[i + 1], value) && value >= 0.0) {
            distances.push_back(value);
            ++i;
        } else if (!arg.empty() && arg[0] != '-' && input.empty())
            input = arg;
        else
            return usage();
    }
    if (input.empty()) return usage();

    auto doc = wfmat::io::read_region_json(input);
    if (!doc) {
        std::cerr << input << ": " << wfmat::to_string(doc.error()) << "\n";
        return 1;
    }
    auto region = wfmat::prepare(doc->region);
    if (!region) {
        std::cerr << input << ": " << wfmat::to_string(region.error()) << "\n";
        return 1;
    }

    std::size_t arcs = 0, convex = 0, tangent = 0, reflex = 0;
    for (const auto& s : region->segments) arcs += s.is_arc() ? 1 : 0;
    for (const auto& j : region->joins) {
        convex += j.kind == wfmat::JoinKind::convex;
        tangent += j.kind == wfmat::JoinKind::tangent;
        reflex += j.kind == wfmat::JoinKind::reflex;
    }
    const double scale = region->xf.scale;
    std::printf("%s: valid; %zu segments (%zu lines, %zu arcs) after merging %zu input segments\n", input.c_str(),
                region->segments.size(), region->segments.size() - arcs, arcs, doc->region.outer.vertices.size());
    std::printf("joins: %zu convex, %zu tangent, %zu reflex; area %.10g%s%s%s; input orientation %s\n", convex, tangent,
                reflex, region->area * scale * scale, doc->units.empty() ? "" : " ", doc->units.c_str(), doc->units.empty() ? "" : "^2",
                region->reversed ? "clockwise (reversed)" : "counter-clockwise");

    auto mat = wfmat::compute_mat(*region);
    if (mat) {
        const auto& st = mat->stats();
        std::printf("MAT: %zu vertices, %zu edges; events: %zu collapses, %zu transitions, %zu curvature ends, "
                    "%zu splits, %zu contacts, %zu clusters, %zu loops vanished\n",
                    mat->vertices().size(), mat->edges().size(), st.collapses, st.transitions, st.curvature_ends,
                    st.splits, st.contacts, st.clusters, st.annihilations);
    } else if (mat.error().code == wfmat::ErrorCode::unsupported) {
        std::printf("MAT: not computed yet: %s\n", wfmat::to_string(mat.error()).c_str());
    } else {
        std::cerr << input << ": " << wfmat::to_string(mat.error()) << "\n";
        return 1;
    }

    std::vector<wfmat::io::OffsetLoops> offsets;
    for (double d : distances) {
        auto loops = wfmat::inward_offset(*region, d);
        if (!loops) {
            std::cerr << input << ": offset " << d << ": " << wfmat::to_string(loops.error()) << "\n";
            return 1;
        }
        std::printf("offset %g: %zu loop%s\n", d, loops->size(), loops->size() == 1 ? "" : "s");
        offsets.push_back({d, std::move(*loops)});
    }

    if (!json_path.empty()) {
        if (!mat) {
            std::cerr << input << ": no MAT to write as JSON\n";
            return 1;
        }
        if (!write_file(json_path, wfmat::io::write_mat_json(*mat, offsets, {doc->units, chord_tol}))) return 1;
    }
    if (!svg_path.empty()) {
        wfmat::io::SvgWriter w(wfmat::io::caller_bbox(*region));
        w.boundary(*region);
        for (const auto& o : offsets) w.loops(o.loops);
        if (mat) w.medial_axis(*mat);
        w.joins(*region);
        if (!write_file(svg_path, w.str())) return 1;
    }
    return 0;
}
