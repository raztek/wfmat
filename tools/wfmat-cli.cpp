// wfmat-cli: validate a region given as bulge JSON, compute its MAT and render both to SVG
// [API-03]. Input the engine does not handle yet (arcs before M3, clusters before M5) is reported
// and rendered without its MAT.

#include <cstdio>
#include <fstream>
#include <iostream>
#include <string>
#include <string_view>

#include "wfmat/io/json.hpp"
#include "wfmat/io/svg.hpp"
#include "wfmat/mat.hpp"
#include "wfmat/prepare.hpp"

namespace {

int usage()
{
    std::cerr << "usage: wfmat-cli <input.json> [--svg <output.svg>]\n";
    return 2;
}

} // namespace

int main(int argc, char** argv)
{
    std::string input, svg_path;
    for (int i = 1; i < argc; ++i) {
        const std::string_view arg = argv[i];
        if (arg == "--svg" && i + 1 < argc)
            svg_path = argv[++i];
        else if (!arg.empty() && arg[0] != '-' && input.empty())
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
        std::printf("MAT: %zu vertices, %zu edges; events: %zu collapses, %zu transitions, %zu splits, %zu contacts, "
                    "%zu loops vanished\n",
                    mat->vertices().size(), mat->edges().size(), st.collapses, st.transitions, st.splits, st.contacts,
                    st.annihilations);
    } else if (mat.error().code == wfmat::ErrorCode::unsupported) {
        std::printf("MAT: not computed yet: %s\n", wfmat::to_string(mat.error()).c_str());
    } else {
        std::cerr << input << ": " << wfmat::to_string(mat.error()) << "\n";
        return 1;
    }

    if (!svg_path.empty()) {
        std::ofstream out(svg_path, std::ios::binary);
        out << (mat ? wfmat::io::render_svg(*region, *mat) : wfmat::io::render_svg(*region));
        if (!out) {
            std::cerr << "cannot write " << svg_path << "\n";
            return 1;
        }
    }
    return 0;
}
