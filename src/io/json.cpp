#include "wfmat/io/json.hpp"

#include <fstream>
#include <sstream>

#include <nlohmann/json.hpp>

namespace wfmat::io {

namespace {

using nlohmann::json;

// One loop: an array of [x, y] or [x, y, bulge] vertices. what names it in diagnostics.
Result<Loop> parse_loop(const json& value, const std::string& what, std::size_t loop)
{
    if (!value.is_array()) return make_error(ErrorCode::invalid_input, "IN-02", what + " must be an array of vertices", loop);
    Loop out;
    for (std::size_t i = 0; i < value.size(); ++i) {
        const json& v = value[i];
        if (!v.is_array() || v.size() < 2 || v.size() > 3)
            return make_error(ErrorCode::invalid_input, "IN-02",
                              what + ": vertex " + std::to_string(i) + " must be [x, y] or [x, y, bulge]", loop, i);
        for (const json& c : v)
            if (!c.is_number())
                return make_error(ErrorCode::invalid_input, "IN-02",
                                  what + ": vertex " + std::to_string(i) + " has a non-numeric entry", loop, i);
        out.vertices.push_back({{v[0].get<double>(), v[1].get<double>()}, v.size() == 3 ? v[2].get<double>() : 0.0});
    }
    return out;
}

// A list of loops: an array of vertex arrays.
Result<void> parse_loops(const json& value, const std::string& key, std::vector<Loop>& out)
{
    if (!value.is_array()) return make_error(ErrorCode::invalid_input, "IN-02", "\"" + key + "\" must be an array of loops");
    for (std::size_t k = 0; k < value.size(); ++k) {
        auto loop = parse_loop(value[k], "\"" + key + "\" loop " + std::to_string(k), out.size());
        if (!loop) return tl::unexpected(loop.error());
        out.push_back(std::move(*loop));
    }
    return {};
}

} // namespace

Result<RegionDocument> parse_region_json(std::string_view text)
{
    const json doc = json::parse(text.begin(), text.end(), nullptr, false);
    if (doc.is_discarded()) return make_error(ErrorCode::invalid_input, "IN-02", "not valid JSON");
    if (!doc.is_object()) return make_error(ErrorCode::invalid_input, "IN-02", "the document must be a JSON object");

    RegionDocument out;
    for (const auto& [key, value] : doc.items()) {
        if (key == "outer" || key == "holes" || key == "loops") continue;
        if (key == "units") {
            if (!value.is_string()) return make_error(ErrorCode::invalid_input, "IN-02", "\"units\" must be a string");
            out.units = value.get<std::string>();
        } else {
            return make_error(ErrorCode::invalid_input, "IN-02", "unknown key \"" + key + "\"");
        }
    }

    const auto loops = doc.find("loops"), outer = doc.find("outer"), holes = doc.find("holes");
    std::vector<Loop>& list = out.region.loops;
    if (loops != doc.end()) {
        if (outer != doc.end() || holes != doc.end())
            return make_error(ErrorCode::invalid_input, "IN-02", "use either \"loops\" or \"outer\" with \"holes\", not both");
        if (auto r = parse_loops(*loops, "loops", list); !r) return tl::unexpected(r.error());
        if (list.empty()) return make_error(ErrorCode::invalid_input, "IN-02", "\"loops\" must hold at least one loop");
        return out;
    }
    if (outer == doc.end())
        return make_error(ErrorCode::invalid_input, "IN-02", "the document needs \"loops\" or \"outer\"");
    if (!outer->empty() && outer->front().is_array() && !outer->front().empty() && outer->front().front().is_array())
        return make_error(ErrorCode::invalid_input, "IN-02", "\"outer\" must be a single loop; use \"loops\" for several");
    auto first = parse_loop(*outer, "\"outer\"", 0);
    if (!first) return tl::unexpected(first.error());
    list.push_back(std::move(*first));
    if (holes != doc.end())
        if (auto r = parse_loops(*holes, "holes", list); !r) return tl::unexpected(r.error());
    return out;
}

Result<RegionDocument> read_region_json(const std::filesystem::path& path)
{
    std::ifstream in(path, std::ios::binary);
    if (!in) return make_error(ErrorCode::invalid_input, "IN-02", "cannot open " + path.string());
    std::ostringstream text;
    text << in.rdbuf();
    return parse_region_json(text.str());
}

std::string write_region_json(const RegionDocument& doc)
{
    json loops = json::array();
    for (const Loop& loop : doc.region.loops) {
        json vertices = json::array();
        for (const BulgeVertex& v : loop.vertices) vertices.push_back({v.p.x, v.p.y, v.bulge});
        loops.push_back(std::move(vertices));
    }
    json j = json::object();
    if (!doc.units.empty()) j["units"] = doc.units;
    j["loops"] = std::move(loops);
    return j.dump() + "\n";
}

namespace {

using ojson = nlohmann::ordered_json;

ojson point(Vec2 p)
{
    return ojson::array({p.x, p.y});
}

const char* kind_name(VertexKind k)
{
    switch (k) {
    case VertexKind::corner: return "corner";
    case VertexKind::curvature_end: return "curvature_end";
    case VertexKind::junction: return "junction";
    case VertexKind::transition: return "transition";
    case VertexKind::extremum_min: return "extremum_min";
    case VertexKind::extremum_max: return "extremum_max";
    }
    return "";
}

const char* kind_name(ConicKind k)
{
    switch (k) {
    case ConicKind::line: return "line";
    case ConicKind::parabola: return "parabola";
    case ConicKind::ellipse: return "ellipse";
    case ConicKind::hyperbola: return "hyperbola";
    case ConicKind::plateau_line: return "plateau_line";
    case ConicKind::plateau_arc: return "plateau_arc";
    }
    return "";
}

} // namespace

std::string write_mat_json(const MedialAxis& mat, std::span<const OffsetLoops> offsets, const MatJsonOptions& options)
{
    ojson sites = ojson::array();
    for (const MatSite& s : mat.sites()) {
        ojson j = {{"loop", s.loop}, {"segment", s.segment}, {"corner", s.corner}};
        const Site& g = s.geometry;
        if (g.is_line()) {
            j["kind"] = "line";
            j["normal"] = point(g.n);
            j["offset"] = g.c;
        } else {
            j["kind"] = "circle";
            j["centre"] = point(g.centre);
            j["radius"] = g.R;
            j["sigma"] = g.sigma;
        }
        sites.push_back(std::move(j));
    }

    ojson vertices = ojson::array();
    for (const MatVertex& v : mat.vertices()) {
        ojson contacts = ojson::array();
        for (const MatContact& c : v.contacts) contacts.push_back({{"site", c.site}, {"foot", point(c.foot)}});
        vertices.push_back(
            {{"p", point(v.p)}, {"r", v.r}, {"kind", kind_name(v.kind)}, {"component", v.component}, {"edges", v.edges},
             {"contacts", contacts}});
    }

    ojson edges = ojson::array();
    std::vector<Vec2> poly;
    for (const MatEdge& e : mat.edges()) {
        ojson bezier = ojson::array();
        for (const RationalQuadBezier& b : e.bezier())
            bezier.push_back(
                {{"p0", point(b.p0)}, {"p1", point(b.p1)}, {"p2", point(b.p2)}, {"w", b.w}, {"u", {b.u0, b.u1}}});
        ojson j = {{"v0", e.v0()},   {"v1", e.v1()},   {"left", e.left()}, {"right", e.right()}, {"component", e.component()},
                  {"kind", kind_name(e.kind())}, {"r0", e.r0()}, {"r1", e.r1()}, {"bezier", std::move(bezier)}};
        if (options.chord_tol > 0.0) {
            poly.clear();
            e.flatten(options.chord_tol, poly);
            ojson pts = ojson::array();
            for (Vec2 p : poly) pts.push_back(point(p));
            j["polyline"] = std::move(pts);
        }
        edges.push_back(std::move(j));
    }

    ojson offs = ojson::array();
    for (const OffsetLoops& o : offsets) {
        ojson comps = ojson::array();
        for (const std::vector<Loop>& component : o.components) {
            ojson loops = ojson::array();
            for (const Loop& loop : component) {
                ojson vs = ojson::array();
                for (const BulgeVertex& v : loop.vertices) vs.push_back({v.p.x, v.p.y, v.bulge});
                loops.push_back(std::move(vs));
            }
            comps.push_back(std::move(loops));
        }
        offs.push_back({{"distance", o.distance}, {"components", std::move(comps)}});
    }

    const RunStats& st = mat.stats();
    const ojson stats = {{"collapses", st.collapses},       {"transitions", st.transitions},
                        {"curvature_ends", st.curvature_ends}, {"splits", st.splits},
                        {"contacts", st.contacts},         {"merges", st.merges},         {"annihilations", st.annihilations},
                        {"clusters", st.clusters},         {"stale_events", st.stale_events},
                        {"windows", st.windows},           {"max_residual", st.max_residual}};

    ojson comps = ojson::array();
    for (const MatComponent& c : mat.components())
        comps.push_back({{"outer", c.outer},
                         {"holes", c.holes},
                         {"vertices", {c.first_vertex, c.vertex_count}},
                         {"edges", {c.first_edge, c.edge_count}}});

    ojson doc = ojson::object();
    if (!options.units.empty()) doc["units"] = options.units;
    doc["components"] = std::move(comps);
    doc["sites"] = std::move(sites);
    doc["vertices"] = std::move(vertices);
    doc["edges"] = std::move(edges);
    if (!offsets.empty()) doc["offsets"] = std::move(offs);
    doc["stats"] = stats;
    return doc.dump() + "\n";
}

} // namespace wfmat::io
