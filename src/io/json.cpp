#include "wfmat/io/json.hpp"

#include <fstream>
#include <sstream>

#include <nlohmann/json.hpp>

namespace wfmat::io {

namespace {

using nlohmann::json;

auto schema_error(std::string message, std::optional<std::size_t> segment = std::nullopt)
{
    return make_error(ErrorCode::invalid_input, "IN-02", std::move(message), 0, segment);
}

} // namespace

Result<RegionDocument> parse_region_json(std::string_view text)
{
    const json doc = json::parse(text.begin(), text.end(), nullptr, false);
    if (doc.is_discarded()) return make_error(ErrorCode::invalid_input, "IN-02", "not valid JSON");
    if (!doc.is_object()) return make_error(ErrorCode::invalid_input, "IN-02", "the document must be a JSON object");

    RegionDocument out;
    for (const auto& [key, value] : doc.items()) {
        if (key == "outer") continue;
        if (key == "units") {
            if (!value.is_string()) return make_error(ErrorCode::invalid_input, "IN-02", "\"units\" must be a string");
            out.units = value.get<std::string>();
        } else if (key == "holes") {
            if (!value.is_array() || !value.empty())
                return make_error(ErrorCode::invalid_input, "IN-05",
                                  "regions with holes are not supported yet; holes arrive in v2", 1);
        } else {
            return make_error(ErrorCode::invalid_input, "IN-02", "unknown key \"" + key + "\"");
        }
    }

    const auto outer = doc.find("outer");
    if (outer == doc.end() || !outer->is_array())
        return make_error(ErrorCode::invalid_input, "IN-02", "\"outer\" must be an array of vertices");
    if (!outer->empty() && outer->front().is_array() && !outer->front().empty() && outer->front().front().is_array())
        return make_error(ErrorCode::invalid_input, "IN-05",
                          "\"outer\" must be a single loop; regions with holes arrive in v2");

    for (std::size_t i = 0; i < outer->size(); ++i) {
        const json& v = (*outer)[i];
        if (!v.is_array() || v.size() < 2 || v.size() > 3)
            return schema_error("vertex " + std::to_string(i) + " must be [x, y] or [x, y, bulge]", i);
        for (const json& c : v)
            if (!c.is_number()) return schema_error("vertex " + std::to_string(i) + " has a non-numeric entry", i);
        out.region.outer.vertices.push_back(
            {{v[0].get<double>(), v[1].get<double>()}, v.size() == 3 ? v[2].get<double>() : 0.0});
    }
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
    json outer = json::array();
    for (const BulgeVertex& v : doc.region.outer.vertices) outer.push_back({v.p.x, v.p.y, v.bulge});
    json j = json::object();
    if (!doc.units.empty()) j["units"] = doc.units;
    j["outer"] = std::move(outer);
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
        ojson j = {{"segment", s.segment}, {"corner", s.corner}};
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
            {{"p", point(v.p)}, {"r", v.r}, {"kind", kind_name(v.kind)}, {"edges", v.edges}, {"contacts", contacts}});
    }

    ojson edges = ojson::array();
    std::vector<Vec2> poly;
    for (const MatEdge& e : mat.edges()) {
        ojson bezier = ojson::array();
        for (const RationalQuadBezier& b : e.bezier())
            bezier.push_back(
                {{"p0", point(b.p0)}, {"p1", point(b.p1)}, {"p2", point(b.p2)}, {"w", b.w}, {"u", {b.u0, b.u1}}});
        ojson j = {{"v0", e.v0()},   {"v1", e.v1()},   {"left", e.left()}, {"right", e.right()},
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
        ojson loops = ojson::array();
        for (const Loop& loop : o.loops) {
            ojson vs = ojson::array();
            for (const BulgeVertex& v : loop.vertices) vs.push_back({v.p.x, v.p.y, v.bulge});
            loops.push_back(std::move(vs));
        }
        offs.push_back({{"distance", o.distance}, {"loops", std::move(loops)}});
    }

    const RunStats& st = mat.stats();
    const ojson stats = {{"collapses", st.collapses},       {"transitions", st.transitions},
                        {"curvature_ends", st.curvature_ends}, {"splits", st.splits},
                        {"contacts", st.contacts},         {"annihilations", st.annihilations},
                        {"clusters", st.clusters},         {"stale_events", st.stale_events},
                        {"windows", st.windows},           {"max_residual", st.max_residual}};

    ojson doc = ojson::object();
    if (!options.units.empty()) doc["units"] = options.units;
    doc["sites"] = std::move(sites);
    doc["vertices"] = std::move(vertices);
    doc["edges"] = std::move(edges);
    if (!offsets.empty()) doc["offsets"] = std::move(offs);
    doc["stats"] = stats;
    return doc.dump() + "\n";
}

} // namespace wfmat::io
