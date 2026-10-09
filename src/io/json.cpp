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

} // namespace wfmat::io
