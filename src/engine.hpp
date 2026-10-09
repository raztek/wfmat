// The propagation engine [ALG-02..ALG-05], [EV-01..EV-13]: the front, the event queue and the
// handlers. Internal to the library; compute_mat in src/mat.cpp is the public entry point.
#pragma once

#include <array>
#include <cstdint>
#include <optional>
#include <queue>
#include <vector>

#include "wfmat/mat.hpp"

namespace wfmat::detail {

// [ALG-02] Front vertices and elements, addressed by stable 32-bit ids that are never reused.
enum class VertexType : std::uint8_t { shock, regular };

struct FrontVertex {
    VertexType type = VertexType::shock;
    std::uint32_t left = no_id, right = no_id;  // the elements before and after it
    double t0 = 0.0;                            // birth time
    Vec2 p0;                                    // birth point; regular vertices: the corner
    Vec2 m;                                     // regular vertices: q(t) = p0 + t m [KN-04]
    int branch = 0;                             // shocks: [KN-03] branch selector
    EdgeId edge = no_id;                        // shocks: the MAT edge being traced
    bool alive = true;
};

struct FrontElement {
    SiteId site = 0;
    std::uint32_t prev = no_id, next = no_id;  // vertices
    std::uint32_t loop = 0;
    std::uint32_t version = 0;                 // bumped whenever either vertex changes [EV-10]
    bool alive = true;
};

enum class EventKind : std::uint8_t { collapse, split, contact, plateau };

struct Event {
    double t = 0.0;
    EventKind kind = EventKind::collapse;
    std::uint32_t key = 0;   // smallest involved element id [RB-01]
    std::uint64_t seq = 0;   // insertion order: the final, deterministic tie-break
    Vec2 p;
    std::array<std::uint32_t, 2> elements{no_id, no_id};
    std::array<std::uint32_t, 2> versions{0, 0};
    std::uint32_t vertex = no_id;  // split: the shock that reaches elements[0]
};

// [RB-01] Earliest time first, then kind, then the smallest id.
struct Later {
    bool operator()(const Event& a, const Event& b) const
    {
        if (a.t != b.t) return a.t > b.t;
        if (a.kind != b.kind) return a.kind > b.kind;
        if (a.key != b.key) return a.key > b.key;
        return a.seq > b.seq;
    }
};

class Engine {
public:
    Engine(const PreparedRegion& region, const Options& options);

    Result<void> run();

    // Results, in the unit frame.
    struct OutVertex {
        Vec2 p;
        double t = 0.0;
        VertexKind kind = VertexKind::junction;
        std::vector<MatContact> contacts;
    };
    std::vector<Site> sites;
    std::vector<MatSite> site_info;  // provenance; geometry is filled in by the caller
    std::vector<OutVertex> mat_vertices;
    std::vector<MatEdge::Data> mat_edges;
    RunStats stats;

private:
    Result<void> initialise();
    Result<void> finish() const;

    // Front navigation and geometry.
    const Site& site(std::uint32_t element) const { return sites[elements_[element].site]; }
    std::uint32_t next_element(std::uint32_t e) const { return vertices_[elements_[e].next].right; }
    std::uint32_t prev_element(std::uint32_t e) const { return vertices_[elements_[e].prev].left; }
    std::optional<Vec2> position(std::uint32_t vertex, double t) const;
    bool live(std::uint32_t element, Vec2 p, double t) const;
    bool branch_ok(std::uint32_t vertex, Vec2 p) const;
    std::size_t loop_size(std::uint32_t element, std::size_t limit) const;

    // Construction.
    SiteId add_site(const Site& s, std::uint32_t segment, bool corner);
    std::uint32_t add_element(SiteId site, std::uint32_t loop);
    std::uint32_t add_mat_vertex(Vec2 p, double t, VertexKind kind, const std::vector<std::uint32_t>& elements);
    Result<std::uint32_t> add_shock(std::uint32_t left, std::uint32_t right, Vec2 p, double t, std::uint32_t mat_vertex,
                                    std::optional<Vec2> direction = std::nullopt);
    void close_edge(std::uint32_t vertex, std::uint32_t mat_vertex, Vec2 p, double t);
    void relabel_loop(std::uint32_t start);

    // [EV-10], [EV-11] Scheduling; all pairs within a loop [EV-12].
    void schedule_collapse(std::uint32_t e);
    void schedule_split(std::uint32_t v, std::uint32_t c);
    void schedule_regular_split(std::uint32_t g, std::uint32_t c);
    void schedule_contact(std::uint32_t x, std::uint32_t y);
    void reschedule_element(std::uint32_t e);
    void reschedule_vertex(std::uint32_t v);
    void push(Event ev);

    // [EV-13] Pop-time validity, [RB-02] cluster detection, and the handlers.
    bool valid(const Event& ev) const;
    Result<void> check_cluster(const Event& ev);
    std::optional<Result<void>> annihilate(const Event& ev);
    Result<void> handle_collapse(const Event& ev);
    Result<void> handle_split(const Event& ev);
    Result<void> handle_contact(const Event& ev);
    Result<void> check_invariants() const;

    const PreparedRegion& region_;
    Tolerances tol_;
    bool debug_ = false;
    double t_now_ = 0.0;
    std::uint32_t next_loop_ = 1;
    std::uint64_t seq_ = 0;
    std::vector<FrontElement> elements_;
    std::vector<FrontVertex> vertices_;
    std::priority_queue<Event, std::vector<Event>, Later> queue_;
};

} // namespace wfmat::detail
