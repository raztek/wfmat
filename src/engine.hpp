// The propagation engine [ALG-02..ALG-05], [EV-01..EV-13]: the front, the event queue and the
// handlers. Internal to the library; compute_mat in src/mat.cpp is the public entry point.
#pragma once

#include <algorithm>
#include <array>
#include <cstdint>
#include <limits>
#include <memory>
#include <optional>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "wfmat/mat.hpp"

namespace wfmat::detail {

inline bool overlap(const Box& a, const Box& b)
{
    return a.lo.x <= b.hi.x && b.lo.x <= a.hi.x && a.lo.y <= b.hi.y && b.lo.y <= a.hi.y;
}

// "(x, y) at t = ..." for diagnostics, in the unit frame.
std::string point_text(Vec2 p, double t);

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
    Vec2 q;  // plateau: the end of the overlap; p is its start along the front of elements[0]
    std::array<std::uint32_t, 2> elements{no_id, no_id};
    std::array<std::uint32_t, 2> versions{0, 0};
    std::uint32_t vertex = no_id;  // split: the shock that reaches elements[0]
};

// [RB-01] Earliest time first, then kind, then the smallest id; then the rest of the event's
// identity, so that the order does not depend on when an event was found [EV-12].
struct Later {
    bool operator()(const Event& a, const Event& b) const
    {
        if (a.t != b.t) return a.t > b.t;
        if (a.kind != b.kind) return a.kind > b.kind;
        if (a.key != b.key) return a.key > b.key;
        const auto id = [](const Event& e) {
            return std::tuple{e.elements[0], e.elements[1], e.vertex, e.versions[0], e.versions[1]};
        };
        if (id(a) != id(b)) return id(a) > id(b);
        return a.seq > b.seq;
    }
};

// The event heap [EV-10]: a binary heap on a vector, so that stale events can be dropped in bulk.
class EventQueue {
public:
    bool empty() const { return heap_.empty(); }
    std::size_t size() const { return heap_.size(); }
    const Event& top() const { return heap_.front(); }
    void push(const Event& ev)
    {
        heap_.push_back(ev);
        std::push_heap(heap_.begin(), heap_.end(), Later{});
    }
    void pop()
    {
        std::pop_heap(heap_.begin(), heap_.end(), Later{});
        heap_.pop_back();
    }
    // Drops the events matching pred and restores the heap; returns how many were dropped.
    template <class Pred>
    std::size_t remove_if(Pred pred)
    {
        const std::size_t before = heap_.size();
        std::erase_if(heap_, pred);
        std::make_heap(heap_.begin(), heap_.end(), Later{});
        return before - heap_.size();
    }

private:
    std::vector<Event> heap_;
};

struct BroadIndex;  // the R-tree of the windowed broad phase, in engine.cpp

class Engine {
public:
    Engine(const PreparedRegion& region, const Options& options);
    ~Engine();

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
    std::uint32_t add_mat_vertex_at_sites(Vec2 p, double t, VertexKind kind, const std::vector<SiteId>& contacts);
    Result<std::uint32_t> add_shock(std::uint32_t left, std::uint32_t right, Vec2 p, double t, std::uint32_t mat_vertex,
                                    std::optional<Vec2> direction = std::nullopt);
    void close_edge(std::uint32_t vertex, std::uint32_t mat_vertex, Vec2 p, double t);
    void relabel_loops(const std::vector<std::uint32_t>& starts);

    // [EV-11] The windowed broad phase; with all_pairs [EV-12], every element of the loop is a
    // candidate and the window is unbounded.
    Box element_box(std::uint32_t e, double t) const;
    double initial_window() const;
    bool open_window();
    void index_element(std::uint32_t e);
    void near_elements(const Box& b, std::uint32_t element, std::vector<std::uint32_t>& out);
    void compact_queue();

    // [EV-10], [EV-11] Scheduling. Non-local candidates are searched from t_from() to the end of
    // the window; earlier roots belong to earlier windows.
    double t_from() const { return std::max(t_now_, push_floor_); }
    void schedule_collapse(std::uint32_t e);
    void schedule_split(std::uint32_t v, std::uint32_t c);
    void schedule_regular_split(std::uint32_t g, std::uint32_t c);
    void schedule_contact(std::uint32_t x, std::uint32_t y);
    std::optional<std::pair<Vec2, Vec2>> plateau_ends(std::uint32_t x, std::uint32_t y, const Contact& k) const;
    void reschedule_element(std::uint32_t e);
    void reschedule_vertex(std::uint32_t v);
    void push(Event ev);

    // [EV-13] Pop-time validity, [RB-02] cluster formation, and the handlers.
    bool valid(const Event& ev) const;
    double cluster_radius() const { return 10.0 * tol_.geom; }
    double plateau_distance(const Event& plateau, Vec2 x) const;
    double event_gap(const Event& a, const Event& b) const;
    std::vector<Event> gather_cluster(const Event& ev);
    Result<void> resolve_cluster(const std::vector<Event>& cluster);
    std::optional<Result<void>> annihilate(const Event& ev);
    Result<void> handle_collapse(const Event& ev);
    Result<void> handle_split(const Event& ev);
    Result<void> handle_contact(const Event& ev);
    Result<void> check_invariants() const;

    const PreparedRegion& region_;
    Tolerances tol_;
    bool debug_ = false;
    bool resolve_all_ = false;
    double t_now_ = 0.0;
    std::uint32_t next_loop_ = 1;
    std::uint64_t seq_ = 0;
    std::vector<FrontElement> elements_;
    std::vector<FrontVertex> vertices_;
    EventQueue queue_;
    std::size_t compact_at_ = 4096;

    // Broad phase state. Non-collapse events are kept only when push_floor_ < t <= win_hi_.
    static constexpr double unbounded = std::numeric_limits<double>::infinity();
    bool windowed_ = true;
    double win_hi_ = unbounded;
    double push_floor_ = -unbounded;
    double delta_ = 0.0;
    std::vector<std::uint32_t> live_;  // elements alive at the last window start, plus newer ones
    std::size_t listed_ = 0;           // elements_ below this index have been through live_
    std::vector<Box> window_box_;      // the indexed box of each element; empty when not indexed
    std::unique_ptr<BroadIndex> index_;
    std::size_t window_pairs_ = 0;      // candidate pairs found when the window opened
    std::size_t window_events_ = 0;     // events handled in the window
    std::vector<std::uint32_t> near_;   // scratch for near_elements

    // relabel_loops: per-element marks, valid for the current epoch only.
    std::vector<std::uint64_t> mark_;
    std::uint32_t epoch_ = 0;
};

} // namespace wfmat::detail
