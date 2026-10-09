#include "engine.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <numbers>
#include <string>

namespace wfmat::detail {

namespace {

constexpr double two_pi = 2.0 * std::numbers::pi;

// Step along a known direction of travel used to pick the branch of a shock born on its axis.
constexpr double branch_probe = 1e-6;

// Upper bound on any event time in the unit frame: an inscribed radius never exceeds the
// half-diagonal of the bounding box, which is 1 [IN-07].
constexpr double t_max = 1.5;

// Distance within which the third vertex of a vanishing three-element loop must meet the other two.
constexpr double annihilation_slack = 1e-7;

std::string point_text(Vec2 p, double t)
{
    char buf[96];
    std::snprintf(buf, sizeof buf, "(%.12g, %.12g) at t = %.12g", p.x, p.y, t);
    return buf;
}

tl::unexpected<Error> unsupported(const char* requirement, const std::string& what, Vec2 p, double t,
                                  const char* milestone)
{
    return make_error(ErrorCode::unsupported, requirement,
                      what + " near " + point_text(p, t) + " (unit frame); handled from milestone " + milestone, 0);
}

// Unit tangent of the front of s at p in front order (the region on the left).
Vec2 front_tangent(const Site& s, Vec2 p)
{
    const Vec2 g = s.gradient(p);
    return {g.y, -g.x};
}

} // namespace

Engine::Engine(const PreparedRegion& region, const Options& options)
    : region_(region), tol_(options.tol), debug_(options.debug_checks)
{
}

Result<void> Engine::run()
{
    if (auto r = initialise(); !r) return r;
    while (!queue_.empty()) {
        const Event ev = queue_.top();
        queue_.pop();
        if (!valid(ev)) {
            ++stats.stale_events;
            continue;
        }
        if (auto r = check_cluster(ev); !r) return r;
        t_now_ = std::max(t_now_, ev.t);
        Result<void> r;
        switch (ev.kind) {
        case EventKind::collapse: r = handle_collapse(ev); break;
        case EventKind::split: r = handle_split(ev); break;
        case EventKind::contact: r = handle_contact(ev); break;
        case EventKind::plateau: r = unsupported("RB-04", "a plateau edge", ev.p, ev.t, "M5"); break;
        }
        if (!r) return r;
        if (debug_)
            if (auto c = check_invariants(); !c) return c;
    }
    return finish();
}

// [ALG-03] One element per site in boundary order, a point element at every reflex corner, a
// shock at every convex corner and a regular vertex on each side of a reflex corner.
Result<void> Engine::initialise()
{
    const auto& segs = region_.segments;
    const std::size_t n = segs.size();
    for (std::size_t i = 0; i < n; ++i) {
        if (segs[i].is_arc())
            return make_error(ErrorCode::unsupported, "M3", "circular arcs are handled from milestone M3", 0,
                              region_.sources[i].front());
        if (region_.joins[i].kind == JoinKind::tangent)
            return make_error(ErrorCode::unsupported, "M3", "tangent joins are handled from milestone M3", 0,
                              region_.sources[i].back());
    }

    std::vector<std::uint32_t> line(n), corner(n, no_id);
    for (std::uint32_t i = 0; i < n; ++i) {
        line[i] = add_element(add_site(site_of(segs[i]), region_.sources[i].front(), false), 0);
        if (region_.joins[i].kind == JoinKind::reflex)
            corner[i] = add_element(add_site(Site::point(region_.joins[i].p), region_.sources[i].back(), true), 0);
    }

    auto link = [&](std::uint32_t v) {
        elements_[vertices_[v].left].next = v;
        elements_[vertices_[v].right].prev = v;
    };
    auto regular = [&](std::uint32_t left, std::uint32_t right, Vec2 corner_point, Vec2 m) {
        FrontVertex v;
        v.type = VertexType::regular;
        v.left = left;
        v.right = right;
        v.p0 = corner_point;
        v.m = m;
        vertices_.push_back(v);
        link(static_cast<std::uint32_t>(vertices_.size() - 1));
    };
    for (std::uint32_t i = 0; i < n; ++i) {
        const Join& join = region_.joins[i];
        const std::uint32_t next = line[(i + 1) % n];
        if (join.kind == JoinKind::reflex) {
            regular(line[i], corner[i], join.p, site(line[i]).n);
            regular(corner[i], next, join.p, site(next).n);
        } else {
            const std::uint32_t mv = add_mat_vertex(join.p, 0.0, VertexKind::corner, {line[i], next});
            auto v = add_shock(line[i], next, join.p, 0.0, mv);
            if (!v) return tl::unexpected(v.error());
            link(*v);
        }
    }

    const auto m = static_cast<std::uint32_t>(elements_.size());
    for (std::uint32_t e = 0; e < m; ++e) schedule_collapse(e);
    for (std::uint32_t v = 0; v < vertices_.size(); ++v)
        if (vertices_[v].type == VertexType::shock)
            for (std::uint32_t e = 0; e < m; ++e) schedule_split(v, e);
    for (std::uint32_t x = 0; x < m; ++x)
        for (std::uint32_t y = x + 1; y < m; ++y) schedule_contact(x, y);
    return {};
}

// [ALG-05] Every loop has vanished and every MAT edge has both ends.
Result<void> Engine::finish() const
{
    const auto left = std::count_if(elements_.begin(), elements_.end(), [](const FrontElement& e) { return e.alive; });
    if (left > 0)
        return make_error(ErrorCode::invariant_violation, "ALG-05",
                          "the event queue ran dry with " + std::to_string(left) + " front elements left", 0);
    for (const MatEdge::Data& d : mat_edges)
        if (d.v1 == no_id)
            return make_error(ErrorCode::invariant_violation, "EV-07", "a MAT edge was never closed", 0);
    return {};
}

std::optional<Vec2> Engine::position(std::uint32_t vertex, double t) const
{
    const FrontVertex& v = vertices_[vertex];
    if (v.type == VertexType::regular) return v.p0 + t * v.m;
    return shock_position(site(v.left), site(v.right), t, v.branch);
}

// [EV-13] p, a point of the element's front at time t, lies on its live part: between the
// positions of its two vertices at t, with eps_geom slack.
bool Engine::live(std::uint32_t element, Vec2 p, double t) const
{
    const FrontElement& e = elements_[element];
    const auto pl = position(e.prev, t), pr = position(e.next, t);
    if (!pl || !pr) return false;
    const Site& s = sites[e.site];
    const double slack = tol_.geom;
    if (s.is_line()) {
        const Vec2 tan{s.n.y, -s.n.x};
        const double sp = dot(tan, p - *pl), sr = dot(tan, *pr - *pl);
        return sp >= -slack && sp <= sr + slack;
    }
    // A point element's front runs clockwise about the site, over less than half a turn.
    const Vec2 c = s.centre;
    const double r = dist(p, c);
    if (r <= slack) return true;
    const double ang = slack / r;
    auto angle = [&](Vec2 q) { return std::atan2(q.y - c.y, q.x - c.x); };
    const double al = angle(*pl);
    double sweep = wrap_2pi(al - angle(*pr));
    double sp = wrap_2pi(al - angle(p));
    if (sweep > two_pi - ang) sweep = 0.0;
    if (sp > two_pi - ang) sp = 0.0;
    return sp <= sweep + ang;
}

// A shock's trajectory never crosses its pair's axis [KN-03], so a candidate point must lie on the
// side fixed at birth. Points on the axis itself (tangencies) are let through.
bool Engine::branch_ok(std::uint32_t vertex, Vec2 p) const
{
    const FrontVertex& v = vertices_[vertex];
    if (v.type == VertexType::regular || v.branch == 0) return true;
    const int b = branch_of(site(v.left), site(v.right), p);
    return b == 0 || b == v.branch;
}

std::size_t Engine::loop_size(std::uint32_t element, std::size_t limit) const
{
    std::size_t k = 0;
    std::uint32_t x = element;
    do {
        ++k;
        x = next_element(x);
    } while (x != element && k < limit);
    return k;
}

SiteId Engine::add_site(const Site& s, std::uint32_t segment, bool corner)
{
    sites.push_back(s);
    site_info.push_back({s, segment, corner});
    return static_cast<SiteId>(sites.size() - 1);
}

std::uint32_t Engine::add_element(SiteId s, std::uint32_t loop)
{
    FrontElement e;
    e.site = s;
    e.loop = loop;
    elements_.push_back(e);
    return static_cast<std::uint32_t>(elements_.size() - 1);
}

std::uint32_t Engine::add_mat_vertex(Vec2 p, double t, VertexKind kind, std::initializer_list<std::uint32_t> elements)
{
    OutVertex v{p, t, kind, {}};
    for (std::uint32_t e : elements) {
        const Site& s = site(e);
        v.contacts.push_back({elements_[e].site, s.foot(p, t)});
        stats.max_residual = std::max(stats.max_residual, std::abs(s.distance(p) - t));
    }
    mat_vertices.push_back(std::move(v));
    return static_cast<std::uint32_t>(mat_vertices.size() - 1);
}

// [EV-06], [EV-07] A new shock between two elements meeting at a corner, and the MAT edge it opens.
// The branch comes from the birth point, or from the direction of travel when the shock is born on
// its axis (a contact, E3).
Result<std::uint32_t> Engine::add_shock(std::uint32_t left, std::uint32_t right, Vec2 p, double t,
                                        std::uint32_t mat_vertex, std::optional<Vec2> direction)
{
    const Site &a = site(left), &b = site(right);
    if (a.is_line() && b.is_line() && std::abs(cross(a.n, b.n)) <= tol_.ang)
        return unsupported("RB-04", "a shock between parallel lines (a plateau)", p, t, "M5");
    int branch = 0;
    if (!(a.is_line() && b.is_line())) {
        branch = direction ? branch_of(a, b, p + branch_probe * *direction) : branch_of(a, b, p);
        if (branch == 0) return unsupported("RB-02", "a shock born on its axis (simultaneous events)", p, t, "M5");
    }

    MatEdge::Data d;
    d.v0 = mat_vertex;
    d.left = elements_[left].site;
    d.right = elements_[right].site;
    d.kind = bisector_kind(a, b, tol_.ang);
    d.branch = branch;
    d.a = a;
    d.b = b;
    d.t0 = t;
    d.p0 = p;
    d.xf = region_.xf;
    mat_edges.push_back(d);

    FrontVertex v;
    v.type = VertexType::shock;
    v.left = left;
    v.right = right;
    v.t0 = t;
    v.p0 = p;
    v.branch = branch;
    v.edge = static_cast<EdgeId>(mat_edges.size() - 1);
    vertices_.push_back(v);
    return static_cast<std::uint32_t>(vertices_.size() - 1);
}

void Engine::close_edge(std::uint32_t vertex, std::uint32_t mat_vertex, Vec2 p, double t)
{
    const FrontVertex& v = vertices_[vertex];
    if (v.type != VertexType::shock) return;
    MatEdge::Data& d = mat_edges[v.edge];
    d.v1 = mat_vertex;
    d.t1 = std::max(t, d.t0);
    d.p1 = p;
}

void Engine::relabel_loop(std::uint32_t start)
{
    const std::uint32_t id = next_loop_++;
    std::uint32_t x = start;
    do {
        elements_[x].loop = id;
        x = next_element(x);
    } while (x != start);
}

void Engine::push(Event ev)
{
    ev.seq = seq_++;
    ev.key = std::min(ev.elements[0], ev.elements[1]);
    if (ev.vertex != no_id)
        ev.key = std::min({ev.key, vertices_[ev.vertex].left, vertices_[ev.vertex].right});
    queue_.push(ev);
}

// [EV-01], [EV-02] The element shrinks to a point when its two vertices meet.
void Engine::schedule_collapse(std::uint32_t e)
{
    const FrontElement& el = elements_[e];
    const std::uint32_t lv = el.prev, rv = el.next;
    const FrontVertex &l = vertices_[lv], &r = vertices_[rv];
    const std::uint32_t a = l.left, c = r.right;
    if (a == e || c == e || a == c) return;

    auto enqueue = [&](double t, Vec2 p) {
        Event ev;
        ev.t = t;
        ev.kind = EventKind::collapse;
        ev.p = p;
        ev.elements = {e, no_id};
        ev.versions = {el.version, 0};
        push(ev);
    };

    if (l.type == VertexType::shock && r.type == VertexType::shock) {
        const SiteId sa = elements_[a].site, sc = elements_[c].site;
        if (sa == el.site || sc == el.site || sa == sc) return;
        const Solve3 sol = solve_three(sites[sa], sites[el.site], sites[sc], t_now_, tol_, t_max);
        for (const Root3& root : sol.roots) {
            if (root.t < std::max(l.t0, r.t0) - tol_.time) continue;
            if (!branch_ok(lv, root.p) || !branch_ok(rv, root.p)) continue;
            enqueue(root.t, root.p);
            return;
        }
        return;
    }
    if (l.type == VertexType::regular && r.type == VertexType::regular) return;

    // E1b: the shock reaches the regular vertex q(t) = q0 + t m when the shock's other site is at
    // distance t from q(t).
    const bool right_regular = r.type == VertexType::regular;
    const FrontVertex& g = right_regular ? r : l;
    const std::uint32_t shock = right_regular ? lv : rv;
    const std::uint32_t other = right_regular ? a : c;
    const auto t = regular_hit_time(g.p0, g.m, site(other), t_now_, tol_);
    if (!t || *t < vertices_[shock].t0 - tol_.time) return;
    const Vec2 q = g.p0 + *t * g.m;
    if (!branch_ok(shock, q)) return;
    enqueue(*t, q);
}

// [EV-04] Shock v reaches the interior of the non-adjacent element c.
void Engine::schedule_split(std::uint32_t v, std::uint32_t c)
{
    const FrontVertex& vx = vertices_[v];
    if (vx.type != VertexType::shock || !vx.alive || !elements_[c].alive) return;
    const std::uint32_t a = vx.left, b = vx.right;
    if (c == a || c == b || c == prev_element(a) || c == next_element(b)) return;
    if (elements_[c].loop != elements_[a].loop) return;
    const SiteId sa = elements_[a].site, sb = elements_[b].site, sc = elements_[c].site;
    if (sc == sa || sc == sb) return;

    const Solve3 sol = solve_three(sites[sa], sites[sb], sites[sc], t_now_, tol_, t_max);
    for (const Root3& root : sol.roots) {
        if (root.t < vx.t0 - tol_.time) continue;
        if (!branch_ok(v, root.p) || !live(c, root.p, root.t)) continue;
        Event ev;
        ev.t = root.t;
        ev.kind = EventKind::split;
        ev.p = root.p;
        ev.elements = {c, no_id};
        ev.versions = {elements_[c].version, 0};
        ev.vertex = v;
        push(ev);
        return;
    }
}

// [EV-05], [KN-05] Two non-adjacent elements touch along a common normal.
void Engine::schedule_contact(std::uint32_t x, std::uint32_t y)
{
    if (!elements_[x].alive || !elements_[y].alive) return;
    if (y == x || y == prev_element(x) || y == next_element(x)) return;
    if (elements_[x].loop != elements_[y].loop) return;
    const SiteId sx = elements_[x].site, sy = elements_[y].site;
    if (sx == sy) return;

    for (const Contact& k : contacts(sites[sx], sites[sy], t_now_, tol_)) {
        Event ev;
        ev.t = k.t;
        ev.p = k.p;
        ev.elements = {x, y};
        ev.versions = {elements_[x].version, elements_[y].version};
        if (k.plateau) {
            // Antiparallel lines: report the plateau only when the live parts overlap at t*.
            const Site& s = sites[sx];
            if (!s.is_line()) continue;
            const Vec2 tan{s.n.y, -s.n.x};
            const auto x0 = position(elements_[x].prev, k.t), x1 = position(elements_[x].next, k.t);
            const auto y0 = position(elements_[y].next, k.t), y1 = position(elements_[y].prev, k.t);
            if (!x0 || !x1 || !y0 || !y1) continue;
            const double lo = std::max(dot(tan, *x0), dot(tan, *y0));
            const double hi = std::min(dot(tan, *x1), dot(tan, *y1));
            if (hi - lo <= tol_.geom) continue;
            ev.kind = EventKind::plateau;
            ev.p = k.p + (0.5 * (lo + hi) - dot(tan, k.p)) * tan;
            push(ev);
            return;
        }
        if (!live(x, k.p, k.t) || !live(y, k.p, k.t)) continue;
        ev.kind = EventKind::contact;
        push(ev);
        return;
    }
}

// The element's own collapse, every shock of its loop against it, and every contact with it.
void Engine::reschedule_element(std::uint32_t e)
{
    schedule_collapse(e);
    std::uint32_t x = e;
    do {
        if (x != e) schedule_contact(e, x);
        schedule_split(elements_[x].next, e);
        x = next_element(x);
    } while (x != e);
}

void Engine::reschedule_vertex(std::uint32_t v)
{
    const std::uint32_t start = vertices_[v].left;
    std::uint32_t x = start;
    do {
        schedule_split(v, x);
        x = next_element(x);
    } while (x != start);
}

bool Engine::valid(const Event& ev) const
{
    const int n = ev.kind == EventKind::contact || ev.kind == EventKind::plateau ? 2 : 1;
    for (int i = 0; i < n; ++i) {
        const FrontElement& e = elements_[ev.elements[i]];
        if (!e.alive || e.version != ev.versions[i]) return false;
    }
    if (n == 2 && elements_[ev.elements[0]].loop != elements_[ev.elements[1]].loop) return false;
    if (ev.kind == EventKind::split) {
        const FrontVertex& v = vertices_[ev.vertex];
        if (!v.alive || elements_[v.left].loop != elements_[ev.elements[0]].loop) return false;
    }
    return true;
}

// [RB-02] Other valid events at the same time and place form a cluster with ev. Cluster
// resolution arrives with M5; until then a cluster is reported, except the three collapses of a
// vanishing three-element loop, which handle_collapse treats as one event.
Result<void> Engine::check_cluster(const Event& ev)
{
    auto same = [](const Event& a, const Event& b) {
        if (a.kind != b.kind || a.vertex != b.vertex) return false;
        return (a.elements[0] == b.elements[0] && a.elements[1] == b.elements[1]) ||
               (a.elements[0] == b.elements[1] && a.elements[1] == b.elements[0]);
    };
    const bool vanishing = ev.kind == EventKind::collapse && loop_size(ev.elements[0], 4) == 3;

    std::vector<Event> held;
    bool cluster = false;
    while (!queue_.empty() && queue_.top().t <= ev.t + tol_.time) {
        const Event e = queue_.top();
        queue_.pop();
        if (!valid(e)) {
            ++stats.stale_events;
            continue;
        }
        held.push_back(e);
        if (dist(e.p, ev.p) > tol_.geom || same(e, ev)) continue;
        if (vanishing && e.kind == EventKind::collapse &&
            elements_[e.elements[0]].loop == elements_[ev.elements[0]].loop)
            continue;
        cluster = true;
    }
    for (const Event& e : held) queue_.push(e);
    if (cluster) return unsupported("RB-03", "simultaneous events (a cluster)", ev.p, ev.t, "M5");
    return {};
}

// [EV-01] E1a and [EV-02] E1b; [EV-08] a three-element loop vanishes at its junction.
Result<void> Engine::handle_collapse(const Event& ev)
{
    const std::uint32_t e = ev.elements[0];
    const std::uint32_t lv = elements_[e].prev, rv = elements_[e].next;
    const std::uint32_t a = vertices_[lv].left, c = vertices_[rv].right;
    const bool lshock = vertices_[lv].type == VertexType::shock;
    const bool rshock = vertices_[rv].type == VertexType::shock;

    const std::uint32_t wv = elements_[c].next;
    if (vertices_[wv].right == a) {
        if (!lshock || !rshock)
            return unsupported("EV-08", "a three-element loop collapsing through a regular vertex", ev.p, ev.t, "M5");
        const auto pw = position(wv, ev.t);
        if (!pw || dist(*pw, ev.p) > annihilation_slack)
            return make_error(ErrorCode::numerical_failure, "EV-08",
                              "the vertices of a vanishing loop do not meet near " + point_text(ev.p, ev.t), 0);
        const bool wshock = vertices_[wv].type == VertexType::shock;
        const std::uint32_t mv =
            add_mat_vertex(ev.p, ev.t, wshock ? VertexKind::junction : VertexKind::extremum_max, {a, e, c});
        for (std::uint32_t v : {lv, rv, wv}) {
            close_edge(v, mv, ev.p, ev.t);
            vertices_[v].alive = false;
        }
        for (std::uint32_t x : {a, e, c}) {
            elements_[x].alive = false;
            ++elements_[x].version;
        }
        ++stats.annihilations;
        return {};
    }

    const std::uint32_t mv =
        add_mat_vertex(ev.p, ev.t, lshock && rshock ? VertexKind::junction : VertexKind::transition, {a, e, c});
    close_edge(lv, mv, ev.p, ev.t);
    close_edge(rv, mv, ev.p, ev.t);
    vertices_[lv].alive = vertices_[rv].alive = false;
    elements_[e].alive = false;
    ++elements_[e].version;

    auto v = add_shock(a, c, ev.p, ev.t, mv);
    if (!v) return tl::unexpected(v.error());
    elements_[a].next = *v;
    elements_[c].prev = *v;
    ++elements_[a].version;
    ++elements_[c].version;
    ++(lshock && rshock ? stats.collapses : stats.transitions);

    reschedule_element(a);
    reschedule_element(c);
    reschedule_vertex(*v);
    return {};
}

// [EV-04] Shock A|B reaches element C: C splits into C1, C2; new shocks A|C2 and C1|B; the loop
// splits in two [EV-09].
Result<void> Engine::handle_split(const Event& ev)
{
    const std::uint32_t v = ev.vertex, c = ev.elements[0];
    const std::uint32_t a = vertices_[v].left, b = vertices_[v].right;
    const std::uint32_t mv = add_mat_vertex(ev.p, ev.t, VertexKind::junction, {a, b, c});
    close_edge(v, mv, ev.p, ev.t);
    vertices_[v].alive = false;

    const FrontElement old = elements_[c];
    elements_[c].alive = false;
    ++elements_[c].version;
    const std::uint32_t c1 = add_element(old.site, old.loop), c2 = add_element(old.site, old.loop);
    elements_[c1].prev = old.prev;
    vertices_[old.prev].right = c1;
    elements_[c2].next = old.next;
    vertices_[old.next].left = c2;

    auto x = add_shock(a, c2, ev.p, ev.t, mv);
    if (!x) return tl::unexpected(x.error());
    auto y = add_shock(c1, b, ev.p, ev.t, mv);
    if (!y) return tl::unexpected(y.error());
    elements_[a].next = *x;
    elements_[c2].prev = *x;
    elements_[c1].next = *y;
    elements_[b].prev = *y;
    ++elements_[a].version;
    ++elements_[b].version;
    relabel_loop(c1);
    ++stats.splits;

    for (std::uint32_t e : {a, c2, b, c1}) reschedule_element(e);
    reschedule_vertex(*x);
    reschedule_vertex(*y);
    return {};
}

// [EV-05] X and Y touch at p: both split, and two shocks X1|Y2 and Y1|X2 move apart along the
// bisector, one on each side of the common normal; the loop splits in two.
Result<void> Engine::handle_contact(const Event& ev)
{
    const std::uint32_t x = ev.elements[0], y = ev.elements[1];
    const std::uint32_t mv = add_mat_vertex(ev.p, ev.t, VertexKind::extremum_min, {x, y});
    const Vec2 tan = front_tangent(site(x), ev.p);

    const FrontElement ox = elements_[x], oy = elements_[y];
    for (std::uint32_t e : {x, y}) {
        elements_[e].alive = false;
        ++elements_[e].version;
    }
    const std::uint32_t x1 = add_element(ox.site, ox.loop), x2 = add_element(ox.site, ox.loop);
    const std::uint32_t y1 = add_element(oy.site, ox.loop), y2 = add_element(oy.site, ox.loop);
    elements_[x1].prev = ox.prev;
    vertices_[ox.prev].right = x1;
    elements_[x2].next = ox.next;
    vertices_[ox.next].left = x2;
    elements_[y1].prev = oy.prev;
    vertices_[oy.prev].right = y1;
    elements_[y2].next = oy.next;
    vertices_[oy.next].left = y2;

    auto u = add_shock(x1, y2, ev.p, ev.t, mv, -tan);
    if (!u) return tl::unexpected(u.error());
    auto w = add_shock(y1, x2, ev.p, ev.t, mv, tan);
    if (!w) return tl::unexpected(w.error());
    elements_[x1].next = *u;
    elements_[y2].prev = *u;
    elements_[y1].next = *w;
    elements_[x2].prev = *w;
    relabel_loop(y1);
    ++stats.contacts;

    for (std::uint32_t e : {x1, y2, y1, x2}) reschedule_element(e);
    reschedule_vertex(*u);
    reschedule_vertex(*w);
    return {};
}

// [RB-06] Loops are closed and consistently linked, live parts are not inverted.
Result<void> Engine::check_invariants() const
{
    auto fail = [&](const std::string& what) {
        return make_error(ErrorCode::invariant_violation, "RB-06", what + " at t = " + std::to_string(t_now_), 0);
    };
    for (std::uint32_t e = 0; e < elements_.size(); ++e) {
        const FrontElement& el = elements_[e];
        if (!el.alive) continue;
        const FrontVertex &l = vertices_[el.prev], &r = vertices_[el.next];
        if (!l.alive || !r.alive || l.right != e || r.left != e) return fail("broken links at element " + std::to_string(e));
        if (elements_[l.left].loop != el.loop || elements_[r.right].loop != el.loop)
            return fail("inconsistent loop ids at element " + std::to_string(e));
        const Site& s = sites[el.site];
        if (s.is_line()) {
            const auto pl = position(el.prev, t_now_), pr = position(el.next, t_now_);
            if (!pl || !pr) return fail("a vertex has no position at element " + std::to_string(e));
            const Vec2 tan{s.n.y, -s.n.x};
            if (dot(tan, *pr - *pl) < -1e-8) return fail("negative live length at element " + std::to_string(e));
        }
    }
    return {};
}

} // namespace wfmat::detail
