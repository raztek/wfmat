// [RB-02], [RB-03] Cluster formation and the generic resolution of simultaneous events.
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <map>
#include <set>
#include <tuple>
#include <numbers>
#include <utility>

#include "engine.hpp"

namespace wfmat::detail {

namespace {

constexpr double two_pi = 2.0 * std::numbers::pi;

// A convex arc centred at a cluster point touches the disk there along its whole length; its ends
// are among the other contacts, as in [EV-08].
constexpr double centre_slack = 1e-7;

// How far a shock that ends at a tangency of its two fronts can be from it a moment earlier.
constexpr double tangency_radius = 1e-4;

double segment_distance(Vec2 x, Vec2 a, Vec2 b)
{
    const Vec2 d = b - a;
    const double len2 = dot(d, d);
    const double s = len2 > 0.0 ? std::clamp(dot(x - a, d) / len2, 0.0, 1.0) : 0.0;
    return dist(x, a + s * d);
}

double angle_of(Vec2 u) { return std::atan2(u.y, u.x); }

// The point of the front of s at time t nearest to x.
Vec2 front_point(const Site& s, double t, Vec2 x)
{
    if (s.is_line()) return x - (s.distance(x) - t) * s.n;
    const Vec2 d = x - s.centre;
    const double r = norm(d);
    const double rho = std::max(s.offset_radius(t), 0.0);
    return r > 0.0 ? s.centre + (rho / r) * d : s.centre + Vec2{rho, 0.0};
}

} // namespace

// Distance from x to the plateau of the event: the segment from p to q, or the arc from p to q
// along the front of elements[0] about the common centre.
double Engine::plateau_distance(const Event& pl, Vec2 x) const
{
    const Site& s = site(pl.elements[0]);
    if (s.is_line()) return segment_distance(x, pl.p, pl.q);
    const Vec2 c = s.centre;
    auto angle = [&](Vec2 q) { return s.sigma * std::atan2(q.y - c.y, q.x - c.x); };
    const double a0 = angle(pl.p);
    if (wrap_2pi(angle(x) - a0) <= wrap_2pi(angle(pl.q) - a0)) return std::abs(dist(x, c) - dist(pl.p, c));
    return std::min(dist(x, pl.p), dist(x, pl.q));
}

double Engine::event_gap(const Event& a, const Event& b) const
{
    const bool pa = a.kind == EventKind::plateau, pb = b.kind == EventKind::plateau;
    if (pa && pb)
        return std::min({plateau_distance(a, b.p), plateau_distance(a, b.q), plateau_distance(b, a.p),
                         plateau_distance(b, a.q)});
    if (pa) return plateau_distance(a, b.p);
    if (pb) return plateau_distance(b, a.p);
    return dist(a.p, b.p);
}

// [RB-02] Every other valid event within eps_t of ev, and within the cluster radius of an event
// already in the cluster, joins it (a plateau counts with its whole overlap). Repeats of one event
// are dropped; the remaining events go back to the queue.
std::vector<Event> Engine::gather_cluster(const Event& ev)
{
    std::vector<Event> held;
    while (!queue_.empty() && queue_.top().t <= ev.t + tol_.time) {
        const Event e = queue_.top();
        queue_.pop();
        if (!valid(e)) {
            ++stats.stale_events;
            continue;
        }
        held.push_back(e);
    }
    std::vector<Event> cluster{ev};
    if (held.empty()) return cluster;

    // Breadth-first growth over a grid of cells of the cluster radius, so that a cluster of many
    // events at one point (all the split candidates at the centre of a regular polygon) costs
    // linear time. Point events sit in the grid; plateaus are few and are checked directly.
    const double r = cluster_radius();
    using Cell = std::pair<std::int64_t, std::int64_t>;
    auto cell_of = [&](Vec2 p) { return Cell{std::llround(std::floor(p.x / r)), std::llround(std::floor(p.y / r))}; };
    std::map<Cell, std::vector<std::size_t>> grid;
    std::vector<std::size_t> flat;
    for (std::size_t i = 0; i < held.size(); ++i) {
        if (held[i].kind == EventKind::plateau)
            flat.push_back(i);
        else
            grid[cell_of(held[i].p)].push_back(i);
    }
    std::vector<char> taken(held.size(), 0);
    std::vector<const Event*> frontier{&ev};
    auto take = [&](std::size_t i, const Event& from) {
        if (taken[i] || event_gap(from, held[i]) > r) return;
        taken[i] = 1;
        frontier.push_back(&held[i]);
    };
    while (!frontier.empty()) {
        const Event& m = *frontier.back();
        frontier.pop_back();
        for (std::size_t i : flat) take(i, m);
        if (m.kind == EventKind::plateau) {
            for (auto& [cell, members] : grid)
                for (std::size_t i : members) take(i, m);
        } else {
            const Cell c = cell_of(m.p);
            for (std::int64_t dx = -1; dx <= 1; ++dx)
                for (std::int64_t dy = -1; dy <= 1; ++dy) {
                    const auto it = grid.find({c.first + dx, c.second + dy});
                    if (it == grid.end()) continue;
                    for (std::size_t i : it->second) take(i, m);
                    // Members taken once never need checking again.
                    std::erase_if(it->second, [&](std::size_t i) { return taken[i] != 0; });
                }
        }
    }

    // Repeats of one event (the same kind on the same elements and vertex) are dropped.
    auto key = [](const Event& e) {
        return std::tuple{e.kind, std::min(e.elements[0], e.elements[1]), std::max(e.elements[0], e.elements[1]), e.vertex};
    };
    std::set<decltype(key(ev))> seen{key(ev)};
    for (std::size_t i = 0; i < held.size(); ++i) {
        if (!taken[i])
            queue_.push(held[i]);
        else if (seen.insert(key(held[i])).second)
            cluster.push_back(held[i]);
    }
    return cluster;
}

// [RB-03] Generic resolution at time t*. The region D is a small disk about the event points, or
// the plateau with a disk at each end; each disk has an anchor that becomes a MAT vertex.
//  1. Vertices whose position lies in D die; shocks close their edges at their anchor.
//  2. Elements whose live part lies in D die. The others are cut where they enter and leave D:
//     the piece arriving at an anchor is an in-port, the piece leaving it an out-port.
//  3. Ports are sorted counter-clockwise around D, by direction and then by curvature, so that
//     tangent fronts come out in their true order. The region near D lies counter-clockwise of an
//     out-port and clockwise of an in-port, so each in-port links to the out-port just before it,
//     with a new shock moving into the sector between them [EV-06], [EV-07].
//  4. Each anchor's vertex takes its kind from its degree; a plateau adds a constant-r edge
//     between its two anchors. With no ports left, the loop has vanished [EV-08].
// The simple handlers are the cases with one anchor and one or two links [RB-04].
Result<void> Engine::resolve_cluster(const std::vector<Event>& cluster)
{
    const Event& first = cluster.front();
    const double t = first.t;
    const std::uint32_t start = first.elements[0];
    t_now_ = std::max(t_now_, t);

    // [EV-08] A loop of two or three elements whose vertices all meet: handled as one event, with
    // the tangency slack a two-shock loop needs.
    if (first.kind != EventKind::plateau && loop_size(start, 4) <= 3)
        if (auto r = annihilate(first)) return *r;

    auto failure = [&](const std::string& what) {
        return make_error(ErrorCode::numerical_failure, "RB-03", what + " near " + point_text(first.p, t), 0);
    };

    const std::uint32_t loop = elements_[start].loop;
    const Event* plateau = nullptr;
    for (const Event& e : cluster) {
        for (std::uint32_t x : e.elements)
            if (x != no_id && elements_[x].loop != loop) return failure("a cluster on two front loops");
        if (e.kind != EventKind::plateau) continue;
        if (plateau) return failure("a cluster with two plateaus");
        plateau = &e;
    }

    // 1. Anchors and the region.
    struct Anchor {
        Vec2 p;
        double cut = 0.0;  // plateau: the angle towards the other end; sorting starts there
        std::uint32_t mv = no_id;
        int shocks_in = 0, links = 0;
        std::vector<SiteId> contacts;
    };
    std::vector<Anchor> anchors;
    auto add_anchor = [&](Vec2 p, double cut) {
        Anchor a;
        a.p = p;
        a.cut = cut;
        anchors.push_back(a);
    };
    std::vector<Vec2> points;
    std::uint32_t px = no_id, py = no_id;
    if (plateau) {
        px = plateau->elements[0];
        py = plateau->elements[1];
        const Site& s = site(px);
        auto tangent = [&](Vec2 q) { const Vec2 g = s.gradient(q); return Vec2{g.y, -g.x}; };
        add_anchor(plateau->p, angle_of(tangent(plateau->p)));
        add_anchor(plateau->q, angle_of(-tangent(plateau->q)));
    } else {
        Vec2 sum;
        for (const Event& e : cluster) sum = sum + e.p;
        add_anchor(sum / static_cast<double>(cluster.size()), 0.0);
    }
    // One point per cell of half the cluster radius is enough to describe the region.
    {
        std::set<std::pair<std::int64_t, std::int64_t>> cells;
        const double h = 0.5 * cluster_radius();
        for (const Event& e : cluster)
            if (e.kind != EventKind::plateau &&
                cells.insert({std::llround(std::floor(e.p.x / h)), std::llround(std::floor(e.p.y / h))}).second)
                points.push_back(e.p);
    }
    for (const Anchor& a : anchors) points.push_back(a.p);

    const double rho = cluster_radius();
    // How far the points of a point cluster are from its anchor; a plateau's anchors are its ends.
    double spread = rho;
    if (!plateau)
        for (Vec2 q : points) spread = std::max(spread, dist(q, anchors[0].p));
    auto gap = [&](Vec2 x) {
        double d = std::numeric_limits<double>::infinity();
        for (Vec2 q : points) d = std::min(d, dist(x, q));
        if (plateau) d = std::min(d, plateau_distance(*plateau, x));
        return d;
    };
    auto nearest = [&](Vec2 x) {
        int best = 0;
        for (int k = 1; k < static_cast<int>(anchors.size()); ++k)
            if (dist(x, anchors[k].p) < dist(x, anchors[best].p)) best = k;
        return best;
    };

    // The elements whose live part comes near D [EV-11], in id order so that both broad phases
    // see them in the same sequence, and their vertices.
    Box reach;
    for (Vec2 q : points) reach.add(q);
    if (plateau) {
        reach.add(element_box(px, t));
        reach.add(element_box(py, t));
    }
    reach = reach.inflated(std::max(rho, tangency_radius));
    std::vector<std::uint32_t> elems;
    near_elements(reach, start, elems);
    std::erase_if(elems, [&](std::uint32_t e) { return !overlap(element_box(e, t), reach); });
    std::sort(elems.begin(), elems.end());
    const std::size_t m = elems.size();
    std::vector<std::uint32_t> verts;
    std::map<std::uint32_t, std::size_t> vindex;
    for (std::uint32_t e : elems)
        for (std::uint32_t v : {elements_[e].prev, elements_[e].next})
            if (vindex.emplace(v, verts.size()).second) verts.push_back(v);
    std::vector<char> dying(verts.size(), 0);
    std::vector<int> vanchor(verts.size(), 0);
    std::vector<std::optional<Vec2>> at(verts.size());  // where a dying vertex meets D
    for (std::size_t i = 0; i < verts.size(); ++i) {
        const auto q = position(verts[i], t);
        if (q && gap(*q) <= rho) {
            dying[i] = 1;
            vanchor[i] = nearest(*q);
            at[i] = q;
        } else if (q) {
            // A shock between nearly parallel fronts moves fast, so rounding in t* puts it far from
            // D at t*; it is dying if it passes D within the cluster's time window [RB-02].
            const auto a = position(verts[i], std::max(vertices_[verts[i]].t0, t - tol_.time));
            const auto b = position(verts[i], t + tol_.time);
            if (a && b && !plateau) {
                double d = std::numeric_limits<double>::infinity();
                Vec2 meet;
                for (Vec2 x : points)
                    if (segment_distance(x, *a, *b) < d) {
                        d = segment_distance(x, *a, *b);
                        meet = x;
                    }
                if (d <= rho) {
                    dying[i] = 1;
                    vanchor[i] = nearest(meet);
                    at[i] = meet;
                }
            }
        } else {
            // Rounding has taken t* just past the tangency where this shock ends; its position
            // a little earlier is near the tangency point, within sqrt(dt R) [EV-08].
            const auto e = position(verts[i], std::max(vertices_[verts[i]].t0, t - 1e-9));
            if (e && gap(*e) <= std::max(rho, tangency_radius)) {
                dying[i] = 1;
                vanchor[i] = nearest(*e);
            }
        }
    }
    // 2. Pieces of the surviving elements.
    struct Piece {
        std::uint32_t old;
        int start = -1, end = -1;  // anchor of an out-port start or an in-port end; -1 keeps the old vertex
    };
    std::vector<Piece> pieces;
    std::vector<char> touched(m, 0);
    auto contact = [&](int k, std::uint32_t e) { anchors[k].contacts.push_back(elements_[e].site); };
    auto midpoint = [&](std::uint32_t e) -> std::optional<Vec2> {
        const FrontElement& el = elements_[e];
        auto where = [&](std::uint32_t v) {
            const auto it = vindex.find(v);
            return it != vindex.end() && at[it->second] ? at[it->second] : position(v, t);
        };
        const auto pl = where(el.prev), pr = where(el.next);
        if (!pl || !pr) return std::nullopt;
        const Site& s = sites[el.site];
        if (s.is_line()) return 0.5 * (*pl + *pr);
        const Vec2 c = s.centre;
        auto angle = [&](Vec2 q) { return s.sigma * std::atan2(q.y - c.y, q.x - c.x); };
        const double al = angle(*pl);
        const double r = std::max(s.offset_radius(t), 0.0);
        double sweep = wrap_2pi(angle(*pr) - al);
        if (sweep > two_pi - (tol_.ang + (r > rho ? spread / r : 0.0))) sweep = 0.0;  // ends crossed within D
        return c + r * polar(s.sigma * (al + 0.5 * sweep));
    };

    for (std::size_t i = 0; i < m; ++i) {
        const std::uint32_t e = elems[i];
        const std::size_t ip = vindex.at(elements_[e].prev), in = vindex.at(elements_[e].next);
        const bool lp = dying[ip], ln = dying[in];
        if (e == px || e == py) {
            // A plateau element runs along the plateau from one anchor to the other.
            const int enter = e == px ? 0 : 1, leave = 1 - enter;
            touched[i] = 1;
            contact(0, e);
            contact(1, e);
            if (!lp) pieces.push_back({e, -1, enter});
            if (!ln) pieces.push_back({e, leave, -1});
            continue;
        }
        if (lp && ln) {
            touched[i] = 1;
            const auto mid = midpoint(e);
            if (!mid || gap(*mid) <= rho) {
                contact(mid ? nearest(*mid) : vanchor[in], e);
                continue;
            }
            pieces.push_back({e, vanchor[ip], vanchor[in]});
            contact(vanchor[ip], e);
            contact(vanchor[in], e);
        } else if (lp) {
            touched[i] = 1;
            pieces.push_back({e, vanchor[ip], -1});
            contact(vanchor[ip], e);
        } else if (ln) {
            touched[i] = 1;
            pieces.push_back({e, -1, vanchor[in]});
            contact(vanchor[in], e);
        } else {
            // The front passes through D between its two surviving vertices.
            const Site& s = site(e);
            for (int k = 0; k < static_cast<int>(anchors.size()); ++k) {
                const Vec2 f = front_point(s, t, anchors[k].p);
                if (gap(f) > rho || !live(e, f, t)) continue;
                touched[i] = 1;
                pieces.push_back({e, -1, k});
                pieces.push_back({e, k, -1});
                contact(k, e);
                break;
            }
        }
    }

    // Dying vertices.
    for (std::size_t i = 0; i < verts.size(); ++i)
        if (dying[i] && vertices_[verts[i]].type == VertexType::shock) ++anchors[vanchor[i]].shocks_in;

    // Kill the touched elements and create the pieces; a piece keeps the vertex it has not lost.
    struct Port {
        std::uint32_t element;
        bool in;
        int anchor;
        double angle;  // direction away from the anchor along the front
        double k;      // lateral curvature: larger is further counter-clockwise
        double slack;  // how far the direction can turn as the anchor moves within the cluster radius
        double rel = 0.0;
    };
    std::vector<Port> ports;
    std::vector<std::uint32_t> fresh;
    auto port = [&](std::uint32_t id, bool in, int k) {
        const Site& s = site(id);
        const Vec2 f = front_point(s, t, anchors[k].p);
        const Vec2 g = s.gradient(f);
        const Vec2 u = in ? Vec2{-g.y, g.x} : Vec2{g.y, -g.x};
        double lateral = 0.0, slack = 0.0;
        if (!s.is_line()) {
            const double r = s.offset_radius(t);
            if (r > rho) {
                lateral = dot(s.centre - f, perp(u)) / (r * r);
                slack = spread / r;
            }
        }
        ports.push_back({id, in, k, angle_of(u), lateral, slack, 0.0});
    };
    for (std::size_t i = 0; i < m; ++i) {
        if (!touched[i]) continue;
        elements_[elems[i]].alive = false;
        ++elements_[elems[i]].version;
    }
    for (const Piece& pc : pieces) {
        const FrontElement old = elements_[pc.old];
        const std::uint32_t id = add_element(old.site, old.loop);
        fresh.push_back(id);
        if (pc.start < 0) {
            elements_[id].prev = old.prev;
            vertices_[old.prev].right = id;
        } else {
            port(id, false, pc.start);
        }
        if (pc.end < 0) {
            elements_[id].next = old.next;
            vertices_[old.next].left = id;
        } else {
            port(id, true, pc.end);
        }
    }

    // 3. Sort the ports counter-clockwise around D. A point cluster starts in the widest gap
    // between port directions; a plateau end starts at the direction of the other end.
    if (!plateau && !ports.empty()) {
        std::vector<double> a;
        for (const Port& p : ports) a.push_back(wrap_2pi(p.angle));
        std::sort(a.begin(), a.end());
        double widest = -1.0, cut = 0.0;
        for (std::size_t i = 0; i < a.size(); ++i) {
            const double next = i + 1 < a.size() ? a[i + 1] : a[0] + two_pi;
            if (next - a[i] > widest) {
                widest = next - a[i];
                cut = a[i] + 0.5 * widest;
            }
        }
        anchors[0].cut = cut;
    }
    for (Port& p : ports) p.rel = wrap_2pi(p.angle - anchors[p.anchor].cut);
    // Anchor 1 (the end of a plateau) comes first: the boundary of D runs counter-clockwise from
    // the right side of the plateau round its end, back along the left side and round its start.
    std::sort(ports.begin(), ports.end(), [](const Port& a, const Port& b) {
        if (a.anchor != b.anchor) return a.anchor > b.anchor;
        return a.rel != b.rel ? a.rel < b.rel : a.element < b.element;
    });
    // Directions equal within eps_ang (tangent fronts) are ordered by curvature. The direction of a
    // curved front of radius r is only known to within spread / r, since the anchor stands for every
    // point of the cluster: a tangency a cluster radius from a vertex is still a tangency [RB-03].
    for (std::size_t i = 0; i < ports.size();) {
        std::size_t j = i + 1;
        while (j < ports.size() && ports[j].anchor == ports[i].anchor &&
               ports[j].rel - ports[j - 1].rel <= tol_.ang + ports[j].slack + ports[j - 1].slack)
            ++j;
        std::sort(ports.begin() + static_cast<std::ptrdiff_t>(i), ports.begin() + static_cast<std::ptrdiff_t>(j),
                  [](const Port& a, const Port& b) { return a.k != b.k ? a.k < b.k : a.element < b.element; });
        i = j;
    }

    // Pair each in-port with the out-port before it.
    const std::size_t np = ports.size();
    std::vector<std::pair<std::size_t, std::size_t>> links;
    for (std::size_t i = 0; i < np; ++i) {
        if (!ports[i].in) continue;
        const std::size_t j = (i + np - 1) % np;
        if (ports[j].in || ports[j].anchor != ports[i].anchor) return failure("front ports that do not alternate");
        if (elements_[ports[i].element].site == elements_[ports[j].element].site)
            return failure("a front element relinked to itself");
        links.emplace_back(i, j);
        ++anchors[ports[i].anchor].links;
    }
    if (2 * links.size() != np) return failure("unpaired front ports");
    if (np == 0)
        for (std::uint32_t x = start;;) {
            if (elements_[x].alive) return failure("a vanishing loop with surviving elements");
            x = next_element(x);
            if (x == start) break;
        }

    // 4. Emit the MAT vertices, close the dying shocks, add the plateau edge and the new shocks.
    for (Anchor& a : anchors) {
        const int degree = a.shocks_in + a.links + (plateau ? 1 : 0);
        VertexKind kind = VertexKind::curvature_end;
        if (degree >= 3)
            kind = VertexKind::junction;
        else if (degree == 2)
            kind = a.shocks_in == 2 ? VertexKind::extremum_max
                   : a.links == 2   ? VertexKind::extremum_min
                                    : VertexKind::transition;
        std::vector<SiteId> contacts;
        for (SiteId id : a.contacts) {
            const Site& s = sites[id];
            const bool centred = !s.is_line() && s.sigma > 0 && dist(s.centre, a.p) <= centre_slack;
            if (!centred && std::find(contacts.begin(), contacts.end(), id) == contacts.end()) contacts.push_back(id);
        }
        a.mv = add_mat_vertex_at_sites(a.p, t, kind, contacts);
    }
    for (std::size_t i = 0; i < verts.size(); ++i) {
        if (!dying[i]) continue;
        close_edge(verts[i], anchors[vanchor[i]].mv, anchors[vanchor[i]].p, t);
        vertices_[verts[i]].alive = false;
    }
    if (plateau) {
        const Site &sx = site(px), &sy = site(py);
        MatEdge::Data d;
        d.v0 = anchors[0].mv;
        d.v1 = anchors[1].mv;
        // Along the plateau, in the direction of the front of x, x lies on the right.
        d.left = elements_[py].site;
        d.right = elements_[px].site;
        d.kind = sx.is_line() ? ConicKind::plateau_line : ConicKind::plateau_arc;
        d.branch = sx.is_line() ? 0 : sx.sigma;
        d.a = sy;
        d.b = sx;
        d.t0 = d.t1 = t;
        d.p0 = anchors[0].p;
        d.p1 = anchors[1].p;
        d.xf = region_.xf;
        mat_edges.push_back(d);
    }

    std::vector<std::uint32_t> born;
    for (const auto& [i, j] : links) {
        const Port &in = ports[i], &out = ports[j];
        const Anchor& a = anchors[in.anchor];
        // Tangent ports ordered by curvature can be measured the other way round.
        double sweep = wrap_2pi(in.angle - out.angle);
        if (sweep > two_pi - (tol_.ang + in.slack + out.slack)) sweep = 0.0;
        auto v = add_shock(in.element, out.element, a.p, t, a.mv, polar(out.angle + 0.5 * sweep));
        if (!v) return tl::unexpected(v.error());
        elements_[in.element].next = *v;
        elements_[out.element].prev = *v;
        born.push_back(*v);
    }
    std::vector<std::uint32_t> starts;
    for (std::uint32_t v : born) starts.push_back(vertices_[v].left);
    relabel_loops(starts);

    ++stats.clusters;
    if (np == 0) ++stats.annihilations;
    for (std::uint32_t e : fresh) reschedule_element(e);
    for (std::uint32_t v : born) reschedule_vertex(v);
    return {};
}

} // namespace wfmat::detail
