// The example of docs/api.md: the medial axis of a filleted slot-like plate, its edges as Bezier
// arcs and polylines, and an inward offset. Built with the tools (WFMAT_BUILD_TOOLS).
#include <cstdio>
#include <vector>

#include "wfmat/wfmat.hpp"

int main()
{
    using namespace wfmat;

    // 100 x 40 with a half-disk on the right end: vertex i starts segment i; bulge = tan(sweep / 4),
    // so 1 is a counter-clockwise half circle.
    Region plate;
    plate.outer.vertices = {{{0, 0}, 0.0}, {{100, 0}, 1.0}, {{100, 40}, 0.0}, {{0, 40}, 0.0}};

    auto mat = compute_mat(plate);
    if (!mat) {
        std::fprintf(stderr, "%s\n", to_string(mat.error()).c_str());
        return 1;
    }
    std::printf("wfmat %s: %zu vertices, %zu edges\n", version, mat->vertices().size(), mat->edges().size());

    for (const MatEdge& e : mat->edges()) {
        const MatVertex& a = mat->vertices()[e.v0()];
        const MatVertex& b = mat->vertices()[e.v1()];
        std::printf("edge (%.3f, %.3f) r=%.3f -> (%.3f, %.3f) r=%.3f\n", a.p.x, a.p.y, e.r0(), b.p.x, b.p.y, e.r1());

        for (const RationalQuadBezier& arc : e.bezier())
            std::printf("  bezier p1 (%.3f, %.3f) w %.6f over u in [%.3f, %.3f]\n", arc.p1.x, arc.p1.y, arc.w, arc.u0,
                        arc.u1);

        std::vector<Vec2> polyline;
        e.flatten(0.01, polyline);  // chord tolerance in the input's units
        std::printf("  %zu polyline points; radius at the middle %.3f\n", polyline.size(), e.radius_at(0.5));
    }

    auto offset = inward_offset(plate, 5.0);
    if (!offset) {
        std::fprintf(stderr, "%s\n", to_string(offset.error()).c_str());
        return 1;
    }
    for (const Loop& loop : *offset) {
        std::printf("offset loop:");
        for (const BulgeVertex& v : loop.vertices) std::printf(" (%.3f, %.3f; %.3f)", v.p.x, v.p.y, v.bulge);
        std::printf("\n");
    }
    return 0;
}
