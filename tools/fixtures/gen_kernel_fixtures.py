#!/usr/bin/env python3
"""Generate tests/data/kernel/three_site.json: reference roots of the three-site solve [VER-01].

Every site parameter is a double; the oracle converts it to an exact rational and solves the
squared site equations exactly with SymPy, then keeps the real roots that satisfy the unsquared
equations (t >= 0, non-negative circle radii), evaluated to 40 digits. kernel_fixtures.mpl does
the same in Maple for cross-checking.

Usage: python3 tools/fixtures/gen_kernel_fixtures.py > tests/data/kernel/three_site.json
"""
import itertools
import json
import math
import random
import sys
from fractions import Fraction

import sympy as sp

X, Y, T = sp.symbols("x y t", real=True)
DIGITS = 40


def rat(v):
    return sp.Rational(Fraction(v))


def line(nx, ny, c):
    return {"kind": "line", "n": [nx, ny], "c": c}


def circle(cx, cy, R, sigma):
    return {"kind": "circle", "centre": [cx, cy], "R": R, "sigma": sigma}


def equation(s):
    if s["kind"] == "line":
        nx, ny = map(rat, s["n"])
        return nx * X + ny * Y - rat(s["c"]) - T
    cx, cy = map(rat, s["centre"])
    R, sg = rat(s["R"]), s["sigma"]
    return X**2 + Y**2 - T**2 - 2 * cx * X - 2 * cy * Y + 2 * sg * R * T + cx**2 + cy**2 - R**2


def unsquared_ok(s, x, y, t):
    if s["kind"] == "line":
        return True
    cx, cy = map(rat, s["centre"])
    R, sg = rat(s["R"]), s["sigma"]
    radius = R - sg * t
    if radius < -sp.Rational(1, 10**30):
        return False
    return abs(sp.sqrt((x - cx) ** 2 + (y - cy) ** 2) - radius) < sp.Rational(1, 10**30)


def reference_roots(sites):
    eqs = [equation(s) for s in sites]
    sols = sp.solve(eqs, [X, Y, T], dict=True)
    roots = []
    for sol in sols:
        if X not in sol or Y not in sol or T not in sol:
            return None  # a continuum of solutions: degenerate
        vals = [sp.N(sol[v], DIGITS + 10) for v in (X, Y, T)]
        if any(abs(sp.im(v)) > 10 ** (-DIGITS) for v in vals):
            continue
        x, y, t = (sp.re(v) for v in vals)
        if t < 0:
            continue
        if all(unsquared_ok(s, x, y, t) for s in sites):
            roots.append((x, y, t))
    roots.sort(key=lambda r: (r[2], r[0], r[1]))
    return roots


def random_site(kind, rng, x0, y0, t0):
    """A site of the given kind at distance t0 from (x0, y0), so a root exists by construction."""
    a = rng.uniform(0, 2 * math.pi)
    u = (math.cos(a), math.sin(a))
    if kind == "L":  # normal u points from the line towards (x0, y0)
        return line(u[0], u[1], u[0] * x0 + u[1] * y0 - t0)
    if kind == "X":  # convex: (x0, y0) inside, tangent internally
        rho = rng.uniform(0.05, 0.5)
        return circle(x0 + rho * u[0], y0 + rho * u[1], rho + t0, 1)
    if kind == "V":  # concave: outside, at distance t0 from the arc
        R = rng.uniform(0.05, 0.3)
        return circle(x0 + (R + t0) * u[0], y0 + (R + t0) * u[1], R, -1)
    return circle(x0 + t0 * u[0], y0 + t0 * u[1], 0.0, -1)  # point


def as_json_roots(roots):
    return [{"x": float(x), "y": float(y), "t": float(t)} for x, y, t in roots]


def main():
    rng = random.Random(20261009)
    cases = []
    for combo in itertools.combinations_with_replacement("LXVP", 3):
        for k in range(3):
            x0, y0, t0 = rng.uniform(-0.5, 0.5), rng.uniform(-0.5, 0.5), rng.uniform(0.05, 0.4)
            sites = [random_site(kind, rng, x0, y0, t0) for kind in combo]
            roots = reference_roots(sites)
            if roots is None:
                continue
            cases.append({"name": f"{''.join(combo)}-{k}", "sites": sites, "t_now": 0.0,
                          "degenerate": False, "roots": as_json_roots(roots)})

    # Rank-deficient systems [KN-02].
    cases.append({"name": "equal-normals", "t_now": 0.0, "degenerate": True, "roots": [],
                  "sites": [line(0.0, 1.0, 0.0), line(0.0, 1.0, 1.0), line(1.0, 0.0, -3.0)]})
    # Facing parallel lines are not rank-deficient with a third site: the root is on the midline.
    facing = [line(0.0, 1.0, 0.0), line(0.0, -1.0, -2.0), line(1.0, 0.0, -3.0)]
    cases.append({"name": "facing-lines", "t_now": 0.0, "degenerate": False, "sites": facing,
                  "roots": as_json_roots(reference_roots(facing))})
    cases.append({"name": "three-concentric-circles", "t_now": 0.0, "degenerate": True, "roots": [],
                  "sites": [circle(0.0, 0.0, 2.0, 1), circle(0.0, 0.0, 1.0, -1), circle(0.0, 0.0, 0.5, -1)]})
    ring = [circle(0.0, 0.0, 2.0, 1), circle(0.0, 0.0, 1.0, -1), line(1.0, 0.0, -1.25)]
    cases.append({"name": "concentric-circles-and-line", "t_now": 0.0, "degenerate": False, "sites": ring,
                  "roots": as_json_roots(reference_roots(ring))})

    # Near-tangent family [KN-06]: line y = 0, concave circle about (0, 2) of radius 1/2, and a point
    # (3/10, u). At u = u* the two roots merge; offsets of 1e-6 and 1e-10 give ill-conditioned pairs.
    U = sp.symbols("u", real=True)
    fam = [line(0.0, 1.0, 0.0), circle(0.0, 2.0, 0.5, -1)]
    eqs = [equation(s) for s in fam] + [(X - sp.Rational(3, 10)) ** 2 + (Y - U) ** 2 - T**2]
    # Eliminate: the first two give y and t in terms of x; the third is then a quadratic in x.
    sub = sp.solve(eqs[:2], [Y, T], dict=True)
    for branch in sub:
        quad = sp.expand(eqs[2].subs(branch))
        disc = sp.discriminant(sp.Poly(quad, X).as_expr(), X)
        for ustar in sp.solve(disc, U):
            if not ustar.is_real:
                continue
            for off in (1e-6, 1e-10):
                for sgn in (1, -1):
                    u = float(sp.N(ustar)) + sgn * off
                    sites = fam + [circle(0.3, u, 0.0, -1)]
                    roots = reference_roots(sites)
                    if roots and len(roots) >= 1:
                        cases.append({"name": f"near-tangent-u{float(sp.N(ustar)):+.3f}{sgn * off:+.0e}",
                                      "sites": sites, "t_now": 0.0, "degenerate": False,
                                      "roots": as_json_roots(roots)})
    json.dump({"generator": "tools/fixtures/gen_kernel_fixtures.py", "digits": DIGITS, "cases": cases},
              sys.stdout, indent=1)
    sys.stdout.write("\n")


if __name__ == "__main__":
    main()
