# Maple cross-check for tests/data/kernel/three_site.json [VER-01].
#
# Solves the squared site equations of section 7 exactly and keeps the real roots that satisfy
# the unsquared equations (t >= 0, non-negative circle radii), as gen_kernel_fixtures.py does
# with SymPy. Copy a case's sites from the JSON file into a call of ThreeSite and compare the
# printed roots with the case's "roots".
#
# A site is [line, nx, ny, c] or [circle, cx, cy, R, sigma]; numbers are converted exactly
# from their decimal form, as the C++ test reads them as doubles.
#
# This script was written without access to Maple and has not been run.

SiteEquation := proc(s)
    local nx, ny, c, cx, cy, R, sg;
    if s[1] = line then
        nx, ny, c := op(map(convert, s[2..4], rational, exact));
        return nx*x + ny*y - c - t;
    else
        cx, cy, R := op(map(convert, s[2..4], rational, exact));
        sg := s[5];
        return x^2 + y^2 - t^2 - 2*cx*x - 2*cy*y + 2*sg*R*t + cx^2 + cy^2 - R^2;
    end if;
end proc:

UnsquaredOK := proc(s, X, Y, Tm)
    local cx, cy, R, sg, radius;
    if s[1] = line then return true; end if;
    cx, cy, R := op(map(convert, s[2..4], rational, exact));
    sg := s[5];
    radius := R - sg*Tm;
    return evalb(radius >= -1e-30 and abs(sqrt((X - cx)^2 + (Y - cy)^2) - radius) < 1e-30);
end proc:

ThreeSite := proc(sites::list)
    local eqs, sols, sol, X, Y, Tm, roots;
    Digits := 50;
    eqs := map(SiteEquation, sites);
    sols := [solve(eqs, {x, y, t}, explicit)];
    roots := [];
    for sol in sols do
        X := evalf(eval(x, sol)); Y := evalf(eval(y, sol)); Tm := evalf(eval(t, sol));
        if abs(Im(X)) + abs(Im(Y)) + abs(Im(Tm)) > 1e-40 then next; end if;
        X, Y, Tm := Re(X), Re(Y), Re(Tm);
        if Tm < 0 then next; end if;
        if andmap(s -> UnsquaredOK(s, X, Y, Tm), sites) then
            roots := [op(roots), [X, Y, Tm]];
        end if;
    end do;
    return sort(roots, (a, b) -> evalb(a[3] < b[3]));
end proc:

# Example: the "facing-lines" case; expected root x = -2, y = 1, t = 1.
ThreeSite([[line, 0.0, 1.0, 0.0], [line, 0.0, -1.0, -2.0], [line, 1.0, 0.0, -3.0]]);
