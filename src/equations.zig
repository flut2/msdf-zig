const std = @import("std");

pub fn solveQuadratic(roots: *[2]f64, a: f64, b: f64, c: f64) u8 {
    if (a == 0 or @abs(b) > 1e12 * @abs(a)) {
        if (b == 0) return 0;
        roots[0] = -c / b;
        return 1;
    }
    const dscr = b * b - 4.0 * a * c;
    if (dscr > 0) {
        const dscr_sqrt = @sqrt(dscr);
        roots[0] = (-b + dscr_sqrt) / (2 * a);
        roots[1] = (-b - dscr_sqrt) / (2 * a);
        return 2;
    } else if (dscr == 0) {
        roots[0] = -b / (2 * a);
        return 1;
    } else return 0;
}

fn solveCubicNormed(roots: *[3]f64, a: f64, b: f64, c: f64) u8 {
    const a2 = a * a;
    var q = 1.0 / 9.0 * (a2 - 3 * b);
    const r = 1.0 / 54.0 * (a * (2 * a2 - 9 * b) + 27 * c);
    const r2 = r * r;
    const q3 = q * q * q;
    const one_third = 1.0 / 3.0;
    const mod_a = a * one_third;
    if (r2 < q3) {
        var t = r / @sqrt(q3);
        if (t < -1) t = -1;
        if (t > 1) t = 1;
        t = std.math.acos(t);
        q = -2 * @sqrt(q);
        roots[0] = q * @cos(one_third * t) - mod_a;
        roots[1] = q * @cos(one_third * (t + 2 * std.math.pi)) - mod_a;
        roots[2] = q * @cos(one_third * (t - 2 * std.math.pi)) - mod_a;
        return 3;
    } else {
        const u = @as(f64, (if (r < 0) 1.0 else -1.0)) * std.math.pow(f64, @abs(r) + @sqrt(r2 - q3), one_third);
        const v = if (u == 0) 0 else q / u;
        roots[0] = (u + v) - mod_a;
        if (u == v or @abs(u - v) < 1e-12 * @abs(u + v)) {
            roots[1] = -0.5 * (u + v) - mod_a;
            return 2;
        }
        return 1;
    }
}

pub fn solveCubic(x: *[3]f64, a: f64, b: f64, c: f64, d: f64) u8 {
    if (a != 0) {
        const bn = b / a;
        if (@abs(bn) < 1e6) return solveCubicNormed(x, bn, c / a, d / a);
    }

    return solveQuadratic(@ptrCast(x), b, c, d);
}
