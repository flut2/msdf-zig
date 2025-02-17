const std = @import("std");

const EdgeColor = @import("edge_color.zig").EdgeColor;
const equations = @import("equations.zig");
const f64i = @import("Generator.zig").f64i;
const math = @import("math.zig");
const mix = math.mix;
const median = math.median;
const Shape = @import("Shape.zig");
const Vec2 = @import("Vec2.zig");

const ErrorCorrection = @This();

const red: u32 = @intFromEnum(EdgeColor.red);
const green: u32 = @intFromEnum(EdgeColor.green);
const blue: u32 = @intFromEnum(EdgeColor.blue);

const artifact_t_epsilon = 0.01;
const protection_radius_tolerance = 1.001;

pub const Mode = enum { indiscriminate, edge_priority, edge_only };
pub const DistanceMode = enum { none, at_edge, always };
pub const Options = struct {
    mode: Mode = .edge_priority,
    /// The distance checking will be forcefully turned off (set to ``.none``) if the scanline pass is enabled
    distance_mode: DistanceMode = .none,
    min_deviation_ratio: f64 = 1.11111111111111111,
    min_improve_ratio: f64 = 1.11111111111111111,
};

const StencilFlag = packed struct {
    err: bool = false,
    protected: bool = false,
};

const ClassifierFlag = packed struct {
    candidate: bool = false,
    artifact: bool = false,
};

options: Options = .{},
stencil: []StencilFlag = &.{},
stencil_w: u16 = 0,
stencil_h: u16 = 0,

pub fn create(allocator: std.mem.Allocator, w: u16, h: u16, options: Options) !ErrorCorrection {
    if (options.distance_mode != .none) @panic("The given distance mode is currently not implemented");
    return .{
        .options = options,
        .stencil = try allocator.alloc(StencilFlag, w * h),
        .stencil_w = w,
        .stencil_h = h,
    };
}

pub fn destroy(self: *ErrorCorrection, allocator: std.mem.Allocator) void {
    allocator.free(self.stencil);
}

pub fn applyProtections(
    self: *ErrorCorrection,
    shape: *Shape,
    scale: f64,
    px_range: f64,
    tx: f64,
    ty: f64,
    sdf_px: []f64,
    sdf_w: u16,
    sdf_h: u16,
    channels: u8,
) void {
    switch (self.options.mode) {
        .edge_priority => {
            self.protectCorners(shape, scale, tx, ty);
            self.protectEdges(scale, px_range, sdf_px, sdf_w, sdf_h, channels);
        },
        .edge_only => self.protectAll(),
        .indiscriminate => {},
    }

    self.findErrors(scale, px_range, sdf_px, sdf_w, sdf_h, channels);
    self.apply(sdf_px, sdf_w, sdf_h, channels);
}

pub fn protectCorners(self: *ErrorCorrection, shape: *Shape, scale: f64, tx: f64, ty: f64) void {
    for (shape.contours.items) |contour| {
        if (contour.edges.items.len == 0) continue;
        var prev_edge = contour.edges.getLast();
        for (contour.edges.items) |edge| {
            const common_color = @intFromEnum(prev_edge.color) & @intFromEnum(edge.color);
            prev_edge = edge;
            if ((common_color & (common_color - 1)) == 0) continue;
            const base_point = edge.point(0).mul(scale).add(Vec2{ .x = tx, .y = ty });
            const left: i32 = @intFromFloat(base_point.x - 0.5);
            const bottom: i32 = @intFromFloat(f64i(self.stencil_h) - base_point.y - 0.5 - 2.0);
            const right = left + 1;
            const top = bottom + 1;
            if (left < self.stencil_w and bottom < self.stencil_h and right >= 0 and top >= 0) {
                if (left >= 0 and bottom >= 0)
                    self.stencil[index(@intCast(left), @intCast(bottom), self.stencil_w, 1)].protected = true;
                if (right < self.stencil_w and bottom >= 0)
                    self.stencil[index(@intCast(right), @intCast(bottom), self.stencil_w, 1)].protected = true;
                if (left >= 0 and top < self.stencil_h)
                    self.stencil[index(@intCast(left), @intCast(top), self.stencil_w, 1)].protected = true;
                if (right < self.stencil_w and top < self.stencil_h)
                    self.stencil[index(@intCast(right), @intCast(top), self.stencil_w, 1)].protected = true;
            }
        }
    }
}

fn index(x: usize, y: usize, w: usize, channels: usize) usize {
    return y * w * channels + x * channels;
}

fn edgeBetweenTexelsChannel(a: *const [3]f64, b: *const [3]f64, channel: u8) bool {
    const delta = a[channel] - b[channel];
    if (delta == 0.0) return false;
    const t = (a[channel] - 0.5) / delta;
    if (t > 0 and t < 1) {
        const c: [3]f64 = .{
            mix(a[0], b[0], t),
            mix(a[1], b[1], t),
            mix(a[2], b[2], t),
        };
        return median(c[0], c[1], c[2]) == c[channel];
    }
    return false;
}

fn edgeBetweenTexels(a: *const [3]f64, b: *const [3]f64) u32 {
    return red * @intFromBool(edgeBetweenTexelsChannel(a, b, 0)) +
        green * @intFromBool(edgeBetweenTexelsChannel(a, b, 1)) +
        blue * @intFromBool(edgeBetweenTexelsChannel(a, b, 2));
}

fn protectExtremeChannels(stencil_point: *StencilFlag, msd: *const [3]f64, m: f64, mask: u32) void {
    if (mask & red != 0 and msd[0] != m or
        mask & green != 0 and msd[1] != m or
        mask & blue != 0 and msd[2] != m)
        stencil_point.protected = true;
}

fn msdfMedian(msdf: *const [3]f64) f64 {
    return median(msdf[0], msdf[1], msdf[2]);
}

fn protectEdges(
    self: *ErrorCorrection,
    scale: f64,
    px_range: f64,
    sdf_px: []const f64,
    sdf_w: u16,
    sdf_h: u16,
    channels: u8,
) void {
    const hori_radius = (Vec2{ .x = px_range, .y = 0 }).div(scale).length() * protection_radius_tolerance;
    for (0..sdf_h) |y| for (0..sdf_w - 1) |x| {
        const left = sdf_px[index(0, y, sdf_w, channels)..][0..3];
        const right = sdf_px[index(1, y, sdf_w, channels)..][0..3];
        const median_left = msdfMedian(left);
        const median_right = msdfMedian(right);
        if (@abs(median_left - 0.5) + @abs(median_right - 0.5) < hori_radius) {
            const mask = edgeBetweenTexels(left, right);
            protectExtremeChannels(&self.stencil[index(x, y, self.stencil_w, 1)], left, median_left, mask);
            protectExtremeChannels(&self.stencil[index(x + 1, y, self.stencil_w, 1)], right, median_right, mask);
        }
    };

    const vert_radius = (Vec2{ .x = 0, .y = px_range }).div(scale).length() * protection_radius_tolerance;
    for (0..sdf_h - 1) |y| for (0..sdf_w) |x| {
        const bottom = sdf_px[index(0, y, sdf_w, channels)..][0..3];
        const top = sdf_px[index(0, y + 1, sdf_w, channels)..][0..3];
        const median_bottom = msdfMedian(bottom);
        const median_top = msdfMedian(top);
        if (@abs(median_bottom - 0.5) + @abs(median_top - 0.5) < vert_radius) {
            const mask = edgeBetweenTexels(bottom, top);
            protectExtremeChannels(&self.stencil[index(x, y, self.stencil_w, 1)], bottom, median_bottom, mask);
            protectExtremeChannels(&self.stencil[index(x, y + 1, self.stencil_w, 1)], top, median_top, mask);
        }
    };

    const diag_radius = (Vec2{ .x = px_range, .y = px_range }).div(scale).length() * protection_radius_tolerance;
    for (0..sdf_h - 1) |y| for (0..sdf_w - 1) |x| {
        const bottom_left = sdf_px[index(0, y, sdf_w, channels)..][0..3];
        const bottom_right = sdf_px[index(1, y, sdf_w, channels)..][0..3];
        const top_left = sdf_px[index(0, y + 1, sdf_w, channels)..][0..3];
        const top_right = sdf_px[index(1, y + 1, sdf_w, channels)..][0..3];
        const median_bottom_left = msdfMedian(bottom_left);
        const median_bottom_right = msdfMedian(bottom_right);
        const median_top_left = msdfMedian(top_left);
        const median_top_right = msdfMedian(top_right);
        if (@abs(median_bottom_left - 0.5) + @abs(median_top_right - 0.5) < diag_radius) {
            const mask = edgeBetweenTexels(bottom_left, top_right);
            protectExtremeChannels(&self.stencil[index(x, y, self.stencil_w, 1)], bottom_left, median_bottom_left, mask);
            protectExtremeChannels(&self.stencil[index(x + 1, y + 1, self.stencil_w, 1)], top_right, median_top_right, mask);
        }
        if (@abs(median_bottom_right - 0.5) + @abs(median_top_left - 0.5) < diag_radius) {
            const mask = edgeBetweenTexels(bottom_right, top_left);
            protectExtremeChannels(&self.stencil[index(x + 1, y, self.stencil_w, 1)], bottom_right, median_bottom_right, mask);
            protectExtremeChannels(&self.stencil[index(x, y + 1, self.stencil_w, 1)], top_left, median_top_left, mask);
        }
    };
}

fn protectAll(self: *ErrorCorrection) void {
    for (self.stencil) |*mask| mask.protected = true;
}

fn interpolatedMedianLinear(a: *const [3]f64, b: *const [3]f64, t: f64) f64 {
    return median(
        mix(a[0], b[0], t),
        mix(a[1], b[1], t),
        mix(a[2], b[2], t),
    );
}

fn interpolatedMedianBilinear(a: *const [3]f64, l: *const [3]f64, q: *const [3]f64, t: f64) f64 {
    return median(
        t * (t * q[0] + l[0]) + a[0],
        t * (t * q[1] + l[1]) + a[1],
        t * (t * q[2] + l[2]) + a[2],
    );
}

fn rangeTest(span: f64, protected: bool, at: f64, bt: f64, xt: f64, am: f64, bm: f64, xm: f64) ClassifierFlag {
    if ((am > 0.5 and bm > 0.5 and xm <= 0.5) or
        (am < 0.5 and bm < 0.5 and xm >= 0.5) or
        (!protected and median(am, bm, xm) != xm))
    {
        const ax_span = (xt - at) * span;
        const bx_span = (bt - xt) * span;
        if (!(xm >= am - ax_span and xm <= am + ax_span and xm >= bm - bx_span and xm <= bm + bx_span))
            return .{ .candidate = true, .artifact = true };
        return .{ .candidate = true };
    }
    return .{};
}

fn isArtifact(is_protected: bool, ax_span: f64, bx_span: f64, median_a: f64, median_b: f64, median_x: f64) bool {
    // zig fmt: off
    return 
        (!is_protected or
        (median_a > 0.5 and median_b > 0.5 and median_x <= 0.5) or
        (median_a < 0.5 and median_b < 0.5 and median_x >= 0.5)) and
        !(median_x >= median_a - ax_span and
        median_x <= median_a + ax_span and
        median_x >= median_b - bx_span and
        median_x <= median_b + bx_span);
    // zig fmt: on
}

fn hasLinearArtifactInner(span: f64, protected: bool, am: f64, bm: f64, a: *const [3]f64, b: *const [3]f64, da: f64, db: f64) bool {
    const delta = da - db;
    if (delta == 0) return false;
    const t = da / (da - db);
    if (t > artifact_t_epsilon and t < 1 - artifact_t_epsilon) {
        const xm = interpolatedMedianLinear(a, b, t);
        return rangeTest(span, protected, 0, 1, t, am, bm, xm).artifact;
    }
    return false;
}

fn hasDiagonalArtifactInner(
    span: f64,
    protected: bool,
    am: f64,
    dm: f64,
    a: *const [3]f64,
    l: *const [3]f64,
    q: *const [3]f64,
    d_a: f64,
    d_bc: f64,
    d_d: f64,
    t_ex_0: f64,
    t_ex_1: f64,
) bool {
    var t: [2]f64 = undefined;
    const solutions = equations.solveQuadratic(&t, d_d - d_bc + d_a, d_bc - d_a - d_a, d_a);
    for (0..solutions) |i| if (t[i] > artifact_t_epsilon and t[i] < 1 - artifact_t_epsilon) {
        const xm = interpolatedMedianBilinear(a, l, q, t[i]);
        if (rangeTest(span, protected, 0, 1, t[i], am, dm, xm).artifact) return true;
        var t_end: [2]f64 = undefined;
        var em: [2]f64 = undefined;

        if (t_ex_0 > 0 and t_ex_0 < 1) {
            t_end[0] = 0;
            t_end[1] = 1;
            em[0] = am;
            em[1] = dm;
            t_end[@intFromBool(t_ex_0 > t[i])] = t_ex_0;
            em[@intFromBool(t_ex_0 > t[i])] = interpolatedMedianBilinear(a, l, q, t_ex_0);
            if (rangeTest(span, protected, t_end[0], t_end[1], t[i], em[0], em[1], xm).artifact) return true;
        }

        if (t_ex_1 > 0 and t_ex_1 < 1) {
            t_end[0] = 0;
            t_end[1] = 1;
            em[0] = am;
            em[1] = dm;
            t_end[@intFromBool(t_ex_1 > t[i])] = t_ex_1;
            em[@intFromBool(t_ex_1 > t[i])] = interpolatedMedianBilinear(a, l, q, t_ex_1);
            if (rangeTest(span, protected, t_end[0], t_end[1], t[i], em[0], em[1], xm).artifact) return true;
        }
    };
    return false;
}

fn hasLinearArtifact(span: f64, protected: bool, am: f64, a: *const [3]f64, b: *const [3]f64) bool {
    const bm = median(b[0], b[1], b[2]);
    return (@abs(am - 0.5) >= @abs(bm - 0.5) and (hasLinearArtifactInner(span, protected, am, bm, a, b, a[1] - a[0], b[1] - b[0]) or
        hasLinearArtifactInner(span, protected, am, bm, a, b, a[2] - a[1], b[2] - b[1]) or
        hasLinearArtifactInner(span, protected, am, bm, a, b, a[0] - a[2], b[0] - b[2])));
}

fn hasDiagonalArtifact(span: f64, protected: bool, am: f64, a: *const [3]f64, b: *const [3]f64, c: *const [3]f64, d: *const [3]f64) bool {
    const dm = median(d[0], d[1], d[2]);
    if (@abs(am - 0.5) >= @abs(dm - 0.5)) {
        const abc: [3]f64 = .{
            a[0] - b[0] - c[0],
            a[1] - b[1] - c[1],
            a[2] - b[2] - c[2],
        };
        const l: [3]f64 = .{
            -a[0] - abc[0],
            -a[1] - abc[1],
            -a[2] - abc[2],
        };
        const q: [3]f64 = .{
            d[0] + abc[0],
            d[1] + abc[1],
            d[2] + abc[2],
        };

        if (q[0] == 0.0 or q[1] == 0.0 or q[2] == 0.0) return false;

        const t_ex: [3]f64 = .{
            -0.5 * l[0] / q[0],
            -0.5 * l[1] / q[1],
            -0.5 * l[2] / q[2],
        };
        return hasDiagonalArtifactInner(
            span,
            protected,
            am,
            dm,
            a,
            &l,
            &q,
            a[1] - a[0],
            b[1] - b[0] + c[1] - c[0],
            d[1] - d[0],
            t_ex[0],
            t_ex[1],
        ) or hasDiagonalArtifactInner(
            span,
            protected,
            am,
            dm,
            a,
            &l,
            &q,
            a[2] - a[1],
            b[2] - b[1] + c[2] - c[1],
            d[2] - d[1],
            t_ex[1],
            t_ex[2],
        ) or hasDiagonalArtifactInner(
            span,
            protected,
            am,
            dm,
            a,
            &l,
            &q,
            a[0] - a[2],
            b[0] - b[2] + c[0] - c[2],
            d[0] - d[2],
            t_ex[2],
            t_ex[0],
        );
    }
    return false;
}

fn findErrors(
    self: *ErrorCorrection,
    scale: f64,
    px_range: f64,
    sdf_px: []const f64,
    sdf_w: u16,
    sdf_h: u16,
    channels: u8,
) void {
    const min_deviation_ratio = self.options.min_deviation_ratio;
    const hori_span = (Vec2{ .x = px_range, .y = 0 }).div(scale).length() * min_deviation_ratio;
    const vert_span = (Vec2{ .x = 0, .y = px_range }).div(scale).length() * min_deviation_ratio;
    const diag_span = (Vec2{ .x = px_range, .y = px_range }).div(scale).length() * min_deviation_ratio;

    for (0..sdf_h) |y| for (0..sdf_w) |x| {
        const current = sdf_px[index(x, y, sdf_w, channels)..][0..3];
        const median_current = msdfMedian(current);
        var current_stencil = &self.stencil[index(x, y, self.stencil_w, 1)];
        const is_protected = current_stencil.protected;

        if (x > 0) {
            const left = sdf_px[index(x - 1, y, sdf_w, channels)..][0..3];
            if (hasLinearArtifact(hori_span, is_protected, median_current, current, left)) {
                current_stencil.err = true;
                continue;
            }

            if (y > 0) {
                const top = sdf_px[index(x, y - 1, sdf_w, channels)..][0..3];
                const top_left = sdf_px[index(x - 1, y - 1, sdf_w, channels)..][0..3];
                if (hasDiagonalArtifact(diag_span, is_protected, median_current, current, left, top, top_left)) {
                    current_stencil.err = true;
                    continue;
                }
            }

            if (y < sdf_h - 1) {
                const bottom = sdf_px[index(x, y + 1, sdf_w, channels)..][0..3];
                const bottom_left = sdf_px[index(x - 1, y + 1, sdf_w, channels)..][0..3];
                if (hasDiagonalArtifact(diag_span, is_protected, median_current, current, left, bottom, bottom_left)) {
                    current_stencil.err = true;
                    continue;
                }
            }
        }

        if (y > 0) {
            const top = sdf_px[index(x, y - 1, sdf_w, channels)..][0..3];
            if (hasLinearArtifact(vert_span, is_protected, median_current, current, top)) {
                current_stencil.err = true;
                continue;
            }
        }

        if (x < sdf_w - 1) {
            const right = sdf_px[index(x + 1, y, sdf_w, channels)..][0..3];
            if (hasLinearArtifact(hori_span, is_protected, median_current, current, right)) {
                current_stencil.err = true;
                continue;
            }

            if (y > 0) {
                const top = sdf_px[index(x, y - 1, sdf_w, channels)..][0..3];
                const top_right = sdf_px[index(x + 1, y - 1, sdf_w, channels)..][0..3];
                if (hasDiagonalArtifact(diag_span, is_protected, median_current, current, right, top, top_right)) {
                    current_stencil.err = true;
                    continue;
                }
            }

            if (y < sdf_h - 1) {
                const bottom = sdf_px[index(x, y + 1, sdf_w, channels)..][0..3];
                const bottom_right = sdf_px[index(x + 1, y + 1, sdf_w, channels)..][0..3];
                if (hasDiagonalArtifact(diag_span, is_protected, median_current, current, right, bottom, bottom_right)) {
                    current_stencil.err = true;
                    continue;
                }
            }
        }

        if (y < sdf_h - 1) {
            const bottom = sdf_px[index(x, y + 1, sdf_w, channels)..][0..3];
            if (hasLinearArtifact(vert_span, is_protected, median_current, current, bottom)) {
                current_stencil.err = true;
                continue;
            }
        }
    };
}

fn apply(self: *ErrorCorrection, sdf_px: []f64, sdf_w: u16, sdf_h: u16, channels: u8) void {
    for (0..sdf_w * sdf_h) |i| if (self.stencil[i].err) {
        const msdf = sdf_px[i * channels ..][0..3];
        const msdf_median = msdfMedian(msdf);
        @memset(msdf, msdf_median);
    };
}
