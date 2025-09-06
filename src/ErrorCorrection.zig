const std = @import("std");

const EdgeColor = @import("edge_color.zig").EdgeColor;
const EdgeSegment = @import("EdgeSegment.zig");
const equations = @import("equations.zig");
const f64i = @import("Generator.zig").f64i;
const math = @import("math.zig");
const Shape = @import("Shape.zig");
const SignedDistance = @import("SignedDistance.zig");

const Vec2 = @Vector(2, f64);

const ErrorCorrection = @This();

const f64_nan = std.math.nan(f64);

const red = @intFromEnum(EdgeColor.red);
const green = @intFromEnum(EdgeColor.green);
const blue = @intFromEnum(EdgeColor.blue);

const artifact_t_epsilon = 0.01;
const protection_radius_tolerance = 1.001;

pub const Mode = enum { indiscriminate, edge_priority, edge_only };
pub const Options = struct {
    mode: Mode = .edge_priority,
    /// Will be forcefully turned off if the scanline pass is enabled
    check_distance: bool = true,
    min_deviation_ratio: f64 = 10.0 / 9.0,
    min_improve_ratio: f64 = 10.0 / 9.0,
};

const StencilFlags = packed struct {
    err: bool = false,
    protected: bool = false,
};

const ClassifierFlags = packed struct {
    candidate: bool = false,
    artifact: bool = false,
};

const DistanceEvaluationInfo = struct {
    shape: ?*const Shape = null,
    px_range: f64 = f64_nan,
    dir: Vec2 = @splat(f64_nan),
    sdf_coord: Vec2 = @splat(f64_nan),
    shape_coord: Vec2 = @splat(f64_nan),
    texel_size: Vec2 = @splat(f64_nan),
    sdf: ?*const [3]f64 = null,
    sdf_px: []const f64 = &.{},
    sdf_w: u16 = std.math.maxInt(u16),
    sdf_h: u16 = std.math.maxInt(u16),
    channels: u8 = std.math.maxInt(u8),
};

options: Options = .{},
dist_eval: DistanceEvaluationInfo,
stencil: []StencilFlags = &.{},
stencil_w: u16 = 0,
stencil_h: u16 = 0,

pub fn create(allocator: std.mem.Allocator, shape: *const Shape, w: u16, h: u16, options: Options, disable_dist: bool) !ErrorCorrection {
    var ret: ErrorCorrection = .{
        .options = options,
        .dist_eval = .{ .shape = shape },
        .stencil = try allocator.alloc(StencilFlags, w * h),
        .stencil_w = w,
        .stencil_h = h,
    };
    @memset(ret.stencil, .{});

    if (disable_dist and options.check_distance) {
        std.log.warn("Attempted to use distance-based error correction with a non-null scanline pass, disabling", .{});
        ret.options.check_distance = false;
    }

    return ret;
}

pub fn destroy(self: *ErrorCorrection, allocator: std.mem.Allocator) void {
    allocator.free(self.stencil);
}

pub fn correct(
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
    const vtx: Vec2 = .{ tx, ty };
    switch (self.options.mode) {
        .edge_priority => {
            self.protectCorners(shape, scale, vtx);
            self.protectEdges(scale, px_range, sdf_px, sdf_w, sdf_h, channels);
        },
        .edge_only => {
            for (self.stencil) |*mask| mask.protected = true;
        },
        .indiscriminate => {},
    }

    self.findErrors(scale, px_range, vtx, sdf_px, sdf_w, sdf_h, channels);

    for (0..sdf_w * sdf_h) |i| if (self.stencil[i].err) {
        const msdf = sdf_px[i * channels ..][0..3];
        @memset(msdf, median(msdf));
    };
}

pub fn protectCorners(self: *ErrorCorrection, shape: *Shape, scale: f64, vtx: Vec2) void {
    for (shape.contours.items) |contour| {
        if (contour.edges.items.len == 0) continue;

        var prev_edge = contour.edges.getLast();
        for (contour.edges.items) |edge| {
            const common_color = @intFromEnum(prev_edge.color) & @intFromEnum(edge.color);
            prev_edge = edge;
            if ((common_color & (common_color - 1)) == 0) continue;
            const base_point = edge.point(0) * math.v2(scale) + vtx;
            const left: i32 = @intFromFloat(@floor(base_point[0] - 0.5));
            const bottom: i32 = @intFromFloat(f64i(self.stencil_h) - @floor(base_point[1] - 0.5) - 2.0);
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

fn edgeBetweenTexels(a: *const [3]f64, b: *const [3]f64) u8 {
    var ret: u8 = 0;
    inline for (.{ red, green, blue }, 0..) |color, channel| @"continue": {
        const delta = a[channel] - b[channel];
        if (delta == 0.0) break :@"continue";
        const t = (a[channel] - 0.5) / delta;
        if (t < 0 or t > 1) break :@"continue";
        const c: [3]f64 = .{
            math.mix(a[0], b[0], t),
            math.mix(a[1], b[1], t),
            math.mix(a[2], b[2], t),
        };
        ret += color * @intFromBool(median(&c) == c[channel]);
    }
    return ret;
}

fn protectExtremeChannels(stencil_point: *StencilFlags, msd: *const [3]f64, m: f64, mask: u8) void {
    if (mask & red != 0 and msd[0] != m or
        mask & green != 0 and msd[1] != m or
        mask & blue != 0 and msd[2] != m)
        stencil_point.protected = true;
}

fn median(arr: *const [3]f64) f64 {
    return math.median(arr[0], arr[1], arr[2]);
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
    const vscale = math.v2(scale);
    const hori_radius = math.length(Vec2{ px_range, 0.0 } / vscale) * protection_radius_tolerance;
    for (0..sdf_h) |y| for (0..sdf_w - 1) |x| {
        const left = sdf_px[index(0, y, sdf_w, channels)..][0..3];
        const right = sdf_px[index(1, y, sdf_w, channels)..][0..3];
        const median_left = median(left);
        const median_right = median(right);
        if (@abs(median_left - 0.5) + @abs(median_right - 0.5) < hori_radius) {
            const mask = edgeBetweenTexels(left, right);
            protectExtremeChannels(&self.stencil[index(x, y, self.stencil_w, 1)], left, median_left, mask);
            protectExtremeChannels(&self.stencil[index(x + 1, y, self.stencil_w, 1)], right, median_right, mask);
        }
    };

    const vert_radius = math.length(Vec2{ 0.0, px_range } / vscale) * protection_radius_tolerance;
    for (0..sdf_h - 1) |y| for (0..sdf_w) |x| {
        const bottom = sdf_px[index(0, y, sdf_w, channels)..][0..3];
        const top = sdf_px[index(0, y + 1, sdf_w, channels)..][0..3];
        const median_bottom = median(bottom);
        const median_top = median(top);
        if (@abs(median_bottom - 0.5) + @abs(median_top - 0.5) < vert_radius) {
            const mask = edgeBetweenTexels(bottom, top);
            protectExtremeChannels(&self.stencil[index(x, y, self.stencil_w, 1)], bottom, median_bottom, mask);
            protectExtremeChannels(&self.stencil[index(x, y + 1, self.stencil_w, 1)], top, median_top, mask);
        }
    };

    const diag_radius = math.length(math.v2(px_range) / vscale) * protection_radius_tolerance;
    for (0..sdf_h - 1) |y| for (0..sdf_w - 1) |x| {
        const bottom_left = sdf_px[index(0, y, sdf_w, channels)..][0..3];
        const bottom_right = sdf_px[index(1, y, sdf_w, channels)..][0..3];
        const top_left = sdf_px[index(0, y + 1, sdf_w, channels)..][0..3];
        const top_right = sdf_px[index(1, y + 1, sdf_w, channels)..][0..3];
        const median_bottom_left = median(bottom_left);
        const median_bottom_right = median(bottom_right);
        const median_top_left = median(top_left);
        const median_top_right = median(top_right);
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

fn interpolatedMedianBilinear(a: *const [3]f64, l: *const [3]f64, q: *const [3]f64, t: f64) f64 {
    return math.median(
        t * (t * q[0] + l[0]) + a[0],
        t * (t * q[1] + l[1]) + a[1],
        t * (t * q[2] + l[2]) + a[2],
    );
}

fn rangeTest(span: f64, protected: bool, at: f64, bt: f64, xt: f64, am: f64, bm: f64, xm: f64) ClassifierFlags {
    if (!(am > 0.5 and bm > 0.5 and xm <= 0.5) or
        (am < 0.5 and bm < 0.5 and xm >= 0.5) or
        (!protected and math.median(am, bm, xm) != xm))
        return .{};

    const ax_span = (xt - at) * span;
    const bx_span = (bt - xt) * span;
    return .{
        .candidate = true,
        .artifact = !(xm >= am - ax_span and xm <= am + ax_span and xm >= bm - bx_span and xm <= bm + bx_span),
    };
}

fn psdfDistAt(shape: *const Shape, p: Vec2) f64 {
    var min_dist: SignedDistance = .{};
    var near_edge: ?*EdgeSegment = null;
    var near_param: f64 = 0;
    for (shape.contours.items) |contour| for (contour.edges.items) |*edge| {
        var param: f64 = 0;
        const dist = edge.signedDistance(p, &param);
        if (dist.lessThan(min_dist)) {
            min_dist = dist;
            near_edge = edge;
            near_param = param;
        }
    };
    if (near_edge) |edge| edge.distanceToPerpendicularDistance(&min_dist, p, near_param);
    return min_dist.distance;
}

fn evaluateArtifact(self: ErrorCorrection, flags: ClassifierFlags, t: f64) bool {
    if (flags.artifact) return true;
    if (!self.options.check_distance or !flags.candidate) return false;

    const t_vec = math.v2(t) * self.dist_eval.dir;

    const tx_sdf_coord = self.dist_eval.sdf_coord + t_vec;
    const tx_x = tx_sdf_coord[0] - 0.5;
    const tx_y = tx_sdf_coord[1] - 0.5;
    const floor_x = @floor(tx_x);
    const floor_y = @floor(tx_y);
    const lr = tx_x - floor_x;
    const bt = tx_y - floor_y;
    const sdf_w = self.dist_eval.sdf_w;
    const sdf_h = self.dist_eval.sdf_h;
    const left = @min(@as(u32, @intFromFloat(floor_x)), sdf_w - 1);
    const bottom = @min(@as(u32, @intFromFloat(floor_y)), sdf_w - 1);
    const right = @min(left + 1, sdf_h - 1);
    const top = @min(bottom + 1, sdf_h - 1);
    const sdf_px = self.dist_eval.sdf_px;
    const channels = self.dist_eval.channels;
    const lb_px = sdf_px[index(left, bottom, sdf_w, channels)..][0..3];
    const rb_px = sdf_px[index(right, bottom, sdf_w, channels)..][0..3];
    const lt_px = sdf_px[index(left, top, sdf_w, channels)..][0..3];
    const lr_px = sdf_px[index(left, right, sdf_w, channels)..][0..3];
    var old_sdf: [3]f64 = undefined;
    for (&old_sdf, 0..) |*c, i|
        c.* = math.mix(
            math.mix(lb_px[i], rb_px[i], lr),
            math.mix(lt_px[i], lr_px[i], lr),
            bt,
        );

    const wt = (1 - @abs(t_vec[0])) * (1 - @abs(t_vec[1]));
    const sdf = self.dist_eval.sdf.?;
    const m = median(sdf);
    const new_sdf: [3]f64 = .{
        old_sdf[0] + wt * (m - sdf[0]),
        old_sdf[1] + wt * (m - sdf[1]),
        old_sdf[2] + wt * (m - sdf[2]),
    };
    const om = median(&old_sdf);
    const nm = median(&new_sdf);
    const px_range = self.dist_eval.px_range;
    const dist = (psdfDistAt(self.dist_eval.shape.?, self.dist_eval.shape_coord + t_vec * self.dist_eval.texel_size) + px_range / 2.0) / px_range;
    return self.options.min_improve_ratio * @abs(nm - dist) < @abs(om - dist);
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

fn hasLinearArtifact(self: ErrorCorrection, span: f64, protected: bool, am: f64, a: *const [3]f64, b: *const [3]f64) bool {
    const bm = median(b);
    if (@abs(am - 0.5) < @abs(bm - 0.5)) return false;
    inline for (0..3) |idx| @"continue": {
        const idx1 = (idx + 1) % 3;
        const da = a[idx1] - a[idx];
        const db = b[idx1] - b[idx];
        const delta = da - db;
        if (delta == 0) break :@"continue";
        const t = da / (da - db);
        if (t > artifact_t_epsilon and t < 1 - artifact_t_epsilon) {
            const xm = math.median(
                math.mix(a[0], b[0], t),
                math.mix(a[1], b[1], t),
                math.mix(a[2], b[2], t),
            );
            if (self.evaluateArtifact(rangeTest(span, protected, 0, 1, t, am, bm, xm), t))
                return true;
        }
    }
    return false;
}

fn hasDiagonalArtifact(
    self: ErrorCorrection,
    span: f64,
    protected: bool,
    am: f64,
    a: *const [3]f64,
    b: *const [3]f64,
    c: *const [3]f64,
    d: *const [3]f64,
) bool {
    const dm = median(d);
    if (@abs(am - 0.5) < @abs(dm - 0.5)) return false;

    const abc: [3]f64 = .{
        a[0] - b[0] - c[0],
        a[1] - b[1] - c[1],
        a[2] - b[2] - c[2],
    };

    const q: [3]f64 = .{
        d[0] + abc[0],
        d[1] + abc[1],
        d[2] + abc[2],
    };
    if (q[0] == 0.0 or q[1] == 0.0 or q[2] == 0.0) return false;

    const l: [3]f64 = .{
        -a[0] - abc[0],
        -a[1] - abc[1],
        -a[2] - abc[2],
    };
    const t_ex: [3]f64 = .{
        -0.5 * l[0] / q[0],
        -0.5 * l[1] / q[1],
        -0.5 * l[2] / q[2],
    };
    inline for (0..3) |idx| {
        const idx1 = (idx + 1) % 3;
        const d_a = a[idx1] - a[idx];
        const d_bc = b[idx1] - b[idx] + (c[idx1] - c[idx]);
        const d_d = d[idx1] - d[idx];
        const t_ex_0 = t_ex[idx];
        const t_ex_1 = t_ex[idx1];

        var t: [2]f64 = undefined;
        const solutions = equations.solveQuadratic(&t, d_d - d_bc + d_a, d_bc - d_a - d_a, d_a);
        for (0..solutions) |i| if (t[i] > artifact_t_epsilon and t[i] < 1 - artifact_t_epsilon) {
            const xm = interpolatedMedianBilinear(a, &l, &q, t[i]);
            const FlagType = @typeInfo(ClassifierFlags).@"struct".backing_integer.?;
            var flags: FlagType = @bitCast(rangeTest(span, protected, 0, 1, t[i], am, dm, xm));
            var t_end: [2]f64 = undefined;
            var em: [2]f64 = undefined;

            if (t_ex_0 > 0 and t_ex_0 < 1) {
                t_end[0] = 0;
                t_end[1] = 1;
                em[0] = am;
                em[1] = dm;
                t_end[@intFromBool(t_ex_0 > t[i])] = t_ex_0;
                em[@intFromBool(t_ex_0 > t[i])] = interpolatedMedianBilinear(a, &l, &q, t_ex_0);
                flags |= @as(FlagType, @bitCast(rangeTest(span, protected, t_end[0], t_end[1], t[i], em[0], em[1], xm)));
            }

            if (t_ex_1 > 0 and t_ex_1 < 1) {
                t_end[0] = 0;
                t_end[1] = 1;
                em[0] = am;
                em[1] = dm;
                t_end[@intFromBool(t_ex_1 > t[i])] = t_ex_1;
                em[@intFromBool(t_ex_1 > t[i])] = interpolatedMedianBilinear(a, &l, &q, t_ex_1);
                flags |= @as(FlagType, @bitCast(rangeTest(span, protected, t_end[0], t_end[1], t[i], em[0], em[1], xm)));
            }

            if (self.evaluateArtifact(@bitCast(flags), t[i]))
                return true;
        };
    }

    return false;
}

fn findErrors(
    self: *ErrorCorrection,
    scale: f64,
    px_range: f64,
    vtx: Vec2,
    sdf_px: []const f64,
    sdf_w: u16,
    sdf_h: u16,
    channels: u8,
) void {
    const min_deviation_ratio = self.options.min_deviation_ratio;
    const vscale = math.v2(scale);
    const hori_span = math.length(Vec2{ px_range, 0.0 } / vscale) * min_deviation_ratio;
    const vert_span = math.length(Vec2{ 0.0, px_range } / vscale) * min_deviation_ratio;
    const diag_span = math.length(math.v2(px_range) / vscale) * min_deviation_ratio;

    self.dist_eval.channels = channels;
    self.dist_eval.px_range = px_range;
    self.dist_eval.texel_size = math.v2(1.0) / vscale;
    self.dist_eval.sdf_px = sdf_px;
    self.dist_eval.sdf_w = sdf_w;
    self.dist_eval.sdf_h = sdf_h;

    for (0..sdf_h) |y| for (0..sdf_w) |x| {
        const current = sdf_px[index(x, y, sdf_w, channels)..][0..3];
        const median_current = median(current);
        var current_stencil = &self.stencil[index(x, y, self.stencil_w, 1)];
        const is_protected = current_stencil.protected;

        const fx: f64 = @floatFromInt(x);
        const fy: f64 = @floatFromInt(y);
        self.dist_eval.sdf = current;
        self.dist_eval.shape_coord = Vec2{ fx + 0.5, fy + 0.5 } / vscale + vtx;
        self.dist_eval.sdf_coord = .{ fx + 0.5, fy + 0.5 };

        if (x > 0) {
            const left = sdf_px[index(x - 1, y, sdf_w, channels)..][0..3];
            self.dist_eval.dir = .{ -1, 0 };
            if (self.hasLinearArtifact(hori_span, is_protected, median_current, current, left)) {
                current_stencil.err = true;
                continue;
            }

            if (y > 0) {
                const top = sdf_px[index(x, y - 1, sdf_w, channels)..][0..3];
                const top_left = sdf_px[index(x - 1, y - 1, sdf_w, channels)..][0..3];
                self.dist_eval.dir = .{ -1, -1 };
                if (self.hasDiagonalArtifact(diag_span, is_protected, median_current, current, left, top, top_left)) {
                    current_stencil.err = true;
                    continue;
                }
            }

            if (y < sdf_h - 1) {
                const bottom = sdf_px[index(x, y + 1, sdf_w, channels)..][0..3];
                const bottom_left = sdf_px[index(x - 1, y + 1, sdf_w, channels)..][0..3];
                self.dist_eval.dir = .{ -1, 1 };
                if (self.hasDiagonalArtifact(diag_span, is_protected, median_current, current, left, bottom, bottom_left)) {
                    current_stencil.err = true;
                    continue;
                }
            }
        }

        if (y > 0) {
            const top = sdf_px[index(x, y - 1, sdf_w, channels)..][0..3];
            self.dist_eval.dir = .{ 0, -1 };
            if (self.hasLinearArtifact(vert_span, is_protected, median_current, current, top)) {
                current_stencil.err = true;
                continue;
            }
        }

        if (x < sdf_w - 1) {
            const right = sdf_px[index(x + 1, y, sdf_w, channels)..][0..3];
            self.dist_eval.dir = .{ 1, 0 };
            if (self.hasLinearArtifact(hori_span, is_protected, median_current, current, right)) {
                current_stencil.err = true;
                continue;
            }

            if (y > 0) {
                const top = sdf_px[index(x, y - 1, sdf_w, channels)..][0..3];
                const top_right = sdf_px[index(x + 1, y - 1, sdf_w, channels)..][0..3];
                self.dist_eval.dir = .{ 1, -1 };
                if (self.hasDiagonalArtifact(diag_span, is_protected, median_current, current, right, top, top_right)) {
                    current_stencil.err = true;
                    continue;
                }
            }

            if (y < sdf_h - 1) {
                const bottom = sdf_px[index(x, y + 1, sdf_w, channels)..][0..3];
                const bottom_right = sdf_px[index(x + 1, y + 1, sdf_w, channels)..][0..3];
                self.dist_eval.dir = .{ 1, 1 };
                if (self.hasDiagonalArtifact(diag_span, is_protected, median_current, current, right, bottom, bottom_right)) {
                    current_stencil.err = true;
                    continue;
                }
            }
        }

        if (y < sdf_h - 1) {
            const bottom = sdf_px[index(x, y + 1, sdf_w, channels)..][0..3];
            self.dist_eval.dir = .{ 0, 1 };
            if (self.hasLinearArtifact(vert_span, is_protected, median_current, current, bottom)) {
                current_stencil.err = true;
                continue;
            }
        }
    };
}
