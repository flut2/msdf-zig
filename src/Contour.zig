const std = @import("std");

const Contour = @import("Contour.zig");
const EdgeSegment = @import("EdgeSegment.zig");
const math = @import("math.zig");
const Vec2 = @Vector(2, f64);

edges: std.ArrayList(EdgeSegment) = .empty,

pub fn boundPoint(l: *f64, b: *f64, r: *f64, t: *f64, p: Vec2) void {
    const x = p[0];
    const y = p[1];
    if (x < l.*) l.* = x;
    if (y < b.*) b.* = y;
    if (x > r.*) r.* = x;
    if (y > t.*) t.* = y;
}

pub fn bound(self: Contour, l: *f64, b: *f64, r: *f64, t: *f64) void {
    for (self.edges.items) |edge| edge.bound(l, b, r, t);
}

pub fn boundMiters(self: Contour, l: *f64, b: *f64, r: *f64, t: *f64, border: f64, miter_limit: f64, polarity: i32) void {
    if (self.edges.items.len == 0) return;
    var prev_dir = math.normal(self.edges.getLast().direction(1), false);
    for (self.edges.items) |edge| {
        const dir = math.normal(edge.direction(0), false) * math.v2(-1.0);
        if (@as(f64, @floatFromInt(polarity)) * math.cross(prev_dir, dir) >= 0) {
            var miter_length = miter_limit;
            const q = 0.5 * (1 - math.dot(prev_dir, dir));
            if (q > 0) miter_length = @min(1 / @sqrt(q), miter_limit);
            const miter = edge.point(0) + math.normal(prev_dir + dir, false) * math.v2(border * miter_length);
            boundPoint(l, b, r, t, miter);
        }
        prev_dir = math.normal(edge.direction(1), false);
    }
}

pub fn reverse(self: *Contour) void {
    std.mem.reverse(EdgeSegment, self.edges.items);
    for (self.edges.items) |*edge| edge.reverse();
}
