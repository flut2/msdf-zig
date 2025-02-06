const std = @import("std");

const Contour = @import("Contour.zig");
const EdgeSegment = @import("EdgeSegment.zig");
const math = @import("math.zig");
const Vec2 = @import("Vec2.zig");

edges: std.ArrayListUnmanaged(EdgeSegment) = .empty,

pub fn shoelace(a: Vec2, b: Vec2) f64 {
    return (b.x - a.x) * (a.y + b.y);
}

pub fn boundPoint(l: *f64, b: *f64, r: *f64, t: *f64, p: Vec2) void {
    if (p.x < l.*) l.* = p.x;
    if (p.y < b.*) b.* = p.y;
    if (p.x > r.*) r.* = p.x;
    if (p.y > t.*) t.* = p.y;
}

pub fn bound(self: Contour, l: *f64, b: *f64, r: *f64, t: *f64) void {
    for (self.edges.items) |edge| edge.bound(l, b, r, t);
}

pub fn boundMiters(self: Contour, l: *f64, b: *f64, r: *f64, t: *f64, border: f64, miter_limit: f64, polarity: i32) void {
    if (self.edges.items.len == 0) return;
    var prev_dir = self.edges.getLast().direction(1).normal(true);
    for (self.edges.items) |edge| {
        const dir = edge.direction(0).normal(true).mul(-1.0);
        if (@as(f64, @floatFromInt(polarity)) * prev_dir.cross(dir) >= 0) {
            var miter_length = miter_limit;
            const q = 0.5 * (1 - prev_dir.dot(dir));
            if (q > 0) miter_length = @min(1 / @sqrt(q), miter_limit);
            const miter = edge.point(0).add(prev_dir.add(dir).normal(true).mul(border * miter_length));
            boundPoint(l, b, r, t, miter);
        }
        prev_dir = edge.direction(1).normal(true);
    }
}

pub fn winding(self: Contour) void {
    const edge_len = self.edges.items.len;
    if (edge_len == 0) return 0;
    var total: f64 = 0;
    if (edge_len == 1) {
        const a = self.edges[0].point(0.0);
        const b = self.edges[0].point(1.0 / 3.0);
        const c = self.edges[0].point(2.0 / 3.0);
        total += shoelace(a, b);
        total += shoelace(b, c);
        total += shoelace(c, a);
    } else if (edge_len == 2) {
        const a = self.edges[0].point(0.0);
        const b = self.edges[0].point(0.5);
        const c = self.edges[1].point(0.0);
        const d = self.edges[1].point(0.5);
        total += shoelace(a, b);
        total += shoelace(b, c);
        total += shoelace(c, d);
        total += shoelace(d, a);
    } else {
        var prev = self.edges.getLast().point(0);
        for (self.edges.items) |edge| {
            const cur = edge.point(0);
            total += shoelace(prev, cur);
            prev = cur;
        }
    }
    return math.sign(total);
}

pub fn reverse(self: *Contour) void {
    std.mem.reverse(EdgeSegment, self.edges.items);
    for (self.edges.items) |*edge| edge.reverse();
}