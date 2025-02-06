const std = @import("std");

const Contour = @import("Contour.zig");
const EdgeSegment = @import("EdgeSegment.zig");
const math = @import("math.zig");
const Scanline = @import("Scanline.zig");
const Vec2 = @import("Vec2.zig");

const deconverge_overshoot = 1.11111111111111111;
const corner_dot_epsilon = 0.000001;

const Shape = @This();
pub const Bounds = struct {
    left: f64 = 0.0,
    right: f64 = 0.0,
    bottom: f64 = 0.0,
    top: f64 = 0.0,
};

contours: std.ArrayListUnmanaged(Contour) = .empty,

pub fn print(self: Shape) void {
    for (self.contours.items, 0..) |contour, i| {
        std.log.info("Contour {}: [", .{i});
        for (contour.edges.items, 0..) |edge, j| switch (edge.segment) {
            .linear => |p| std.log.info(
                "Edge {}: linear [x1={d:.2},y1={d:.2},x2={d:.2},y2={d:.2}]",
                .{ j, p[0].x, p[0].y, p[1].x, p[1].y },
            ),
            .quadratic_bezier => |p| std.log.info(
                "Edge {}: quadratic [x1={d:.2},y1={d:.2},x2={d:.2},y2={d:.2},x3={d:.2},y3={d:.2}]",
                .{ j, p[0].x, p[0].y, p[1].x, p[1].y, p[2].x, p[2].y },
            ),
            .cubic_bezier => |p| std.log.info(
                "Edge {}: cubic [x1={d:.2},y1={d:.2},x2={d:.2},y2={d:.2},x3={d:.2},y3={d:.2},x4={d:.2},y4={d:.2}]",
                .{ j, p[0].x, p[0].y, p[1].x, p[1].y, p[2].x, p[2].y, p[3].x, p[3].y },
            ),
        };
        std.log.info("];", .{});
    }
}

pub fn validate(self: Shape) bool {
    for (self.contours.items) |contour| if (contour.edges.items.len > 0) {
        var corner = contour.edges.getLast().point(1);
        for (contour.edges.items) |edge| {
            const p0 = edge.point(0);
            if (!p0.eql(corner))
                return false;
            corner = edge.point(1);
        }
    };

    return true;
}

pub fn normalize(self: *Shape, allocator: std.mem.Allocator) !void {
    for (self.contours.items) |*contour| {
        if (contour.edges.items.len == 1) {
            var parts: [3]EdgeSegment = @splat(.{});
            contour.edges.items[0].splitInThirds(&parts[0], &parts[1], &parts[2]);
            contour.edges.clearRetainingCapacity();
            try contour.edges.appendSlice(allocator, &parts);
        } else {
            var prev_edge = &contour.edges.items[contour.edges.items.len - 1];
            for (contour.edges.items) |*edge| {
                var prev_dir = prev_edge.direction(1).normal(false);
                const cur_dir = edge.direction(0).normal(false);
                if (prev_dir.dot(cur_dir) < corner_dot_epsilon - 1) {
                    const factor = deconverge_overshoot *
                        @sqrt(1 - (corner_dot_epsilon - 1) * (corner_dot_epsilon - 1)) / (corner_dot_epsilon - 1);
                    var axis = cur_dir.sub(prev_dir).normal(false).mul(factor);
                    if (prev_edge.directionChange(1).cross(edge.direction(0)) + (edge.directionChange(0).cross(prev_edge.direction(1))) < 0)
                        axis = axis.mul(-1.0);
                    prev_edge.deconverge(1, axis.ortho(true));
                    edge.deconverge(0, axis.ortho(false));
                }
                prev_edge = edge;
            }
        }
    }
}

pub fn bound(self: Shape, l: *f64, b: *f64, r: *f64, t: *f64) void {
    for (self.contours.items) |contour| contour.bound(l, b, r, t);
}

pub fn boundMiters(self: Shape, l: *f64, b: *f64, r: *f64, t: *f64, border: f64, miter_limit: f64, polarity: i32) void {
    for (self.contours.items) |contour| contour.boundMiters(l, b, r, t, border, miter_limit, polarity);
}

pub fn getBounds(self: Shape, border: f64, miter_limit: f64, polarity: i32) Bounds {
    const large_value = 1e240;
    var bounds: Bounds = .{ .left = large_value, .bottom = large_value, .right = -large_value, .top = -large_value };
    self.bound(&bounds.left, &bounds.bottom, &bounds.right, &bounds.top);
    if (border > 0) {
        bounds.left -= border;
        bounds.bottom -= border;
        bounds.right += border;
        bounds.top += border;
        if (miter_limit > 0)
            self.boundMiters(&bounds.left, &bounds.bottom, &bounds.right, &bounds.top, border, miter_limit, polarity);
    }
    return bounds;
}

pub fn scanline(self: Shape, line: *Scanline, y: f64, allocator: std.mem.Allocator) !void {
    line.intersections.clearRetainingCapacity();
    defer line.preprocess();

    var x: [3]f64 = undefined;
    var dy: [3]i32 = undefined;
    for (self.contours.items) |contour| for (contour.edges.items) |edge| {
        for (0..edge.scanlineIntersections(&x, &dy, y)) |i|
            try line.intersections.append(allocator, .{ .x = x[i], .dir = dy[i] });
    };
}

pub fn edgeCount(self: Shape) u32 {
    var total: u32 = 0;
    for (self.contours.items) |contour| total += contour.edges.items.len;
    return total;
}

pub fn orientContours(self: *Shape, allocator: std.mem.Allocator) !void {
    const Intersection = struct {
        x: f64,
        direction: i32,
        contour_index: i32,

        pub fn lessThan(_: void, a: @This(), b: @This()) bool {
            return a.x < b.x;
        }
    };

    const ratio = 0.5 * (@sqrt(5.0) - 1);
    var intersections: std.ArrayListUnmanaged(Intersection) = .empty;
    defer intersections.deinit(allocator);
    var orientations: std.ArrayListUnmanaged(i32) = .empty;
    defer orientations.deinit(allocator);
    const contours_len = self.contours.items.len;
    try orientations.ensureTotalCapacity(allocator, contours_len);
    try orientations.appendNTimes(allocator, 0, contours_len);
    for (0..contours_len) |i| {
        if (orientations.items[i] == 0 or self.contours.items[i].edges.items.len == 0) continue;
        const y0 = self.contours.items[i].edges.items[0].point(0).y;
        var y1 = y0;
        for (self.contours.items[i].edges.items) |edge| y1 = edge.point(1).y;
        for (self.contours.items[i].edges.items) |edge| y1 = edge.point(ratio).y;
        const y = math.mix(y0, y1, ratio);
        var x: [3]f64 = undefined;
        var dy: [3]i32 = undefined;
        for (0..self.contours.items.len) |j|
            for (self.contours.items[j].edges.items) |edge|
                for (0..edge.scanlineIntersections(&x, &dy, y)) |k|
                    try intersections.append(allocator, .{ .x = x[k], .direction = dy[k], .contour_index = @intCast(j) });

        if (intersections.items.len == 0) continue;
        std.sort.pdq(Intersection, intersections.items, {}, Intersection.lessThan);
        for (1..intersections.items.len) |j| if (intersections.items[j].x == intersections.items[j - 1].x) {
            intersections.items[j].direction = 0;
            intersections.items[j - 1].direction = 0;
        };
        for (0..intersections.items.len) |j| if (intersections.items[j].direction != 0) {
            orientations.items[@intCast(intersections.items[j].contour_index)] +=
                2 * ((@as(i32, @intCast(j)) & 1) ^ @intFromBool(intersections.items[j].direction > 0)) - 1;
        };
        intersections.clearRetainingCapacity();
    }

    for (self.contours.items, orientations.items) |*contour, orientation| if (orientation < 0) contour.reverse();
}
