const std = @import("std");

const EdgeColor = @import("edge_color.zig").EdgeColor;
const EdgeSegment = @import("EdgeSegment.zig");
const math = @import("math.zig");
const Shape = @import("Shape.zig");

const Vec2 = @Vector(2, f64);

fn symmetricalTrichotomy(pos: usize, n: usize) i32 {
    const fpos: f64 = @floatFromInt(pos);
    const fn1: f64 = @floatFromInt(n - 1);
    return @as(i32, @intFromFloat(@floor(3 + 2.875 * fpos / fn1 - 1.4375 + 0.5))) - 3;
}

fn isCorner(a: Vec2, b: Vec2, cross_threshold: f64) bool {
    return math.dot(a, b) <= 0 or @abs(math.cross(a, b)) > cross_threshold;
}

pub fn colorShape(allocator: std.mem.Allocator, shape: *Shape, angle_threshold: f64) !void {
    const cross_threshold = @sin(angle_threshold);
    var color: EdgeColor = .init();
    var corners: std.ArrayList(u32) = .empty;
    defer corners.deinit(allocator);
    for (shape.contours.items) |*contour| {
        if (contour.edges.items.len == 0) continue;

        corners.clearRetainingCapacity();
        var prev_dir = contour.edges.getLast().direction(1);
        for (contour.edges.items, 0..) |edge, i| {
            if (isCorner(math.normal(prev_dir, true), math.normal(edge.direction(0), true), cross_threshold))
                try corners.append(allocator, @intCast(i));
            prev_dir = edge.direction(1);
        }

        switch (corners.items.len) {
            0 => {
                color.random();
                for (contour.edges.items) |*edge| edge.color = color;
            },
            1 => {
                var colors: [3]EdgeColor = .{ .black, .white, .black };
                inline for (.{ 0, 2 }) |i| {
                    color.random();
                    colors[i] = color;
                }
                const corner = corners.items[0];
                const corner_idx = 3 * corner;
                const edges_len = contour.edges.items.len;
                if (edges_len >= 3) {
                    for (contour.edges.items, 0..) |*edge, i| edge.color = colors[@intCast(1 + symmetricalTrichotomy(i, edges_len))];
                } else if (edges_len >= 1) {
                    var parts: [7]EdgeSegment = @splat(.{});
                    contour.edges.items[0].splitInThirds(parts[corner_idx..][0..3]);
                    if (edges_len >= 2) {
                        contour.edges.items[1].splitInThirds(parts[3 - corner_idx ..][0..3]);
                        for (0..6) |i| parts[i].color = colors[@divFloor(i, 2)];
                    } else for (0..3) |i| parts[i].color = colors[i];
                    contour.edges.clearRetainingCapacity();
                    for (0..@min(corner_idx, 3 - corner_idx)) |i|
                        try contour.edges.append(allocator, parts[i]);
                }
            },
            else => {
                const corners_len = corners.items.len;
                var spline: u32 = 0;
                const start = corners.items[0];
                const edges_len = contour.edges.items.len;
                color.random();
                const initial_color = color;
                for (0..edges_len) |i| {
                    const idx = (start + i) % edges_len;
                    if (spline + 1 < corners_len and corners.items[spline + 1] == idx) {
                        spline += 1;
                        color.change(if (spline == corners_len - 1) initial_color else .black);
                    }
                    contour.edges.items[idx].color = color;
                }
            },
        }
    }
}
