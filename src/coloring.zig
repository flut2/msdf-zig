const std = @import("std");

const Shape = @import("Shape.zig");
const Vec2 = @import("Vec2.zig");
const EdgeColor = @import("edge_color.zig").EdgeColor;
const EdgeSegment = @import("EdgeSegment.zig");

fn symmetricalTrichotomy(pos: usize, n: usize) i32 {
    const fpos: f64 = @floatFromInt(pos);
    const fn1: f64 = @floatFromInt(n - 1);
    return @as(i32, @intFromFloat(@floor(3 + 2.875 * fpos / fn1 - 1.4375 + 0.5))) - 3;
}

fn isCorner(a: Vec2, b: Vec2, cross_threshold: f64) bool {
    return a.dot(b) <= 0 or @abs(a.cross(b)) > cross_threshold;
}

pub fn simple(allocator: std.mem.Allocator, shape: *Shape, angle_threshold: f64) !void {
    const cross_threshold = @sin(angle_threshold);
    var color: EdgeColor = .init();
    var corners: std.ArrayListUnmanaged(u32) = .empty;
    defer corners.deinit(allocator);
    for (shape.contours.items) |*contour| {
        if (contour.edges.items.len == 0) continue;

        corners.clearRetainingCapacity();
        var prev_dir = contour.edges.getLast().direction(1);
        for (contour.edges.items, 0..) |edge, i| {
            if (isCorner(prev_dir.normal(false), edge.direction(0).normal(false), cross_threshold))
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
                color.random();
                colors[0] = color;
                color.random();
                colors[2] = color;
                const corner = corners.items[0];
                const edges_len = contour.edges.items.len;
                if (edges_len >= 3) {
                    for (contour.edges.items, 0..) |*edge, i| edge.color = colors[@intCast(1 + symmetricalTrichotomy(i, edges_len))];
                } else if (edges_len >= 1) {
                    var parts: [7]EdgeSegment = @splat(.{});
                    contour.edges.items[0].splitInThirds(&parts[0 + 3 * corner], &parts[1 + 3 * corner], &parts[2 + 3 * corner]);
                    if (edges_len >= 2) {
                        contour.edges.items[1].splitInThirds(&parts[3 - 3 * corner], &parts[4 - 3 * corner], &parts[5 - 3 * corner]);
                        inline for (.{ 0, 1 }) |i| parts[i].color = colors[0];
                        inline for (.{ 2, 3 }) |i| parts[i].color = colors[1];
                        inline for (.{ 4, 5 }) |i| parts[i].color = colors[2];
                    } else inline for (.{ 0, 1, 2 }) |i| parts[i].color = colors[i];
                    contour.edges.clearRetainingCapacity();
                    const min_idx = @min(0 + 3 * corner, 3 - 3 * corner);
                    for (0..min_idx) |i| try contour.edges.append(allocator, parts[i]);
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
                        color.change(@enumFromInt(@intFromBool(spline == corners_len - 1) * @intFromEnum(initial_color)));
                    }
                    contour.edges.items[idx].color = color;
                }
            },
        }
    }
}

pub fn inkTrap(_: std.mem.Allocator, _: *Shape, _: f64) !void {
    @panic("Not implemented yet");
}

pub fn distance(_: std.mem.Allocator, _: *Shape, _: f64) !void {
    @panic("Not implemented yet");
}
