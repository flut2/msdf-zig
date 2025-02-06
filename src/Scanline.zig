const std = @import("std");

const math = @import("math.zig");

const Scanline = @This();

pub const Intersection = struct {
    x: f64,
    dir: i32,
};

const FillRule = enum {
    non_zero,
    odd,
    positive,
    negative,
};

intersections: std.ArrayListUnmanaged(Intersection) = .empty,
last_index: i32 = 0,

pub fn interpretFillRule(intersections: i32, fill_rule: FillRule) bool {
    return switch (fill_rule) {
        .non_zero => intersections != 0,
        .odd => intersections % 2 == 1,
        .positive => intersections > 0,
        .negative => intersections < 0,
    };
}

pub fn overlap(a: Scanline, b: Scanline, x_from: f64, x_to: f64, fill_rule: FillRule) f64 {
    var total: f64 = 0.0;
    var a_inside = false;
    var b_inside = false;
    var ai: i32 = 0;
    var bi: i32 = 0;
    var ax = if (!a.intersections.empty()) a.intersections[ai].x else x_to;
    var bx = if (!b.intersections.empty()) b.intersections[bi].x else x_to;
    while (ax < x_from or bx < x_from) {
        const x_next = @min(ax, bx);
        if (ax == x_next and ai < a.intersections.items.len) {
            a_inside = interpretFillRule(a.intersections.items[ai].direction, fill_rule);
            ai += 1;
            ax = if (ai < a.intersections.items.len) a.intersections[ai].x else x_to;
        }
        if (bx == x_next and bi < b.intersections.items.len) {
            b_inside = interpretFillRule(b.intersections[bi].direction, fill_rule);
            bi += 1;
            bx = if (bi < b.intersections.items.len) b.intersections[bi].x else x_to;
        }
    }
    var x = x_from;
    while (ax < x_to or bx < x_to) {
        const x_next = @min(ax, bx);
        if (a_inside == b_inside) total += x_next - x;
        if (ax == x_next and ai < a.intersections.items.len) {
            a_inside = interpretFillRule(a.intersections[ai].direction, fill_rule);
            ai += 1;
            ax = if (ai < a.intersections.items.len) a.intersections[ai].x else x_to;
        }
        if (bx == x_next and bi < b.intersections.items.len) {
            b_inside = interpretFillRule(b.intersections[bi].direction, fill_rule);
            bi += 1;
            bx = if (bi < b.intersections.items.len) b.intersections[bi].x else x_to;
        }
        x = x_next;
    }
    if (a_inside == b_inside) total += x_to - x;
    return total;
}

pub fn countIntersections(self: *Scanline, x: f64) i32 {
    return self.moveTo(x) + 1;
}

pub fn sumIntersections(self: *Scanline, x: f64) i32 {
    const index = self.moveTo(x);
    if (index >= 0) return self.intersections[@intCast(index)].dir;
    return 0;
}

pub fn filled(self: *Scanline, x: f64, fill_rule: FillRule) bool {
    return interpretFillRule(self.sumIntersections(x), fill_rule);
}

fn lessThan(_: void, a: Intersection, b: Intersection) bool {
    return a.x < b.x;
}

pub fn preprocess(self: *Scanline) void {
    self.last_index = 0;
    if (self.intersections.items.len == 0) return;
    std.sort.pdq(void, self.intersections.items, {}, lessThan);
    var total_direction: i32 = 0;
    for (self.intersections.items) |*intersection| {
        total_direction += intersection.direction;
        intersection.direction = total_direction;
    }
}

pub fn moveTo(self: *Scanline, x: f64) i32 {
    if (self.intersections.items.len == 0) return -1;

    var index = self.last_index;
    if (x < self.intersections.items[index].x) {
        while (x < self.intersections.items[index].x) {
            if (index == 0) {
                self.last_index = 0;
                return -1;
            }
            index -= 1;
        }
    } else while (index < self.intersections.items.len - 1 and x >= self.intersections.items[index + 1].x)
        index += 1;

    self.last_index = index;
    return index;
}
