const std = @import("std");

const EdgeColor = @import("edge_color.zig").EdgeColor;
const equations = @import("equations.zig");
const math = @import("math.zig");
const mix = math.mix;
const dot = math.dot;
const cross = math.cross;
const normal = math.normal;
const v2 = math.v2;
const SignedDistance = @import("SignedDistance.zig");

const Vec2 = @Vector(2, f64);
const EdgeSegment = @This();

const cubic_starts = 4;
const cubic_steps = 4;

color: EdgeColor = .white,
segment: union(enum) {
    linear: [2]Vec2,
    quadratic_bezier: [3]Vec2,
    cubic_bezier: [4]Vec2,
} = .{ .linear = @splat(@splat(0.0)) },

pub fn format(self: EdgeSegment, writer: *std.Io.Writer) std.Io.Writer.Error!void {
    switch (self.segment) {
        inline else => |vals, tag| try writer.print("type={t}, color={t}, vals={}", .{ tag, self.color, vals }),
    }
}

fn createLinear(p0: Vec2, p1: Vec2, color: EdgeColor) EdgeSegment {
    return .{ .color = color, .segment = .{ .linear = .{ p0, p1 } } };
}

fn createQuadratic(p0: Vec2, p1: Vec2, p2: Vec2, color: EdgeColor) EdgeSegment {
    return .{ .color = color, .segment = .{ .quadratic_bezier = .{ p0, p1, p2 } } };
}

fn createCubic(p0: Vec2, p1: Vec2, p2: Vec2, p3: Vec2, color: EdgeColor) EdgeSegment {
    return .{ .color = color, .segment = .{ .cubic_bezier = .{ p0, p1, p2, p3 } } };
}

pub fn create(p0: Vec2, p1: Vec2, point_2: ?Vec2, point_3: ?Vec2, color: EdgeColor) EdgeSegment {
    if (point_3) |p3| {
        if (point_2 == null)
            @panic("Invalid parameters, you need to specify `point_2` if you specify `point_3`.");

        const p2 = point_2.?;
        var p12: Vec2 = p2 - p1;
        if (cross(p1 - p0, p12) == 0.0 and cross(p12, p3 - p2) == 0.0)
            return createLinear(p0, p3, color);

        p12 = p1 * v2(1.5) - p0 * v2(0.5);
        if (std.meta.eql(p12, p2 * v2(1.5) - p3 * v2(0.5)))
            return createQuadratic(p0, p12, p3, color);

        return createCubic(p0, p1, p2, p3, color);
    }

    if (point_2) |p2| {
        if (cross(p1 - p0, p2 - p1) == 0.0)
            return createLinear(p0, p2, color);
        return createQuadratic(p0, p1, p2, color);
    }

    return createLinear(p0, p1, color);
}

pub fn distanceToPerpendicularDistance(self: EdgeSegment, distance: *SignedDistance, origin: Vec2, param: f64) void {
    if (param < 0) {
        const dir = normal(self.direction(0), true);
        const aq = origin - self.point(0);
        if (dot(aq, dir) < 0) {
            const perp_dist = cross(aq, dir);
            if (@abs(perp_dist) <= @abs(distance.distance))
                distance.* = .{ .distance = perp_dist };
        }
    } else if (param > 1) {
        const dir = normal(self.direction(1), true);
        const bq = origin - self.point(1);
        if (dot(bq, dir) > 0) {
            const perp_dist = cross(bq, dir);
            if (@abs(perp_dist) <= @abs(distance.distance))
                distance.* = .{ .distance = perp_dist };
        }
    }
}

pub fn point(self: EdgeSegment, param: f64) Vec2 {
    switch (self.segment) {
        .linear => |p| return mix(p[0], p[1], param),
        .quadratic_bezier => |p| return mix(
            mix(p[0], p[1], param),
            mix(p[1], p[2], param),
            param,
        ),
        .cubic_bezier => |p| {
            const p12 = mix(p[1], p[2], param);
            return mix(
                mix(mix(p[0], p[1], param), p12, param),
                mix(p12, mix(p[2], p[3], param), param),
                param,
            );
        },
    }
}

pub fn direction(self: EdgeSegment, comptime index: u1) Vec2 {
    switch (self.segment) {
        .linear => |p| return p[1] - p[0],
        .quadratic_bezier => |p| {
            const tangent = switch (index) {
                0 => p[1] - p[0],
                1 => p[2] - p[1],
            };
            if (std.meta.eql(tangent, v2(0.0))) return p[2] - p[0];
            return tangent;
        },
        .cubic_bezier => |p| {
            const tangent = switch (index) {
                0 => p[1] - p[0],
                1 => p[3] - p[2],
            };
            if (std.meta.eql(tangent, v2(0.0))) switch (index) {
                0 => return p[2] - p[0],
                1 => return p[3] - p[1],
            };
            return tangent;
        },
    }
}

pub fn directionChange(self: EdgeSegment, comptime index: u1) Vec2 {
    switch (self.segment) {
        .linear => return .{ 0.0, 0.0 },
        .quadratic_bezier => |p| return (p[2] - p[1]) - (p[1] - p[0]),
        .cubic_bezier => |p| return switch (index) {
            0 => (p[2] - p[1]) - (p[1] - p[0]),
            1 => (p[3] - p[2]) - (p[2] - p[1]),
        },
    }
}

pub fn signedDistance(self: EdgeSegment, origin: Vec2, param: *f64) SignedDistance {
    switch (self.segment) {
        .linear => |p| {
            const aq = origin - p[0];
            const ab = p[1] - p[0];
            const new_param = dot(aq, ab) / dot(ab, ab);
            param.* = new_param;
            const eq = p[@intFromBool(new_param > 0.5)] - origin;
            const endpoint_dist = math.length(eq);
            if (new_param > 0.0 and new_param < 1.0) {
                const ortho_dist = dot(math.orthonormal(ab, false, true), aq);
                if (@abs(ortho_dist) < endpoint_dist) return .{ .distance = ortho_dist };
            }
            return .{
                .distance = math.nonZeroSign(cross(aq, ab)) * endpoint_dist,
                .dot = @abs(dot(normal(ab, true), normal(eq, true))),
            };
        },
        .quadratic_bezier => |p| {
            const qa = p[0] - origin;
            const ab = p[1] - p[0];
            const br = p[2] - p[1] - ab;
            const a = dot(br, br);
            const b = 3.0 * dot(ab, br);
            const c = 2.0 * dot(ab, ab) + dot(qa, br);
            const d = dot(qa, ab);
            var roots: [3]f64 = undefined;
            const num_solutions = equations.solveCubic(&roots, a, b, c, d);

            var ep_dir = self.direction(0);
            var min_dist = math.nonZeroSign(cross(ep_dir, qa)) * math.length(qa);
            param.* = -dot(qa, ep_dir) / dot(ep_dir, ep_dir);
            ep_dir = self.direction(1);
            var dist = math.length(p[2] - origin);
            if (dist < @abs(min_dist)) {
                min_dist = math.nonZeroSign(cross(ep_dir, p[2] - origin)) * dist;
                param.* = dot(origin - p[1], ep_dir) / dot(ep_dir, ep_dir);
            }

            for (roots[0..num_solutions]) |root| if (root > 0 and root < 1) {
                const qe = qa + ab * v2(root * 2.0) + br * v2(root * root);
                dist = math.length(qe);
                if (dist < @abs(min_dist)) {
                    min_dist = math.nonZeroSign(cross(ab + br * v2(root), qe)) * dist;
                    param.* = root;
                }
            };

            if (param.* < 0.0)
                return .{
                    .distance = min_dist,
                    .dot = @abs(dot(
                        normal(self.direction(0), true),
                        normal(qa, true),
                    )),
                }
            else if (param.* > 1.0)
                return .{
                    .distance = min_dist,
                    .dot = @abs(dot(
                        normal(self.direction(1), true),
                        normal(p[2] - origin, true),
                    )),
                };
            return .{ .distance = min_dist };
        },
        .cubic_bezier => |p| {
            const qa = p[0] - origin;
            const ab = p[1] - p[0];
            const br = p[2] - p[1] - ab;
            const as = (p[3] - p[2]) - (p[2] - p[1]) - br;

            var ep_dir = self.direction(0);
            var min_dist = math.nonZeroSign(cross(ep_dir, qa)) * math.length(qa);
            param.* = -dot(qa, ep_dir) / dot(ep_dir, ep_dir);

            ep_dir = self.direction(1);
            var dist = math.length(p[3] - origin);
            if (dist < @abs(min_dist)) {
                min_dist = math.nonZeroSign(cross(ep_dir, p[3] - origin)) * dist;
                param.* = dot(ep_dir - (p[3] - origin), ep_dir) / dot(ep_dir, ep_dir);
            }

            for (0..cubic_starts) |i| {
                const fi: f64 = @floatFromInt(i);
                var t = fi / cubic_starts;
                const t_cubed = v2(t * t * t);
                const t_sqr_3 = v2(t * t * 3.0);
                const t_6 = v2(t * 6.0);
                const t_3 = v2(t * 3.0);
                var qe = qa + ab * t_3 + br * t_sqr_3 + as * t_cubed;
                for (0..cubic_steps) |_| {
                    const d1 = ab * v2(3.0) + br * t_6 + as * t_sqr_3;
                    const d2 = br * v2(6.0) + as * t_6;
                    t -= dot(qe, d1) / (dot(d1, d1) + dot(qe, d2));
                    if (t <= 0 or t >= 1) break;
                    qe = qa + ab * t_3 + br * t_sqr_3 + as * t_cubed;
                    dist = math.length(qe);
                    if (dist < @abs(min_dist)) {
                        min_dist = math.nonZeroSign(cross(d1, qe)) * dist;
                        param.* = t;
                    }
                }
            }

            if (param.* < 0.0)
                return .{
                    .distance = min_dist,
                    .dot = @abs(dot(
                        normal(self.direction(0), true),
                        normal(qa, true),
                    )),
                }
            else if (param.* > 1.0)
                return .{
                    .distance = min_dist,
                    .dot = @abs(dot(
                        normal(self.direction(1), true),
                        normal(p[3] - origin, true),
                    )),
                };
            return .{ .distance = min_dist };
        },
    }
}

pub fn scanlineIntersections(self: EdgeSegment, x: *[3]f64, dy: *[3]i32, y: f64) u32 {
    switch (self.segment) {
        .linear => |p| {
            if (y >= p[0][1] and y < p[1][1] or y >= p[1][1] and y < p[0][1]) {
                const param = (y - p[0][1]) / (p[1][1] - p[0][1]);
                x[0] = mix(p[0][0], p[1][0], param);
                dy[0] = @intFromFloat(std.math.sign(p[1][1] - p[0][1]));
                return 1;
            }
            return 0;
        },
        .quadratic_bezier => |p| {
            var total: u32 = 0;
            var next_dy: i32 = math.boolSign(y > p[0][1]);
            x[total] = p[0][0];
            if (p[0][1] == y) {
                if (p[0][1] < p[1][1] or p[0][1] == p[1][1] and p[0][1] < p[2][1]) {
                    dy[total] = 1;
                    total += 1;
                } else next_dy = 1;
            }

            const ab = p[1] - p[0];
            const br = p[2] - p[1] - ab;
            var roots: [2]f64 = undefined;
            const num_solutions = equations.solveQuadratic(&roots, br[1], 2 * ab[1], p[0][1] - y);
            if (num_solutions >= 2 and roots[0] > roots[1]) std.mem.swap(f64, &roots[0], &roots[1]);
            for (roots[0..num_solutions]) |root| if (root >= 0 and root <= 1) {
                x[total] = p[0][0] + 2 * root * ab[0] + root * root * br[0];
                if (@as(f64, @floatFromInt(next_dy)) * (ab[1] + root * br[1]) >= 0) {
                    dy[total] = next_dy;
                    total += 1;
                    next_dy = -next_dy;
                }
            };

            if (p[2][1] == y) {
                if (next_dy > 0 and total > 0) {
                    total -= 1;
                    next_dy = -1;
                }
                if ((p[2][1] < p[1][1] or p[2][1] == p[1][1] and p[2][1] < p[0][1]) and total < 2) {
                    x[total] = p[2][0];
                    if (next_dy < 0) {
                        dy[total] = -1;
                        total += 1;
                        next_dy = 1;
                    }
                }
            }

            if (next_dy != math.boolSign(y >= p[2][1])) {
                if (total > 0)
                    total -= 1
                else {
                    if (@abs(p[2][1] - y) < @abs(p[0][1] - y)) x[total] = p[2][0];
                    dy[total] = next_dy;
                    total += 1;
                }
            }

            return total;
        },
        .cubic_bezier => |p| {
            var total: u32 = 0;
            var next_dy: i32 = math.boolSign(y > p[0][1]);
            x[total] = p[0][0];
            if (p[0][1] == y) {
                if (p[0][1] < p[1][1] or (p[0][1] == p[1][1] and (p[0][1] < p[2][1] or (p[0][1] == p[2][1] and p[0][1] < p[3][1])))) {
                    dy[total] = 1;
                    total += 1;
                } else next_dy = 1;
            }

            const ab = p[1] - p[0];
            const br = p[2] - p[1] - ab;
            const as = (p[3] - p[2]) - (p[2] - p[1]) - br;
            var roots: [3]f64 = undefined;
            const num_solutions = equations.solveCubic(&roots, as[1], 3 * br[1], 3 * ab[1], p[0][1] - y);
            if (num_solutions >= 2) {
                if (roots[0] > roots[1]) std.mem.swap(f64, &roots[0], &roots[1]);

                if (num_solutions >= 3 and roots[1] > roots[2]) {
                    std.mem.swap(f64, &roots[1], &roots[2]);
                    if (roots[0] > roots[1]) std.mem.swap(f64, &roots[0], &roots[1]);
                }
            }

            for (roots[0..num_solutions]) |root| if (root >= 0 and root <= 1) {
                x[total] = p[0][0] + 3 * root * ab[0] + 3 * root * root * br[0] + root * root * root * as[0];
                if (@as(f64, @floatFromInt(next_dy)) * (ab[1] + 2 * root * br[1] + root * root * as[1]) >= 0) {
                    dy[total] = next_dy;
                    total += 1;
                    next_dy = -next_dy;
                }
            };

            if (p[3][1] == y) {
                if (next_dy > 0 and total > 0) {
                    total -= 1;
                    next_dy = -1;
                }

                if ((p[3][1] < p[2][1] or (p[3][1] == p[2][1] and (p[3][1] < p[1][1] or (p[3][1] == p[1][1] and p[3][1] < p[0][1])))) and total < 3) {
                    x[total] = p[3][0];
                    if (next_dy < 0) {
                        dy[total] = -1;
                        total += 1;
                        next_dy = 1;
                    }
                }
            }

            if (next_dy != math.boolSign(y >= p[3][1])) {
                if (total > 0)
                    total -= 1
                else {
                    if (@abs(p[3][1] - y) < @abs(p[0][1] - y)) x[total] = p[3][0];
                    dy[total] = next_dy;
                    total += 1;
                }
            }

            return total;
        },
    }
}

fn pointBounds(p: Vec2, l: *f64, b: *f64, r: *f64, t: *f64) void {
    const x = p[0];
    const y = p[1];
    if (x < l.*) l.* = x;
    if (y < b.*) b.* = y;
    if (x > r.*) r.* = x;
    if (y > t.*) t.* = y;
}

pub fn bound(self: EdgeSegment, l: *f64, b: *f64, r: *f64, t: *f64) void {
    switch (self.segment) {
        .linear => |p| {
            pointBounds(p[0], l, b, r, t);
            pointBounds(p[1], l, b, r, t);
        },
        .quadratic_bezier => |p| {
            pointBounds(p[0], l, b, r, t);
            pointBounds(p[2], l, b, r, t);
            const bot = (p[1] - p[0]) - (p[2] - p[1]);
            if (bot[0] != 0.0) {
                const param = (p[1][0] - p[0][0]) / bot[0];
                if (param > 0 and param < 1) pointBounds(self.point(param), l, b, r, t);
            }
            if (bot[1] != 0.0) {
                const param = (p[1][1] - p[0][1]) / bot[1];
                if (param > 0 and param < 1) pointBounds(self.point(param), l, b, r, t);
            }
        },
        .cubic_bezier => |p| {
            pointBounds(p[0], l, b, r, t);
            pointBounds(p[3], l, b, r, t);
            const a0 = p[1] - p[0];
            const a1 = (p[2] - p[1] - a0) * v2(2.0);
            const a2 = p[3] - (p[2] * v2(3.0) + p[1] * v2(3.0)) - p[0];
            var roots: [2]f64 = undefined;
            var roots_len = equations.solveQuadratic(&roots, a2[0], a1[0], a0[0]);
            for (roots[0..roots_len]) |root| if (root > 0 and root < 1) pointBounds(self.point(root), l, b, r, t);
            roots_len = equations.solveQuadratic(&roots, a2[1], a1[1], a0[1]);
            for (roots[0..roots_len]) |root| if (root > 0 and root < 1) pointBounds(self.point(root), l, b, r, t);
        },
    }
}

pub fn reverse(self: *EdgeSegment) void {
    switch (self.segment) {
        .linear => |*p| std.mem.swap(Vec2, &p[0], &p[1]),
        .quadratic_bezier => |*p| std.mem.swap(Vec2, &p[0], &p[2]),
        .cubic_bezier => |*p| {
            std.mem.swap(Vec2, &p[0], &p[3]);
            std.mem.swap(Vec2, &p[1], &p[2]);
        },
    }
}

pub fn splitInThirds(self: EdgeSegment, out_p: *[3]EdgeSegment) void {
    switch (self.segment) {
        .linear => |p| {
            out_p[0] = createLinear(p[0], self.point(1.0 / 3.0), self.color);
            out_p[1] = createLinear(self.point(1.0 / 3.0), self.point(2.0 / 3.0), self.color);
            out_p[2] = createLinear(self.point(2.0 / 3.0), p[1], self.color);
        },
        .quadratic_bezier => |p| {
            out_p[0] = createQuadratic(p[0], mix(p[0], p[1], 1.0 / 3.0), self.point(1.0 / 3.0), self.color);
            out_p[1] = createQuadratic(
                self.point(1.0 / 3.0),
                mix(
                    mix(p[0], p[1], 5.0 / 9.0),
                    mix(p[1], p[2], 4.0 / 9.0),
                    0.5,
                ),
                self.point(2.0 / 3.0),
                self.color,
            );
            out_p[2] = createQuadratic(self.point(2.0 / 3.0), mix(p[1], p[2], 2.0 / 3.0), p[2], self.color);
        },
        .cubic_bezier => |p| {
            out_p[0] = createCubic(
                p[0],
                if (std.meta.eql(p[0], p[1])) p[0] else mix(p[0], p[1], 1.0 / 3.0),
                mix(
                    mix(p[0], p[1], 1.0 / 3.0),
                    mix(p[1], p[2], 1.0 / 3.0),
                    1.0 / 3.0,
                ),
                self.point(1.0 / 3.0),
                self.color,
            );
            out_p[1] = createCubic(
                self.point(1.0 / 3.0),
                mix(
                    mix(
                        mix(p[0], p[1], 1.0 / 3.0),
                        mix(p[1], p[2], 1.0 / 3.0),
                        1.0 / 3.0,
                    ),
                    mix(
                        mix(p[1], p[2], 1.0 / 3.0),
                        mix(p[2], p[3], 1.0 / 3.0),
                        1.0 / 3.0,
                    ),
                    2.0 / 3.0,
                ),
                mix(
                    mix(
                        mix(p[0], p[1], 2.0 / 3.0),
                        mix(p[1], p[2], 2.0 / 3.0),
                        2.0 / 3.0,
                    ),
                    mix(
                        mix(p[1], p[2], 2.0 / 3.0),
                        mix(p[2], p[3], 2.0 / 3.0),
                        2.0 / 3.0,
                    ),
                    1.0 / 3.0,
                ),
                self.point(2.0 / 3.0),
                self.color,
            );
            out_p[2] = createCubic(
                self.point(2.0 / 3.0),
                mix(
                    mix(p[1], p[2], 2.0 / 3.0),
                    mix(p[2], p[3], 2.0 / 3.0),
                    2.0 / 3.0,
                ),
                if (std.meta.eql(p[2], p[3])) p[3] else mix(p[2], p[3], 2.0 / 3.0),
                p[3],
                self.color,
            );
        },
    }
}

pub fn convertToCubic(self: *EdgeSegment) void {
    if (self.segment != .quadratic_bezier) @panic("This function is only supported on quadratic beziers");
    const p = self.segment.quadratic_bezier;
    self.* = createCubic(
        p[0],
        mix(p[0], p[1], 2.0 / 3.0),
        mix(p[1], p[2], 1.0 / 3.0),
        p[2],
        self.color,
    );
}

pub fn deconverge(self: *EdgeSegment, param: u32, vector: Vec2) void {
    switch (self.segment) {
        .linear => @panic("Deconverging an edge is only supported on quadratic and cubic beziers"),
        .cubic_bezier => {},
        .quadratic_bezier => self.convertToCubic(),
    }

    const p = self.segment.cubic_bezier;
    switch (param) {
        0 => self.segment.cubic_bezier[1] = p[1] + v2(math.length(vector * (p[1] - p[0]))),
        1 => self.segment.cubic_bezier[2] = p[2] + v2(math.length(vector * (p[2] - p[3]))),
        else => @panic("Unsupported operation"),
    }
}
