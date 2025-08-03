const std = @import("std");

const EdgeColor = @import("edge_color.zig").EdgeColor;
const equations = @import("equations.zig");
const math = @import("math.zig");
const mix = math.mix;
const nonZeroSign = math.nonZeroSign;
const SignedDistance = @import("SignedDistance.zig");
const Vec2 = @import("Vec2.zig");

const EdgeSegment = @This();

const cubic_starts = 4;
const cubic_steps = 4;

color: EdgeColor = .white,
segment: union(enum) {
    linear: [2]Vec2,
    quadratic_bezier: [3]Vec2,
    cubic_bezier: [4]Vec2,
} = .{ .linear = @splat(.{ .x = 0.0, .y = 0.0 }) },

pub fn format(self: EdgeSegment, writer: *std.io.Writer) std.io.Writer.Error!void {
    switch (self.segment) {
        .linear => |vals| try writer.print("linear, color={t}, p1={f}, p2={f}", .{ self.color, vals[0], vals[1] }),
        .quadratic_bezier => |vals| try writer.print("quadratic, color={t}, p1={f}, p2={f}, p3={f}", .{ self.color, vals[0], vals[1], vals[2] }),
        .cubic_bezier => |vals| try writer.print("cubic, color={t}, p1={f}, p2={f}, p3={f}, p4={f}", .{ self.color, vals[0], vals[1], vals[2], vals[3] }),
    }
}

pub fn create(p0: Vec2, p1: Vec2, p2: ?Vec2, p3: ?Vec2, color: EdgeColor) EdgeSegment {
    if (p3 != null) {
        if (p2 == null) @panic("Invalid parameters, you need to specify `p2` if you specify `p3`.");

        var p12: Vec2 = p2.?.sub(p1);
        if (p1.sub(p0).cross(p12) == 0.0 and p12.cross(p3.?.sub(p2.?)) == 0.0)
            return .{
                .color = color,
                .segment = .{ .linear = .{ p0, p3.? } },
            };
        p12 = p1.mul(1.5).sub(p0.mul(0.5));
        if (p12.eql(p2.?.mul(1.5).sub(p3.?.mul(0.5))))
            return .{
                .color = color,
                .segment = .{ .quadratic_bezier = .{ p0, p12, p3.? } },
            };

        return .{
            .color = color,
            .segment = .{ .cubic_bezier = .{ p0, p1, p2.?, p3.? } },
        };
    }

    if (p2 != null) {
        if (p1.sub(p0).cross(p2.?.sub(p1)) == 0.0)
            return .{
                .color = color,
                .segment = .{ .linear = .{ p0, p2.? } },
            };

        return .{
            .color = color,
            .segment = .{ .quadratic_bezier = .{ p0, p1, p2.? } },
        };
    }

    return .{
        .color = color,
        .segment = .{ .linear = .{ p0, p1 } },
    };
}

pub fn distanceToPerpendicularDistance(self: EdgeSegment, distance: *SignedDistance, origin: Vec2, param: f64) void {
    if (param < 0) {
        const dir = self.direction(0).normal(false);
        const aq = origin.sub(self.point(0));
        const ts = aq.dot(dir);
        if (ts < 0) {
            const perp_dist = aq.cross(dir);
            if (@abs(perp_dist) <= @abs(distance.distance)) {
                distance.distance = perp_dist;
                distance.dot = 0;
            }
        }
    } else if (param > 1) {
        const dir = self.direction(1).normal(false);
        const bq = origin.sub(self.point(1));
        const ts = bq.dot(dir);
        if (ts > 0) {
            const perp_dist = bq.cross(dir);
            if (@abs(perp_dist) <= @abs(distance.distance)) {
                distance.distance = perp_dist;
                distance.dot = 0;
            }
        }
    }
}

pub fn point(self: EdgeSegment, param: f64) Vec2 {
    return switch (self.segment) {
        .linear => |p| mix(p[0], p[1], param),
        .quadratic_bezier => |p| mix(mix(p[0], p[1], param), mix(p[1], p[2], param), param),
        .cubic_bezier => |p| blk: {
            const p12 = mix(p[1], p[2], param);
            break :blk mix(mix(mix(p[0], p[1], param), p12, param), mix(p12, mix(p[2], p[3], param), param), param);
        },
    };
}

pub fn direction(self: EdgeSegment, param: f64) Vec2 {
    switch (self.segment) {
        .linear => |p| return p[1].sub(p[0]),
        .quadratic_bezier => |p| {
            const tangent = mix(p[1].sub(p[0]), p[2].sub(p[1]), param);
            if (tangent.zero()) return p[2].sub(p[0]);
            return tangent;
        },
        .cubic_bezier => |p| {
            const tangent = mix(mix(p[1].sub(p[0]), p[2].sub(p[1]), param), mix(p[2].sub(p[1]), p[3].sub(p[2]), param), param);
            if (tangent.zero()) {
                if (param == 0) return p[2].sub(p[0]);
                if (param == 1) return p[3].sub(p[1]);
            }
            return tangent;
        },
    }
}

pub fn directionChange(self: EdgeSegment, param: f64) Vec2 {
    switch (self.segment) {
        .linear => return .{},
        .quadratic_bezier => |p| return p[2].sub(p[1]).sub(p[1].sub(p[0])),
        .cubic_bezier => |p| return mix(p[2].sub(p[1]).sub(p[1].sub(p[0])), p[3].sub(p[2]).sub(p[2].sub(p[1])), param),
    }
}

pub fn length(self: EdgeSegment) f64 {
    switch (self.segment) {
        .linear => |p| p[1].sub(p[0]).length(),
        .quadratic_bezier => |p| {
            const ab = p[1].sub(p[0]);
            const br = p[2].sub(p[1]).sub(ab);
            const abab = ab.dot(ab);
            const abbr = ab.dot(br);
            const brbr = br.dot(br);
            const ab_len = @sqrt(abab);
            const br_len = @sqrt(brbr);
            const crs = ab.cross(br);
            const h = @sqrt(abab + abbr + abbr + brbr);
            // zig fmt: off
            return (
                br_len * ((abbr + brbr) * h - abbr * ab_len) +
                crs * crs * @log((br_len * h + abbr + brbr) / (br_len * ab_len + abbr))
            ) / (brbr * br_len);
            // zig fmt: on
        },
        .cubic_bezier => @panic("Length calculation for cubic beziers is not implemented"),
    }
}

pub fn signedDistance(self: EdgeSegment, origin: Vec2, param: *f64) SignedDistance {
    switch (self.segment) {
        .linear => |p| {
            const aq = origin.sub(p[0]);
            const ab = p[1].sub(p[0]);
            const new_param = aq.dot(ab) / ab.dot(ab);
            param.* = new_param;
            const eq = p[@intFromBool(new_param > 0.5)].sub(origin);
            const endpoint_dist = eq.length();
            if (new_param > 0.0 and new_param < 1.0) {
                const ortho_dist = ab.orthonormal(false, false).dot(aq);
                if (@abs(ortho_dist) < endpoint_dist) return .{ .distance = ortho_dist, .dot = 0 };
            }
            return .{ .distance = nonZeroSign(aq.cross(ab)) * endpoint_dist, .dot = @abs(ab.normal(false).dot(eq.normal(false))) };
        },
        .quadratic_bezier => |p| {
            const qa = p[0].sub(origin);
            const ab = p[1].sub(p[0]);
            const br = p[2].sub(p[1]).sub(ab);
            const a = br.dot(br);
            const b = 3.0 * ab.dot(br);
            const c = 2.0 * ab.dot(ab) + qa.dot(br);
            const d = qa.dot(ab);
            var roots: [3]f64 = undefined;
            const num_solutions = equations.solveCubic(&roots, a, b, c, d);

            var ep_dir = self.direction(0);
            var min_dist = nonZeroSign(ep_dir.cross(qa)) * qa.length();
            param.* = -qa.dot(ep_dir) / ep_dir.dot(ep_dir);
            {
                ep_dir = self.direction(1);
                const dist = p[2].sub(origin).length();
                if (dist < @abs(min_dist)) {
                    min_dist = nonZeroSign(ep_dir.cross(p[2].sub(origin))) * dist;
                    param.* = origin.sub(p[1]).dot(ep_dir) / ep_dir.dot(ep_dir);
                }
            }

            for (roots[0..num_solutions]) |root| if (root > 0 and root < 1) {
                const qe = qa.add(ab.mul(root * 2.0)).add(br.mul(root * root));
                const dist = qe.length();
                if (dist < @abs(min_dist)) {
                    min_dist = nonZeroSign(ab.add(br.mul(root)).cross(qe)) * dist;
                    param.* = root;
                }
            };

            if (param.* >= 0 and param.* <= 1) return .{ .distance = min_dist, .dot = 0 };
            if (param.* < 0.5)
                return .{ .distance = min_dist, .dot = @abs(self.direction(0).normal(false).dot(qa.normal(false))) }
            else
                return .{ .distance = min_dist, .dot = @abs(self.direction(1).normal(false).dot(p[2].sub(origin).normal(false))) };
        },
        .cubic_bezier => |p| {
            const qa = p[0].sub(origin);
            const ab = p[1].sub(p[0]);
            const br = p[2].sub(p[1]).sub(ab);
            const as = p[3].sub(p[2]).sub(p[2].sub(p[1])).sub(br);

            var ep_dir = self.direction(0);
            var min_dist = nonZeroSign(ep_dir.cross(qa)) * qa.length();
            param.* = -qa.dot(ep_dir) / ep_dir.dot(ep_dir);
            {
                ep_dir = self.direction(1);
                const dist = p[3].sub(origin).length();
                if (dist < @abs(min_dist)) {
                    min_dist = nonZeroSign(ep_dir.cross(p[3].sub(origin))) * dist;
                    param.* = ep_dir.sub(p[3].sub(origin)).dot(ep_dir) / ep_dir.dot(ep_dir);
                }
            }
            for (0..cubic_starts) |i| {
                const fi: f64 = @floatFromInt(i);
                var t = fi / cubic_starts;
                var qe = qa.add(ab.mul(t * 3.0)).add(br.mul(t * t * 3.0)).add(as.mul(t * t * t));
                for (0..cubic_steps) |_| {
                    const d1 = ab.mul(3.0).add(br.mul(t * 6.0)).add(as.mul(t * t * 3.0));
                    const d2 = br.mul(6.0).add(as.mul(t * 6.0));
                    t -= qe.dot(d1) / (d1.dot(d1) + qe.dot(d2));
                    if (t <= 0 or t >= 1) break;
                    qe = qa.add(ab.mul(t * 3.0)).add(br.mul(t * t * 3.0)).add(as.mul(t * t * t));
                    const dist = qe.length();
                    if (dist < @abs(min_dist)) {
                        min_dist = nonZeroSign(d1.cross(qe)) * dist;
                        param.* = t;
                    }
                }
            }

            if (param.* >= 0 and param.* <= 1) return .{ .distance = min_dist, .dot = 0 };
            if (param.* < 0.5)
                return .{ .distance = min_dist, .dot = @abs(self.direction(0).normal(false).dot(qa.normal(false))) }
            else
                return .{ .distance = min_dist, .dot = @abs(self.direction(1).normal(false).dot(p[3].sub(origin).normal(false))) };
        },
    }
}

pub fn scanlineIntersections(self: EdgeSegment, x: *[3]f64, dy: *[3]i32, y: f64) u32 {
    switch (self.segment) {
        .linear => |p| {
            if (y >= p[0].y and y < p[1].y or y >= p[1].y and y < p[0].y) {
                const param = (y - p[0].y) / (p[1].y - p[0].y);
                x[0] = mix(p[0].x, p[1].x, param);
                dy[0] = @intFromFloat(math.sign(p[1].y - p[0].y));
                return 1;
            }
            return 0;
        },
        .quadratic_bezier => |p| {
            var total: u32 = 0;
            var next_dy: i32 = if (y > p[0].y) 1 else -1;
            x[total] = p[0].x;
            if (p[0].y == y) {
                if (p[0].y < p[1].y or p[0].y == p[1].y and p[0].y < p[2].y) {
                    dy[total] = 1;
                    total += 1;
                } else next_dy = 1;
            }
            {
                const ab = p[1].sub(p[0]);
                const br = p[2].sub(p[1]).sub(ab);
                var roots: [2]f64 = undefined;
                const num_solutions = equations.solveQuadratic(&roots, br.y, 2 * ab.y, p[0].y - y);
                if (num_solutions >= 2 and roots[0] > roots[1]) std.mem.swap(f64, &roots[0], &roots[1]);
                for (roots[0..@min(2, num_solutions)]) |root| {
                    if (root >= 0 and root <= 1) {
                        x[total] = p[0].x + 2 * root * ab.x + root * root * br.x;
                        if (@as(f64, @floatFromInt(next_dy)) * (ab.y + root * br.y) >= 0) {
                            dy[total] = next_dy;
                            total += 1;
                            next_dy = -next_dy;
                        }
                    }
                }
            }

            if (p[2].y == y) {
                if (next_dy > 0 and total > 0) {
                    total -= 1;
                    next_dy = -1;
                }
                if ((p[2].y < p[1].y or p[2].y == p[1].y and p[2].y < p[0].y) and total < 2) {
                    x[total] = p[2].x;
                    if (next_dy < 0) {
                        dy[total] = -1;
                        total += 1;
                        next_dy = 1;
                    }
                }
            }

            if (next_dy != @as(i32, if (y >= p[2].y) 1 else -1)) {
                if (total > 0)
                    total -= 1
                else {
                    if (@abs(p[2].y - y) < @abs(p[0].y - y)) x[total] = p[2].x;
                    dy[total] = next_dy;
                    total += 1;
                }
            }

            return total;
        },
        .cubic_bezier => |p| {
            var total: u32 = 0;
            var next_dy: i32 = if (y > p[0].y) 1 else -1;
            x[total] = p[0].x;
            if (p[0].y == y) {
                if (p[0].y < p[1].y or (p[0].y == p[1].y and (p[0].y < p[2].y or (p[0].y == p[2].y and p[0].y < p[3].y)))) {
                    dy[total] = 1;
                    total += 1;
                } else next_dy = 1;
            }
            {
                const ab = p[1].sub(p[0]);
                const br = p[2].sub(p[1]).sub(ab);
                const as = p[3].sub(p[2]).sub(p[2].sub(p[1])).sub(br);
                var roots: [3]f64 = undefined;
                const num_solutions = equations.solveCubic(&roots, as.y, 3 * br.y, 3 * ab.y, p[0].y - y);
                if (num_solutions >= 2) {
                    if (roots[0] > roots[1]) std.mem.swap(f64, &roots[0], &roots[1]);

                    if (num_solutions >= 3 and roots[1] > roots[2]) {
                        std.mem.swap(f64, &roots[1], &roots[2]);
                        if (roots[0] > roots[1]) std.mem.swap(f64, &roots[0], &roots[1]);
                    }
                }

                for (roots[0..@min(3, num_solutions)]) |root| {
                    if (root >= 0 and root <= 1) {
                        x[total] = p[0].x + 3 * root * ab.x + 3 * root * root * br.x + root * root * root * as.x;
                        if (@as(f64, @floatFromInt(next_dy)) * (ab.y + 2 * root * br.y + root * root * as.y) >= 0) {
                            dy[total] = next_dy;
                            total += 1;
                            next_dy = -next_dy;
                        }
                    }
                }
            }

            if (p[3].y == y) {
                if (next_dy > 0 and total > 0) {
                    total -= 1;
                    next_dy = -1;
                }

                if ((p[3].y < p[2].y or (p[3].y == p[2].y and (p[3].y < p[1].y or (p[3].y == p[1].y and p[3].y < p[0].y)))) and total < 3) {
                    x[total] = p[3].x;
                    if (next_dy < 0) {
                        dy[total] = -1;
                        total += 1;
                        next_dy = 1;
                    }
                }
            }

            if (next_dy != @as(i32, if (y >= p[3].y) 1 else -1)) {
                if (total > 0)
                    total -= 1
                else {
                    if (@abs(p[3].y - y) < @abs(p[0].y - y)) x[total] = p[3].x;
                    dy[total] = next_dy;
                    total += 1;
                }
            }

            return total;
        },
    }
}

fn pointBounds(p: Vec2, l: *f64, b: *f64, r: *f64, t: *f64) void {
    if (p.x < l.*) l.* = p.x;
    if (p.y < b.*) b.* = p.y;
    if (p.x > r.*) r.* = p.x;
    if (p.y > t.*) t.* = p.y;
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
            const bot = p[1].sub(p[0]).sub(p[2].sub(p[1]));
            if (bot.x != 0.0) {
                const param = (p[1].x - p[0].x) / bot.x;
                if (param > 0 and param < 1) pointBounds(self.point(param), l, b, r, t);
            }
            if (bot.y != 0.0) {
                const param = (p[1].y - p[0].y) / bot.y;
                if (param > 0 and param < 1) pointBounds(self.point(param), l, b, r, t);
            }
        },
        .cubic_bezier => |p| {
            pointBounds(p[0], l, b, r, t);
            pointBounds(p[3], l, b, r, t);
            const a0 = p[1].sub(p[0]);
            const a1 = p[2].sub(p[1]).sub(a0).mul(2.0);
            const a2 = p[3].sub(p[2].mul(3.0).add(p[1].mul(3.0)).sub(p[0]));
            var roots: [2]f64 = undefined;
            _ = equations.solveQuadratic(&roots, a2.x, a1.x, a0.x);
            for (roots) |root| if (root > 0 and root < 1) pointBounds(self.point(root), l, b, r, t);
            _ = equations.solveQuadratic(&roots, a2.y, a1.y, a0.y);
            for (roots) |root| if (root > 0 and root < 1) pointBounds(self.point(root), l, b, r, t);
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

pub fn moveStartPoint(self: *EdgeSegment, to: Vec2) void {
    switch (self.segment) {
        .linear => |*p| p[0] = to,
        .quadratic_bezier => |*p| {
            const orig_s_dir = p[0].sub(p[1]);
            const orig_p1 = p[1];
            p[1] = p[1].add(p[2].sub(p[1]).mul(p[0].sub(p[1]).cross(to.sub(p[0])) / p[0].sub(p[1]).cross(p[2].sub(p[1]))));
            p[0] = to;
            if (orig_s_dir.dot(p[0].sub(p[1])) < 0) p[1] = orig_p1;
        },
        .cubic_bezier => |*p| {
            p[1] = p[1].add(to.sub(p[0]));
            p[0] = to;
        },
    }
}

pub fn moveEndPoint(self: *EdgeSegment, to: Vec2) void {
    switch (self.segment) {
        .linear => |*p| p[1] = to,
        .quadratic_bezier => |*p| {
            const orig_e_dir = p[2].sub(p[1]);
            const orig_p1 = p[1];
            p[1] = p[1].add(p[0].sub(p[1]).mul(p[2].sub(p[1]).cross(to.sub(p[2])) / p[2].sub(p[1]).cross(p[0].sub(p[1]))));
            p[2] = to;
            if (orig_e_dir.dot(p[2].sub(p[1])) < 0) p[1] = orig_p1;
        },
        .cubic_bezier => |*p| {
            p[2] = p[2].add(to.sub(p[3]));
            p[3] = to;
        },
    }
}

pub fn splitInThirds(self: EdgeSegment, p0: *EdgeSegment, p1: *EdgeSegment, p2: *EdgeSegment) void {
    const one_third = 1.0 / 3.0;
    const two_thirds = 2.0 / 3.0;
    switch (self.segment) {
        .linear => |p| {
            p0.* = .{
                .color = self.color,
                .segment = .{ .linear = .{ p[0], self.point(one_third) } },
            };
            p1.* = .{
                .color = self.color,
                .segment = .{ .linear = .{ self.point(one_third), self.point(two_thirds) } },
            };
            p2.* = .{
                .color = self.color,
                .segment = .{ .linear = .{ self.point(two_thirds), p[1] } },
            };
        },
        .quadratic_bezier => |p| {
            p0.* = .{
                .color = self.color,
                .segment = .{ .quadratic_bezier = .{ p[0], mix(p[0], p[1], one_third), self.point(one_third) } },
            };
            p1.* = .{
                .color = self.color,
                .segment = .{ .quadratic_bezier = .{
                    self.point(one_third),
                    mix(mix(p[0], p[1], 5.0 / 9.0), mix(p[1], p[2], 4.0 / 9.0), 0.5),
                    self.point(two_thirds),
                } },
            };
            p2.* = .{
                .color = self.color,
                .segment = .{ .quadratic_bezier = .{ self.point(two_thirds), mix(p[1], p[2], 2.0 / 3.0), p[2] } },
            };
        },
        .cubic_bezier => |p| {
            p0.* = .{
                .color = self.color,
                .segment = .{ .cubic_bezier = .{
                    p[0],
                    if (p[0].eql(p[1])) p[0] else mix(p[0], p[1], one_third),
                    mix(mix(p[0], p[1], one_third), mix(p[1], p[2], one_third), one_third),
                    self.point(one_third),
                } },
            };
            p1.* = .{
                .color = self.color,
                .segment = .{ .cubic_bezier = .{
                    self.point(one_third),
                    mix(
                        mix(mix(p[0], p[1], one_third), mix(p[1], p[2], one_third), one_third),
                        mix(mix(p[1], p[2], one_third), mix(p[2], p[3], one_third), one_third),
                        two_thirds,
                    ),
                    mix(
                        mix(mix(p[0], p[1], two_thirds), mix(p[1], p[2], two_thirds), two_thirds),
                        mix(mix(p[1], p[2], two_thirds), mix(p[2], p[3], two_thirds), two_thirds),
                        one_third,
                    ),
                    self.point(two_thirds),
                } },
            };
            p2.* = .{
                .color = self.color,
                .segment = .{ .cubic_bezier = .{
                    self.point(two_thirds),
                    mix(mix(p[1], p[2], two_thirds), mix(p[2], p[3], two_thirds), two_thirds),
                    if (p[2].eql(p[3])) p[3] else mix(p[2], p[3], two_thirds),
                    p[3],
                } },
            };
        },
    }
}

pub fn convertToCubic(self: *EdgeSegment) void {
    if (self.segment != .quadratic_bezier) @panic("This function is only supported on quadratic beziers");
    const p = self.segment.quadratic_bezier;
    self.* = .{
        .color = self.color,
        .segment = .{ .cubic_bezier = .{ p[0], mix(p[0], p[1], 2.0 / 3.0), mix(p[1], p[2], 1.0 / 3.0), p[2] } },
    };
}

pub fn deconverge(self: *EdgeSegment, param: u32, vector: Vec2) void {
    switch (self.segment) {
        .linear => @panic("Deconverging an edge is only supported on quadratic and cubic beziers"),
        .cubic_bezier => {},
        .quadratic_bezier => self.convertToCubic(),
    }

    const p = self.segment.cubic_bezier;
    switch (param) {
        0 => self.segment.cubic_bezier[1] = p[1].add(vector.mul(p[1].sub(p[0]).length())),
        1 => self.segment.cubic_bezier[2] = p[2].add(vector.mul(p[2].sub(p[3]).length())),
        else => @panic("Unsupported operation"),
    }
}
