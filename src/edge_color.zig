const std = @import("std");

pub var rng: std.Random.DefaultPrng = .init(0);

pub const EdgeColor = enum(u8) {
    const len = @typeInfo(@This()).@"enum".fields.len;

    black = 0,
    red = 1,
    green = 2,
    yellow = 3,
    blue = 4,
    magenta = 5,
    cyan = 6,
    white = 7,

    pub fn init() EdgeColor {
        const two_channel: [3]EdgeColor = .{ .cyan, .magenta, .yellow };
        return two_channel[rng.next() % two_channel.len];
    }

    fn cmyWithExcls(excl_1: EdgeColor, excl_2: EdgeColor) EdgeColor {
        const T = EdgeColor;
        const cmy_tuple = .{ T.cyan, T.magenta, T.yellow };

        var validation: [2]bool = @splat(false);
        for (cmy_tuple) |c| {
            if (c == excl_1) validation[0] = true;
            if (c == excl_2) validation[1] = true;
        }
        if (!std.meta.eql(validation, @splat(true)))
            @compileError("You must provide CMY exclusion values");

        for (cmy_tuple) |c| if (c != excl_1 and c != excl_2)
            return c;
    }

    // Manually specifying each with the current values
    // seems to yield a better coloring than the old ordered one
    fn colorSwitchTargets(self: EdgeColor) [2]EdgeColor {
        return switch (self) {
            .black, .white => @compileError("Switching black/white colors is not supported"),
            .red => .{ .green, .blue },
            .green => .{ .blue, .red },
            .blue => .{ .red, .green },
            .cyan => .{ .magenta, .yellow },
            .magenta => .{ .yellow, .cyan },
            .yellow => .{ .cyan, .magenta },
        };
    }

    pub fn change(self: *EdgeColor, banned: EdgeColor) void {
        switch (self.*) {
            inline .cyan, .magenta, .yellow => |c| switch (banned) {
                inline .cyan, .magenta, .yellow => |bc| if (comptime c != bc) {
                    self.* = comptime cmyWithExcls(c, bc);
                    return;
                },
                else => {},
            },
            else => {},
        }

        self.random();
    }

    pub fn random(self: *EdgeColor) void {
        switch (self.*) {
            .black, .white => return,
            inline else => |c| {
                const switch_targets = comptime colorSwitchTargets(c);
                self.* = switch_targets[rng.next() % switch_targets.len];
            },
        }
    }
};
