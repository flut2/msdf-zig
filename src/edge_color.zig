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

    fn channelValuesWithExcl(excl: EdgeColor, base_vals: anytype) [2]EdgeColor {
        var ret: [2]EdgeColor = undefined;
        var i: usize = 0;
        for (base_vals) |c| if (c != excl) {
            ret[i] = c;
            i += 1;
        };
        return ret;
    }

    pub fn random(self: *EdgeColor) void {
        const T = EdgeColor;
        switch (self.*) {
            .black, .white => return,
            inline else => |c| {
                const vals_with_excl = comptime channelValuesWithExcl(c, switch (c) {
                    .red, .green, .blue => .{ T.red, T.green, T.blue },
                    .cyan, .magenta, .yellow => .{ T.cyan, T.magenta, T.yellow },
                    else => unreachable,
                });
                self.* = vals_with_excl[rng.next() % vals_with_excl.len];
            },
        }
    }

    pub fn change(self: *EdgeColor, banned: EdgeColor) void {
        const combined = @intFromEnum(self.*) & @intFromEnum(banned);
        switch (@as(EdgeColor, @enumFromInt(combined))) {
            .red, .green, .blue => self.* = @enumFromInt(EdgeColor.len - 1 - combined),
            else => self.random(),
        }
    }
};
