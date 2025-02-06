const std = @import("std");

pub var rng: std.Random.DefaultPrng = .init(0);

pub const EdgeColor = enum(u8) {
    black = 0,
    red = 1,
    green = 2,
    yellow = 3,
    blue = 4,
    magenta = 5,
    cyan = 6,
    white = 7,

    pub fn init() EdgeColor {
        const default_colors: [3]EdgeColor = .{ .cyan, .magenta, .yellow };
        return default_colors[rng.next() % 3];
    }

    pub fn random(self: *EdgeColor) void {
        const shifted = @intFromEnum(self.*) << @as(u3, @intCast(1 + rng.next() % 2));
        self.* = @enumFromInt((shifted | shifted >> 3) & @intFromEnum(EdgeColor.white));
    }

    pub fn change(self: *EdgeColor, banned: EdgeColor) void {
        const combined: EdgeColor = @enumFromInt(@intFromEnum(self.*) & @intFromEnum(banned));
        switch (combined) {
            .red, .green, .blue => self.* = @enumFromInt(@intFromEnum(combined) ^ @intFromEnum(EdgeColor.white)),
            else => self.random(),
        }
    }
};
