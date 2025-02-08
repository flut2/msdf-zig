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

    pub fn random(self: *EdgeColor) void {
        switch (self.*) {
            .black, .white => return,
            .red => {
                const one_channel_minus_red: [2]EdgeColor = .{ .green, .blue };
                self.* = one_channel_minus_red[rng.next() % one_channel_minus_red.len];
            },
            .green => {
                const one_channel_minus_green: [2]EdgeColor = .{ .blue, .red };
                self.* = one_channel_minus_green[rng.next() % one_channel_minus_green.len];
            },
            .blue => {
                const one_channel_minus_blue: [2]EdgeColor = .{ .red, .green };
                self.* = one_channel_minus_blue[rng.next() % one_channel_minus_blue.len];
            },
            .yellow => {
                const two_channel_minus_yellow: [2]EdgeColor = .{ .cyan, .magenta };
                self.* = two_channel_minus_yellow[rng.next() % two_channel_minus_yellow.len];
            },
            .magenta => {
                const two_channel_minus_magenta: [2]EdgeColor = .{ .yellow, .cyan };
                self.* = two_channel_minus_magenta[rng.next() % two_channel_minus_magenta.len];
            },
            .cyan => {
                const two_channel_minus_cyan: [2]EdgeColor = .{ .magenta, .yellow };
                self.* = two_channel_minus_cyan[rng.next() % two_channel_minus_cyan.len];
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
