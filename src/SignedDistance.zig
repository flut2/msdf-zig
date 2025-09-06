const std = @import("std");

const SignedDistance = @This();

distance: f64 = std.math.floatMax(f64),
dot: f64 = 0,

pub fn lessThan(self: SignedDistance, other: SignedDistance) bool {
    return @abs(self.distance) < @abs(other.distance) or
        (@abs(self.distance) == @abs(other.distance) and self.dot < other.dot);
}

pub fn format(self: SignedDistance, writer: *std.Io.Writer) std.Io.Writer.Error!void {
    try writer.print("dist={d:.2}, dot={d:.2}", .{ self.distance, self.dot });
}
