const std = @import("std");

const SignedDistance = @This();

distance: f64 = std.math.floatMax(f64),
dot: f64 = 0,

pub fn lt(self: SignedDistance, other: SignedDistance) bool {
    return @abs(self.distance) < @abs(other.distance) or (@abs(self.distance) == @abs(other.distance) and self.dot < other.dot);
}

pub fn gt(self: SignedDistance, other: SignedDistance) bool {
    return @abs(self.distance) > @abs(other.distance) or (@abs(self.distance) == @abs(other.distance) and self.dot > other.dot);
}

pub fn le(self: SignedDistance, other: SignedDistance) bool {
    return @abs(self.distance) <= @abs(other.distance) or (@abs(self.distance) == @abs(other.distance) and self.dot <= other.dot);
}

pub fn ge(self: SignedDistance, other: SignedDistance) bool {
    return @abs(self.distance) >= @abs(other.distance) or (@abs(self.distance) == @abs(other.distance) and self.dot >= other.dot);
}
