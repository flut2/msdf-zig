const Vec2 = @import("Vec2.zig");

pub fn median(a: anytype, b: anytype, c: anytype) @TypeOf(a) {
    return @max(@min(a, b), @min(@max(a, b), c));
}

pub fn mix(a: anytype, b: anytype, weight: anytype) @TypeOf(a) {
    if (@TypeOf(a) == Vec2) {
        const one_vec: Vec2 = .{ .x = 1.0, .y = 1.0 };
        return one_vec.sub(weight).mul(a).add(b.mul(weight));
    }
    return (1 - weight) * a + weight * b;
}

pub fn clampNorm(n: anytype) @TypeOf(n) {
    if (n >= 0 and n <= 1) return n;
    return if (n > 0) 1 else 0;
}

pub fn clampMax(n: anytype, b: anytype) @TypeOf(n) {
    if (n >= 0 and n <= b) return n;
    return if (n > 0) b else 0;
}

pub fn clampRange(n: anytype, a: anytype, b: anytype) @TypeOf(n) {
    if (n >= a and n <= b) return n;
    return if (n < a) a else b;
}

pub fn sign(n: anytype) @TypeOf(n) {
    if (n > 0) return 1;
    return if (n < 0) -1 else 0;
}

pub fn nonZeroSign(n: anytype) @TypeOf(n) {
    return if (n > 0) 1 else -1;
}
