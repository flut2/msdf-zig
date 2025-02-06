const Vec2 = @This();

x: f64 = 0.0,
y: f64 = 0.0,

pub fn lengthSqr(self: Vec2) f64 {
    return self.x * self.x + self.y * self.y;
}

pub fn length(self: Vec2) f64 {
    return @sqrt(self.lengthSqr());
}

pub fn normal(self: Vec2, allow_zero: bool) Vec2 {
    const len = self.length();
    if (len != 0.0) return .{ .x = self.x / len, .y = self.y / len };
    return .{ .y = @floatFromInt(@intFromBool(!allow_zero)) };
}

pub fn ortho(self: Vec2, polarity: bool) Vec2 {
    return if (polarity) .{ .x = -self.y, .y = self.x } else .{ .x = self.y, .y = -self.x };
}

pub fn orthonormal(self: Vec2, polarity: bool, allow_zero: bool) Vec2 {
    const len = self.length();
    if (len != 0.0) return if (polarity) .{ .x = -self.y / len, .y = self.x / len } else .{ .x = self.y / len, .y = -self.x / len };
    return if (polarity) 
        .{ .y = @floatFromInt(@intFromBool(!allow_zero)) } 
    else 
        .{ .y = -@as(f64, @floatFromInt(@intFromBool(!allow_zero))) };
}

pub fn zero(self: Vec2) bool {
    return self.x == 0.0 and self.y == 0.0;
}

pub fn eql(self: Vec2, other: Vec2) bool {
    return self.x == other.x and self.y == other.y;
}

pub fn add(self: Vec2, other: anytype) Vec2 {
    switch (@TypeOf(other)) {
        Vec2, *Vec2 => return .{ .x = self.x + other.x, .y = self.y + other.y },
        comptime_float, f64, f32 => return .{ .x = self.x + other, .y = self.y + other },
        else => @compileError("Unsupported type"),
    }
}

pub fn sub(self: Vec2, other: anytype) Vec2 {
    switch (@TypeOf(other)) {
        Vec2, *Vec2 => return .{ .x = self.x - other.x, .y = self.y - other.y },
        comptime_float, f64, f32 => return .{ .x = self.x - other, .y = self.y - other },
        else => @compileError("Unsupported type"),
    }
}

pub fn mul(self: Vec2, other: anytype) Vec2 {
    switch (@TypeOf(other)) {
        Vec2, *Vec2 => return .{ .x = self.x * other.x, .y = self.y * other.y },
        comptime_float, f64, f32 => return .{ .x = self.x * other, .y = self.y * other },
        else => @compileError("Unsupported type"),
    }
}

pub fn div(self: Vec2, other: anytype) Vec2 {
    switch (@TypeOf(other)) {
        Vec2, *Vec2 => return .{ .x = self.x / other.x, .y = self.y / other.y },
        comptime_float, f64, f32 => return .{ .x = self.x / other, .y = self.y / other },
        else => @compileError("Unsupported type"),
    }
}

pub fn dot(self: Vec2, other: Vec2) f64 {
    return self.x * other.x + self.y * other.y;
}

pub fn cross(self: Vec2, other: Vec2) f64 {
    return self.x * other.y - self.y * other.x;
}
