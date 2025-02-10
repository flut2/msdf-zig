const std = @import("std");

const ft = @import("mach-freetype");

const coloring = @import("coloring.zig");
const Contour = @import("Contour.zig");
const edge_color = @import("edge_color.zig");
const EdgeColor = edge_color.EdgeColor;
const EdgeSegment = @import("EdgeSegment.zig");
const math = @import("math.zig");
const clampNorm = math.clampNorm;
const median = math.median;
const Scanline = @import("Scanline.zig");
const Shape = @import("Shape.zig");
const SignedDistance = @import("SignedDistance.zig");
const Vec2 = @import("Vec2.zig");

const Generator = @This();

pub const FontMetrics = struct {
    line_height: f64,
    ascender: f64,
    descender: f64,
    underline_y: f64,
    underline_thickness: f64,
};

pub const GlyphData = struct {
    advance: f64,
    pixels: []const u8,

    pub fn deinit(self: GlyphData, allocator: std.mem.Allocator) void {
        allocator.free(self.pixels);
    }
};

pub const SdfType = enum {
    sdf,
    psdf,
    msdf,
    mtsdf,
};

pub const GenerationOptions = struct {
    sdf_type: SdfType,
    px_size: u16,
    px_range: u16,
    corner_angle_threshold: f64 = 3.0,
    orientation: enum { guess, keep, reverse } = .reverse,
    geometry_preprocess: bool = false,
    /// Requires geometry preprocessing to be disabled. A ``null`` value also disables this
    scanline_fill_rule: ?Scanline.FillRule = .non_zero,
};

const FreetypeContext = struct {
    allocator: std.mem.Allocator,
    scale: f64,
    shape: *Shape,
    pos: Vec2 = .{ .x = 0.0, .y = 0.0 },
    contour: ?*Contour = null,
};

library: ft.Library = undefined,
face: ft.Face = undefined,

/// ``font_memory`` is the raw font file data
pub fn create(font_memory: []const u8) !Generator {
    var library: ft.Library = try .init();
    return .{
        .library = library,
        .face = try library.createFaceMemory(font_memory, 0),
    };
}

pub fn destroy(self: *Generator) void {
    self.library.deinit();
}

pub fn fontMetrics(self: *Generator) !FontMetrics {
    const scale = 1.0 / f64i(self.face.unitsPerEM());
    return .{
        .line_height = scale * f64i(self.face.height()),
        .ascender = scale * f64i(self.face.ascender()),
        .descender = scale * f64i(self.face.descender()),
        .underline_y = scale * f64i(self.face.underlinePosition()),
        .underline_thickness = scale * f64i(self.face.underlineThickness()),
    };
}

/// The resulting GlyphData is under the caller's ownership
pub fn generate(
    self: *Generator,
    allocator: std.mem.Allocator,
    codepoint: u21,
    opts: GenerationOptions,
) !GlyphData {
    const scale = 1.0 / f64i(self.face.unitsPerEM());
    const glyph_index = self.face.getCharIndex(codepoint) orelse return error.InvalidCodepoint;
    try self.face.loadGlyph(glyph_index, .{ .no_scale = true, .no_bitmap = true });

    var shape: Shape = .{};
    defer {
        for (shape.contours.items) |*contour| contour.edges.deinit(allocator);
        shape.contours.deinit(allocator);
    }

    var context: FreetypeContext = .{
        .allocator = allocator,
        .scale = scale,
        .shape = &shape,
    };

    const outline = self.face.glyph().outline().?;
    try ft.intToError(ft.c.FT_Outline_Decompose(
        outline.handle,
        &ft.c.FT_Outline_Funcs{
            .move_to = ftMoveTo,
            .line_to = ftLineTo,
            .conic_to = ftConicTo,
            .cubic_to = ftCubicTo,
            .shift = 0,
            .delta = 0,
        },
        &context,
    ));

    if (shape.contours.items.len != 0 and shape.contours.getLast().edges.items.len == 0)
        _ = shape.contours.orderedRemove(shape.contours.items.len - 1);

    if (!shape.validate()) return error.InvalidShape;
    if (opts.geometry_preprocess) try shape.orientContours(allocator);
    try shape.normalize(allocator);

    const f_px_size = f64i(opts.px_size);
    const px_range = f64i(opts.px_range) / f_px_size;

    var bounds = shape.getBounds(0, 0, 0);
    if (bounds.left >= bounds.right or bounds.bottom >= bounds.top)
        bounds = .{ .left = 0, .bottom = 0, .right = 1, .top = 1 };

    const translate_x = 0.5 * (1 - (bounds.right - bounds.left)) - bounds.left;
    const translate_y = 0.5 * (1 - (bounds.top - bounds.bottom)) - bounds.bottom;

    const oob_point: Vec2 = if (opts.orientation == .guess)
        .{ .x = bounds.left - (bounds.right - bounds.left) - 1, .y = bounds.bottom - (bounds.top - bounds.bottom) - 1 }
    else
        undefined;

    return .{
        .advance = scale * f64i(self.face.glyph().advance().x),
        .pixels = blk: switch (opts.sdf_type) {
            .sdf => {
                const channels = 1;
                const float_pixels = try allocator.alloc(f64, opts.px_size * opts.px_size * channels);
                defer allocator.free(float_pixels);
                const pixels = try allocator.alloc(u8, opts.px_size * opts.px_size * channels);
                generateSdf(float_pixels, opts.px_size, f_px_size, shape, px_range, translate_x, translate_y);
                if (opts.orientation == .reverse or opts.orientation == .guess and findDistanceAt(shape, oob_point, px_range) > 0)
                    invertPixels(float_pixels);
                if (!opts.geometry_preprocess) if (opts.scanline_fill_rule) |fill_rule| try sdfSignCorrection(
                    allocator,
                    float_pixels,
                    opts.px_size,
                    opts.px_size,
                    shape,
                    f_px_size,
                    translate_x,
                    translate_y,
                    fill_rule,
                );
                convertPixels(float_pixels, pixels, opts.px_size, opts.px_size, channels);
                break :blk pixels;
            },
            .psdf => {
                const channels = 1;
                const float_pixels = try allocator.alloc(f64, opts.px_size * opts.px_size * channels);
                defer allocator.free(float_pixels);
                const pixels = try allocator.alloc(u8, opts.px_size * opts.px_size * channels);
                generatePsdf(float_pixels, opts.px_size, f_px_size, shape, px_range, translate_x, translate_y);
                if (opts.orientation == .reverse or opts.orientation == .guess and findDistanceAt(shape, oob_point, px_range) > 0)
                    invertPixels(float_pixels);
                if (!opts.geometry_preprocess) if (opts.scanline_fill_rule) |fill_rule| try sdfSignCorrection(
                    allocator,
                    float_pixels,
                    opts.px_size,
                    opts.px_size,
                    shape,
                    f_px_size,
                    translate_x,
                    translate_y,
                    fill_rule,
                );
                convertPixels(float_pixels, pixels, opts.px_size, opts.px_size, channels);
                break :blk pixels;
            },
            .msdf => {
                const channels = 3;
                const float_pixels = try allocator.alloc(f64, opts.px_size * opts.px_size * channels);
                defer allocator.free(float_pixels);
                const pixels = try allocator.alloc(u8, opts.px_size * opts.px_size * channels);
                try coloring.simple(allocator, &shape, opts.corner_angle_threshold);
                generateMsdf(float_pixels, opts.px_size, f_px_size, shape, px_range, translate_x, translate_y);
                if (opts.orientation == .reverse or opts.orientation == .guess and findDistanceAt(shape, oob_point, px_range) > 0)
                    invertPixels(float_pixels);
                if (!opts.geometry_preprocess) if (opts.scanline_fill_rule) |fill_rule| try msdfSignCorrection(
                    allocator,
                    float_pixels,
                    opts.px_size,
                    opts.px_size,
                    shape,
                    f_px_size,
                    translate_x,
                    translate_y,
                    fill_rule,
                    channels,
                );
                convertPixels(float_pixels, pixels, opts.px_size, opts.px_size, channels);
                break :blk pixels;
            },
            .mtsdf => {
                const channels = 4;
                const float_pixels = try allocator.alloc(f64, opts.px_size * opts.px_size * channels);
                defer allocator.free(float_pixels);
                const pixels = try allocator.alloc(u8, opts.px_size * opts.px_size * channels);
                try coloring.simple(allocator, &shape, opts.corner_angle_threshold);
                generateMtsdf(float_pixels, opts.px_size, f_px_size, shape, px_range, translate_x, translate_y);
                if (opts.orientation == .reverse or opts.orientation == .guess and findDistanceAt(shape, oob_point, px_range) > 0)
                    invertPixels(float_pixels);
                if (!opts.geometry_preprocess) if (opts.scanline_fill_rule) |fill_rule| try msdfSignCorrection(
                    allocator,
                    float_pixels,
                    opts.px_size,
                    opts.px_size,
                    shape,
                    f_px_size,
                    translate_x,
                    translate_y,
                    fill_rule,
                    channels,
                );
                convertPixels(float_pixels, pixels, opts.px_size, opts.px_size, channels);
                break :blk pixels;
            },
        },
    };
}

fn convertPixels(in_pixels: []f64, out_pixels: []u8, w: u16, h: u16, channels: u8) void {
    for (0..h) |y| for (0..w) |x| {
        const idx = y * w * channels + x * channels;
        for (0..channels) |i| out_pixels[idx + i] = u8f(in_pixels[idx + i]);
    };
}

fn invertPixels(in_pixels: []f64) void {
    for (in_pixels) |*pixel| pixel.* = 1.0 - pixel.*;
}

fn sdfSignCorrection(
    allocator: std.mem.Allocator,
    out_pixels: []f64,
    w: u16,
    h: u16,
    shape: Shape,
    scale: f64,
    tx: f64,
    ty: f64,
    fill_rule: Scanline.FillRule,
) !void {
    var scanline: Scanline = .{};
    defer scanline.intersections.deinit(allocator);
    for (0..h) |y| {
        const row = h - y - 1;
        try shape.scanline(&scanline, (f64i(y) + 0.5) / scale - ty, allocator);
        for (0..w) |x| {
            const idx = row * w + x;
            const distance = out_pixels[idx];
            if ((distance > 0.5) != scanline.filled((f64i(x) + 0.5) / scale - tx, fill_rule))
                out_pixels[idx] = 1.0 - distance;
        }
    }
}

fn msdfSignCorrection(
    allocator: std.mem.Allocator,
    out_pixels: []f64,
    w: u16,
    h: u16,
    shape: Shape,
    scale: f64,
    tx: f64,
    ty: f64,
    fill_rule: Scanline.FillRule,
    channels: u8,
) !void {
    var scanline: Scanline = .{};
    defer scanline.intersections.deinit(allocator);
    var ambiguous = false;
    var match_map: []i32 = try allocator.alloc(i32, w * h);
    defer allocator.free(match_map);
    var match_idx: usize = 0;
    const scaled_w = w * channels;
    for (0..h) |y| {
        const row = h - y - 1;
        try shape.scanline(&scanline, (f64i(y) + 0.5) / scale - ty, allocator);
        for (0..w) |x| {
            const filled = scanline.filled((f64i(x) + 0.5) / scale - tx, fill_rule);
            const idx = row * scaled_w + x * channels;
            const distance = median(out_pixels[idx], out_pixels[idx + 1], out_pixels[idx + 2]);
            if (distance == 0.5) {
                ambiguous = true;
            } else if ((distance > 0.5) != filled) {
                for (0..3) |i| out_pixels[idx + i] = 1.0 - out_pixels[idx + i];
                match_map[match_idx] = -1;
            } else match_map[match_idx] = 1;
            if (channels >= 4 and (out_pixels[idx + 3] > 0.5) != filled)
                out_pixels[idx + 3] = 1.0 - out_pixels[idx + 3];
            match_idx += 1;
        }
    }

    if (ambiguous) {
        match_idx = 0;
        for (0..h) |y| {
            const row = h - y - 1;
            for (0..w) |x| {
                const match = match_map[match_idx];
                if (match == 0) {
                    var neighbor_match: i32 = 0;
                    if (x > 0) neighbor_match += match - 1;
                    if (x < w - 1) neighbor_match += match + 1;
                    if (y > 0) neighbor_match += match - w;
                    if (y < h - 1) neighbor_match += match + w;
                    if (neighbor_match < 0) {
                        const idx = row * scaled_w + x * channels;
                        for (0..3) |i| out_pixels[idx + i] = 1.0 - out_pixels[idx + i];
                    }
                }
                match_idx += 1;
            }
        }
    }
}

fn findDistanceAt(shape: Shape, p: Vec2, px_range: f64) f64 {
    var dummy: f64 = 0;
    var min_dist: SignedDistance = .{};
    for (shape.contours.items) |contour| for (contour.edges.items) |*edge| {
        const dist = edge.signedDistance(p, &dummy);
        if (dist.lt(min_dist)) min_dist = dist;
    };
    return (min_dist.distance + px_range / 2.0) / px_range;
}

fn generateSdf(out_pixels: []f64, size: u16, scale: f64, shape: Shape, px_range: f64, tx: f64, ty: f64) void {
    for (0..size) |y| {
        const row = size - y - 1;
        for (0..size) |x| {
            var dummy: f64 = 0;
            const p: Vec2 = .{
                .x = (f64i(x) + 0.5) / scale - tx,
                .y = (f64i(y) + 0.5) / scale - ty,
            };
            var min_dist: SignedDistance = .{};
            for (shape.contours.items) |contour| for (contour.edges.items) |*edge| {
                const dist = edge.signedDistance(p, &dummy);
                if (dist.lt(min_dist)) min_dist = dist;
            };
            out_pixels[row * size + x] = (min_dist.distance + px_range / 2.0) / px_range;
        }
    }
}

fn generatePsdf(out_pixels: []f64, size: u16, scale: f64, shape: Shape, px_range: f64, tx: f64, ty: f64) void {
    for (0..size) |y| {
        const row = size - y - 1;
        for (0..size) |x| {
            const p: Vec2 = .{
                .x = (f64i(x) + 0.5) / scale - tx,
                .y = (f64i(y) + 0.5) / scale - ty,
            };
            var min_dist: SignedDistance = .{};
            var near_edge: ?*EdgeSegment = null;
            var near_param: f64 = 0;
            for (shape.contours.items) |contour| for (contour.edges.items) |*edge| {
                var param: f64 = 0;
                const dist = edge.signedDistance(p, &param);
                if (dist.lt(min_dist)) {
                    min_dist = dist;
                    near_edge = edge;
                    near_param = param;
                }
            };
            if (near_edge) |edge| edge.distanceToPerpendicularDistance(&min_dist, p, near_param);
            out_pixels[row * size + x] = (min_dist.distance + px_range / 2.0) / px_range;
        }
    }
}

fn generateMsdf(out_pixels: []f64, size: u16, scale: f64, shape: Shape, px_range: f64, tx: f64, ty: f64) void {
    for (0..size) |y| {
        const row = size - y - 1;
        for (0..size) |x| {
            const p: Vec2 = .{
                .x = (f64i(x) + 0.5) / scale - tx,
                .y = (f64i(y) + 0.5) / scale - ty,
            };
            const PsdfData = struct {
                min_dist: SignedDistance = .{},
                near_edge: ?*EdgeSegment = null,
                near_param: f64 = 0,
            };
            var r: PsdfData = .{};
            var g: PsdfData = .{};
            var b: PsdfData = .{};
            for (shape.contours.items) |contour| for (contour.edges.items) |*edge| {
                var param: f64 = 0;
                const dist = edge.signedDistance(p, &param);
                if ((@intFromEnum(edge.color) & @intFromEnum(EdgeColor.red)) != 0 and dist.lt(r.min_dist)) {
                    r.min_dist = dist;
                    r.near_edge = edge;
                    r.near_param = param;
                }
                if ((@intFromEnum(edge.color) & @intFromEnum(EdgeColor.green)) != 0 and dist.lt(g.min_dist)) {
                    g.min_dist = dist;
                    g.near_edge = edge;
                    g.near_param = param;
                }
                if ((@intFromEnum(edge.color) & @intFromEnum(EdgeColor.blue)) != 0 and dist.lt(b.min_dist)) {
                    b.min_dist = dist;
                    b.near_edge = edge;
                    b.near_param = param;
                }
            };
            if (r.near_edge) |edge| edge.distanceToPerpendicularDistance(&r.min_dist, p, r.near_param);
            if (g.near_edge) |edge| edge.distanceToPerpendicularDistance(&g.min_dist, p, g.near_param);
            if (b.near_edge) |edge| edge.distanceToPerpendicularDistance(&b.min_dist, p, b.near_param);

            const channels = 3;
            const scaled_size = size * channels;
            const scaled_x = x * channels;
            out_pixels[row * scaled_size + scaled_x] = (r.min_dist.distance + px_range / 2.0) / px_range;
            out_pixels[row * scaled_size + scaled_x + 1] = (g.min_dist.distance + px_range / 2.0) / px_range;
            out_pixels[row * scaled_size + scaled_x + 2] = (b.min_dist.distance + px_range / 2.0) / px_range;
        }
    }
}

fn generateMtsdf(out_pixels: []f64, size: u16, scale: f64, shape: Shape, px_range: f64, tx: f64, ty: f64) void {
    for (0..size) |y| {
        const row = size - y - 1;
        for (0..size) |x| {
            const p: Vec2 = .{
                .x = (f64i(x) + 0.5) / scale - tx,
                .y = (f64i(y) + 0.5) / scale - ty,
            };
            const PsdfData = struct {
                min_dist: SignedDistance = .{},
                near_edge: ?*EdgeSegment = null,
                near_param: f64 = 0,
            };
            var r: PsdfData = .{};
            var g: PsdfData = .{};
            var b: PsdfData = .{};
            var min_dist: SignedDistance = .{};
            for (shape.contours.items) |contour| for (contour.edges.items) |*edge| {
                var param: f64 = 0;
                const dist = edge.signedDistance(p, &param);
                if (dist.lt(min_dist)) min_dist = dist;
                if ((@intFromEnum(edge.color) & @intFromEnum(EdgeColor.red)) != 0 and dist.lt(r.min_dist)) {
                    r.min_dist = dist;
                    r.near_edge = edge;
                    r.near_param = param;
                }
                if ((@intFromEnum(edge.color) & @intFromEnum(EdgeColor.green)) != 0 and dist.lt(g.min_dist)) {
                    g.min_dist = dist;
                    g.near_edge = edge;
                    g.near_param = param;
                }
                if ((@intFromEnum(edge.color) & @intFromEnum(EdgeColor.blue)) != 0 and dist.lt(b.min_dist)) {
                    b.min_dist = dist;
                    b.near_edge = edge;
                    b.near_param = param;
                }
            };
            if (r.near_edge) |edge| edge.distanceToPerpendicularDistance(&r.min_dist, p, r.near_param);
            if (g.near_edge) |edge| edge.distanceToPerpendicularDistance(&g.min_dist, p, g.near_param);
            if (b.near_edge) |edge| edge.distanceToPerpendicularDistance(&b.min_dist, p, b.near_param);

            const channels = 4;
            const scaled_size = size * channels;
            const scaled_x = x * channels;
            out_pixels[row * scaled_size + scaled_x] = (r.min_dist.distance + px_range / 2.0) / px_range;
            out_pixels[row * scaled_size + scaled_x + 1] = (g.min_dist.distance + px_range / 2.0) / px_range;
            out_pixels[row * scaled_size + scaled_x + 2] = (b.min_dist.distance + px_range / 2.0) / px_range;
            out_pixels[row * scaled_size + scaled_x + 3] = (min_dist.distance + px_range / 2.0) / px_range;
        }
    }
}

fn ftMoveTo(to: [*c]const ft.Vector, ud: ?*anyopaque) callconv(.C) i32 {
    var context: *FreetypeContext = @alignCast(@ptrCast(ud));
    if (!(context.contour != null and context.contour.?.edges.items.len == 0)) {
        context.contour = context.shape.contours.addOne(context.allocator) catch return ft.c.FT_Err_Out_Of_Memory;
        context.contour.?.* = .{};
    }
    context.pos = .{ .x = f64i(to.*.x) * context.scale, .y = f64i(to.*.y) * context.scale };
    return 0;
}

fn ftLineTo(to: [*c]const ft.Vector, ud: ?*anyopaque) callconv(.C) i32 {
    var context: *FreetypeContext = @alignCast(@ptrCast(ud));
    const endpoint: Vec2 = .{ .x = f64i(to.*.x) * context.scale, .y = f64i(to.*.y) * context.scale };
    if (!endpoint.eql(context.pos)) {
        context.contour.?.edges.append(
            context.allocator,
            .create(context.pos, endpoint, null, null, .white),
        ) catch return ft.c.FT_Err_Out_Of_Memory;
        context.pos = endpoint;
    }
    return 0;
}

fn ftConicTo(control: [*c]const ft.Vector, to: [*c]const ft.Vector, ud: ?*anyopaque) callconv(.C) i32 {
    var context: *FreetypeContext = @alignCast(@ptrCast(ud));
    const endpoint: Vec2 = .{ .x = f64i(to.*.x) * context.scale, .y = f64i(to.*.y) * context.scale };
    if (!endpoint.eql(context.pos)) {
        context.contour.?.edges.append(context.allocator, .create(
            context.pos,
            .{ .x = f64i(control.*.x) * context.scale, .y = f64i(control.*.y) * context.scale },
            endpoint,
            null,
            .white,
        )) catch return ft.c.FT_Err_Out_Of_Memory;
        context.pos = endpoint;
    }
    return 0;
}

fn ftCubicTo(control1: [*c]const ft.Vector, control2: [*c]const ft.Vector, to: [*c]const ft.Vector, ud: ?*anyopaque) callconv(.C) i32 {
    var context: *FreetypeContext = @alignCast(@ptrCast(ud));
    const endpoint: Vec2 = .{ .x = f64i(to.*.x) * context.scale, .y = f64i(to.*.y) * context.scale };
    const scaled_c1: Vec2 = .{ .x = f64i(control1.*.x) * context.scale, .y = f64i(control1.*.y) * context.scale };
    const scaled_c2: Vec2 = .{ .x = f64i(control2.*.x) * context.scale, .y = f64i(control2.*.y) * context.scale };
    if (!endpoint.eql(context.pos) or scaled_c1.sub(endpoint).cross(scaled_c2.sub(endpoint)) != 0.0) {
        context.contour.?.edges.append(
            context.allocator,
            .create(context.pos, scaled_c1, scaled_c2, endpoint, .white),
        ) catch return ft.c.FT_Err_Out_Of_Memory;
        context.pos = endpoint;
    }
    return 0;
}

fn f64i(int: anytype) f32 {
    return @floatFromInt(int);
}

fn u8f(float: anytype) u8 {
    const int_form: u8 = @intFromFloat(255.5 - 255.0 * clampNorm(float));
    return ~int_form;
}
