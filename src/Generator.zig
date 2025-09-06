const std = @import("std");

const ft = @import("mach-freetype");
const pack = @import("turbopack");

const coloring = @import("coloring.zig");
const Contour = @import("Contour.zig");
const edge_color = @import("edge_color.zig");
const EdgeColor = edge_color.EdgeColor;
const EdgeSegment = @import("EdgeSegment.zig");
const ErrorCorrection = @import("ErrorCorrection.zig");
const math = @import("math.zig");
const Scanline = @import("Scanline.zig");
const Shape = @import("Shape.zig");
const SignedDistance = @import("SignedDistance.zig");

const Vec2 = @Vector(2, f64);

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
    bearing_x: f64,
    bearing_y: f64,
    width: u16,
    height: u16,
};

pub const KerningPair = struct {
    codepoint_1: u21,
    codepoint_2: u21,
    x: f64,
    y: f64,
};

pub const AtlasGlyphData = struct {
    glyph_data: GlyphData,
    codepoint: u21,
    tex_u: f64,
    tex_v: f64,
    tex_w: f64,
    tex_h: f64,
};

pub const SingleGlyphData = struct {
    glyph_data: GlyphData,
    pixels: []const u8,

    pub fn deinit(self: SingleGlyphData, allocator: std.mem.Allocator) void {
        allocator.free(self.pixels);
    }
};

pub const AtlasData = struct {
    glyphs: []const AtlasGlyphData,
    kernings: []const KerningPair,
    pixels: []const u8,

    pub fn deinit(self: AtlasData, allocator: std.mem.Allocator) void {
        allocator.free(self.glyphs);
        allocator.free(self.kernings);
        allocator.free(self.pixels);
    }
};

pub const SdfType = enum {
    sdf,
    psdf,
    msdf,
    mtsdf,

    pub fn numChannels(self: SdfType) u8 {
        return switch (self) {
            .sdf, .psdf => 1,
            .msdf => 3,
            .mtsdf => 4,
        };
    }
};

pub const OrientationType = enum {
    guess,
    keep,
    reverse,
};

pub const VarFontArgument = struct {
    name: []const u8,
    value: f64,
};

pub const GenerationOptions = struct {
    sdf_type: SdfType,
    px_size: u16,
    px_range: u16,
    coloring_rng_seed: u64 = 0,
    corner_angle_threshold: f64 = 3.0,
    orientation: OrientationType = .guess,
    geometry_preprocess: bool = false,
    /// Requires geometry preprocessing to be disabled
    scanline_fill_rule: ?Scanline.FillRule = null,
    /// Only MSDFs and MTSDFs can be error corrected
    error_correction_opts: ?ErrorCorrection.Options = .{},
    var_font_args: []const VarFontArgument = &.{},
};

const FreetypeContext = struct {
    allocator: std.mem.Allocator,
    scale: f64,
    shape: *Shape,
    pos: Vec2 = @splat(0.0),
    contour: ?*Contour = null,
};

library: ft.Library = undefined,
face: ft.Face = undefined,

/// `font_memory` is the raw font file data
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

fn handleVarFont(self: *Generator, allocator: std.mem.Allocator, var_args: []const VarFontArgument, face_flags: ft.FaceFlags) !void {
    if (var_args.len == 0) return;

    if (face_flags.multiple_masters) {
        const vf = try self.face.createVarFontInfo();
        if (vf) |var_font| if (var_font.num_axis > 0) {
            var coords = try allocator.alloc(ft.c.FT_Fixed, var_font.num_axis);
            defer allocator.free(coords);
            try self.face.getVarDesignCoords(coords);

            for (var_args) |args|
                for (var_font.axis[0..var_font.num_axis], 0..) |axis, i|
                    if (std.mem.eql(u8, std.mem.span(axis.name), args.name)) {
                        coords[i] = @intFromFloat(std.math.maxInt(u16) * args.value);
                    };
            try self.face.setVarDesignCoords(coords);
        };
        try self.library.destroyVarFontInfo(vf);
    } else std.log.warn("Var font args supplied, but the face only has a single master", .{});
}

/// The result is under the caller's ownership (call `deinit()` or deallocate fields manually)
pub fn generateSingle(
    self: *Generator,
    allocator: std.mem.Allocator,
    codepoint: u21,
    gen_opts: GenerationOptions,
) !SingleGlyphData {
    edge_color.rng.seed(gen_opts.coloring_rng_seed);

    try self.handleVarFont(allocator, gen_opts.var_font_args, self.face.faceFlags());

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
    if (gen_opts.geometry_preprocess) try shape.orientContours(allocator);
    try shape.normalize(allocator);

    const f_px_size = f64i(gen_opts.px_size);
    const px_range = f64i(gen_opts.px_range) / f_px_size;

    var bounds = shape.getBounds(0, 0, 0);
    if (bounds.left >= bounds.right or bounds.bottom >= bounds.top)
        bounds = .{ .left = 0, .bottom = 0, .right = 1, .top = 1 };

    const translate_x = -bounds.left + px_range / 2.0;
    const translate_y = -bounds.bottom + px_range / 2.0;
    const w: u16 = @intFromFloat((bounds.right - bounds.left + px_range) * f_px_size);
    const h: u16 = @intFromFloat((bounds.top - bounds.bottom + px_range) * f_px_size);

    const oob_point: Vec2 = if (gen_opts.orientation == .guess)
        .{ bounds.left - (bounds.right - bounds.left) - 1, bounds.bottom - (bounds.top - bounds.bottom) - 1 }
    else
        undefined;

    const metrics = self.face.glyph().metrics();
    return .{
        .glyph_data = .{
            .advance = scale * f64i(self.face.glyph().advance().x),
            .bearing_x = scale * f64i(metrics.horiBearingX),
            .bearing_y = scale * f64i(metrics.horiBearingY),
            .width = w,
            .height = h,
        },
        .pixels = try getSdfPixels(allocator, gen_opts, w, h, &shape, translate_x, translate_y, oob_point),
    };
}

/// The result is under the caller's ownership (call `deinit()` or deallocate fields manually)
pub fn generateAtlas(
    self: *Generator,
    allocator: std.mem.Allocator,
    codepoints: []const u21,
    w: u16,
    h: u16,
    padding: u8,
    use_kerning: bool,
    gen_opts: GenerationOptions,
) !AtlasData {
    edge_color.rng.seed(gen_opts.coloring_rng_seed);

    const face_flags = self.face.faceFlags();
    try self.handleVarFont(allocator, gen_opts.var_font_args, face_flags);

    const channels = gen_opts.sdf_type.numChannels();
    const glyphs = try allocator.alloc(AtlasGlyphData, codepoints.len);
    errdefer allocator.free(glyphs);
    const pixels = try allocator.alloc(u8, @as(u32, w) * @as(u32, h) * channels);
    errdefer allocator.free(pixels);
    @memset(pixels, 0);

    var pack_ctx: pack.Context = try .create(allocator, w, h, .{});
    defer pack_ctx.deinit();

    const char_indices = try allocator.alloc(u32, codepoints.len);
    defer allocator.free(char_indices);
    for (codepoints, char_indices) |c, *i|
        i.* = self.face.getCharIndex(c) orelse return error.InvalidCodepoint;

    var kernings: std.ArrayList(KerningPair) = .empty;
    errdefer kernings.deinit(allocator);

    if (use_kerning) {
        if (face_flags.kerning) {
            for (char_indices, codepoints) |i, ci| for (char_indices, codepoints) |j, cj|
                if (i != j) {
                    const kern = try self.face.getKerning(i, j, .default);
                    if (kern.x != 0 or kern.y != 0)
                        try kernings.append(allocator, .{
                            .codepoint_1 = ci,
                            .codepoint_2 = cj,
                            .x = f64i(kern.x),
                            .y = f64i(kern.y),
                        });
                };
        } else std.log.warn("Kerning requested, but none were found in the font file. Note: FreeType doesn't support GPOS kerning", .{});
    }

    for (codepoints, 0..) |codepoint, i| {
        const scale = 1.0 / f64i(self.face.unitsPerEM());
        const idx = char_indices[i];
        try self.face.loadGlyph(idx, .{ .no_scale = true, .no_bitmap = true });

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
        if (gen_opts.geometry_preprocess) try shape.orientContours(allocator);
        try shape.normalize(allocator);

        const f_px_size = f64i(gen_opts.px_size);
        const px_range = f64i(gen_opts.px_range) / f_px_size;

        var bounds = shape.getBounds(0, 0, 0);
        if (bounds.left >= bounds.right or bounds.bottom >= bounds.top)
            bounds = .{ .left = 0, .bottom = 0, .right = 1, .top = 1 };

        const translate_x = -bounds.left + px_range / 2.0;
        const translate_y = -bounds.bottom + px_range / 2.0;
        const glyph_w: u16 = @intFromFloat((bounds.right - bounds.left + px_range) * f_px_size);
        const glyph_h: u16 = @intFromFloat((bounds.top - bounds.bottom + px_range) * f_px_size);

        const oob_point: Vec2 = if (gen_opts.orientation == .guess)
            .{ bounds.left - (bounds.right - bounds.left) - 1, bounds.bottom - (bounds.top - bounds.bottom) - 1 }
        else
            undefined;

        if (glyph_w <= 0 or glyph_h <= 0) {
            glyphs[i] = .{
                .glyph_data = .{
                    .advance = scale * f64i(self.face.glyph().advance().x),
                    .bearing_x = 0,
                    .bearing_y = 0,
                    .width = glyph_w,
                    .height = glyph_h,
                },
                .codepoint = codepoint,
                .tex_u = 1.0,
                .tex_v = 1.0,
                .tex_w = 0.0,
                .tex_h = 0.0,
            };
            continue;
        }

        const sdf_pixels = try getSdfPixels(allocator, gen_opts, glyph_w, glyph_h, &shape, translate_x, translate_y, oob_point);
        defer allocator.free(sdf_pixels);

        var rect: [1]pack.Rect = .{.{
            .w = glyph_w + padding * 2,
            .h = glyph_h + padding * 2,
        }};
        try pack.pack(pack.Rect, &pack_ctx, &rect, .{});

        const cur_atlas_x: u32 = @intCast(rect[0].x + padding);
        const cur_atlas_y: u32 = @intCast(rect[0].y + padding);
        for (0..glyph_h) |j| {
            const atlas_idx = ((cur_atlas_y + j) * w + cur_atlas_x) * channels;
            const src_idx = (j * glyph_w) * channels;
            @memcpy(pixels[atlas_idx .. atlas_idx + glyph_w * channels], sdf_pixels[src_idx .. src_idx + glyph_w * channels]);
        }

        const metrics = self.face.glyph().metrics();
        glyphs[i] = .{
            .glyph_data = .{
                .advance = scale * f64i(self.face.glyph().advance().x),
                .bearing_x = scale * f64i(metrics.horiBearingX),
                .bearing_y = scale * f64i(metrics.horiBearingY),
                .width = @intCast(rect[0].w),
                .height = @intCast(rect[0].h),
            },
            .codepoint = codepoint,
            .tex_u = f64i(rect[0].x) / f64i(w),
            .tex_v = f64i(rect[0].y) / f64i(h),
            .tex_w = f64i(rect[0].w) / f64i(w),
            .tex_h = f64i(rect[0].h) / f64i(h),
        };
    }

    return .{
        .glyphs = glyphs,
        .pixels = pixels,
        .kernings = if (use_kerning and kernings.items.len > 0)
            try kernings.toOwnedSlice(allocator)
        else
            &.{},
    };
}

fn getSdfPixels(
    allocator: std.mem.Allocator,
    opts: GenerationOptions,
    w: u16,
    h: u16,
    shape: *Shape,
    translate_x: f64,
    translate_y: f64,
    oob_point: Vec2,
) ![]const u8 {
    const f_px_size = f64i(opts.px_size);
    const px_range = f64i(opts.px_range) / f_px_size;

    var error_correction: ?ErrorCorrection =
        if (opts.error_correction_opts) |ec_opts| b: {
            break :b if (opts.sdf_type == .msdf or opts.sdf_type == .mtsdf)
                try .create(allocator, shape, w, h, ec_opts, opts.scanline_fill_rule != null)
            else
                null;
        } else null;
    defer if (error_correction) |*ec| ec.destroy(allocator);

    const channels = opts.sdf_type.numChannels();

    const float_pixels = try allocator.alloc(f64, w * h * channels);
    defer allocator.free(float_pixels);
    const pixels = try allocator.alloc(u8, w * h * channels);
    const invert_pixels = opts.orientation == .reverse or
        (opts.orientation == .guess and findDistanceAt(shape.*, oob_point, px_range) > 0);
    switch (opts.sdf_type) {
        .sdf => generateSdf(float_pixels, w, h, f_px_size, shape.*, px_range, translate_x, translate_y, invert_pixels),
        .psdf => generatePsdf(float_pixels, w, h, f_px_size, shape.*, px_range, translate_x, translate_y, invert_pixels),
        .msdf => {
            try coloring.colorShape(allocator, shape, opts.corner_angle_threshold);
            generateMsdf(float_pixels, w, h, f_px_size, shape.*, px_range, translate_x, translate_y, invert_pixels);
        },
        .mtsdf => {
            try coloring.colorShape(allocator, shape, opts.corner_angle_threshold);
            generateMtsdf(float_pixels, w, h, f_px_size, shape.*, px_range, translate_x, translate_y, invert_pixels);
        },
    }

    if (!opts.geometry_preprocess) if (opts.scanline_fill_rule) |fill_rule|
        switch (opts.sdf_type) {
            .sdf, .psdf => try sdfSignCorrection(
                allocator,
                float_pixels,
                w,
                h,
                f_px_size,
                shape.*,
                translate_x,
                translate_y,
                fill_rule,
            ),
            .msdf, .mtsdf => try msdfSignCorrection(
                allocator,
                float_pixels,
                w,
                h,
                f_px_size,
                shape.*,
                translate_x,
                translate_y,
                fill_rule,
                channels,
            ),
        };

    if (error_correction) |*ec| ec.correct(shape, f_px_size, px_range, translate_x, translate_y, float_pixels, w, h, channels);
    convertPixels(float_pixels, pixels, w, h, channels);
    return pixels;
}

fn convertPixels(in_pixels: []f64, out_pixels: []u8, w: u16, h: u16, channels: u8) void {
    for (0..h) |y| for (0..w) |x| {
        const idx = y * w * channels + x * channels;
        for (0..channels) |i|
            out_pixels[idx + i] = @intFromFloat(255.0 * std.math.clamp(in_pixels[idx + i], 0.0, 1.0));
    };
}

fn sdfSignCorrection(
    allocator: std.mem.Allocator,
    out_pixels: []f64,
    w: u16,
    h: u16,
    scale: f64,
    shape: Shape,
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
    scale: f64,
    shape: Shape,
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
            const distance = math.median(out_pixels[idx], out_pixels[idx + 1], out_pixels[idx + 2]);
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
                        for (out_pixels[row * scaled_w + x * channels ..][0..3]) |*px|
                            px.* = 1.0 - px.*;
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
        if (dist.lessThan(min_dist)) min_dist = dist;
    };
    return (min_dist.distance + px_range / 2.0) / px_range;
}

fn generateSdf(out_pixels: []f64, w: u16, h: u16, scale: f64, shape: Shape, px_range: f64, tx: f64, ty: f64, invert_pixels: bool) void {
    for (0..h) |y| {
        const row = h - y - 1;
        for (0..w) |x| {
            var dummy: f64 = 0;
            const p: Vec2 = .{
                (f64i(x) + 0.5) / scale - tx,
                (f64i(y) + 0.5) / scale - ty,
            };
            var min_dist: SignedDistance = .{};
            for (shape.contours.items) |contour| for (contour.edges.items) |*edge| {
                const dist = edge.signedDistance(p, &dummy);
                if (dist.lessThan(min_dist)) min_dist = dist;
            };
            const out = &out_pixels[row * w + x];
            out.* = (min_dist.distance + px_range / 2.0) / px_range;
            if (invert_pixels) out.* = 1.0 - out.*;
        }
    }
}

fn generatePsdf(out_pixels: []f64, w: u16, h: u16, scale: f64, shape: Shape, px_range: f64, tx: f64, ty: f64, invert_pixels: bool) void {
    for (0..h) |y| {
        const row = h - y - 1;
        for (0..w) |x| {
            const p: Vec2 = .{
                (f64i(x) + 0.5) / scale - tx,
                (f64i(y) + 0.5) / scale - ty,
            };
            var min_dist: SignedDistance = .{};
            var near_edge: ?*EdgeSegment = null;
            var near_param: f64 = 0;
            for (shape.contours.items) |contour| for (contour.edges.items) |*edge| {
                var param: f64 = 0;
                const dist = edge.signedDistance(p, &param);
                if (dist.lessThan(min_dist)) {
                    min_dist = dist;
                    near_edge = edge;
                    near_param = param;
                }
            };
            if (near_edge) |edge| edge.distanceToPerpendicularDistance(&min_dist, p, near_param);
            const out = &out_pixels[row * w + x];
            out.* = (min_dist.distance + px_range / 2.0) / px_range;
            if (invert_pixels) out.* = 1.0 - out.*;
        }
    }
}

fn generateMsdf(out_pixels: []f64, w: u16, h: u16, scale: f64, shape: Shape, px_range: f64, tx: f64, ty: f64, invert_pixels: bool) void {
    for (0..h) |y| {
        const row = h - y - 1;
        for (0..w) |x| {
            const p: Vec2 = .{
                (f64i(x) + 0.5) / scale - tx,
                (f64i(y) + 0.5) / scale - ty,
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
                if ((@intFromEnum(edge.color) & @intFromEnum(EdgeColor.red)) != 0 and dist.lessThan(r.min_dist)) {
                    r.min_dist = dist;
                    r.near_edge = edge;
                    r.near_param = param;
                }
                if ((@intFromEnum(edge.color) & @intFromEnum(EdgeColor.green)) != 0 and dist.lessThan(g.min_dist)) {
                    g.min_dist = dist;
                    g.near_edge = edge;
                    g.near_param = param;
                }
                if ((@intFromEnum(edge.color) & @intFromEnum(EdgeColor.blue)) != 0 and dist.lessThan(b.min_dist)) {
                    b.min_dist = dist;
                    b.near_edge = edge;
                    b.near_param = param;
                }
            };
            if (r.near_edge) |edge| edge.distanceToPerpendicularDistance(&r.min_dist, p, r.near_param);
            if (g.near_edge) |edge| edge.distanceToPerpendicularDistance(&g.min_dist, p, g.near_param);
            if (b.near_edge) |edge| edge.distanceToPerpendicularDistance(&b.min_dist, p, b.near_param);

            const channels = 3;
            const sc_w = w * channels;
            const sc_x = x * channels;
            const out = out_pixels[row * sc_w + sc_x ..];
            out[0] = (r.min_dist.distance + px_range / 2.0) / px_range;
            out[1] = (g.min_dist.distance + px_range / 2.0) / px_range;
            out[2] = (b.min_dist.distance + px_range / 2.0) / px_range;
            if (invert_pixels) {
                for (out[0..channels]) |*v| v.* = 1.0 - v.*;
            }
        }
    }
}

fn generateMtsdf(out_pixels: []f64, w: u16, h: u16, scale: f64, shape: Shape, px_range: f64, tx: f64, ty: f64, invert_pixels: bool) void {
    for (0..h) |y| {
        const row = h - y - 1;
        for (0..w) |x| {
            const p: Vec2 = .{
                (f64i(x) + 0.5) / scale - tx,
                (f64i(y) + 0.5) / scale - ty,
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
                if (dist.lessThan(min_dist)) min_dist = dist;
                if ((@intFromEnum(edge.color) & @intFromEnum(EdgeColor.red)) != 0 and dist.lessThan(r.min_dist)) {
                    r.min_dist = dist;
                    r.near_edge = edge;
                    r.near_param = param;
                }
                if ((@intFromEnum(edge.color) & @intFromEnum(EdgeColor.green)) != 0 and dist.lessThan(g.min_dist)) {
                    g.min_dist = dist;
                    g.near_edge = edge;
                    g.near_param = param;
                }
                if ((@intFromEnum(edge.color) & @intFromEnum(EdgeColor.blue)) != 0 and dist.lessThan(b.min_dist)) {
                    b.min_dist = dist;
                    b.near_edge = edge;
                    b.near_param = param;
                }
            };
            if (r.near_edge) |edge| edge.distanceToPerpendicularDistance(&r.min_dist, p, r.near_param);
            if (g.near_edge) |edge| edge.distanceToPerpendicularDistance(&g.min_dist, p, g.near_param);
            if (b.near_edge) |edge| edge.distanceToPerpendicularDistance(&b.min_dist, p, b.near_param);

            const channels = 4;
            const sc_w = w * channels;
            const sc_x = x * channels;
            const out = out_pixels[row * sc_w + sc_x ..];
            out[0] = (r.min_dist.distance + px_range / 2.0) / px_range;
            out[1] = (g.min_dist.distance + px_range / 2.0) / px_range;
            out[2] = (b.min_dist.distance + px_range / 2.0) / px_range;
            out[3] = (min_dist.distance + px_range / 2.0) / px_range;
            if (invert_pixels) {
                for (out[0..channels]) |*v| v.* = 1.0 - v.*;
            }
        }
    }
}

fn ftMoveTo(to: [*c]const ft.Vector, ud: ?*anyopaque) callconv(.c) i32 {
    var context: *FreetypeContext = @ptrCast(@alignCast(ud));
    if (!(context.contour != null and context.contour.?.edges.items.len == 0)) {
        context.contour = context.shape.contours.addOne(context.allocator) catch return ft.c.FT_Err_Out_Of_Memory;
        context.contour.?.* = .{};
    }
    context.pos = .{ f64i(to.*.x) * context.scale, f64i(to.*.y) * context.scale };
    return 0;
}

fn ftLineTo(to: [*c]const ft.Vector, ud: ?*anyopaque) callconv(.c) i32 {
    var context: *FreetypeContext = @ptrCast(@alignCast(ud));
    const endpoint: Vec2 = .{ f64i(to.*.x) * context.scale, f64i(to.*.y) * context.scale };
    if (!std.meta.eql(endpoint, context.pos)) {
        context.contour.?.edges.append(
            context.allocator,
            .create(context.pos, endpoint, null, null, .white),
        ) catch return ft.c.FT_Err_Out_Of_Memory;
        context.pos = endpoint;
    }
    return 0;
}

fn ftConicTo(control: [*c]const ft.Vector, to: [*c]const ft.Vector, ud: ?*anyopaque) callconv(.c) i32 {
    var context: *FreetypeContext = @ptrCast(@alignCast(ud));
    const endpoint: Vec2 = .{ f64i(to.*.x) * context.scale, f64i(to.*.y) * context.scale };
    if (!std.meta.eql(endpoint, context.pos)) {
        context.contour.?.edges.append(context.allocator, .create(
            context.pos,
            .{ f64i(control.*.x) * context.scale, f64i(control.*.y) * context.scale },
            endpoint,
            null,
            .white,
        )) catch return ft.c.FT_Err_Out_Of_Memory;
        context.pos = endpoint;
    }
    return 0;
}

fn ftCubicTo(control1: [*c]const ft.Vector, control2: [*c]const ft.Vector, to: [*c]const ft.Vector, ud: ?*anyopaque) callconv(.c) i32 {
    var context: *FreetypeContext = @ptrCast(@alignCast(ud));
    const endpoint: Vec2 = .{ f64i(to.*.x) * context.scale, f64i(to.*.y) * context.scale };
    const scaled_c1: Vec2 = .{ f64i(control1.*.x) * context.scale, f64i(control1.*.y) * context.scale };
    const scaled_c2: Vec2 = .{ f64i(control2.*.x) * context.scale, f64i(control2.*.y) * context.scale };
    if (!std.meta.eql(endpoint, context.pos) or math.cross(scaled_c1 - endpoint, scaled_c2 - endpoint) != 0.0) {
        context.contour.?.edges.append(
            context.allocator,
            .create(context.pos, scaled_c1, scaled_c2, endpoint, .white),
        ) catch return ft.c.FT_Err_Out_Of_Memory;
        context.pos = endpoint;
    }
    return 0;
}

pub fn f64i(int: anytype) f32 {
    return @floatFromInt(int);
}
