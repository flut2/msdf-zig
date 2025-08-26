const std = @import("std");

const Generator = @import("msdf-zig");
const zstbi = @import("zstbi");

pub fn main() !void {
    var dbg_alloc: std.heap.DebugAllocator(.{ .stack_trace_frames = 10 }) = .init;
    defer _ = dbg_alloc.deinit();

    const allocator = dbg_alloc.allocator();

    zstbi.init(allocator);
    defer zstbi.deinit();

    var file = try std.fs.cwd().openFile("assets/OpenSans-Bold.ttf", .{});
    defer file.close();

    var read_buf: [4096]u8 = undefined;
    var reader = file.reader(&read_buf);

    const font_memory = try reader.interface.allocRemaining(allocator, .unlimited);
    defer allocator.free(font_memory);

    var gen: Generator = try .create(font_memory);
    defer gen.destroy();

    const metrics = try gen.fontMetrics();
    std.log.info(
        \\Font Metrics:
        \\Ascender: {d:.2}
        \\Descender: {d:.2}
        \\Underline Y: {d:.2}
        \\Underline Thickness: {d:.2}
        \\Line Height: {d:.2}
    , .{
        metrics.ascender,
        metrics.descender,
        metrics.underline_y,
        metrics.underline_thickness,
        metrics.line_height,
    });

    const gen_opts: Generator.GenerationOptions = .{ .sdf_type = .mtsdf, .px_size = 64, .px_range = 8 };

    inline for (.{ 'A', 'B', 'C' }) |codepoint| {
        var timer: std.time.Timer = try .start();
        const data = try gen.generateSingle(allocator, codepoint, gen_opts);
        defer data.deinit(allocator);
        std.log.info("SDF for codepoint {u} generated in: {d}us", .{ codepoint, @divFloor(timer.read(), std.time.ns_per_us) });

        std.log.info(
            \\Single Glyph Data for "{u}":
            \\Advance: {d:.2}
            \\X Bearing: {d:.2}
            \\Y Bearing: {d:.2}
            \\Width: {d:.2}
            \\Height: {d:.2}
            \\
        , .{
            codepoint,
            data.glyph_data.advance,
            data.glyph_data.bearing_x,
            data.glyph_data.bearing_y,
            data.glyph_data.width,
            data.glyph_data.height,
        });

        var image: zstbi.Image = try .createEmpty(data.glyph_data.width, data.glyph_data.height, gen_opts.sdf_type.numChannels(), .{});
        defer image.deinit();
        @memcpy(image.data, data.pixels);

        const path = std.fmt.comptimePrint("{u}_sdf.png", .{codepoint});
        try image.writeToFile(path, .png);
    }

    const atlas_w = 512;
    const atlas_h = 512;
    var timer: std.time.Timer = try .start();
    const data = try gen.generateAtlas(
        allocator,
        &.{ 'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L' },
        atlas_w,
        atlas_h,
        2,
        gen_opts,
    );
    defer data.deinit(allocator);
    std.log.info("SDF for atlas generated in: {d}us", .{@divFloor(timer.read(), std.time.ns_per_us)});
    for (data.glyphs) |atlas_glyph| std.log.info(
        \\Atlas Glyph Data for "{u}":
        \\Advance: {d:.2}
        \\X Bearing: {d:.2}
        \\Y Bearing: {d:.2}
        \\Width: {d:.2}
        \\Height: {d:.2}
        \\Texture U: {d:.2}
        \\Texture V: {d:.2}
        \\Texture W: {d:.2}
        \\Texture H: {d:.2}
        \\
    , .{
        atlas_glyph.codepoint,
        atlas_glyph.glyph_data.advance,
        atlas_glyph.glyph_data.bearing_x,
        atlas_glyph.glyph_data.bearing_y,
        atlas_glyph.glyph_data.width,
        atlas_glyph.glyph_data.height,
        atlas_glyph.tex_u,
        atlas_glyph.tex_v,
        atlas_glyph.tex_w,
        atlas_glyph.tex_h,
    });

    var image: zstbi.Image = try .createEmpty(atlas_w, atlas_h, gen_opts.sdf_type.numChannels(), .{});
    defer image.deinit();
    @memcpy(image.data, data.pixels);

    try image.writeToFile("atlas_sdf.png", .png);
}
