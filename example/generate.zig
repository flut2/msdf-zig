const std = @import("std");

const Generator = @import("msdf-zig");
const zstbi = @import("zstbi");

fn printableAscii() []const u21 {
    var ret: []const u21 = &.{};
    for (32..127) |i| ret = ret ++ [_]u21{i};
    return ret;
}

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

        var image: zstbi.Image = try .createEmpty(data.glyph_data.width, data.glyph_data.height, gen_opts.sdf_type.numChannels(), .{});
        defer image.deinit();
        @memcpy(image.data, data.pixels);

        const path = std.fmt.comptimePrint("{u}_sdf.png", .{codepoint});
        try image.writeToFile(path, .png);
    }

    const atlas_w = 1024;
    const atlas_h = 512;
    var timer: std.time.Timer = try .start();
    const data = try gen.generateAtlas(
        allocator,
        comptime printableAscii(),
        atlas_w,
        atlas_h,
        2,
        true,
        gen_opts,
    );
    defer data.deinit(allocator);
    std.log.info("SDF for atlas generated in: {d}us", .{@divFloor(timer.read(), std.time.ns_per_us)});

    var image: zstbi.Image = try .createEmpty(atlas_w, atlas_h, gen_opts.sdf_type.numChannels(), .{});
    defer image.deinit();
    @memcpy(image.data, data.pixels);

    try image.writeToFile("atlas_sdf.png", .png);
}
