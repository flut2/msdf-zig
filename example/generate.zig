const std = @import("std");

const Generator = @import("msdf-zig");
const zstbi = @import("zstbi");

const font_data = @embedFile("OpenSans-Bold.ttf");

pub fn main() !void {
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer _ = gpa.deinit();
    const allocator = gpa.allocator();

    zstbi.init(allocator);
    defer zstbi.deinit();

    // In a real world scenario you'd open the font from disk rather than embedding it
    // var file = try std.fs.cwd().openFile("OpenSans-Bold.ttf", .{});
    // defer file.close();

    // const font_memory = try file.readToEndAlloc(allocator, std.math.maxInt(u32));
    // defer allocator.free(font_memory);

    // var gen: Generator = try .create(font_memory);
    // defer gen.destroy();

    var gen: Generator = try .create(font_data);
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

    const px_size = 64;
    inline for (.{ 'A', 'S' }) |codepoint| {
        const data = try gen.generate(allocator, codepoint, .{ .sdf_type = .msdf, .px_size = px_size, .px_range = 8 });
        defer data.deinit(allocator);
        std.log.info("Glyph Advance for \'{c}\': {d:.2}", .{codepoint, data.advance});

        var image: zstbi.Image = try .createEmpty(px_size, px_size, 3, .{});
        defer image.deinit();
        @memcpy(image.data, data.pixels);

        const path = std.fmt.comptimePrint("{c}_sdf.png", .{codepoint});
        try image.writeToFile(path, .png);
    }
}
