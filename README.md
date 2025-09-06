# msdf-zig
A Zig implementation of [Viktor Chlumský's signed distance field generator](https://github.com/Chlumsky/msdfgen).

## Usage
```zig
const Generator = @import("msdf-zig");
const font_data = @embedFile("OpenSans-Bold.ttf");

var gen: Generator = try .create(font_data);
defer gen.destroy();

inline for (.{ 'A', 'B', 'C' }) |codepoint| {
    const data = try gen.generateSingle(allocator, codepoint, .{ .sdf_type = .mtsdf, .px_size = 64, .px_range = 8 });
    defer data.deinit(allocator);
    
    var image: zstbi.Image = try .createEmpty(data.glyph_data.width, data.glyph_data.height, Generator.SdfType.numChannels(.mtsdf), .{});
    defer image.deinit();
    @memcpy(image.data, data.pixels);

    const path = std.fmt.comptimePrint("{u}_sdf.png", .{codepoint});
    try image.writeToFile(path, .png);
}
```

A more in-depth example can be found in `example/generate.zig`.

## Disclaimer
This library might provide an option for it later, but you currently need to preprocess your fonts manually to resolve overlapping contours (if the font has them).