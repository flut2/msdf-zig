# msdf-zig
A Zig implementation of [Viktor Chlumský's signed distance field generator](https://github.com/Chlumsky/msdfgen).

The code is still very much a work in progress, current TODOs:
- Multi-channel error correction
- Output glyph framing lacks border trimming
- Atlas packing option
- Support for overlapping contours
- SVG support
- Variable fonts
- Option for Skia geometry preprocessing