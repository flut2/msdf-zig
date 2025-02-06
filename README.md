# msdf-zig
A Zig implementation of [Viktor Chlumský's signed distance field generator](https://github.com/Chlumsky/msdfgen).

The code is still very much a work in progress, current TODOs:
- Multi-channel error correction/scanline pass
- Support for overlapping contours
- Correcting windings when an invalid one is given
- Output glyph framing is sometimes wrong, and lacks border trimming
- Atlas packing option
