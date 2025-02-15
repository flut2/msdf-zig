# msdf-zig
A Zig implementation of [Viktor Chlumský's signed distance field generator](https://github.com/Chlumsky/msdfgen).

The code is still very much a work in progress, current TODOs:
- Rest of the error correction distance modes (other than ``.none``)
- Support for overlapping contours
- Different coloring strategies
- SVG support
- Option for kerning pairs when generating an atlas
- Variable fonts
- Option for Skia geometry preprocessing