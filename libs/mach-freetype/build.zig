const std = @import("std");

pub fn build(b: *std.Build) !void {
    const optimize = b.standardOptimizeOption(.{});
    const target = b.standardTargetOptions(.{});
    const use_system_zlib = b.option(bool, "use_system_zlib", "Use system zlib") orelse false;
    const enable_brotli = b.option(bool, "enable_brotli", "Build brotli") orelse true;

    const freetype_module = b.addModule("mach-freetype", .{
        .root_source_file = b.path("src/freetype.zig"),
    });
    const harfbuzz_module = b.addModule("mach-harfbuzz", .{
        .root_source_file = b.path("src/harfbuzz.zig"),
        .imports = &.{.{ .name = "freetype", .module = freetype_module }},
    });

    const freetype_tests = b.addTest(.{
        .name = "freetype-tests",
        .root_module = b.createModule(.{
            .root_source_file = b.path("src/freetype.zig"),
            .target = target,
            .optimize = optimize,
        }),
    });
    freetype_tests.root_module.addImport("freetype", freetype_module);

    const harfbuzz_tests = b.addTest(.{
        .name = "harfbuzz-tests",
        .root_module = b.createModule(.{
            .root_source_file = b.path("src/harfbuzz.zig"),
            .target = target,
            .optimize = optimize,
        }),
    });
    harfbuzz_tests.root_module.addImport("freetype", freetype_module);
    harfbuzz_tests.root_module.addImport("harfbuzz", harfbuzz_module);

    const freetype_dep = b.dependency("freetype", .{
        .target = target,
        .optimize = optimize,
        .use_system_zlib = use_system_zlib,
        .enable_brotli = enable_brotli,
    });
    freetype_tests.root_module.linkLibrary(freetype_dep.artifact("freetype"));
    freetype_module.linkLibrary(freetype_dep.artifact("freetype"));
    harfbuzz_module.linkLibrary(freetype_dep.artifact("freetype"));
    harfbuzz_tests.root_module.linkLibrary(freetype_dep.artifact("freetype"));

    const harfbuzz_dep = b.dependency("harfbuzz", .{
        .target = target,
        .optimize = optimize,
        .enable_freetype = true,
        .freetype_use_system_zlib = use_system_zlib,
        .freetype_enable_brotli = enable_brotli,
    });
    harfbuzz_module.linkLibrary(harfbuzz_dep.artifact("harfbuzz"));
    harfbuzz_tests.root_module.linkLibrary(harfbuzz_dep.artifact("harfbuzz"));

    const test_step = b.step("test", "Run library tests");
    test_step.dependOn(&b.addRunArtifact(freetype_tests).step);
    test_step.dependOn(&b.addRunArtifact(harfbuzz_tests).step);
}
