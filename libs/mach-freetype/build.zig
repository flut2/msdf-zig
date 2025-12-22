const std = @import("std");

pub fn build(b: *std.Build) !void {
    const optimize = b.standardOptimizeOption(.{});
    const target = b.standardTargetOptions(.{});
    const use_system_zlib = b.option(bool, "use_system_zlib", "Use system zlib") orelse false;

    const freetype_module = b.addModule("mach-freetype", .{
        .root_source_file = b.path("src/freetype.zig"),
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

    const freetype_dep = b.dependency("freetype", .{
        .target = target,
        .optimize = optimize,
        .use_system_zlib = use_system_zlib,
    });
    freetype_tests.root_module.linkLibrary(freetype_dep.artifact("freetype"));
    freetype_module.linkLibrary(freetype_dep.artifact("freetype"));

    const test_step = b.step("test", "Run library tests");
    test_step.dependOn(&b.addRunArtifact(freetype_tests).step);
}
