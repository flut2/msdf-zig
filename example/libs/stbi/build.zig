const std = @import("std");

pub fn build(b: *std.Build) void {
    const optimize = b.standardOptimizeOption(.{});
    const target = b.standardTargetOptions(.{});

    _ = b.addModule("root", .{ .root_source_file = b.path("stbi.zig") });

    const stbi_lib = b.addLibrary(.{
        .name = "stbi",
        .linkage = .static,
        .root_module = b.createModule(.{
            .target = target,
            .optimize = optimize,
            .link_libc = true,
        }),
    });

    const stb_dep = b.dependency("stb", .{ .target = target, .optimize = optimize });
    stbi_lib.root_module.addIncludePath(stb_dep.path("."));
    stbi_lib.root_module.addCSourceFile(.{
        .file = b.path("stbi.c"),
        .flags = &.{ "-std=c99", "-fno-sanitize=undefined" },
    });
    b.installArtifact(stbi_lib);
}
