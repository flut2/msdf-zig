const std = @import("std");

pub fn build(b: *std.Build) void {
    const optimize = b.standardOptimizeOption(.{});
    const target = b.standardTargetOptions(.{});

    const msdf_zig_mod = b.addModule("msdf-zig", .{ .root_source_file = b.path("src/Generator.zig") });
    const freetype_dep = b.dependency("mach_freetype", .{
        .optimize = optimize,
        .target = target,
    });
    msdf_zig_mod.addImport("mach-freetype", freetype_dep.module("mach-freetype"));

    const turbopack_dep = b.dependency("turbopack", .{
        .target = target,
        .optimize = optimize,
    });
    msdf_zig_mod.addImport("turbopack", turbopack_dep.module("turbopack"));
}
