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

    const exe = b.addExecutable(.{
        .name = "Example",
        .root_module = b.createModule(.{
            .root_source_file = b.path("example/generate.zig"),
            .target = target,
            .optimize = optimize,
        }),
    });
    exe.root_module.link_libc = true;
    exe.root_module.addImport("msdf-zig", msdf_zig_mod);

    const zstbi_dep = b.dependency("zstbi", .{
        .target = target,
        .optimize = optimize,
    });
    exe.root_module.addImport("zstbi", zstbi_dep.module("root"));
    exe.linkLibrary(zstbi_dep.artifact("zstbi"));

    const run_step = b.addRunArtifact(exe);
    const run = b.step("run", "Run the example");
    run.dependOn(&run_step.step);

    b.installArtifact(exe);
}
