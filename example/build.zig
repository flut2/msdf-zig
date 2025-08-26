const std = @import("std");

pub fn build(b: *std.Build) void {
    const optimize = b.standardOptimizeOption(.{});
    const target = b.standardTargetOptions(.{});

    const exe = b.addExecutable(.{
        .name = "Example",
        .root_module = b.createModule(.{
            .root_source_file = b.path("generate.zig"),
            .target = target,
            .optimize = optimize,
        }),
    });
    exe.root_module.link_libc = true;

    const msdf_zig_mod = b.dependency("msdf_zig", .{
        .target = target,
        .optimize = optimize,
    });
    exe.root_module.addImport("msdf-zig", msdf_zig_mod.module("msdf-zig"));

    const zstbi_dep = b.dependency("zstbi", .{
        .target = target,
        .optimize = optimize,
    });
    exe.root_module.addImport("zstbi", zstbi_dep.module("root"));
    exe.linkLibrary(zstbi_dep.artifact("zstbi"));

    exe.step.dependOn(&b.addInstallDirectory(.{
        .source_dir = b.path("assets/"),
        .install_dir = .{ .bin = {} },
        .install_subdir = "assets",
    }).step);

    const run_step = b.addRunArtifact(exe);
    const run = b.step("run", "Run the example");
    run.dependOn(&run_step.step);

    b.installArtifact(exe);
}
