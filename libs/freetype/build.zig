const std = @import("std");

const ft_sources = [_][]const u8{
    "autofit/autofit.c",
    "base/ftbase.c",
    "base/ftsystem.c",
    "base/ftdebug.c",
    "base/ftbbox.c",
    "base/ftbdf.c",
    "base/ftbitmap.c",
    "base/ftcid.c",
    "base/ftfstype.c",
    "base/ftgasp.c",
    "base/ftglyph.c",
    "base/ftgxval.c",
    "base/ftinit.c",
    "base/ftmm.c",
    "base/ftotval.c",
    "base/ftpatent.c",
    "base/ftpfr.c",
    "base/ftstroke.c",
    "base/ftsynth.c",
    "base/fttype1.c",
    "base/ftwinfnt.c",
    "bdf/bdf.c",
    "bzip2/ftbzip2.c",
    "cache/ftcache.c",
    "cff/cff.c",
    "cid/type1cid.c",
    "gzip/ftgzip.c",
    "lzw/ftlzw.c",
    "pcf/pcf.c",
    "pfr/pfr.c",
    "psaux/psaux.c",
    "pshinter/pshinter.c",
    "psnames/psnames.c",
    "raster/raster.c",
    "sdf/sdf.c",
    "sfnt/sfnt.c",
    "smooth/smooth.c",
    "svg/svg.c",
    "truetype/truetype.c",
    "type1/type1.c",
    "type42/type42.c",
    "winfonts/winfnt.c",
};

const brotli_sources = [_][]const u8{
    "common/constants.c",
    "common/context.c",
    "common/dictionary.c",
    "common/platform.c",
    "common/shared_dictionary.c",
    "common/transform.c",
    "dec/bit_reader.c",
    "dec/decode.c",
    "dec/huffman.c",
    "dec/state.c",
    "enc/backward_references.c",
    "enc/backward_references_hq.c",
    "enc/bit_cost.c",
    "enc/block_splitter.c",
    "enc/brotli_bit_stream.c",
    "enc/cluster.c",
    "enc/command.c",
    "enc/compound_dictionary.c",
    "enc/compress_fragment.c",
    "enc/compress_fragment_two_pass.c",
    "enc/dictionary_hash.c",
    "enc/encode.c",
    "enc/encoder_dict.c",
    "enc/entropy_encode.c",
    "enc/fast_log.c",
    "enc/histogram.c",
    "enc/literal_cost.c",
    "enc/memory.c",
    "enc/metablock.c",
    "enc/static_dict.c",
    "enc/utf8_util.c",
};

pub fn build(b: *std.Build) void {
    const target = b.standardTargetOptions(.{});
    const optimize = b.standardOptimizeOption(.{});

    const use_system_zlib = b.option(bool, "use_system_zlib", "Use system zlib") orelse false;

    const lib = b.addLibrary(.{
        .name = "freetype",
        .linkage = .static,
        .root_module = b.createModule(.{
            .target = target,
            .optimize = optimize,
            .link_libc = true,
        }),
    });

    const ft_dep = b.dependency("freetype", .{});
    lib.root_module.addIncludePath(ft_dep.path("include"));
    lib.root_module.addCMacro("FT2_BUILD_LIBRARY", "1");
    lib.root_module.addCMacro("TT_CONFIG_OPTION_GPOS_KERNING", "");

    if (use_system_zlib)
        lib.root_module.addCMacro("FT_CONFIG_OPTION_SYSTEM_ZLIB", "1");

    lib.root_module.addCMacro("HAVE_UNISTD_H", "1");
    lib.root_module.addCSourceFiles(.{
        .root = ft_dep.path("src"),
        .files = &ft_sources,
        .flags = &.{},
    });
    if (target.result.os.tag == .macos) lib.root_module.addCSourceFile(.{
        .file = ft_dep.path("src/base/ftmac.c"),
        .flags = &.{},
    });
    lib.installHeadersDirectory(ft_dep.path("include/freetype"), "freetype", .{});
    lib.installHeader(ft_dep.path("include/ft2build.h"), "ft2build.h");
    b.installArtifact(lib);
}
