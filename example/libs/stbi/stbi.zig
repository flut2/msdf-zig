const std = @import("std");

var stbi_allocator: ?std.mem.Allocator = null;
var pointer_size_map: std.AutoHashMapUnmanaged(usize, usize) = .empty;
var alloc_mutex: std.Thread.Mutex = .{};
const alignment: std.mem.Alignment = .of(std.c.max_align_t);

fn allocatorMissing() noreturn {
    @panic("stbi: Allocator is missing, set it through `stbi.init()`");
}

fn outOfMemory() noreturn {
    @panic("stbi: Out of memory");
}

fn stbiMalloc(size: usize) callconv(.c) ?*anyopaque {
    const allocator = stbi_allocator orelse allocatorMissing();

    alloc_mutex.lock();
    defer alloc_mutex.unlock();

    const mem = allocator.alignedAlloc(u8, alignment, size) catch outOfMemory();
    pointer_size_map.put(allocator, @intFromPtr(mem.ptr), size) catch outOfMemory();
    return mem.ptr;
}

fn stbiRealloc(maybe_ptr: ?*anyopaque, new_size: usize) callconv(.c) ?*anyopaque {
    const allocator = stbi_allocator orelse allocatorMissing();

    alloc_mutex.lock();
    defer alloc_mutex.unlock();

    const old_size = if (maybe_ptr) |p| pointer_size_map.fetchRemove(@intFromPtr(p)).?.value else 0;
    const old_mem: [*]align(alignment.toByteUnits()) u8 = if (maybe_ptr) |p| @ptrCast(@alignCast(p)) else &.{};
    const new_mem = allocator.realloc(old_mem[0..old_size], new_size) catch outOfMemory();
    pointer_size_map.put(allocator, @intFromPtr(new_mem.ptr), new_size) catch outOfMemory();
    return new_mem.ptr;
}

fn stbiFree(maybe_ptr: ?*anyopaque) callconv(.c) void {
    const allocator = stbi_allocator orelse allocatorMissing();
    const ptr = maybe_ptr orelse return;

    alloc_mutex.lock();
    defer alloc_mutex.unlock();

    const kv = pointer_size_map.fetchRemove(@intFromPtr(ptr)) orelse {
        std.log.err("stbi: Invalid free attempted on {*}", .{ptr});
        return;
    };
    const mem: [*]align(alignment.toByteUnits()) u8 = @ptrCast(@alignCast(ptr));
    allocator.free(mem[0..kv.value]);
}

fn stbirMalloc(size: usize, _: ?*anyopaque) callconv(.c) ?*anyopaque {
    return stbiMalloc(size);
}

fn stbirFree(maybe_ptr: ?*anyopaque, _: ?*anyopaque) callconv(.c) void {
    stbiFree(maybe_ptr);
}

pub fn init(allocator: std.mem.Allocator) void {
    if (stbi_allocator != null)
        @panic("stbi: Library already initialized");
    stbi_allocator = allocator;

    stbiMallocPtr = stbiMalloc;
    stbiReallocPtr = stbiRealloc;
    stbiFreePtr = stbiFree;
    stbirMallocPtr = stbirMalloc;
    stbirFreePtr = stbirFree;
    stbiwMallocPtr = stbiMalloc;
    stbiwReallocPtr = stbiRealloc;
    stbiwFreePtr = stbiFree;
}

pub fn deinit() void {
    const allocator = stbi_allocator orelse return;
    pointer_size_map.deinit(allocator);
    stbi_allocator = null;
}

pub const JpgWriteSettings = struct {
    quality: u32,
};

pub const ImageWriteFormat = union(enum) {
    png,
    jpg: JpgWriteSettings,
};

pub const ImageWriteError = error{CouldNotWriteImage};

pub const Image = struct {
    data: []u8,
    width: u32,
    height: u32,
    num_components: u32,
    bytes_per_component: u32,
    bytes_per_row: u32,

    pub const invalid: Image = .{
        .data = &.{},
        .width = std.math.maxInt(u32),
        .height = std.math.maxInt(u32),
        .num_components = std.math.maxInt(u32),
        .bytes_per_component = std.math.maxInt(u32),
        .bytes_per_row = std.math.maxInt(u32),
    };

    pub fn info(pathname: [:0]const u8) struct {
        is_supported: bool,
        width: u32,
        height: u32,
        num_components: u32,
    } {
        if (stbi_allocator == null) allocatorMissing();

        var w: c_int = 0;
        var h: c_int = 0;
        var c: c_int = 0;
        const is_supported = stbi_info(pathname, &w, &h, &c);
        return .{
            .is_supported = is_supported == 1,
            .width = @intCast(w),
            .height = @intCast(h),
            .num_components = @intCast(c),
        };
    }

    pub fn loadFromFile(path: [:0]const u8, forced_comps: u32) !Image {
        if (stbi_allocator == null) allocatorMissing();

        var x: c_int = undefined;
        var y: c_int = undefined;
        var c: c_int = undefined;
        const data = stbi_load(
            path,
            &x,
            &y,
            &c,
            @intCast(forced_comps),
        ) orelse {
            std.log.err("stbi: Loading image from path `{s}` failed", .{path});
            return error.ImageInitFailed;
        };

        const num_components: u32 = if (forced_comps == 0) @intCast(c) else forced_comps;
        const width: u32 = @intCast(x);
        const height: u32 = @intCast(y);
        const bytes_per_row = width * num_components;

        return .{
            .data = data[0 .. height * bytes_per_row],
            .width = width,
            .height = height,
            .num_components = num_components,
            .bytes_per_component = 1,
            .bytes_per_row = bytes_per_row,
        };
    }

    pub fn loadFromMemory(data: []const u8, forced_comps: u32) !Image {
        if (stbi_allocator == null) allocatorMissing();

        var x: c_int = undefined;
        var y: c_int = undefined;
        var c: c_int = undefined;
        const image_data = stbi_load_from_memory(
            data.ptr,
            @intCast(data.len),
            &x,
            &y,
            &c,
            @intCast(forced_comps),
        ) orelse {
            std.log.err("stbi: Loading image from data `{*}` failed", .{data.ptr});
            return error.ImageInitFailed;
        };

        const num_components: u32 = if (forced_comps == 0) @intCast(c) else forced_comps;
        const width: u32 = @intCast(x);
        const height: u32 = @intCast(y);
        const bytes_per_row = width * num_components;

        return .{
            .data = image_data[0 .. height * bytes_per_row],
            .width = width,
            .height = height,
            .num_components = num_components,
            .bytes_per_component = 1,
            .bytes_per_row = bytes_per_row,
        };
    }

    pub fn createEmpty(width: u32, height: u32, num_components: u32, opts: struct {
        bytes_per_component: u32 = 0,
        bytes_per_row: u32 = 0,
    }) !Image {
        if (stbi_allocator == null) allocatorMissing();

        const bytes_per_component = if (opts.bytes_per_component == 0) 1 else opts.bytes_per_component;
        const bytes_per_row = if (opts.bytes_per_row == 0)
            width * num_components * bytes_per_component
        else
            opts.bytes_per_row;

        const size = height * bytes_per_row;
        const data: [*]u8 = @ptrCast(stbiMalloc(size));
        const data_slice = data[0..size];
        @memset(data_slice, 0);

        return .{
            .data = data_slice,
            .width = width,
            .height = height,
            .num_components = num_components,
            .bytes_per_component = bytes_per_component,
            .bytes_per_row = bytes_per_row,
        };
    }

    pub fn resize(image: *const Image, new_width: u32, new_height: u32) Image {
        if (stbi_allocator == null) allocatorMissing();

        const new_bytes_per_row = new_width * image.num_components * image.bytes_per_component;
        const new_size = new_height * new_bytes_per_row;
        const new_data: [*]u8 = @ptrCast(stbiMalloc(new_size));
        stbir_resize_uint8(
            image.data.ptr,
            @intCast(image.width),
            @intCast(image.height),
            0,
            new_data,
            @intCast(new_width),
            @intCast(new_height),
            0,
            @intCast(image.num_components),
        );
        return .{
            .data = new_data[0..new_size],
            .width = new_width,
            .height = new_height,
            .num_components = image.num_components,
            .bytes_per_component = image.bytes_per_component,
            .bytes_per_row = new_bytes_per_row,
        };
    }

    pub fn writeToFile(
        image: Image,
        filename: [:0]const u8,
        image_format: ImageWriteFormat,
    ) ImageWriteError!void {
        if (stbi_allocator == null) allocatorMissing();

        const w: c_int = @intCast(image.width);
        const h: c_int = @intCast(image.height);
        const c: c_int = @intCast(image.num_components);
        const result = switch (image_format) {
            .png => stbi_write_png(filename.ptr, w, h, c, image.data.ptr, 0),
            .jpg => |settings| stbi_write_jpg(
                filename.ptr,
                w,
                h,
                c,
                image.data.ptr,
                @intCast(settings.quality),
            ),
        };

        if (result == 0) return ImageWriteError.CouldNotWriteImage;
    }

    pub fn writeToFn(
        image: Image,
        writeFn: *const fn (ctx: ?*anyopaque, data: ?*anyopaque, size: c_int) callconv(.c) void,
        context: ?*anyopaque,
        image_format: ImageWriteFormat,
    ) ImageWriteError!void {
        if (stbi_allocator == null) allocatorMissing();

        const w: c_int = @intCast(image.width);
        const h: c_int = @intCast(image.height);
        const c: c_int = @intCast(image.num_components);
        const result = switch (image_format) {
            .png => stbi_write_png_to_func(writeFn, context, w, h, c, image.data.ptr, 0),
            .jpg => |settings| stbi_write_jpg_to_func(
                writeFn,
                context,
                w,
                h,
                c,
                image.data.ptr,
                @intCast(settings.quality),
            ),
        };

        if (result == 0) return ImageWriteError.CouldNotWriteImage;
    }

    pub fn deinit(image: *Image) void {
        stbi_image_free(image.data.ptr);
        image.* = undefined;
    }
};

pub fn setFlipVerticallyOnLoad(should_flip: bool) void {
    stbi_set_flip_vertically_on_load(@intFromBool(should_flip));
}

pub fn setFlipVerticallyOnWrite(should_flip: bool) void {
    stbi_flip_vertically_on_write(@intFromBool(should_flip));
}

extern var stbiMallocPtr: ?*const fn (size: usize) callconv(.c) ?*anyopaque;
extern var stbiReallocPtr: ?*const fn (ptr: ?*anyopaque, size: usize) callconv(.c) ?*anyopaque;
extern var stbiFreePtr: ?*const fn (maybe_ptr: ?*anyopaque) callconv(.c) void;
extern var stbiwMallocPtr: ?*const fn (size: usize) callconv(.c) ?*anyopaque;
extern var stbiwReallocPtr: ?*const fn (ptr: ?*anyopaque, size: usize) callconv(.c) ?*anyopaque;
extern var stbiwFreePtr: ?*const fn (maybe_ptr: ?*anyopaque) callconv(.c) void;
extern var stbirMallocPtr: ?*const fn (size: usize, maybe_context: ?*anyopaque) callconv(.c) ?*anyopaque;
extern var stbirFreePtr: ?*const fn (maybe_ptr: ?*anyopaque, maybe_context: ?*anyopaque) callconv(.c) void;

extern fn stbi_info(filename: [*:0]const u8, x: *c_int, y: *c_int, comp: *c_int) c_int;

extern fn stbi_load(
    filename: [*:0]const u8,
    x: *c_int,
    y: *c_int,
    channels_in_file: *c_int,
    desired_channels: c_int,
) ?[*]u8;

extern fn stbi_load_16(
    filename: [*:0]const u8,
    x: *c_int,
    y: *c_int,
    channels_in_file: *c_int,
    desired_channels: c_int,
) ?[*]u16;

extern fn stbi_loadf(
    filename: [*:0]const u8,
    x: *c_int,
    y: *c_int,
    channels_in_file: *c_int,
    desired_channels: c_int,
) ?[*]f32;

pub extern fn stbi_load_from_memory(
    buffer: [*]const u8,
    len: c_int,
    x: *c_int,
    y: *c_int,
    channels_in_file: *c_int,
    desired_channels: c_int,
) ?[*]u8;

pub extern fn stbi_loadf_from_memory(
    buffer: [*]const u8,
    len: c_int,
    x: *c_int,
    y: *c_int,
    channels_in_file: *c_int,
    desired_channels: c_int,
) ?[*]f32;

extern fn stbi_image_free(image_data: ?[*]u8) void;

extern fn stbi_set_flip_vertically_on_load(flag_true_if_should_flip: c_int) void;
extern fn stbi_flip_vertically_on_write(flag: c_int) void; // flag is non-zero to flip data vertically

extern fn stbir_resize_uint8(
    input_pixels: [*]const u8,
    input_w: c_int,
    input_h: c_int,
    input_stride_in_bytes: c_int,
    output_pixels: [*]u8,
    output_w: c_int,
    output_h: c_int,
    output_stride_in_bytes: c_int,
    num_channels: c_int,
) void;

extern fn stbi_write_jpg(
    filename: [*:0]const u8,
    w: c_int,
    h: c_int,
    comp: c_int,
    data: [*]const u8,
    quality: c_int,
) c_int;

extern fn stbi_write_png(
    filename: [*:0]const u8,
    w: c_int,
    h: c_int,
    comp: c_int,
    data: [*]const u8,
    stride_in_bytes: c_int,
) c_int;

extern fn stbi_write_png_to_func(
    func: *const fn (?*anyopaque, ?*anyopaque, c_int) callconv(.c) void,
    context: ?*anyopaque,
    w: c_int,
    h: c_int,
    comp: c_int,
    data: [*]const u8,
    stride_in_bytes: c_int,
) c_int;

extern fn stbi_write_jpg_to_func(
    func: *const fn (?*anyopaque, ?*anyopaque, c_int) callconv(.c) void,
    context: ?*anyopaque,
    x: c_int,
    y: c_int,
    comp: c_int,
    data: [*]const u8,
    quality: c_int,
) c_int;
