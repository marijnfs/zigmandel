const std = @import("std");
const sort = std.sort.sort;
const math = std.math;
const complex = std.math.complex;

const c = @cImport({
    @cInclude("webp/encode.h");
    @cInclude("webp/decode.h");
});

var allocator = std.heap.page_allocator;

const ITERATIONS: i64 = 1000;

const Result = struct {
    val: f64,
    x: f64,
    y: f64,

    fn init() Result {
        return Result{ .val = 0, .x = 0, .y = 0 };
    }
};

const C64 = complex.Complex(f64);

fn zeta(cx: f64, cy: f64) Result {
    // std.log.info("{} {}", .{cx, cy});

    var sum = C64{ .re = 0, .im = 0 };
    // Correctly define s = cx + i*cy (not -cx - i*cy)
    const s = C64{ .re = cx, .im = cy };

    // Compute n^{-s} = exp(-s * ln(n))
    for (0..ITERATIONS) |n| { // Start at n=1 (not 0) to avoid undefined terms
        // for (1..1001) |n| {  // Start at n=1 (not 0) to avoid undefined terms
        const n_complex = C64{ .re = @floatFromInt(n), .im = 0 };
        const term = complex.pow(n_complex, s.neg());
        sum = sum.add(term);
    }

    // Return actual zeta value (not s.re) and correct coordinates

    return Result{ .val = sum.re, .x = cx, .y = cy };
}

fn mandel(cx: f64, cy: f64) Result {
    var a: f64 = 0;
    var b: f64 = 0;

    var i: i64 = 0;
    while (i < ITERATIONS) {
        const as = a * a;
        const bs = b * b;
        if (as + bs > 4.0)
            return Result{ .val = math.sqrt(as + bs) + @as(f32, @floatFromInt(ITERATIONS - i)), .x = a, .y = b };

        b = 2.0 * a * b + cy;
        a = as - bs + cx;
        i += 1;
    }
    return Result{ .val = math.sqrt(a * a + b * b), .x = a, .y = b };
}

fn clip(a: f64, mina: f64, maxa: f64) f64 {
    if (a == std.math.inf(f64))
        return maxa;
    if (a == -std.math.inf(f64))
        return mina;
    if (a == std.math.nan(f64))
        return maxa;
    if (a == -std.math.nan(f64))
        return mina;
    if (a > maxa)
        return maxa;
    if (a < mina)
        return mina;
    return a;
}

pub fn quickSort(comptime T: type, items: []T, lessThan: fn (lhs: T, rhs: T) bool) void {
    // std.log.info("{}:{}\n", items.ptr, items.len);
    var tmp: T = undefined;
    if (items.len <= 1)
        return;
    if (items.len == 2) {
        if (!lessThan(items[0], items[1])) {
            tmp = items[0];
            items[0] = items[1];
            items[1] = tmp;
        }
        return;
    }

    // if (items.len < 8) {
    //     std.sort.insertionSort(f64, items, lessThan);
    //     return;
    // }

    const a: T = items[0];
    const b: T = items[items.len / 2];
    const e: T = items[items.len - 1];

    const pivot = blk: {
        if (lessThan(b, a)) {
            if (lessThan(a, e)) {
                break :blk a;
            } else if (lessThan(e, b)) {
                break :blk b;
            } else {
                break :blk e;
            }
        } else {
            if (lessThan(e, a)) {
                break :blk a;
            } else if (lessThan(b, e)) {
                break :blk b;
            } else {
                break :blk e;
            }
        }
    };

    var cur: [*]T = items.ptr;
    var it: [*]T = items.ptr;
    var pivot_len: usize = 0;

    for (items) |v| {
        if (lessThan(v, pivot)) {
            if (cur == it) {
                it += 1;
                continue;
            }
            tmp = cur[0];
            cur[0] = it[0];
            it[0] = tmp;
            cur += 1;
            pivot_len += 1;
        }
        it += 1;
    }
    // std.log.info("{} {}\n", pivot_len, items.len);
    if (pivot_len == 0 or pivot_len == items.len)
        return;
    quickSort(T, items[0 .. pivot_len - 1], lessThan);
    quickSort(T, items[pivot_len..], lessThan);
}

pub fn main() anyerror!void {
    std.log.info("Zig Mandelbrot\n", .{});

    const width: i64 = 5000;
    const height: i64 = 5000;

    const values = try allocator.alloc(f64, width * height);
    defer allocator.free(values);

    for (values) |*v|
        v.* = 0;

    // const cx : f64 = -1.7400623825;
    // const cy : f64 = 0.0281753397;
    // const rx : f64 = .2;
    // const ry : f64 = .2;
    // const cx : f64 = -0.5;
    // const cy : f64 = 0;
    // const rx : f64 = 1.5;
    // const ry : f64 = 1.2;

    //Star
    // const basephase = math.pi * 1.9;
    // const cx : f64 = -0.707132;
    // const cy : f64 = -0.353294;
    // const rx : f64 = 0.015;
    // const ry : f64 = rx * (@intToFloat(f64, height) / @intToFloat(f64, width));

    //Curv
    // const basephase = math.pi * 1.9;
    // const cx : f64 = -0.664092;
    // const cy : f64 = -0.327854;
    // const rx : f64 = 0.0000329;
    // const ry : f64 = rx * (@intToFloat(f64, height) / @intToFloat(f64, width));

    //seahorse

    const basephase = math.pi * 1.9;
    //const cx: f64 = -0.664092;
    //const cy: f64 = -0.327654;
    //const rx: f64 = 0.0000329;
    //const ry: f64 = rx * (@as(f32, @floatFromInt(height)) / @as(f32, @floatFromInt(width)));

    const cx: f64 = -3.0;
    const cy: f64 = -3.0;
    const rx: f64 = 6.0;
    const ry: f64 = 6.0;

    //spiral
    //const basephase = math.pi * 0.0;
    //const cx : f64 = -0.6818117;
    //const cy : f64 = -0.32214989;
    //const rx : f64 = 0.0102;
    //const ry : f64 = rx * (@intToFloat(f64, height) / @intToFloat(f64, width));

    const dx = rx / @as(f64, @floatFromInt(width));
    const dy = ry / @as(f64, @floatFromInt(height));

    // Function for threads

    const F = struct {
        width_: usize,
        cx_: f64,
        cy_: f64,
        dx_: f64,
        dy_: f64,
        fy_: f64,
        i_: usize,

        fn f(self: *@This(), values_: []f64) void {
            const alias_dx = self.dx_ / 4.0;
            const alias_dy = self.dx_ / 4.0;

            var x: usize = 0;
            while (x < self.width_) {
                const fx = self.cx_ + @as(f64, @floatFromInt(x)) * self.dx_;

                var response = Result.init();

                //const f = mandel;
                const eval = zeta;
                const r1 = eval(fx + alias_dx, self.fy_ + alias_dy);
                const r2 = eval(fx + alias_dx, self.fy_ - alias_dy);
                const r3 = eval(fx - alias_dx, self.fy_ + alias_dy);
                const r4 = eval(fx - alias_dx, self.fy_ - alias_dy);
                response.val = r1.val + r2.val + r3.val + r4.val;

                values_[self.i_] = response.val;
                x += 1;
                self.i_ += 1;
            }
        }
    };

    {

        // Thread pool

        var pool: std.Thread.Pool = undefined;
        try pool.init(.{
            .allocator = allocator,
        });
        defer pool.deinit();

        var i: usize = 0;
        var y: i64 = 0;

        while (y < height) {
            const fy = cy + @as(f64, @floatFromInt(y)) * dy;
            const f = try allocator.create(F);
            f.* = .{
                .width_ = width,
                .cx_ = cx,
                .cy_ = cy,
                .dx_ = dx,
                .dy_ = dy,
                .fy_ = fy,
                .i_ = i,
            };
            i += width;
            try pool.spawn(F.f, .{ f, values });

            // var x: i64 = -width / 2;
            // while (x < width / 2) {
            //     const fx = cx + @as(f64, @floatFromInt(x)) * dx;

            //     var response = Result.init();

            //     //const f = mandel;
            //     const f = zeta;
            //     const r1 = f(fx + alias_dx, fy + alias_dy);
            //     const r2 = f(fx + alias_dx, fy - alias_dy);
            //     const r3 = f(fx - alias_dx, fy + alias_dy);
            //     const r4 = f(fx - alias_dx, fy - alias_dy);
            //     response.val = r1.val + r2.val + r3.val + r4.val;

            //     values[i] = response.val;
            //     x += 1;
            //     i += 1;
            // }
            y += 1;
        }
    }
    var image_data = try allocator.alloc(u8, 3 * width * height);
    defer allocator.free(image_data);

    var sorted_values = try allocator.alloc(f64, width * height);
    for (0.., values) |n, v| {
        sorted_values[n] = v;
    }

    std.log.info("sorting {}\n", .{sorted_values.len});

    const start_time = std.time.milliTimestamp();
    const ValueIndex = struct {
        n: usize,
        val: f64,
    };
    var ivalues = try allocator.alloc(ValueIndex, values.len);

    for (0.., values) |n, v| {
        ivalues[n].n = n;
        ivalues[n].val = v;
    }

    const compare = struct {
        fn inner(_: void, v1: ValueIndex, v2: ValueIndex) bool {
            return v1.val < v2.val;
        }
    };

    std.sort.heap(ValueIndex, ivalues, {}, compare.inner);
    // quickSort(ValueIndex, ivalues, compare.inner);

    // quickSort(f64, values, std.sort.asc(f64));
    std.log.info("done in {}\n", .{std.time.milliTimestamp() - start_time});

    for (0.., ivalues) |n, v| {
        const rel = @as(f64, @floatFromInt(n)) / @as(f64, @floatFromInt(ivalues.len));
        // const r_phase = (math.sin(math.pi * 2.0 * rel) + 1.0) / 2.0;
        // const g_phase = (math.sin(math.pi * 2.0 * rel) + 1.0) / 2.0;
        // const b_phase = (math.sin(math.pi * 2.0 * rel) + 1.0) / 2.0;
        // const r_phase = (math.sin(math.pi * 2.0 * rel) + 1.0) / 2.0;
        // const g_phase = (math.sin(math.pi * 2.0 * rel + math.pi * 2.0 / 3.0) + 1.0) / 2.0;
        // const b_phase = (math.sin(math.pi * 2.0 * rel + math.pi * 4.0 / 3.0) + 1.0) / 2.0;
        // const within = v.val < 1.0;

        const diff = math.pi + basephase;

        const phase1 = (math.sin(math.pi * 10.0 * rel + diff) + 1.0) / 2.0;
        const phase2 = (math.sin(math.pi * 3.33 * rel + diff + math.pi / 2.0) + 1.0) / 2.0;

        // std.log.info("{} {} {}\n", r_phase, g_phase, b_phase);

        const r_phase: f64 = @max(0.0, @min(1.0, 0.1 + 0.3 * phase1 + 0.1 * phase2));
        //const g_phase = 0.31 + 0.3  * phase1 - 0.1543 * phase2;
        const g_phase: f64 = @max(0.0, @min(1.0, 0.1 + 0.2 * phase1 - 0.1543 * phase2));
        const b_phase: f64 = @max(0.0, @min(1.0, 0.32 + 0.0 * phase2 + 0.6 * phase2));

        image_data[v.n * 3] = @as(u8, @intFromFloat(r_phase * 255));
        image_data[v.n * 3 + 1] = @as(u8, @intFromFloat(g_phase * 255));
        image_data[v.n * 3 + 2] = @as(u8, @intFromFloat(b_phase * 255));
    }

    var output_ptr: [*c]u8 = undefined;

    const encode_len = c.WebPEncodeLosslessRGB(image_data.ptr, width, height, width * 3, &output_ptr);
    std.log.info("encode len: {} bytes\n", .{encode_len});
    defer c.WebPFree(output_ptr);

    const slice = output_ptr[0..encode_len];
    var file = try std.fs.cwd().createFile("test.webp", std.fs.File.CreateFlags{});
    try file.writeAll(slice);
}
