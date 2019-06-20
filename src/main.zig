const std = @import("std");
const sort = std.sort.sort;
const warn = std.debug.warn;
const math = std.math;


const c = @cImport({
    @cInclude("webp/encode.h");
    @cInclude("webp/decode.h"); 
});

var da = std.heap.DirectAllocator.init();
var allocator = &da.allocator;

const ITERATIONS : i64 = 100;

const Result = struct {
    val : f64,
    x : f64,
    y : f64
};

fn mandel(cx : f64, cy : f64) Result {
    var a : f64 = 0;
    var b : f64 = 0;

    var i : i64 = 0;
    while (i < ITERATIONS) {
        const as = a * a;
        const bs = b * b;
        if (as + bs > 4.0)
            return  Result{.val = as + bs,//math.log(f64, 2.0, as + bs) + @intToFloat(f64, ITERATIONS - i),
                .x = a,
                .y = b};

        b = 2.0 * a * b + cy;
        a = as - bs + cx;
        i += 1;
    }
    return Result{.val = a * a + b * b,
                .x = a,
                .y = b};
}

fn clip(a : f64, mina : f64, maxa : f64) f64 {
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
    // warn("{}:{}\n", items.ptr, items.len);
    var tmp : T = 0;
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

    const a : T = items[0];
    const b : T = items[items.len / 2];
    const e : T = items[items.len - 1];
    
    const pivot = blk: {
        if (a > b) {
            if (a < e) {
                break :blk a;
            }
            else if (b > e) {
                break :blk b;
            } else {
                break :blk e;
            }
        } else {
            if (a > e) {
                break :blk a;
            }
            else if (b < e) {
                break :blk b;
            } else {
                break :blk e;
            }
        }
    };

    var cur : [*]T = items.ptr;
    var it : [*]T = items.ptr;
    var pivot_len : usize = 0;

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
    // warn("{} {}\n", pivot_len, items.len);
    if (pivot_len == 0 or pivot_len == items.len)
        return;
    quickSort(T, items[0..pivot_len-1], lessThan);
    quickSort(T, items[pivot_len..], lessThan);
}

pub fn main() anyerror!void {
    warn("Zig Mandelbrot\n");

    const width : i64 = 4000;
    const height : i64 = 4000;
    
    var values = try allocator.alloc(f64, width * height);
    var values_x = try allocator.alloc(f64, width * height);
    var values_y = try allocator.alloc(f64, width * height);
    defer allocator.free(values);

    for (values) |*v|
       v.* = 0;
    

    // const cx : f64 = -1.7400623825;
    // const cy : f64 = 0.0281753397;
    // const rx : f64 = .2;
    // const ry : f64 = .2; 
    const cx : f64 = -0.5;
    const cy : f64 = 0;
    const rx : f64 = 1.5;
    const ry : f64 = 1.2;

    const stepx = rx / f64(width / 2);
    const stepy = ry / f64(height / 2);

    var i : usize = 0;
    var y : i64 = -height / 2;
    while (y < height / 2) {
        var fy = cy + @intToFloat(f64, y) * stepy;

        var x : i64 = -width / 2;
        while (x < width / 2) {
            var fx = cx + @intToFloat(f64, x) * stepx;
            const response = mandel(fx, fy);
            values[i] = response.val;
            values_x[i] = response.x;
            values_y[i] = response.y;
            x += 1;
            i += 1;
        }
        y += 1;
    }


    var image_data = try allocator.alloc(u8, 3 * width * height);
    defer allocator.free(image_data);

    var sorted_values = try allocator.alloc(f64, width * height);
    for (values) |v, n| 
    {
        sorted_values[n] = v;
    }

    warn("sorting {}\n", sorted_values.len);
    // std.sort.sort(f64, sorted_values, std.sort.asc(f64));
    // std.sort.sort(f64, values, std.sort.asc(f64));
    // std.sort.insertionSort(f64, values_x, std.sort.asc(f64));
    // std.sort.insertionSort(f64, values_y, std.sort.asc(f64));
    // std.sort.insertionSort(f64, values, std.sort.asc(f64));

    const start_time = std.time.milliTimestamp();

    // @newStackCall(somestack, quickSort, f64, values, std.sort.asc(f64));
    quickSort(f64, values, std.sort.asc(f64));
    // std.sort.sort(f64, values, std.sort.asc(f64));
    warn("done in {}\n", std.time.milliTimestamp() - start_time);

    // std.sort.sort(f64, values, std.sort.asc(f64));
    
    for (values) |v, n| {
        // warn("\n{} {}\n", v, n);
        var vx = values_x[n];
        var vy = values_y[n];
        const norm = vx * vx + vy * vy + 0.00001;
        vx /= norm;
        vy /= norm;
        // var r = vx;
        // var g = vx * -0.5 - vy * 0.8660;
        // var b = vx * -0.5 + vy * 0.8660;
        var r = v;
        var b = v;
        var g = v;
        // norm = 1.0 / norm;

        const low_bound = 0.0;//math.exp(-norm) / math.e;
        //warn("{}\n", norm);
        r = clip(r, -1, 1) * 0.5 + 0.5;
        g = clip(g, -1, 1) * 0.5 + 0.5;
        b = clip(b, -1, 1) * 0.5 + 0.5;
        const up_bound = 1.0;//clip(math.exp(1.0 - norm) / math.e, 0, 1);
        r = low_bound + (r * (up_bound - low_bound));
        g = low_bound + (g * (up_bound - low_bound));
        b = low_bound + (b * (up_bound - low_bound));

        image_data[n * 3]     = @floatToInt(u8, r * 255);
        image_data[n * 3 + 1] = @floatToInt(u8, g * 255);
        image_data[n * 3 + 2] = @floatToInt(u8, b * 255);
    }

    // warn("before");
    // quickSort(u8, image_data[0..], std.sort.asc(u8));
    // warn("after");
    var output_ptr : [*c]u8 = undefined;


    const encode_len = c.WebPEncodeLosslessRGB(image_data.ptr, c_int(width), c_int(height), c_int(width * 3), &output_ptr);
    warn("encode len: {} bytes\n", encode_len);
    defer c.WebPFree(output_ptr);

    const slice = output_ptr[0..encode_len];
    try std.io.writeFile("test.webp", slice);
}
