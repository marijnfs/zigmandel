const std = @import("std");
const warn = std.debug.warn;
const math = std.math;


const c = @cImport({
    @cInclude("webp/encode.h");
    @cInclude("webp/decode.h"); 
});

var da = std.heap.DirectAllocator.init();
var allocator = &da.allocator;

const ITERATIONS : i64 = 10;

fn mandel(cx : f64, cy : f64) f64 {
    var a : f64 = 0;
    var b : f64 = 0;

    var i : i64 = 0;
    while (i < ITERATIONS) {
        const as = a * a;
        const bs = b * b;
        if (as + bs > 4.0)
            return  math.log(f64, 2.0, as + bs) + @intToFloat(f64, ITERATIONS - i);

        b = 2.0 * a * b + cy;
        a = as - bs + cx;
        i += 1;
    }
    return  math.log(f64, 2.0, a * a + b * b);
}

pub fn main() anyerror!void {
    std.debug.warn("All your base are belong to us.\n");
    const width : i64 = 2000;
    const height : i64 = 2000;
    
    var values = try allocator.alloc(f64, width * height);
    defer allocator.free(values);

    for (values) |*v|
       v.* = 0;
    

    const cx : f64 = -0.5;
    const cy : f64 = 0;
    const rx : f64 = 1;
    const ry : f64 = 1;

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
            values[i] = response;
            x += 1;
            i += 1;
        }
        y += 1;
    }


    var image_data = try allocator.alloc(u8, 3 * width * height);
    defer allocator.free(image_data);

    for (values) |v, n| {
        // warn("\n{} {}\n", v, n);
        var value : f64 = v;
        value *= 255;
        if (value < 0)
            value = 0;
        if (value > 255.0)
            value = 255.0;

        if (v > 1.0) {
            @memset(image_data.ptr + n * 3, @floatToInt(u8, value), 3);
        } else {
            @memset(image_data.ptr + n * 3, @floatToInt(u8, value), 3);
        }
    }

    var output_ptr : [*c]u8 = undefined;
    warn("{}\n", output_ptr);
    const encode_len = c.WebPEncodeLosslessRGB(image_data.ptr, c_int(width), c_int(height), c_int(width * 3), &output_ptr);
    warn("{}\n", output_ptr);
    warn("encode len {}\n", encode_len);
    defer c.WebPFree(output_ptr);

    const slice = output_ptr[0..encode_len];
    try std.io.writeFile("test.webp", slice);
}
