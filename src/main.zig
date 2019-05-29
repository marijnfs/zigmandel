const std = @import("std");
const warn = std.debug.warn;
const math = std.math;

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
    const width : i64 = 200;
    const height : i64 = 30;
    
    var values = try allocator.alloc(f64, width * height);

    for (values) |*v|
       v.* = 0;
    

    const cx : f64 = 0;
    const cy : f64 = 0;
    const rx : f64 = 3;
    const ry : f64 = 1;

    const stepx = rx / f64(width / 2);
    const stepy = ry / f64(height / 2);

    var y : i64 = -height / 2;
    while (y < height / 2) {
        var fy = cy + @intToFloat(f64, y) * stepy;

        var x : i64 = -width / 2;
        while (x < width / 2) {
            var fx = cx + @intToFloat(f64, x) * stepx;
            const response = mandel(fx, fy);
            if (response > 1) {
                warn("x");
            }
            else {
                warn(" ");
            }
            
            x += 1;
        }
        warn("\n");
        y += 1;
    }
}
