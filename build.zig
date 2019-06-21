const std = @import("std");
const Builder = std.build.Builder;

pub fn build(b: *Builder) void {
    const mode = b.standardReleaseOptions();

    const exe = b.addExecutable("zigmandel", "src/mandel.zig");
    exe.setBuildMode(mode);

    exe.addLibPath("/usr/lib64");
    exe.linkSystemLibrary("c");
    exe.linkSystemLibrary("webp");

    const run_cmd = exe.run();

    const run_step = b.step("run", "Run the app");
    run_step.dependOn(&run_cmd.step);

    b.default_step.dependOn(&exe.step);
    b.installArtifact(exe);
}
