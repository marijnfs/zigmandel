all:
	zig build

release:
	zig build --release=fast
