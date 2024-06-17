const std = @import("std");
const testing = std.testing;

pub const dsyrk = @import("./BLAS/level3/dsyrk.zig").dsyrk;

test {
    testing.refAllDecls(@This());
}
