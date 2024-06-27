const std = @import("std");
const testing = std.testing;

pub const dtrmv = @import("./BLAS/level2/dtrmv.zig").dtrmv;

pub const dgemm = @import("./BLAS/level3/dgemm.zig").dgemm;
pub const dsyrk = @import("./BLAS/level3/dsyrk.zig").dsyrk;

test {
    testing.refAllDecls(@This());
}
