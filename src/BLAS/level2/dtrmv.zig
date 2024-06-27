pub fn dtrmv(
    comptime T: type,
    comptime uplo: u8,
    comptime trans: u8,
    comptime diag: u8,
    n: usize,
    // A is an n by n unit, or non-unit,
    // upper or lower triangular matrix.
    A: [][]T,
    lda: usize,
    // x is an n element vector.
    x: []T,
    incx: ?usize,
) void {

    // Quick return if possible.

    if (n == 0) return;

    const nounit: bool = comptime diag == 'N';

    // Set up the start point in X if the increment is not unity. This
    // will be  ( N - 1 )*INCX  too small for descending loops.

    const kx: usize = incx orelse 1;

    // Start the operations. In this version the elements of A are
    // accessed sequentially with one pass through A.

    var temp: T = undefined;

    switch (trans) {
        'N', 'n' => {
            switch (uplo) {
                'U', 'u' => {
                    for (A, x, 0..) |A_i, *x_i, i| {
                        temp = x_i.*;
                        if (nounit) temp *= A_i[i];
                        for (A_i[i + 1 ..], x[i + 1 ..]) |A_ij, x_j| {
                            temp += A_ij * x_j;
                        }
                        x_i.* = temp;
                    }
                },
                'L', 'l' => @compileError(""),
                else => @compileError(""),
            }
        },
        'T', 't' => @compileError(""),
        else => @compileError(""),
    }

    _ = .{ kx, A, lda };

    return;
}

test "dtrmv('U', 'N', 'N')" {
    const ArrI32 = Array(i32){ .allocator = testing.allocator };

    const A: [][]i32 = try ArrI32.matrix(3, 3);
    defer ArrI32.free(A);

    inline for (.{ 2, 6, 8 }, A[0]) |v, *p| p.* = v;
    inline for (.{ 0, 1, 5 }, A[1]) |v, *p| p.* = v;
    inline for (.{ 0, 0, 3 }, A[2]) |v, *p| p.* = v;

    const x: []i32 = try ArrI32.vector(3);
    defer ArrI32.free(x);

    inline for (.{ 1, 5, 3 }, x) |v, *p| p.* = v;

    const y: [3]i32 = .{ 56, 20, 9 };

    dtrmv(i32, 'U', 'N', 'N', 3, A, 3, x, null);
    try testing.expect(mem.eql(i32, &y, x));

    // return error.SkipZigTest;
}

const std = @import("std");
const mem = std.mem;
const testing = std.testing;

const Array = @import("../../utils/array.zig").Array;
