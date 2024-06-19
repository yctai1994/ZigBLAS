pub fn dsyrk(
    comptime T: type,
    comptime uplo: u8,
    comptime trans: u8,
    n: usize,
    k: usize,
    alpha: T,
    A: [][]T,
    lda: usize,
    beta: T,
    C: [][]T,
    ldc: usize,
) void {
    const zero: T = 0;
    const one: T = 1;

    switch (trans) {
        'N', 'n' => if (lda < @max(1, n)) unreachable, // nrowa = n
        'T', 't' => if (lda < @max(1, k)) unreachable, // nrowa = k
        else => @compileError(""),
    }

    if (ldc < @max(1, n)) unreachable;

    if (n == 0) return;
    if (beta == one and (alpha == zero or k == 0)) return;

    // And when alpha == zero.

    if (alpha == zero) {
        switch (uplo) {
            'U', 'u' => {
                if (beta == zero) {
                    for (C, 0..) |C_i, jmin| {
                        for (C_i[jmin..]) |*C_ij| C_ij.* = zero;
                    }
                } else if (beta != one) {
                    for (C, 0..) |C_i, jmin| {
                        for (C_i[jmin..]) |*C_ij| C_ij.* *= beta;
                    }
                }
            },
            'L', 'l' => {
                if (beta == zero) {
                    for (C, 1..) |C_i, jmax| {
                        for (C_i[0..jmax]) |*C_ij| C_ij.* = zero;
                    }
                } else if (beta != one) {
                    for (C, 1..) |C_i, jmax| {
                        for (C_i[0..jmax]) |*C_ij| C_ij.* *= beta;
                    }
                }
            },
            else => @compileError(""),
        }
        return;
    }

    // Start the operations.
    var temp: T = undefined;

    switch (trans) {
        'N', 'n' => {
            switch (uplo) {

                // Form  C = alpha⋅A⋅Aᵀ + beta⋅C.

                'U', 'u' => {
                    for (A, C, 0..) |A_i, C_i, jmin| {
                        for (A[jmin..], C_i[jmin..]) |A_j, *C_ij| {
                            // temp = A[i,:] ⋅ Aᵀ[:,j]
                            temp = zero;
                            for (A_i, A_j) |A_ik, A_jk| temp += A_ik * A_jk;

                            if (beta == zero) {
                                C_ij.* = alpha * temp;
                            } else if (beta == one) {
                                C_ij.* += alpha * temp;
                            } else {
                                C_ij.* = alpha * temp + beta * C_ij.*;
                            }
                        }
                    }
                },
                'L', 'l' => {
                    for (A, C, 1..) |A_i, C_i, jmax| { // i = 0..
                        for (A[0..jmax], C_i[0..jmax]) |A_j, *C_ij| {
                            // temp = A[i,:] ⋅ Aᵀ[:,j]
                            temp = zero;
                            for (A_i, A_j) |A_ik, A_jk| temp += A_ik * A_jk;

                            if (beta == zero) {
                                C_ij.* = alpha * temp;
                            } else if (beta == one) {
                                C_ij.* += alpha * temp;
                            } else {
                                C_ij.* = alpha * temp + beta * C_ij.*;
                            }
                        }
                    }
                },
                else => @compileError(""),
            }
        },
        'T', 't' => {

            // Form  C := alpha⋅Aᵀ⋅A + beta⋅C.

            switch (uplo) {
                'U', 'u' => {
                    for (C, 0..) |C_i, jmin| {
                        if (beta == zero) {
                            for (C_i[jmin..]) |*C_ij| C_ij.* = zero;
                        } else if (beta != one) {
                            for (C_i[jmin..]) |*C_ij| C_ij.* *= beta;
                        }

                        for (A) |A_k| {
                            temp = alpha * A_k[jmin];
                            for (A_k[jmin..], C_i[jmin..]) |A_kj, *C_ij| {
                                C_ij.* += temp * A_kj;
                            }
                        }
                    }
                },
                'L', 'l' => {
                    for (C, 1..) |C_i, jmax| {
                        if (beta == zero) {
                            for (C_i[0..jmax]) |*C_ij| C_ij.* = zero;
                        } else if (beta != one) {
                            for (C_i[0..jmax]) |*C_ij| C_ij.* *= beta;
                        }

                        for (A) |A_k| {
                            temp = alpha * A_k[jmax - 1];
                            for (A_k[0..jmax], C_i[0..jmax]) |A_kj, *C_ij| {
                                C_ij.* += temp * A_kj;
                            }
                        }
                    }
                },
                else => @compileError(""),
            }
        },
        else => @compileError(""),
    }

    return;
}

test "dsyrk: alpha == 0" {
    const ArrF32 = Array(f32){ .allocator = testing.allocator };

    const A: [][]f32 = try ArrF32.matrix(4, 3);
    defer ArrF32.free(A);

    const C: [][]f32 = try ArrF32.matrix(4, 4);
    defer ArrF32.free(C);

    const alpha: comptime_int = 0;
    {
        const Us: [3][4][4]f32 = .{
            .{ .{ 0, 0, 0, 0 }, .{ 2, 0, 0, 0 }, .{ 1, 0, 0, 0 }, .{ 0, 3, 2, 0 } }, // beta = 0
            .{ .{ 3, 2, 1, 0 }, .{ 2, 1, 0, 3 }, .{ 1, 0, 3, 2 }, .{ 0, 3, 2, 1 } }, // beta = 1
            .{ .{ 6, 4, 2, 0 }, .{ 2, 2, 0, 6 }, .{ 1, 0, 6, 4 }, .{ 0, 3, 2, 2 } }, // beta = 2
        };

        inline for (.{ 0, 1, 2 }, Us) |beta, U| {
            inline for (.{ 3, 2, 1, 0 }, C[0]) |v, *p| p.* = v;
            inline for (.{ 2, 1, 0, 3 }, C[1]) |v, *p| p.* = v;
            inline for (.{ 1, 0, 3, 2 }, C[2]) |v, *p| p.* = v;
            inline for (.{ 0, 3, 2, 1 }, C[3]) |v, *p| p.* = v;
            dsyrk(f32, 'U', 'N', 4, 3, alpha, A, 4, beta, C, 4);
            for (U, C) |rowU, rowC| try testing.expect(mem.eql(f32, &rowU, rowC));
        }
    }
    {
        const Ls: [3][4][4]f32 = .{
            .{ .{ 0, 2, 1, 0 }, .{ 0, 0, 0, 3 }, .{ 0, 0, 0, 2 }, .{ 0, 0, 0, 0 } }, // beta = 0
            .{ .{ 3, 2, 1, 0 }, .{ 2, 1, 0, 3 }, .{ 1, 0, 3, 2 }, .{ 0, 3, 2, 1 } }, // beta = 1
            .{ .{ 6, 2, 1, 0 }, .{ 4, 2, 0, 3 }, .{ 2, 0, 6, 2 }, .{ 0, 6, 4, 2 } }, // beta = 2
        };

        inline for (.{ 0, 1, 2 }, Ls) |beta, L| {
            inline for (.{ 3, 2, 1, 0 }, C[0]) |v, *p| p.* = v;
            inline for (.{ 2, 1, 0, 3 }, C[1]) |v, *p| p.* = v;
            inline for (.{ 1, 0, 3, 2 }, C[2]) |v, *p| p.* = v;
            inline for (.{ 0, 3, 2, 1 }, C[3]) |v, *p| p.* = v;
            dsyrk(f32, 'L', 'N', 4, 3, 0, A, 4, beta, C, 4);
            for (L, C) |rowL, rowC| try testing.expect(mem.eql(f32, &rowL, rowC));
        }
    }
    // return error.SkipZigTest;
}

test "dsyrk('U' / 'L', 'N'), alpha = 3" {
    const ArrI32 = Array(i32){ .allocator = testing.allocator };

    const A: [][]i32 = try ArrI32.matrix(4, 3);
    defer ArrI32.free(A);

    inline for (.{ 1, 5, 3 }, A[0]) |v, *p| p.* = v;
    inline for (.{ 3, 0, 2 }, A[1]) |v, *p| p.* = v;
    inline for (.{ 2, 1, 0 }, A[2]) |v, *p| p.* = v;
    inline for (.{ 0, 3, 1 }, A[3]) |v, *p| p.* = v;

    const C: [][]i32 = try ArrI32.matrix(4, 4);
    defer ArrI32.free(C);

    const alpha: comptime_int = 3;
    {
        const Us: [3][4][4]i32 = .{
            .{ .{ 105, 27, 21, 54 }, .{ 2, 39, 18, 6 }, .{ 1, 0, 15, 9 }, .{ 0, 3, 2, 30 } }, // beta = 0
            .{ .{ 108, 29, 22, 54 }, .{ 2, 40, 18, 9 }, .{ 1, 0, 18, 11 }, .{ 0, 3, 2, 31 } }, // beta = 1
            .{ .{ 111, 31, 23, 54 }, .{ 2, 41, 18, 12 }, .{ 1, 0, 21, 13 }, .{ 0, 3, 2, 32 } }, // beta = 2
        };

        inline for (.{ 0, 1, 2 }, Us) |beta, U| {
            inline for (.{ 3, 2, 1, 0 }, C[0]) |v, *p| p.* = v;
            inline for (.{ 2, 1, 0, 3 }, C[1]) |v, *p| p.* = v;
            inline for (.{ 1, 0, 3, 2 }, C[2]) |v, *p| p.* = v;
            inline for (.{ 0, 3, 2, 1 }, C[3]) |v, *p| p.* = v;
            dsyrk(i32, 'U', 'N', 4, 3, alpha, A, 4, beta, C, 4);
            for (U, C) |rowU, rowC| try testing.expect(mem.eql(i32, &rowU, rowC));
        }
    }
    {
        const Ls: [3][4][4]i32 = .{
            .{ .{ 105, 2, 1, 0 }, .{ 27, 39, 0, 3 }, .{ 21, 18, 15, 2 }, .{ 54, 6, 9, 30 } }, // beta = 0
            .{ .{ 108, 2, 1, 0 }, .{ 29, 40, 0, 3 }, .{ 22, 18, 18, 2 }, .{ 54, 9, 11, 31 } }, // beta = 1
            .{ .{ 111, 2, 1, 0 }, .{ 31, 41, 0, 3 }, .{ 23, 18, 21, 2 }, .{ 54, 12, 13, 32 } }, // beta = 2
        };

        inline for (.{ 0, 1, 2 }, Ls) |beta, L| {
            inline for (.{ 3, 2, 1, 0 }, C[0]) |v, *p| p.* = v;
            inline for (.{ 2, 1, 0, 3 }, C[1]) |v, *p| p.* = v;
            inline for (.{ 1, 0, 3, 2 }, C[2]) |v, *p| p.* = v;
            inline for (.{ 0, 3, 2, 1 }, C[3]) |v, *p| p.* = v;
            dsyrk(i32, 'L', 'N', 4, 3, alpha, A, 4, beta, C, 4);
            for (L, C) |rowL, rowC| try testing.expect(mem.eql(i32, &rowL, rowC));
        }
    }
    return error.SkipZigTest;
}

test "dsyrk('U' / 'L', 'T'), alpha = 3" {
    const ArrI32 = Array(i32){ .allocator = testing.allocator };

    const A: [][]i32 = try ArrI32.matrix(4, 3);
    defer ArrI32.free(A);

    inline for (.{ 1, 5, 3 }, A[0]) |v, *p| p.* = v;
    inline for (.{ 3, 0, 2 }, A[1]) |v, *p| p.* = v;
    inline for (.{ 2, 1, 0 }, A[2]) |v, *p| p.* = v;
    inline for (.{ 0, 3, 1 }, A[3]) |v, *p| p.* = v;

    const C: [][]i32 = try ArrI32.matrix(3, 3);
    defer ArrI32.free(C);

    const alpha: comptime_int = 3;
    {
        const Us: [3][3][3]i32 = .{
            .{ .{ 42, 21, 27 }, .{ 1, 105, 54 }, .{ 0, 2, 42 } }, // beta = 0
            .{ .{ 44, 22, 27 }, .{ 1, 105, 56 }, .{ 0, 2, 43 } }, // beta = 1
            .{ .{ 46, 23, 27 }, .{ 1, 105, 58 }, .{ 0, 2, 44 } }, // beta = 2
        };

        inline for (.{ 0, 1, 2 }, Us) |beta, U| {
            inline for (.{ 2, 1, 0 }, C[0]) |v, *p| p.* = v;
            inline for (.{ 1, 0, 2 }, C[1]) |v, *p| p.* = v;
            inline for (.{ 0, 2, 1 }, C[2]) |v, *p| p.* = v;
            dsyrk(i32, 'U', 'T', 3, 3, alpha, A, 4, beta, C, 3);
            for (U, C) |rowU, rowC| try testing.expect(mem.eql(i32, &rowU, rowC));
        }
    }
    {
        const Ls: [3][3][3]i32 = .{
            .{ .{ 42, 1, 0 }, .{ 21, 105, 2 }, .{ 27, 54, 42 } }, // beta = 0
            .{ .{ 44, 1, 0 }, .{ 22, 105, 2 }, .{ 27, 56, 43 } }, // beta = 1
            .{ .{ 46, 1, 0 }, .{ 23, 105, 2 }, .{ 27, 58, 44 } }, // beta = 2
        };

        inline for (.{ 0, 1, 2 }, Ls) |beta, L| {
            inline for (.{ 2, 1, 0 }, C[0]) |v, *p| p.* = v;
            inline for (.{ 1, 0, 2 }, C[1]) |v, *p| p.* = v;
            inline for (.{ 0, 2, 1 }, C[2]) |v, *p| p.* = v;
            dsyrk(i32, 'L', 'T', 3, 3, alpha, A, 4, beta, C, 3);
            for (L, C) |rowL, rowC| try testing.expect(mem.eql(i32, &rowL, rowC));
        }
    }
    return error.SkipZigTest;
}

const std = @import("std");
const mem = std.mem;
const testing = std.testing;

const Array = @import("../../utils/array.zig").Array;
