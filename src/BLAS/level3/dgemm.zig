pub fn dgemm(
    comptime T: type,
    comptime transa: u8,
    comptime transb: u8,
    m: usize,
    n: usize,
    k: usize,
    alpha: T,
    A: [][]T,
    lda: usize,
    B: [][]T,
    ldb: usize,
    beta: T,
    C: [][]T,
    ldc: usize,
) void {
    const zero: T = 0;
    const one: T = 1;

    switch (transa) {
        'N', 'n' => if (lda < @max(1, m)) unreachable,
        'T', 't' => if (lda < @max(1, k)) unreachable,
        else => @compileError(""),
    }

    switch (transb) {
        'N', 'n' => if (ldb < @max(1, k)) unreachable,
        'T', 't' => if (ldb < @max(1, n)) unreachable,
        else => @compileError(""),
    }

    if (ldc < @max(1, m)) unreachable;

    if (m == 0 or n == 0) return;
    if (beta == one and (alpha == zero or k == 0)) return;

    if (alpha == zero) {
        if (beta == zero) {
            for (C) |C_i| {
                for (C_i) |*C_ij| C_ij.* = zero;
            }
        } else if (beta != one) {
            for (C) |C_i| {
                for (C_i) |*C_ij| C_ij.* *= beta;
            }
        }
        return;
    }

    // Start the operations.
    var temp: T = undefined;

    switch (transb) {
        'N', 'n' => {
            switch (transa) {

                // Form  C := alpha⋅A⋅B + beta⋅C

                'N', 'n' => {
                    for (A, C) |A_i, C_i| {
                        if (beta == zero) {
                            for (C_i) |*C_ij| C_ij.* = zero;
                        } else if (beta != one) {
                            for (C_i) |*C_ij| C_ij.* *= beta;
                        }
                        for (A_i, B) |A_ik, B_k| {
                            if (A_ik != zero) {
                                temp = alpha * A_ik;
                                for (B_k, C_i) |B_kj, *C_ij| C_ij.* += temp * B_kj;
                            }
                        }
                    }
                },

                // Form  C := alpha⋅Aᵀ⋅B + beta⋅C

                'T', 't' => {
                    for (C, 0..) |C_i, i| {
                        if (beta == zero) {
                            for (C_i) |*C_ij| C_ij.* = zero;
                        } else if (beta != one) {
                            for (C_i) |*C_ij| C_ij.* *= beta;
                        }
                        for (A, B) |A_k, B_k| {
                            if (A_k[i] != zero) {
                                temp = alpha * A_k[i];
                                for (B_k, C_i) |B_kj, *C_ij| C_ij.* += temp * B_kj;
                            }
                        }
                    }
                },
                else => @compileError(""),
            }
        },
        'T', 't' => {
            switch (transa) {

                // Form  C := alpha⋅A⋅Bᵀ + beta⋅C

                'N', 'n' => {
                    for (A, C) |A_i, C_i| {
                        for (B, C_i) |B_j, *C_ij| {
                            temp = zero;
                            for (A_i, B_j) |A_ik, B_jk| {
                                temp += A_ik * B_jk;
                            }
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

                // Form  C := alpha⋅Aᵀ⋅Bᵀ + beta⋅C

                'T', 't' => {
                    for (C, 0..) |C_i, i| {
                        for (B, C_i) |B_j, *C_ij| {
                            temp = zero;
                            for (A, B_j) |A_k, B_jk| {
                                temp += A_k[i] * B_jk;
                            }
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
        else => @compileError(""),
    }

    return;
}

test "dgemm: alpha == 0" {
    const ArrI32 = Array(i32){ .allocator = testing.allocator };

    const A: [][]i32 = try ArrI32.matrix(4, 3);
    defer ArrI32.free(A);

    const B: [][]i32 = try ArrI32.matrix(3, 4);
    defer ArrI32.free(B);

    const C: [][]i32 = try ArrI32.matrix(4, 4);
    defer ArrI32.free(C);

    const alpha: comptime_int = 0;
    {
        const Rs: [3][4][4]i32 = .{ // references
            .{ .{ 0, 0, 0, 0 }, .{ 0, 0, 0, 0 }, .{ 0, 0, 0, 0 }, .{ 0, 0, 0, 0 } }, // beta = 0
            .{ .{ 3, 2, 1, 0 }, .{ 2, 1, 0, 3 }, .{ 1, 0, 3, 2 }, .{ 0, 3, 2, 1 } }, // beta = 1
            .{ .{ 6, 4, 2, 0 }, .{ 4, 2, 0, 6 }, .{ 2, 0, 6, 4 }, .{ 0, 6, 4, 2 } }, // beta = 2
        };

        inline for (.{ 0, 1, 2 }, Rs) |beta, R| {
            inline for (.{ 3, 2, 1, 0 }, C[0]) |v, *p| p.* = v;
            inline for (.{ 2, 1, 0, 3 }, C[1]) |v, *p| p.* = v;
            inline for (.{ 1, 0, 3, 2 }, C[2]) |v, *p| p.* = v;
            inline for (.{ 0, 3, 2, 1 }, C[3]) |v, *p| p.* = v;
            dgemm(i32, 'N', 'N', 4, 4, 3, alpha, A, 4, B, 3, beta, C, 4);
            for (R, C) |rowR, rowC| try testing.expect(mem.eql(i32, &rowR, rowC));
        }
    }
}

test "dgemm('N', 'N'), alpha = 3" {
    const ArrI32 = Array(i32){ .allocator = testing.allocator };

    const A: [][]i32 = try ArrI32.matrix(4, 3);
    defer ArrI32.free(A);

    inline for (.{ 1, 5, 3 }, A[0]) |v, *p| p.* = v;
    inline for (.{ 3, 0, 2 }, A[1]) |v, *p| p.* = v;
    inline for (.{ 2, 1, 0 }, A[2]) |v, *p| p.* = v;
    inline for (.{ 0, 3, 1 }, A[3]) |v, *p| p.* = v;

    const B: [][]i32 = try ArrI32.matrix(3, 4);
    defer ArrI32.free(B);

    inline for (.{ 2, 1, 0, 3 }, B[0]) |v, *p| p.* = v;
    inline for (.{ 1, 5, 3, 0 }, B[1]) |v, *p| p.* = v;
    inline for (.{ 0, 3, 1, 2 }, B[2]) |v, *p| p.* = v;

    const C: [][]i32 = try ArrI32.matrix(4, 4);
    defer ArrI32.free(C);

    const alpha: comptime_int = 3;
    {
        const Rs: [3][4][4]i32 = .{
            .{ .{ 21, 105, 54, 27 }, .{ 18, 27, 6, 39 }, .{ 15, 21, 9, 18 }, .{ 9, 54, 30, 6 } }, // beta = 0
            .{ .{ 24, 107, 55, 27 }, .{ 20, 28, 6, 42 }, .{ 16, 21, 12, 20 }, .{ 9, 57, 32, 7 } }, // beta = 1
            .{ .{ 27, 109, 56, 27 }, .{ 22, 29, 6, 45 }, .{ 17, 21, 15, 22 }, .{ 9, 60, 34, 8 } }, // beta = 2
        };

        inline for (.{ 0, 1, 2 }, Rs) |beta, R| {
            inline for (.{ 3, 2, 1, 0 }, C[0]) |v, *p| p.* = v;
            inline for (.{ 2, 1, 0, 3 }, C[1]) |v, *p| p.* = v;
            inline for (.{ 1, 0, 3, 2 }, C[2]) |v, *p| p.* = v;
            inline for (.{ 0, 3, 2, 1 }, C[3]) |v, *p| p.* = v;
            dgemm(i32, 'N', 'N', 4, 4, 3, alpha, A, 4, B, 3, beta, C, 4);
            for (R, C) |rowR, rowC| try testing.expect(mem.eql(i32, &rowR, rowC));
        }
    }
}

test "dgemm('T', 'N'), alpha = 3" {
    const ArrI32 = Array(i32){ .allocator = testing.allocator };

    const A: [][]i32 = try ArrI32.matrix(3, 4);
    defer ArrI32.free(A);

    inline for (.{ 1, 3, 2, 0 }, A[0]) |v, *p| p.* = v;
    inline for (.{ 5, 0, 1, 3 }, A[1]) |v, *p| p.* = v;
    inline for (.{ 3, 2, 0, 1 }, A[2]) |v, *p| p.* = v;

    const B: [][]i32 = try ArrI32.matrix(3, 4);
    defer ArrI32.free(B);

    inline for (.{ 2, 1, 0, 3 }, B[0]) |v, *p| p.* = v;
    inline for (.{ 1, 5, 3, 0 }, B[1]) |v, *p| p.* = v;
    inline for (.{ 0, 3, 1, 2 }, B[2]) |v, *p| p.* = v;

    const C: [][]i32 = try ArrI32.matrix(4, 4);
    defer ArrI32.free(C);

    const alpha: comptime_int = 3;
    {
        const Rs: [3][4][4]i32 = .{
            .{ .{ 21, 105, 54, 27 }, .{ 18, 27, 6, 39 }, .{ 15, 21, 9, 18 }, .{ 9, 54, 30, 6 } }, // beta = 0
            .{ .{ 24, 107, 55, 27 }, .{ 20, 28, 6, 42 }, .{ 16, 21, 12, 20 }, .{ 9, 57, 32, 7 } }, // beta = 1
            .{ .{ 27, 109, 56, 27 }, .{ 22, 29, 6, 45 }, .{ 17, 21, 15, 22 }, .{ 9, 60, 34, 8 } }, // beta = 2
        };

        inline for (.{ 0, 1, 2 }, Rs) |beta, R| {
            inline for (.{ 3, 2, 1, 0 }, C[0]) |v, *p| p.* = v;
            inline for (.{ 2, 1, 0, 3 }, C[1]) |v, *p| p.* = v;
            inline for (.{ 1, 0, 3, 2 }, C[2]) |v, *p| p.* = v;
            inline for (.{ 0, 3, 2, 1 }, C[3]) |v, *p| p.* = v;
            dgemm(i32, 'T', 'N', 4, 4, 3, alpha, A, 3, B, 3, beta, C, 4);
            for (R, C) |rowR, rowC| try testing.expect(mem.eql(i32, &rowR, rowC));
        }
    }
}

test "dgemm('N', 'T'), alpha = 3" {
    const ArrI32 = Array(i32){ .allocator = testing.allocator };

    const A: [][]i32 = try ArrI32.matrix(4, 3);
    defer ArrI32.free(A);

    inline for (.{ 1, 5, 3 }, A[0]) |v, *p| p.* = v;
    inline for (.{ 3, 0, 2 }, A[1]) |v, *p| p.* = v;
    inline for (.{ 2, 1, 0 }, A[2]) |v, *p| p.* = v;
    inline for (.{ 0, 3, 1 }, A[3]) |v, *p| p.* = v;

    const B: [][]i32 = try ArrI32.matrix(4, 3);
    defer ArrI32.free(B);

    inline for (.{ 2, 1, 0 }, B[0]) |v, *p| p.* = v;
    inline for (.{ 1, 5, 3 }, B[1]) |v, *p| p.* = v;
    inline for (.{ 0, 3, 1 }, B[2]) |v, *p| p.* = v;
    inline for (.{ 3, 0, 2 }, B[3]) |v, *p| p.* = v;

    const C: [][]i32 = try ArrI32.matrix(4, 4);
    defer ArrI32.free(C);

    const alpha: comptime_int = 3;
    {
        const Rs: [3][4][4]i32 = .{
            .{ .{ 21, 105, 54, 27 }, .{ 18, 27, 6, 39 }, .{ 15, 21, 9, 18 }, .{ 9, 54, 30, 6 } }, // beta = 0
            .{ .{ 24, 107, 55, 27 }, .{ 20, 28, 6, 42 }, .{ 16, 21, 12, 20 }, .{ 9, 57, 32, 7 } }, // beta = 1
            .{ .{ 27, 109, 56, 27 }, .{ 22, 29, 6, 45 }, .{ 17, 21, 15, 22 }, .{ 9, 60, 34, 8 } }, // beta = 2
        };

        inline for (.{ 0, 1, 2 }, Rs) |beta, R| {
            inline for (.{ 3, 2, 1, 0 }, C[0]) |v, *p| p.* = v;
            inline for (.{ 2, 1, 0, 3 }, C[1]) |v, *p| p.* = v;
            inline for (.{ 1, 0, 3, 2 }, C[2]) |v, *p| p.* = v;
            inline for (.{ 0, 3, 2, 1 }, C[3]) |v, *p| p.* = v;
            dgemm(i32, 'N', 'T', 4, 4, 3, alpha, A, 4, B, 4, beta, C, 4);
            for (R, C) |rowR, rowC| try testing.expect(mem.eql(i32, &rowR, rowC));
        }
    }
}

test "dgemm('T', 'T'), alpha = 3" {
    const ArrI32 = Array(i32){ .allocator = testing.allocator };

    const A: [][]i32 = try ArrI32.matrix(3, 4);
    defer ArrI32.free(A);

    inline for (.{ 1, 3, 2, 0 }, A[0]) |v, *p| p.* = v;
    inline for (.{ 5, 0, 1, 3 }, A[1]) |v, *p| p.* = v;
    inline for (.{ 3, 2, 0, 1 }, A[2]) |v, *p| p.* = v;

    const B: [][]i32 = try ArrI32.matrix(4, 3);
    defer ArrI32.free(B);

    inline for (.{ 2, 1, 0 }, B[0]) |v, *p| p.* = v;
    inline for (.{ 1, 5, 3 }, B[1]) |v, *p| p.* = v;
    inline for (.{ 0, 3, 1 }, B[2]) |v, *p| p.* = v;
    inline for (.{ 3, 0, 2 }, B[3]) |v, *p| p.* = v;

    const C: [][]i32 = try ArrI32.matrix(4, 4);
    defer ArrI32.free(C);

    const alpha: comptime_int = 3;
    {
        const Rs: [3][4][4]i32 = .{
            .{ .{ 21, 105, 54, 27 }, .{ 18, 27, 6, 39 }, .{ 15, 21, 9, 18 }, .{ 9, 54, 30, 6 } }, // beta = 0
            .{ .{ 24, 107, 55, 27 }, .{ 20, 28, 6, 42 }, .{ 16, 21, 12, 20 }, .{ 9, 57, 32, 7 } }, // beta = 1
            .{ .{ 27, 109, 56, 27 }, .{ 22, 29, 6, 45 }, .{ 17, 21, 15, 22 }, .{ 9, 60, 34, 8 } }, // beta = 2
        };

        inline for (.{ 0, 1, 2 }, Rs) |beta, R| {
            inline for (.{ 3, 2, 1, 0 }, C[0]) |v, *p| p.* = v;
            inline for (.{ 2, 1, 0, 3 }, C[1]) |v, *p| p.* = v;
            inline for (.{ 1, 0, 3, 2 }, C[2]) |v, *p| p.* = v;
            inline for (.{ 0, 3, 2, 1 }, C[3]) |v, *p| p.* = v;
            dgemm(i32, 'T', 'T', 4, 4, 3, alpha, A, 3, B, 4, beta, C, 4);
            for (R, C) |rowR, rowC| try testing.expect(mem.eql(i32, &rowR, rowC));
        }
    }
}

const std = @import("std");
const mem = std.mem;
const testing = std.testing;

const Array = @import("../../utils/array.zig").Array;
