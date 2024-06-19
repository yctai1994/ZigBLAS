pub fn dtrsm(
    comptime T: type,
    comptime side: u8,
    comptime uplo: u8,
    comptime transa: u8,
    comptime diag: u8,
    m: usize,
    n: usize,
    alpha: T,
    A: [][]T,
    lda: usize,
    B: [][]T,
    ldb: usize,
) void {
    const zero: T = 0;
    const one: T = 1;

    const lside: bool = side == 'L';
    const upper: bool = uplo == 'U';
    const nounit: bool = diag == 'N';

    switch (side) {
        'L', 'l' => if (lda < @max(1, m)) unreachable,
        'R', 'r' => if (lda < @max(1, n)) unreachable,
        else => @compileError(""),
    }

    switch (uplo) {
        'U', 'u' => {},
        'L', 'l' => {},
        else => @compileError(""),
    }

    switch (transa) {
        'N', 'n' => {},
        'T', 't' => {},
        'C', 'c' => {},
        else => @compileError(""),
    }

    switch (diag) {
        'U', 'u' => {},
        'N', 'n' => {},
        else => @compileError(""),
    }

    if (m < 0) unreachable;
    if (n < 0) unreachable;
    if (ldb < @max(1, m)) unreachable;

    // Quick return if possible.

    if (m == 0 or n == 0) return;

    // And when  alpha == zero.

    if (alpha == zero) {
        for (B) |B_i| {
            for (B_i) |*B_ij| B_ij.* = zero;
        }
        return;
    }

    // Start the operations.
    var temp: T = undefined;

    if (lside) {
        if (transa == 'N') {

            // Form  B := alpha*inv( A )*B.

            if (upper) {
                for (1..n) |j| {
                    if (alpha != one) {
                        for (1..m) |i| {
                            B[i][j] = alpha * B[i][j];
                        }
                    }
                    for (m..1) |k| {
                        if (B[k][j] != zero) {
                            if (nounit) B[k][j] = B[k][j] / A[k][k];
                            for (1..k - 1) |i| {
                                B[i][j] = B[i][j] - B[k][j] * A[i][k];
                            }
                        }
                    }
                }
            } else {
                for (1..n) |j| {
                    if (alpha != one) {
                        for (1..m) |i| {
                            B[i][j] = alpha * B[i][j];
                        }
                    }
                    for (1..m) |k| {
                        if (B[k][j] != zero) {
                            if (nounit) B[k][j] = B[k][j] / A[k][k];
                            for (k + 1..m) |i| {
                                B[i][j] = B[i][j] - B[k][j] * A[i][k];
                            }
                        }
                    }
                }
            }
        } else {

            // Form  B := alpha*inv( A**T )*B.

            if (upper) {
                for (1..n) |j| {
                    for (1..m) |i| {
                        temp = alpha * B[i][j];
                        for (1..i - 1) |k| {
                            temp = temp - A[k][i] * B[k][j];
                        }
                        if (nounit) temp = temp / A[i][i];
                        B[i][j] = temp;
                    }
                }
            } else {
                for (1..n) |j| {
                    for (m..1) |i| {
                        temp = alpha * B[i][j];
                        for (i + 1..m) |k| {
                            temp = temp - A[k][i] * B[k][j];
                        }
                        if (nounit) temp = temp / A[i][i];
                        B[i][j] = temp;
                    }
                }
            }
        }
    } else {
        if (transa == 'N') {

            // Form  B := alpha*B*inv( A ).

            if (upper) {
                for (1..n) |j| {
                    if (alpha != one) {
                        for (1..m) |i| {
                            B[i][j] = alpha * B[i][j];
                        }
                    }
                    for (1..j - 1) |k| {
                        if (A[k][j] != zero) {
                            for (1..m) |i| {
                                B[i][j] = B[i][j] - A[k][j] * B[i][k];
                            }
                        }
                    }
                    if (nounit) {
                        temp = one / A[j][j];
                        for (1..m) |i| {
                            B[i][j] = temp * B[i][j];
                        }
                    }
                }
            } else {
                for (n..1) |j| {
                    if (alpha != one) {
                        for (1..m) |i| {
                            B[i][j] = alpha * B[i][j];
                        }
                    }
                    for (j + 1..n) |k| {
                        if (A[k][j] != zero) {
                            for (1..m) |i| {
                                B[i][j] = B[i][j] - A[k][j] * B[i][k];
                            }
                        }
                    }
                    if (nounit) {
                        temp = one / A[j][j];
                        for (1..m) |i| {
                            B[i][j] = temp * B[i][j];
                        }
                    }
                }
            }
        } else {

            // Form  B := alpha*B*inv( A**T ).

            if (upper) {
                for (n..1) |k| {
                    if (nounit) {
                        temp = one / A[k][k];
                        for (1..m) |i| {
                            B[i][k] = temp * B[i][k];
                        }
                    }
                    for (1..k - 1) |j| {
                        if (A[j][k] != zero) {
                            temp = A[j][k];
                            for (1..m) |i| {
                                B[i][j] = B[i][j] - temp * B[i][k];
                            }
                        }
                    }
                    if (alpha != one) {
                        for (1..m) |i| {
                            B[i][k] = alpha * B[i][k];
                        }
                    }
                }
            } else {
                for (1..n) |k| {
                    if (nounit) {
                        temp = one / A[k][k];
                        for (1..m) |i| {
                            B[i][k] = temp * B[i][k];
                        }
                    }
                    for (k + 1..n) |j| {
                        if (A[j][k] != zero) {
                            temp = A[j][k];
                            for (1..m) |i| {
                                B[i][j] = B[i][j] - temp * B[i][k];
                            }
                        }
                    }
                    if (alpha != one) {
                        for (1..m) |i| {
                            B[i][k] = alpha * B[i][k];
                        }
                    }
                }
            }
        }
    }
    return;
}
