pub fn dtrsv(
    comptime T: type,
    comptime uplo: u8,
    comptime trans: u8,
    comptime diag: u8,
    n: usize,
    A: [][]T,
    lda: usize,
    x: []T,
    incx: ?usize,
) void {

    // Quick return if possible.

    if (n == 0) return;

    // Set up the start point in X if the increment is not unity. This
    // will be  ( N - 1 ) * INCX  too small for descending loops.

    const kx: usize = incx orelse 1;

    // Start the operations. In this version the elements of A are
    // accessed sequentially with one pass through A.
    var t: T = undefined;

    switch (trans) {
        'N' => {

            // Form  x := (A)⁻¹⋅x,
            // with A(LDA ≥ n, n), x(n)

            switch (uplo) {
                'U' => {
                    var i: usize = n - 1;
                    var j: usize = undefined;
                    if (kx == 1) {
                        while (true) : (i -= 1) {
                            t = x[i];
                            j = i + 1;
                            for (A[i][j..n], x[j..n]) |A_ij, x_j| t -= A_ij * x_j;
                            switch (diag) {
                                'U' => x[i] = t,
                                'N' => x[i] = t / A[i][i],
                            }
                            if (i == 0) break;
                        }
                    } else unreachable;
                },
                'L' => {
                    if (kx == 1) {
                        for (A, x, 0..) |A_i, *x_i, i| {
                            t = x_i.*;
                            for (A_i[0..i], x[0..i]) |A_ij, x_j| {
                                t -= A_ij * x_j;
                            }
                            x_i.* = t / A_i[i];
                        }
                    } else unreachable;
                },
                else => @compileError(""),
            }
        },
        'T' => {

            // Form  x := (Aᵀ)⁻¹⋅x.

            switch (uplo) {
                'U' => {
                    if (kx == 1) {
                        for (A, x, 1..) |A_i, *x_i, ip1| {
                            switch (diag) {
                                'U' => {},
                                'N' => x_i.* /= A_i[ip1 - 1],
                                else => @compileError(""),
                            }
                            for (A_i[ip1..n], x[ip1..n]) |A_ij, *x_j| {
                                x_j.* -= A_ij * x_i.*;
                            }
                        }
                    } else unreachable;
                },
                'L' => @compileError("x := (Lᵀ)⁻¹⋅x doesn't support yet."),
                else => @compileError(""),
            }
        },
    }

    _ = .{ diag, A, lda, x, incx, kx };
}
