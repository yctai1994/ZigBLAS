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
    var temp: T = undefined;

    switch (trans) {
        'N', 'n' => {

            // Form  x := (A)⁻¹⋅x,
            // with A(LDA ≥ n, n), x(n)

            switch (uplo) {
                'U', 'u' => {
                    if (kx == 1) {
                        var ix: usize = n - 1; // x.len - 1
                        var jx: usize = undefined;

                        while (true) : (ix -= 1) {
                            temp = x[ix]; // ix = 0
                            jx = ix + 1; // jx = 1
                            for (A[ix][jx..], x[jx..]) |A_ij, x_j| temp -= A_ij * x_j;
                            switch (diag) {
                                'U', 'u' => x[ix] = temp,
                                'N', 'n' => x[ix] = temp / A[ix][ix],
                                else => @compileError(""),
                            }
                            if (ix == 0) break;
                        }
                    } else unreachable;
                },
                'L', 'l' => {
                    if (kx == 1) {
                        for (A, x, 0..) |A_i, *x_i, ix| { // ix = i = 2
                            temp = x_i.*;
                            for (A_i[0..ix], x[0..ix]) |A_ij, x_j| temp -= A_ij * x_j;
                            switch (diag) {
                                'U', 'u' => x_i.* = temp,
                                'N', 'n' => x_i.* = temp / A_i[ix],
                                else => @compileError(""),
                            }
                        }
                    } else unreachable;
                },
                else => @compileError(""),
            }
        },
        'T', 't' => {

            // Form  x := (Aᵀ)⁻¹⋅x.

            switch (uplo) {
                'U', 'u' => {
                    if (kx == 1) {
                        //
                    } else unreachable;
                },
                'L', 'l' => {
                    if (kx == 1) {
                        //
                    } else unreachable;
                },
                else => @compileError(""),
            }
        },
    }

    _ = .{ uplo, trans, diag, A, lda, x, incx, kx };
}
