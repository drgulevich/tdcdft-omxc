#ifndef KERNEL_HEADER
#define KERNEL_HEADER
// -------------------------------------------------------
// Copyright (C) 2019 by Dmitry R. Gulevich
// You may use or modify this code but not distribute it.
// -------------------------------------------------------

#define ALDAM 0 // Number of poles (0 for ALDA, not defined for Hartree)

// void function
cx_vec p1(vec rho) {
	cx_vec result = zeros<cx_vec>(rho.n_elem);
	return result;
}

// void function
cx_mat get_n23Coeffs(cx_mat p, vec n13) {
	cx_mat result = zeros<cx_mat>(p.n_rows,p.n_cols);
	return result;
}

#endif // KERNEL_HEADER
