#ifndef KERNEL_HEADER
#define KERNEL_HEADER
// -------------------------------------------------------
// Copyright (C) 2019 by Dmitry R. Gulevich
// You may use or modify this code but not distribute it.
// -------------------------------------------------------

#define ALDAM 1 // Number of poles (0 for ALDA, not defined for Hartree)
#define EXP_ORDER 3

#define XI 2.0
// XI = infinity for ALDA (i.e. short memory)
//#define XI 1.53 // QV results for sqrt version

cx_vec p1(vec rho) {
	cx_double i(0.,1.);
	return (2.-i*XI)*xc::omega_pl(rho);
}

// Returns: n^(2/3) * Cm
cx_mat get_n23Coeffs(cx_mat p, vec n13) {
	xc::Fxc n23fxc = xc::get_n23fxc(n13);
	cx_mat Coeffs(p.n_rows,p.n_cols);
	Coeffs.col(0) = conj(p.col(0)) % (n23fxc.finf-n23fxc.f0) / real(p.col(0));
	return Coeffs;
}

#endif // KERNEL_HEADER