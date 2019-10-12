#ifndef KERNEL_HEADER
#define KERNEL_HEADER
// -------------------------------------------------------
// Copyright (C) 2019 by Dmitry R. Gulevich
// You may use or modify this code but not distribute it.
// -------------------------------------------------------

#define ALDAM 5 // Number of poles (0 for ALDA, not defined for Hartree)
#define EXP_ORDER 3
#define XI 2.0

cx_vec p1(vec rho) {
	double xi = XI;
	cx_double i(0.,1.);
	vec ompl = xc::omega_pl(rho);
	cx_double p = sqrt(4.-xi*xi) - i*xi;
	return p*ompl;
}

// Returns: n^(2/3) * Cm
cx_mat get_n23Coeffs(cx_mat p, vec n13) {
	xc::Fxc n23fxc = xc::get_n23fxc(n13);
	cx_mat Coeffs(p.n_rows,p.n_cols);
	Coeffs.col(0) = conj(p.col(0)) % (n23fxc.finf-n23fxc.f0) / real(p.col(0));
	for(int m=1;m<p.n_cols;++m)
		Coeffs.col(m) = 0.49*Coeffs.col(m-1); // some test values of parameters
	return Coeffs;
}

#endif // KERNEL_HEADER
