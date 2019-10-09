#ifndef KERNEL_HEADER
#define KERNEL_HEADER
// -------------------------------------------------------
// Copyright (C) 2019 by Dmitry R. Gulevich
// You may use or modify this code but not distribute it.
// -------------------------------------------------------

//#define IMDFXC_ON // ImDfxc from QV

#include "../xc.hpp"

#define ALDAM 10 // Number of poles (0 for ALDA, not defined for Hartree)
#define EXP_ORDER 3

const double omega1 = 2.0;
const double gamma1 = 1.0; // less than omega1

// TO DO: assertion gamma1 < omega1

const double imp1 = -gamma1;
const double rep1 = sqrt(omega1*omega1-gamma1*gamma1);
cx_double i(0.,1.);

vec omega_pl(vec rho) {
	static const double factor = sqrt(4.*M_PI);
	return factor*sqrt(rho);
}

cx_vec p1(vec rho) {
	vec ompl = omega_pl(rho);
	return (rep1+i*imp1)*ompl;
}


// QV Eq.(16)
vec S3L(vec lam) {
	vec result(lam.n_elem);
	double factor = -1./(45.*M_PI);
	for(uword m=0;m<lam.n_elem;++m) {
		double part1 = 5.-(lam[m]+5./lam[m])*atan(lam[m]);
		double part2 = -(2./lam[m])*asin(lam[m]/sqrt(1.+lam[m]*lam[m]));
		double arg = 1./(lam[m]*sqrt(2.+lam[m]*lam[m]));
		double part3 = arg*(M_PI - 2.*atan(arg));
		result[m] = factor*(part1+part2+part3);
	}
	return result;
}

const double factorkf = cbrt(3.*M_PI*M_PI); // kf = factorkf*n13;
    
// Returns: n^(2/3) * (d/domega)Imfxc(n,0)
vec get_n23ImDfxc(vec rho, vec n13) {
	vec kf = factorkf*n13;
	vec lam = sqrt(M_PI*kf); // 2kf/ks    
	return -factorkf*S3L(lam)/(M_PI*M_PI*rho);
}

// Returns: n^(2/3) * Cm
cx_mat get_n23Coeffs(cx_mat p, vec n13) {
	xc::Fxc n23fxc = xc::get_n23fxc(n13);
	cx_mat Coeffs(p.n_rows,p.n_cols);
	#ifdef IMDFXC_ON
		vec rho = n13 % n13 % n13; // inefficient: testing
		vec n23ImDfxc = get_n23ImDfxc(rho, n13);
		Coeffs.col(0) = p.col(0) % (n23fxc.finf - n23fxc.f0 - i*conj(p.col(0))%n23ImDfxc ) / real(p.col(0));
		//vec ompl = real( p.col(0) % conj(p.col(0)) );
		//Coeffs.col(0) = ( (1.+i*imp1/rep1)*(n23fxc.finf-n23fxc.f0) - i*(omega1*omega1/rep1)*n23ImDfxc ) % ompl;
	#else
		Coeffs.col(0) = p.col(0) % (n23fxc.finf - n23fxc.f0) / real(p.col(0));
	#endif
	for(int m=1;m<p.n_cols;++m)
		Coeffs.col(m) = 0.1*Coeffs.col(m-1); // some test values of parameters
	return Coeffs;
}

#endif // KERNEL_HEADER
