#ifndef FXCMODELS_HEADER
#define FXCMODELS_HEADER
// -------------------------------------------------------
// Copyright (C) 2019 by Dmitry R. Gulevich
// You may use or modify this code but not distribute it.
// -------------------------------------------------------
// fxc implementations:
//		Fxc_ALDA
//		Fxc_M1_test
// 		Fxc_M1_D0
//		Fxc_M1_DQV
// -------------------------------------------------------
#include "tdcdft-omxc/xc.hpp"
#include <armadillo>

using namespace arma;

/**
* fxc for ALDA
*/
struct Fxc_ALDA : xc::Omxc {
	Fxc_ALDA() {
		Mosc=0;
    	cout << "# ALDA" << endl;
	}
	cx_mat get_p(vec rho) {}
	cx_mat get_n23Coeffs(cx_mat p, vec n13) {}
};

/**
* fxc with M=1, used in one of the tests
*/
struct Fxc_M1_test : xc::Omxc {

	Fxc_M1_test() {
		Mosc=1;
    	cout << "# Number of oscillators: " << Mosc << endl;
	}

	cx_mat get_p(vec rho) {
		cx_double i(0.,1.);
		cx_mat p(rho.n_elem, Mosc);
		p.col(0) = (2.-2.*i)*xc::omega_pl(rho);
		return p;
	}

	// Returns: n^(2/3) * Cm
	cx_mat get_n23Coeffs(cx_mat p, vec n13) {
		mat n23f0finf = xc::get_n23f0finf(n13);
		cx_mat Coeffs(p.n_rows,p.n_cols);
		Coeffs.col(0) = conj(p.col(0)) % (n23f0finf.col(1)-n23f0finf.col(0)) / real(p.col(0));
		return Coeffs;
	}
};


/**
* fxc with M=1, used to produce results in Fig. 2 of D. R. Gulevich et al. PRB Rapid Communications (2019).
*/
struct Fxc_M1_D0 : xc::Omxc {

	Fxc_M1_D0() {
		Mosc=1;
    	cout << "# Fxc_M1_new kernel with D=0 (used in PRB Rapid paper)" << endl;
    	cout << "# Number of oscillators: " << Mosc << endl;
	}

	const double omega1 = 2.0;
	const double gamma1 = 1.0; // assert gamma1 < omega1
	const double imp1 = -gamma1;
	const double rep1 = sqrt(omega1*omega1-gamma1*gamma1);
	cx_double i=cx_double(0.,1.);

	cx_mat get_p(vec rho) {
		cx_mat p(rho.n_elem, Mosc);
		p.col(0) = (rep1+imp1*i)*xc::omega_pl(rho);
		return p;
	}

	// Returns: n^(2/3) * Cm
	cx_mat get_n23Coeffs(cx_mat p, vec n13) {
		mat n23fxc = xc::get_n23f0finf(n13);
		cx_mat Coeffs(p.n_rows,p.n_cols);
		Coeffs.col(0) = p.col(0) % (n23fxc.col(1) - n23fxc.col(0)) / real(p.col(0));
		return Coeffs;
	}

};

/**
* fxc with M=1 with the tangent of Im(fxc) at omega=0 given by Eq.(16) of Z. Qian and G. Vignale, Phys. Rev. B65, 235121 (2002).
*/
struct Fxc_M1_DQV : xc::Omxc {

	Fxc_M1_DQV() {
		Mosc=1;
    	cout << "# Fxc_M1_new kernel with D from QV" << endl;
    	cout << "# Number of oscillators: " << Mosc << endl;
	}

	const double omega1 = 2.0;
	const double gamma1 = 1.0; // assert gamma1 < omega1
	const double imp1 = -gamma1;
	const double rep1 = sqrt(omega1*omega1-gamma1*gamma1);
	cx_double i=cx_double(0.,1.);

	cx_mat get_p(vec rho) {
		cx_mat p(rho.n_elem, Mosc);
		p.col(0) = (rep1+imp1*i)*xc::omega_pl(rho);
		return p;
	}

	// Returns: n^(2/3) * Cm
	cx_mat get_n23Coeffs(cx_mat p, vec n13) {
		mat n23fxc = xc::get_n23f0finf(n13);
		cx_mat Coeffs(p.n_rows,p.n_cols);
		vec rho = n13 % n13 % n13; // can be imporved by passing rho directly
		vec n23ImDfxc = xc::get_n23ImDfxc(rho, n13);
		Coeffs.col(0) = p.col(0) % (n23fxc.col(1) - n23fxc.col(0) - i*conj(p.col(0))%n23ImDfxc ) / real(p.col(0));
		return Coeffs;
	}

};

#endif // FXCMODELS_HEADER
