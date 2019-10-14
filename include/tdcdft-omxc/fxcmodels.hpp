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

struct Fxc_ALDA : xc::FXC {
	Fxc_ALDA() {
		Mosc=0;
    	cout << "# ALDA" << endl;
	}
	cx_mat get_p(vec rho) {}
	cx_mat get_n23Coeffs(cx_mat p, vec n13) {}
};

struct Fxc_M1_test : xc::FXC {

	Fxc_M1_test() {
		Mosc=1;
    	cout << "# Number of oscillators: " << Mosc << endl;
	}

	cx_mat get_p(vec rho) {
		cx_double i(0.,1.);
		cx_mat p(rho.n_elem, Mosc);
		p.col(0) = (2.-2.*i)*omega_pl(rho);
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



struct Fxc_M1_D0 : xc::FXC {

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
		p.col(0) = (rep1+imp1*i)*omega_pl(rho);
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

struct Fxc_M1_DQV : xc::FXC {

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
		p.col(0) = (rep1+imp1*i)*omega_pl(rho);
		return p;
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
		mat n23fxc = xc::get_n23f0finf(n13);
		cx_mat Coeffs(p.n_rows,p.n_cols);
		vec rho = n13 % n13 % n13; // can be imporved by passing rho directly
		vec n23ImDfxc = get_n23ImDfxc(rho, n13);
		Coeffs.col(0) = p.col(0) % (n23fxc.col(1) - n23fxc.col(0) - i*conj(p.col(0))%n23ImDfxc ) / real(p.col(0));
		return Coeffs;
	}

};

#endif // FXCMODELS_HEADER
