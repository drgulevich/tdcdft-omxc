// -------------------------------------------------------
// Copyright (C) 2019 by Dmitry R. Gulevich
// You may use or modify this code but not distribute it.
// -------------------------------------------------------
/*
	Exchange and correlation energies per particle

	Exchange: 
		standard formula, see, e.g., Eq.(34) in: 
	[1] C. A. Ullrich and Z.-H. Yang, A Brief Compendium of Time-Dependent Density Functional Theory, Braz. J. Phys. 44, 154 (2014).

	Correlation:
		spin-unpolarized case in Chachiyo-Karasiev parametrization
	[2] T. Chachiyo, “Communication: Simple and accurate uniform electron gas correlation energy for the full range of densities”, 
		J. Chem. Phys., vol. 	145, p. 021101, 2016.
	[3] V. V. Karasiev, “Comment on ‘Communication: Simple and accurate uniform electron gas correlation energy for the full range of densities”,
		J. Chem. Phys., vol. 145, p. 157101, 2016.
*/
#include "tdcdft-omxc/xc.hpp"
#include <armadillo>
#include <cmath>

namespace tdcdft { namespace xc {

	const double ax = -(3./4.)*cbrt(3./M_PI);
	const double ac = (log(2.)-1.)/(2.*M_PI*M_PI);
	const double B1 = 21.7392245*cbrt(4.*M_PI/3.);
	const double B2 = 20.4562557*cbrt(16.*M_PI*M_PI/9.);

	// exchange energy
	vec get_ex(vec cr) {
		return ax*cr;
	}

	// correlation energy
	vec get_ec(vec cr) {
		vec result(cr.n_elem);
		for(uword m=0;m<cr.n_elem;++m) {
			double cr1 = cr(m);
			result(m) = ac*log(1. + B1*cr1 + B2*cr1*cr1);
		}
		return result;
	}

	// exchange-correlation energy
	vec get_exc(vec cr) {
		vec result(cr.n_elem);
		for(uword m=0;m<cr.n_elem;++m) {
			double cr1 = cr(m);
			result(m) = ax*cr1 + ac*log(1. + B1*cr1 + B2*cr1*cr1);
		}
		return result;
	}

	// LDA xc potential
	vec get_VxcLDA(vec cr) {
		vec result(cr.n_elem);
		for(uword m=0;m<cr.n_elem;++m) {
			// exchange
			double cr1 = cr(m);
			result(m) = 4*ax*cr1/3.;
	
			// correlation
			double cr2 = cr1*cr1;
			double arg = 1. + B1*cr1 + B2*cr2;
			result(m) += ac*((B1*cr1 + 2*B2*cr2)/(3.*arg) + log(arg));
		}
		return result;
	}


	// Returns: n^(2/3) * f0, 
	//		where f0 = d^2(n*exc(n))/dn^2 from Eq.(i) of GK
	vec get_n23f0(vec cr) {
		vec result(cr.n_elem);
		for(uword m=0;m<cr.n_elem;++m) {
			// exchange
			result(m) = 4.*ax/9.;
	
			// correlation
			double cr1 = cr(m);
			double cr2 = cr1*cr1;
			double cr3 = cr1*cr2;
			double arg = 1. + B1*cr1 + B2*cr2;
			result(m) += ac*(4.*B1 + (3.*B1*B1 + 10.*B2)*cr1 + 10.*B1*B2*cr2 + 6.*B2*B2*cr3)/(9.*arg*arg);
		}
		return result;
	}

	// Returns: n^(2/3) * finf, 
	//		where finf is given by Eq.(2) Iwamoto
	// Warning! Division by zero at n=0: log(arg)/cr1 -> 1 at n->0
	vec get_n23finf(vec cr) {
		vec result(cr.n_elem);	
		for(uword m=0;m<cr.n_elem;++m) {
			// exchange
			result(m) = 4.*ax/15.;
	
			// correlation
			double cr1 = cr(m);
			double cr2 = cr1*cr1;
			double arg = 1. + B1*cr1 + B2*cr2;
			result(m) += 2.*ac*(13.*(B1 + 2.*B2*cr1)/arg - 11.*log(arg)/cr1)/15.;
		}
		return result;
	}


	mat get_n23f0finf(vec cr) {
		mat n23f0finf(cr.n_elem,2);
		for(uword m=0;m<cr.n_elem;++m) {
			// exchange
			n23f0finf(m,0) = 4.*ax/9.; // f0
			n23f0finf(m,1) = 4.*ax/15.; // finf
	
			// correlation
			double cr1 = cr(m);
			double cr2 = cr1*cr1;
			double cr3 = cr1*cr2; // rho(m)
			double arg = 1. + B1*cr1 + B2*cr2;
			//fxc.crho(m) = cr1;
			n23f0finf(m,0) += ac*(4.*B1 + (3.*B1*B1 + 10.*B2)*cr1 + 10.*B1*B2*cr2 + 6.*B2*B2*cr3)/(9.*arg*arg); // f0
			n23f0finf(m,1) += 2.*ac*(13.*(B1 + 2.*B2*cr1)/arg - 11.*log(arg)/cr1)/15.; // finf
		}
		return n23f0finf;
	}

	/** 
	* Eq.(16) of Z. Qian and G. Vignale, Phys. Rev. B65, 235121 (2002).
	*/
	vec S3L(vec lam) {
		vec result(lam.n_elem);
		static const double factor = -1./(45.*M_PI);
		for(uword m=0;m<lam.n_elem;++m) {
			double part1 = 5.-(lam[m]+5./lam[m])*atan(lam[m]);
			double part2 = -(2./lam[m])*asin(lam[m]/sqrt(1.+lam[m]*lam[m]));
			double arg = 1./(lam[m]*sqrt(2.+lam[m]*lam[m]));
			double part3 = arg*(M_PI - 2.*atan(arg));
			result[m] = factor*(part1+part2+part3);
		}
		return result;
	}
    
	// Returns: n^(2/3) * (d/domega)Imfxc(n,0)
	vec get_n23ImDfxc(vec rho, vec n13) {
		static const double factorkf = cbrt(3.*M_PI*M_PI); // kf = factorkf*n13;
		vec kf = factorkf*n13;
		vec lam = sqrt(M_PI*kf); // 2kf/ks    
		return -factorkf*S3L(lam)/(M_PI*M_PI*rho);
	}

	vec omega_pl(vec rho) {
		static const double factor = sqrt(4.*M_PI);
		return factor*sqrt(rho);
	}

} }
