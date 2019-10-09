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

#include <cmath>
#include "xc.hpp"

namespace xc {

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


	/*
	// This function gives no advantage over pow(rho,1/3.). 
	// Actually it turns out even a bit slower. Commented out.
	vec get_cbrt(vec rho) {
		vec result(rho.n_elem);
		for(uword m=0;m<rho.n_elem;++m) 
			result(m) = cbrt(rho(m));
		return result;
	}*/

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
	// Warning: division by zero at n=0: 
	//			log(arg)/cr1 -> 1 at n->0
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


	Fxc get_n23fxc(vec cr) {
		Fxc n23fxc(cr.n_elem);
		for(uword m=0;m<cr.n_elem;++m) {
			// exchange
			n23fxc.f0(m) = 4.*ax/9.;
			n23fxc.finf(m) = 4.*ax/15.;
	
			// correlation
			double cr1 = cr(m);
			double cr2 = cr1*cr1;
			double cr3 = cr1*cr2; // rho(m)
			double arg = 1. + B1*cr1 + B2*cr2;
			//fxc.crho(m) = cr1;
			n23fxc.f0(m) += ac*(4.*B1 + (3.*B1*B1 + 10.*B2)*cr1 + 10.*B1*B2*cr2 + 6.*B2*B2*cr3)/(9.*arg*arg);
			n23fxc.finf(m) += 2.*ac*(13.*(B1 + 2.*B2*cr1)/arg - 11.*log(arg)/cr1)/15.;
		}
		return n23fxc;
	}


	vec omega_pl(vec rho) {
		static const double factor = sqrt(4.*M_PI);
		return factor*sqrt(rho);
	}

/*
class Fxc {
	public:
		Fxc();
		std::vector<double> fxc0(std::vector<double> rho);

		int npoles;
    	void omega_pl(double *rho, double *result, int dim);
};
*/

} // end of xc namespace

