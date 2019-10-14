#ifndef XC_HEADER
#define XC_HEADER
// -------------------------------------------------------
// Copyright (C) 2019 by Dmitry R. Gulevich
// You may use or modify this code but not distribute it.
// -------------------------------------------------------
#include <armadillo>

using namespace arma;

namespace xc {

	struct FXC {
		uword Mosc;
		vec omega_pl(vec rho) {
			static const double factor = sqrt(4.*M_PI);
			return factor*sqrt(rho);
		}
		virtual cx_mat get_p(vec rho) = 0; // pure virtual
		virtual cx_mat get_n23Coeffs(cx_mat p, vec n13) = 0; // pure virtual
		virtual ~FXC() {} // virtual destructor
	};

	// Functions of cr=pow(rho,1/3.):
	vec get_ex(vec cr);
	vec get_ec(vec cr);
	vec get_exc(vec cr);
	vec get_VxcLDA(vec cr); // passing cr to VxcLDA is faster than passing rho
	vec get_n23f0(vec cr);
	vec get_n23finf(vec cr);
	mat get_n23f0finf(vec cr);

	// Functions of rho
	vec omega_pl(vec rho);
}

#endif // XC_HEADER
