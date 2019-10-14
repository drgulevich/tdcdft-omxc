#ifndef XC_HEADER
#define XC_HEADER
// -------------------------------------------------------
// Copyright (C) 2019 by Dmitry R. Gulevich
// You may use or modify this code but not distribute it.
// -------------------------------------------------------
#include <armadillo>

using namespace arma;

namespace xc {

	struct Fxc {
		Fxc(uword size) {
			f0.set_size(size);
			finf.set_size(size);
		}
		vec f0;
		vec finf;
	};

	struct FXC {
		/*Kernel_values get_kernel_values(vec rho, vec n13) {
			kernel.p = get_p(rho);
			kernel.n23Coeffs = get_n23Coeffs(p,n13);
			return kernel;
		}*/
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
	Fxc get_n23fxc(vec cr);

	vec omega_pl(vec rho);
}

#endif // XC_HEADER
