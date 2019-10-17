#ifndef XC_HEADER
#define XC_HEADER
// -------------------------------------------------------
// Copyright (C) 2019 by Dmitry R. Gulevich
// You may use or modify this code but not distribute it.
// -------------------------------------------------------
#include <armadillo>

using namespace arma;

namespace tdcdft { namespace xc {

	/**
	* An abstract class for xc kernels in the OMXC form
	*/
	struct Omxc {
		uword Mosc;
		virtual cx_mat get_p(vec rho) = 0; // pure virtual
		virtual cx_mat get_n23Coeffs(cx_mat p, vec n13) = 0; // pure virtual
		virtual ~Omxc() {} // virtual destructor
	};

	/**
	* Exchange energy per particle.
	* See, e.g. Eq.(34) in C. A. Ullrich and Z.-H. Yang, A Brief Compendium of Time-Dependent Density Functional Theory, Braz. J. Phys. 44, 154 (2014).
	*/
	vec get_ex(vec cr /** Cubic root of the electron density, pow(rho,1/3.) */);

	/**
	* Correlation energy per particle for spin-unpolarized case in Chachiyo-Karasiev parametrization:
	* T. Chachiyo, "Communication: Simple and accurate uniform electron gas correlation energy for the full range of densities", 
	* J. Chem. Phys. 145, 021101 (2016).
	* V. V. Karasiev, "Comment on ‘Communication: Simple and accurate uniform electron gas correlation energy for the full range of densities",
	* J. Chem. Phys. 145, 157101 (2016).
	*/
	vec get_ec(vec cr /** Cubic root of the electron density, pow(rho,1/3.) */);

	/**
	* Exchange-correlation energy per particle. A sum of get_ex(vec cr) and get_ec(vec cr) defined above.
	*/
	vec get_exc(vec cr /** Cubic root of the electron density, pow(rho,1/3.) */);

	/**
	* LDA exchange-correlation potential.
	* Note, the cubic root of the electron density cr=pow(rho,1/3.) is passed as an argument due to performance considerations.
	*/
	vec get_VxcLDA(vec cr /** Cubic root of the electron density, pow(rho,1/3.) */);

	/**
	* Zero-frequency limit of longitudinal fxc.
	* Returns n^(2/3) * f0, where f0 = d^2(n*exc(n))/dn^2 given by Eq.(i) in Ref. E. K. U. Gross and W. Kohn, Phys. Rev. Lett. 55, 2850 (1985).
	*/
	vec get_n23f0(vec cr /** Cubic root of the electron density, pow(rho,1/3.) */);

	/**
	* Infinite-frequency limit of longitudinal fxc.
	* Returns n^(2/3) * finf, where finf is given by Eq.(2) in Ref. N. Iwamoto and E. K. U. Gross, Phys. Rev. B 35, 3003 (1987).
	*/
	vec get_n23finf(vec cr /** Cubic root of the electron density, pow(rho,1/3.) */);

	/**
	* Zero- and infinite-frequency limits of longitudinal fxc calculated simultaneously. 
	* Column (0) of the result coincides with get_n23f0(vec cr) defined above.
	* Column (1) of the result coincides with get_n23finf(vec cr) defined above.
	*/
	mat get_n23f0finf(vec cr /** Cubic root of the electron density, pow(rho,1/3.) */);

	/**
	* Eq.(16) of Z. Qian and G. Vignale, Phys. Rev. B65, 235121 (2002).
	*/
	vec S3L(vec lam);
     
	/**
	* Returns n^(2/3) * (d/domega)Imfxc(n,0), see Ref. Z. Qian and G. Vignale, Phys. Rev. B65, 235121 (2002).
	*/
	vec get_n23ImDfxc(vec rho, vec n13);

	/**
	* Plasma frequency of 3D electron gas
	*/
	vec omega_pl(vec rho /** Electron density */);

} }

#endif // XC_HEADER
