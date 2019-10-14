#ifndef TDDFT_HEADER
#define TDDFT_HEADER
// -------------------------------------------------------
// Copyright (C) 2019 by Dmitry R. Gulevich
// You may use or modify this code but not distribute it.
// -------------------------------------------------------
#include <armadillo>
#include "tdcdft-omxc/systems.hpp"
#include "tdcdft-omxc/xc.hpp"

namespace tddft {

	// Time-dependent Kohn-Sham equation parameters
	struct Args {
		double tmin;
		double tmax;
		double dt;
		int ncorrections;
	};


	struct Result {
		vec t;
		vec dipole;
	};

	// Declarations for public functions
	Result Tdks(QuantumWell &qwell, Mesh<QuantumWell> &mesh, xc::FXC &fxc, dft::KsGs &ks, Args args);

} // end of namespace dft

#endif // TDDFT_HEADER
