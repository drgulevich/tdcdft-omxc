#ifndef TDDFT_HEADER
#define TDDFT_HEADER
// -------------------------------------------------------
// Copyright (C) 2019 by Dmitry R. Gulevich
// You may use or modify this code but not distribute it.
// -------------------------------------------------------
#include "tdcdft-omxc/systems.hpp"
#include "tdcdft-omxc/mesh.hpp"
#include "tdcdft-omxc/xc.hpp"
#include "tdcdft-omxc/dft.hpp"
#include <armadillo>

namespace tddft {

	/** 
	* Parameters for the time-dependent Kohn-Sham equation (TDKS)
	*/
	struct Args {
		double tmin;
		double tmax;
		double dt;
		int ncorrections;
	};

	/** 
	* Desired struct for the results of the TDKS evolution
	*/
	struct Result {
		vec t;
		vec dipole;
	};

	/**
	* Evolve the time-dependent Kohn-Sham equation with memory
	*/
	Result Tdks(QuantumWell &qwell, Mesh<QuantumWell> &mesh, xc::Omxc &fxc, dft::KsGs &ks, Args args);

} // end of namespace dft

#endif // TDDFT_HEADER
