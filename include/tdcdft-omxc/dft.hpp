#ifndef DFT_HEADER
#define DFT_HEADER
// -------------------------------------------------------
// Copyright (C) 2019 by Dmitry R. Gulevich
// You may use or modify this code but not distribute it.
// -------------------------------------------------------
#include <armadillo>
#include "tdcdft-omxc/systems.hpp"

namespace dft {

	struct KsGs {
		vec rho; // rho_gs = ns*psi*psi.t();
		vec Veff; // KS potential corresponding to rho
		vec ens; // KS energies
		mat orbs; // KS orbitals, assumed to be real
		std::vector<vec> rhoarr; // container for densities to check convergence
	};

	// Stationary Kohn-Sham equation parameters
	struct KsArgs {
		int nbands;
		int niters; // number of iterations
		double theta; // parameter for smooth update
	};

	// Declarations for public functions
	vec get_VHartree(Mesh<QuantumWell> &mesh, vec &rho, double ns);
	KsGs Ks(QuantumWell &qwell, Mesh<QuantumWell> &mesh, KsArgs ksargs);

} // end of namespace dft

#endif // DFT_HEADER
