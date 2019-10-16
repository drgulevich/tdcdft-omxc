#ifndef DFT_HEADER
#define DFT_HEADER
// -------------------------------------------------------
// Copyright (C) 2019 by Dmitry R. Gulevich
// You may use or modify this code but not distribute it.
// -------------------------------------------------------
#include "tdcdft-omxc/systems.hpp"
#include "tdcdft-omxc/mesh.hpp"
#include <armadillo>

using namespace arma;

namespace dft {

	/**
	* Kohn-Sham ground state.
	*/
	struct KsGs {
		vec rho; /**< Electron density. */
		vec Veff; /**< Effective Kohn-Sham potential. */
		vec ens; /**< Kohn-Sham energies. */
		mat orbs; /**< Kohn-Sham orbitals, assumed to be real. */
		std::vector<vec> rhoarr; /**< Container for electron density values for convergence checks. */
	};

	/** 
	* Parameters for the stationary Kohn-Sham equation.
	*/
	struct KsArgs {
		int nbands; /**< Number of bands to calculate for. */
		int niters; /**< Number of iterations. */
		double theta; /**< Parameter for smooth update of Veff from one iteration to the next one. */
	};

	/**
	* Calculate Hartree potential for a quantum well.
	* @param mesh Mesh<QuantumWell> object.
	* @param rho Electron density array.
	* @param ns Sheet electron density in a quantum well.
	*/
	vec get_VHartree(Mesh<QuantumWell> &mesh, vec &rho, double ns);

	/**
	* Solve the stationary Kohn-Sham equation.
	* @param qwell QuantumWell object.
	* @param mesh Mesh<QuantumWell> object.
	* @param ksargs Parameters for the stationary Kohn-Sham equation as a KsArgs struct.
	*/
	KsGs Ks(QuantumWell &qwell, Mesh<QuantumWell> &mesh, KsArgs ksargs);

} // end of namespace dft

#endif // DFT_HEADER
