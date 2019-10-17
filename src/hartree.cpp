// -------------------------------------------------------
// Copyright (C) 2019 by Dmitry R. Gulevich
// You may use or modify this code but not distribute it.
// -------------------------------------------------------
#include "tdcdft-omxc/hartree.hpp"
#include "tdcdft-omxc/systems.hpp"
#include "tdcdft-omxc/mesh.hpp"
#include <armadillo>

using namespace arma;

namespace tdcdft {

	/**
	* Calculate Hartree potential for a quantum well discretized by Mesh<QuantumWell>
	*/
	vec get_VHartree(Mesh<QuantumWell> &mesh, vec &rho, double ns) {
		double dzdzfactor = -4.*M_PI*mesh.dz*mesh.dz/12.;
		vec VH(mesh.M);
		VH(0) = 0.;
		VH(1) = 2.*VH(0) + dzdzfactor*(rho(1) + 10.*rho(0)); // assuming zero derivative, i.e. VH(0)=VH(-1)=0
		for(uword m=1; m<mesh.M-1; ++m)
			VH(m+1) = 2.*VH(m) - VH(m-1) + dzdzfactor*(rho(m+1) + 10.*rho(m) + rho(m-1));
		VH += 2.*M_PI*ns*mesh.z;
		return VH;
	}

}
