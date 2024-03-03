#ifndef HARTREE_HEADER
#define HARTREE_HEADER
// -------------------------------------------------------
// Copyright (C) 2019 by Dmitry R. Gulevich
// You may use or modify this code but not distribute it.
// -------------------------------------------------------
#include <armadillo>

#include "tdcdft-omxc/mesh.hpp"
#include "tdcdft-omxc/systems.hpp"

using namespace arma;

namespace tdcdft {

/**
 * Calculate Hartree potential for a quantum well.
 * @param mesh Mesh<QuantumWell> object.
 * @param rho Electron density array.
 * @param ns Sheet electron density in a quantum well.
 */
vec get_VHartree(Mesh<QuantumWell> &mesh, vec &rho, double ns);

}  // namespace tdcdft

#endif  // HARTREE_HEADER
