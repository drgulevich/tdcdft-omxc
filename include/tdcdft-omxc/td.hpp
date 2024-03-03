#ifndef TDDFT_HEADER
#define TDDFT_HEADER
// -------------------------------------------------------
// Copyright (C) 2019 by Dmitry R. Gulevich
// You may use or modify this code but not distribute it.
// -------------------------------------------------------
#include <armadillo>

#include "tdcdft-omxc/gs.hpp"
#include "tdcdft-omxc/mesh.hpp"
#include "tdcdft-omxc/systems.hpp"
#include "tdcdft-omxc/xc.hpp"

/**
 * Order of expansion of exponential in the memory propagator.
 * Implemented values : 1, 2 and 3.
 * EXP_ORDER 3 is recommended.
 * Full exponential is used if not defined.
 */
#define EXP_ORDER 3

using namespace arma;

namespace tdcdft {

/**
 * Parameters for the time-dependent Kohn-Sham equation (TDKS).
 */
struct TdArgs {
    double tmin;
    double tmax;
    double dt;
    int ncorrections;
};

/**
 * Struct for the results of the TDKS calculation : evolution of the dipole moment.
 */
struct TdDipole {
    vec t;
    vec value;
};

/**
 * Evolve the time-dependent Kohn-Sham equation with memory.
 */
TdDipole Tdks(QuantumWell &qwell, Mesh<QuantumWell> &mesh, xc::Omxc &fxc, KsGs &ks, TdArgs args);

}  // namespace tdcdft

#endif  // TDDFT_HEADER
