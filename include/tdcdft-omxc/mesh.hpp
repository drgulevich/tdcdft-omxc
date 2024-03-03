#ifndef MESH_HEADER
#define MESH_HEADER
// -------------------------------------------------------
// Copyright (C) 2019 by Dmitry R. Gulevich
// You may use or modify this code but not distribute it.
// -------------------------------------------------------
#include <armadillo>

#include "tdcdft-omxc/systems.hpp"

using namespace arma;

namespace tdcdft {

/**
 * Mesh templates for various types of systems
 */
template <class T>
struct Mesh {};

/** Mesh for quantum well :
 *		M -- number of nodes
 *		dz -- spatial discretization, given by W/(M+1), where W is the width of the box
 * 	The nodes |.....| are defined as follows:
 * 		The width is W = (M+1)*dz
 * 		0...(M-1) label the internal nodes
 * 		The left and right infinite boundaries are (formally) at z(-1)=-W/2 and z(M)=W/2
 */
template <>
struct Mesh<QuantumWell> {
    Mesh(uword M, double dz) : M(M), dz(dz) {
        z.set_size(M);
        for (uword m = 0; m < M; ++m) z(m) = dz * (static_cast<double>(m) - 0.5 * static_cast<double>(M - 1));
    }
    uword M;
    double dz;
    vec z;
};

/**
 * Blank mesh template for atom (under development)
 */
template <>
struct Mesh<Atom> {};

/**
 * Blank mesh template for atomic cluster (under development)
 */
template <>
struct Mesh<Cluster> {};

/**
 * Blank mesh template for fullerene (under development)
 */
template <>
struct Mesh<Fullerene> {};

}  // namespace tdcdft

#endif  // MESH_HEADER
