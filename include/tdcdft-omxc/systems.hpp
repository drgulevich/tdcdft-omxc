#ifndef SYSTEMS_HEADER
#define SYSTEMS_HEADER
// -------------------------------------------------------
// Copyright (C) 2019 by Dmitry R. Gulevich
// You may use or modify this code but not distribute it.
// -------------------------------------------------------
#include <armadillo>

using namespace arma;

struct ElectronSystem {
	virtual double Vext(double z) = 0; // External potential
	virtual double rho_seed(double z) = 0; // Seed electron density distribution
	virtual ~ElectronSystem() {} // virtual destructor
};

struct QuantumWell : ElectronSystem {
	double ns; 	// Sheet electron density
	double EF; // Fermi energy
	uword Nocc; // Number of occupied orbitals
};

struct Fullerene : ElectronSystem {
	int Ne; // Number of electrons
};

struct Atom : ElectronSystem {
	int Z; // Ion charge
};


// Mesh templates for various types of systems
template <class T>
struct Mesh {};

// Mesh for quantum well
// 	Input:
//		M -- number of nodes
//		dz -- discretization, W/(M+1)
// 	The nodes |.....| are defined as follows:
// 		The width is W = (M+1)*dz
// 		0...(M-1) label the internal nodes 
// 		The left and right infinite boundaries are, formally, at z(-1)=-W/2 and z(M)=W/2
template <>
struct Mesh <QuantumWell> {
	Mesh (uword M, double dz) : M(M), dz(dz) {
		z.set_size(M);
		for(uword m=0; m<M; ++m)
			z(m) = dz*(static_cast<double>(m)-0.5*static_cast<double>(M-1));
	} 
	uword M;
	double dz;
	vec z;
};

template <>
struct Mesh <Fullerene> {
};

#endif // SYSTEMS_HEADER
