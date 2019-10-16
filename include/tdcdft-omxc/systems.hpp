#ifndef SYSTEMS_HEADER
#define SYSTEMS_HEADER
// -------------------------------------------------------
// Copyright (C) 2019 by Dmitry R. Gulevich
// You may use or modify this code but not distribute it.
// -------------------------------------------------------

/**
* Electron system (an abstract class used by derived classes).
*/
struct ElectronSystem {
	virtual double Vext(double z) = 0; // External potential
	virtual double rho_seed(double z) = 0; // Seed electron density distribution
	virtual ~ElectronSystem() {} // virtual destructor
};

/**
* Quantum well
*/
struct QuantumWell : ElectronSystem {
	double ns; 	// Sheet electron density
	double EF; // Fermi energy
	uword Nocc; // Number of occupied orbitals
};

/**
* Atom (under development)
*/
struct Atom : ElectronSystem {
	int Z; // Ion charge
};

/**
* Atomic cluster (under development)
*/
struct Cluster : ElectronSystem {
	int Ne; // Number of electrons
};

/**
* Fullerene (under development)
*/
struct Fullerene : ElectronSystem {
	int Ne; // Number of electrons
};

#endif // SYSTEMS_HEADER
