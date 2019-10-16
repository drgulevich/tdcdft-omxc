#ifndef QWMODELS_HEADER
#define QWMODELS_HEADER
// -------------------------------------------------------
// Copyright (C) 2019 by Dmitry R. Gulevich
// You may use or modify this code but not distribute it.
// -------------------------------------------------------
// Quantum well implementations:
// 		QWell_UV_PRB_1998
//		DQWell_UV_PRB_1998
//		QWell_WU_PRL_2005
//		QWell_WU_PRL_2008
// -------------------------------------------------------
#include "tdcdft-omxc/au.hpp"
#include "tdcdft-omxc/systems.hpp"
#include "tdcdft-omxc/mesh.hpp"

/**
* Quantum well AlGaAs/GaAs/AlGaAs with one occupied orbital.
* Parameters of the well are from Ref. C. A. Ullrich and G. Vignale,
* "Collective intersubband transitions in quantum wells: A comparative density-functional study",
* Phys. Rev. B 58, 15756 (1998).
*/ 
struct QWell_UV_PRB_1998 : QuantumWell {
	const double eps=13.0;  // dielectric constant
	const double meff=0.07; // effective mass in units m0
	au::Effau effau = au::Effau(eps,meff); // set up effective atomic units
	const double W_GaAs = effau.to_au(384.,"Angstrom");
	const double W_AlGaAs = effau.to_au(384.,"Angstrom"); // AlGaAs layer from either side of GaAs
	const double Ewell=effau.to_au(250.,"meV"); // Quantum well depth, normalized
	double Efield = effau.to_au(0.,"mV/nm");
	QWell_UV_PRB_1998() {
		cout << "# Quantum well AlGaAs/GaAs/AlGaAs (Ullrich-Vignale-PRB-1998)" << endl;
		ns=0.97e11*effau.a0cm*effau.a0cm; // Sheet density, normalized
	}

	// Quantum well potential, including electric field
	double Vext(double z) {
		double Vqw = (fabs(z)<0.5*W_GaAs)?0.:Ewell;
		return Vqw - Efield*z;
	}

	// Seed electron density for DFT calculations (initial guess)
	double rho_seed(double z) {
		double value=0.;
		if(fabs(z)<0.5*W_GaAs) {
			double cosz = cos(M_PI*z/W_GaAs);
			value = 2.*(ns/W_GaAs)*cosz*cosz;
		}
		return value;
	}

	// Mesh generator for AlGaAs/GaAs/AlGaAs quantum well
	Mesh<QuantumWell> GaAs_mesh(int M_GaAs) {
		double dz = W_GaAs/M_GaAs;
		int M_AlGaAs = round(W_AlGaAs/dz);
		int M = M_AlGaAs + M_GaAs + M_AlGaAs;
		Mesh<QuantumWell> mesh(M, dz);
		return mesh;
	}
};


/**
* Double quantum well AlGaAs/GaAs/AlGaAs/GaAs/AlGaAs with one occupied orbital.
* Parameters of the well are from Ref. C. A. Ullrich and G. Vignale,
* "Collective intersubband transitions in quantum wells: A comparative density-functional study",
* Phys. Rev. B 58, 15756 (1998).
*/ 
struct DQWell_UV_PRB_1998 : QuantumWell {
	const double eps=13.0;  // dielectric constant
	const double meff=0.07; // effective mass in units m0
	au::Effau effau = au::Effau(eps,meff); // set up effective atomic units
	const double W1_GaAs = effau.to_au(85.,"Angstrom");
	const double Wmid_AlGaAs = effau.to_au(23.,"Angstrom"); // AlGaAs layer in the middle
	const double W2_GaAs = effau.to_au(73.,"Angstrom");
	const double W_AlGaAs = effau.to_au(200.,"Angstrom"); // AlGaAs layer from either side of GaAs
	const double W = W_AlGaAs + W1_GaAs + Wmid_AlGaAs + W2_GaAs + W_AlGaAs;	
	const double Ewell=effau.to_au(250.,"meV"); // Quantum well depth, normalized
	double Efield = effau.to_au(0.,"mV/nm");
	DQWell_UV_PRB_1998() {
		cout << "# Double quantum well lGaAs/GaAs/AlGaAs/GaAs/AlGaAs from Ullrich-Vignale-PRB-1998" << endl;
		ns=2.0e11*effau.a0cm*effau.a0cm; // Sheet density, normalized
	}

	// Double quantum well potential, including electric field
	// W_AlGaAs + W1_GaAs + Wmid_AlGaAs + W2_GaAs + W_AlGaAs
	double Vext(double z) {
		double Vqw;
		double z0 = 0.5*W;
		if( z < W_AlGaAs - z0 )
			Vqw=Ewell;
		else if( z < W_AlGaAs + W1_GaAs - z0 )
			Vqw=0.;
		else if( z < W_AlGaAs + W1_GaAs + Wmid_AlGaAs - z0 )
			Vqw=Ewell;
		else if( z < W_AlGaAs + W1_GaAs + Wmid_AlGaAs + W2_GaAs - z0 )
			Vqw=0.;
		else
			Vqw=Ewell;			
		return Vqw - Efield*z;
	}

	// Seed electron density for DFT calculations (initial guess)
	double rho_seed(double z) {
		double value=0.;
		double W1mid2 = W1_GaAs + Wmid_AlGaAs + W2_GaAs;
		if(fabs(z)<0.5*W1mid2) {
			double cosz = cos(M_PI*z/W1mid2);
			value = 2.*(ns/W1mid2)*cosz*cosz;
		}
		return value;
	}

	// Mesh generator for double quantum well
	Mesh<QuantumWell> DQWmesh(int M) {
		double dz = W/(M+1);
		Mesh<QuantumWell> mesh(M, dz);
		return mesh;
	}
};


/**
* Quantum well AlGaAs/GaAs/AlGaAs with one occupied orbital.
* Parameters of the well are from Ref. H. O. Wijewardane and C. A. Ullrich,
* "Time-Dependent Kohn-Sham Theory with Memory",
* Phys. Rev. Lett. 95, 086401 (2005).
*/ 
struct QWell_WU_PRL_2005 : QuantumWell {
	const double eps=13.0;  // dielectric constant
	const double meff=0.067; // effective mass in units m0
	au::Effau effau = au::Effau(eps,meff); // set up effective atomic units
	const double W_GaAs = effau.to_au(400.,"Angstrom");
	const double W_AlGaAs = 400.0/effau.a0; // AlGaAs layer from either side of GaAs
	const double Ewell=257.6/effau.Eh; // Quantum well depth, normalized
	double Efield = effau.to_au(0.,"mV/nm");
	QWell_WU_PRL_2005() {
		cout << "# Quantum well AlGaAs/GaAs/AlGaAs (WU-PRL-2005)" << endl;
		ns=1.0e11*effau.a0cm*effau.a0cm; // Sheet density, normalized
	}

	// Quantum well potential, including electric field
	double Vext(double z) {
		double Vqw = (fabs(z)<0.5*W_GaAs)?0.:Ewell;
		return Vqw - Efield*z;
	}

	// Seed electron density for DFT calculations (initial guess)
	double rho_seed(double z) {
		double value=0.;
		if(fabs(z)<0.5*W_GaAs) {
			double cosz = cos(M_PI*z/W_GaAs);
			value = 2.*(ns/W_GaAs)*cosz*cosz;
		}
		return value;
	}

	// Mesh generator for AlGaAs/GaAs/AlGaAs quantum well
	Mesh<QuantumWell> GaAs_mesh(int M_GaAs) {
		double dz = W_GaAs/M_GaAs;
		int M_AlGaAs = round(W_AlGaAs/dz);
		int M = M_AlGaAs + M_GaAs + M_AlGaAs;
		Mesh<QuantumWell> mesh(M, dz);
		return mesh;
	}
};


/**
* Quantum well AlGaAs/GaAs/AlGaAs with two occupied orbitals.
* Parameters of the well are from Ref. H. O. Wijewardane and C. A. Ullrich,
* "Real-Time Electron Dynamics with Exact-Exchange Time-Dependent Density-Functional Theory"
* Phys. Rev. Lett. 100, 056404 (2008).
*/ 
struct QWell_WU_PRL_2008 : QuantumWell {
	const double eps=13.0;  // dielectric constant
	const double meff=0.067; // effective mass in units m0
	au::Effau effau = au::Effau(eps,meff); // set up effective atomic units
	const double W_GaAs = effau.to_au(400.,"Angstrom");
	const double W_AlGaAs = 400.0/effau.a0; // AlGaAs layer from either side of GaAs
	const double Ewell=257.6/effau.Eh; // Quantum well depth, normalized
	double Efield = effau.to_au(0.,"mV/nm");
	QWell_WU_PRL_2008() {
		cout << "# Quantum well AlGaAs/GaAs/AlGaAs with Nocc=2 (WU-PRL-2008)" << endl;
		ns=2.2e11*effau.a0cm*effau.a0cm; // Sheet density, normalized
	}

	// Quantum well potential, including electric field
	double Vext(double z) {
		double Vqw = (fabs(z)<0.5*W_GaAs)?0.:Ewell;
		return Vqw - Efield*z;
	}

	// Seed electron density for DFT calculations (initial guess)
	double rho_seed(double z) {
		double value=0.;
		if(fabs(z)<0.5*W_GaAs) {
			double cosz = cos(M_PI*z/W_GaAs);
			value = 2.*(ns/W_GaAs)*cosz*cosz;
		}
		return value;
	}

	// Mesh generator for AlGaAs/GaAs/AlGaAs quantum well
	Mesh<QuantumWell> GaAs_mesh(int M_GaAs) {
		double dz = W_GaAs/M_GaAs;
		int M_AlGaAs = round(W_AlGaAs/dz);
		int M = M_AlGaAs + M_GaAs + M_AlGaAs;
		Mesh<QuantumWell> mesh(M, dz);
		return mesh;
	}
};

#endif // QWMODELS_HEADER
