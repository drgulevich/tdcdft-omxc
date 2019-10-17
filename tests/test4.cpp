// -------------------------------------------------------
// Copyright (C) 2019 by Dmitry R. Gulevich
// You may use or modify this code but not distribute it.
// -------------------------------------------------------
#include <iostream>
#include <chrono>  // timing
#include "catch2/catch.hpp" // testing framework
#include "tdcdft-omxc/xc.hpp"
#include "tdcdft-omxc/qwmodels.hpp" // supplementary QuantumWell structs
#include "tdcdft-omxc/gs.hpp"
#include "tdcdft-omxc/td.hpp"
#include "tdcdft-omxc/tools.hpp"
#include "tdcdft-omxc/fxcmodels.hpp"

using namespace std; 
using namespace arma; 
using namespace tdcdft;

struct InfQW : QuantumWell {

	const double eps=13.0;  // dielectric constant, from WU-PRL-2005
	const double meff=0.067; // effective mass in units m0, from WU-PRL-2005
	au::Effau effau = au::Effau(eps,meff); // set up effective atomic units
	const double W = effau.to_au(400.,"Angstrom"); // from WU-PRL-2005
	const double Ewell=257.6/effau.Eh; // Quantum well depth, from WU-PRL-2005
	double Efield = effau.to_au(0.,"mV/nm");
	InfQW() {
		cout << "# GaAs quantum well with infinite walls" << endl;
		ns=1.0e11*effau.a0cm*effau.a0cm; // Sheet density, from WU-PRL-2005
	}

	// Quantum well potential, including electric field
	double Vext(double z) {
		return -Efield*z;
	}

	// Seed electron density for DFT calculations (initial guess)
	double rho_seed(double z) {
		double cosz = cos(M_PI*z/W);
		return 2.*(ns/W)*cosz*cosz;
	}

	// Mesh generator
	Mesh<QuantumWell> mesh(int M) {
		double dz = W/(M+1);
		Mesh<QuantumWell> mesh(M, dz);
		return mesh;
	}
};



TEST_CASE( "ALDA test", "[tddft]" ) {

	cout << "Test 4: " << endl;
	tools::recommend_num_threads(1);

	InfQW qwell;
	Mesh<QuantumWell> mesh = qwell.mesh(100);
	KsGs ks;
	KsArgs ksargs = {10,100,0.1}; // { number of bands, iterations, smoothness }

	qwell.Efield = qwell.effau.to_au(0.5,"mV/nm");
	ks = Ks(qwell, mesh, ksargs); // Solve KS equation 

	qwell.Efield = qwell.effau.to_au(0.0,"mV/nm");
	TdArgs args = {0.,500.,0.05,3};

	xc::Fxc_ALDA fxc;

	auto start = std::chrono::high_resolution_clock::now();

	TdDipole dipole = Tdks(qwell, mesh, fxc, ks, args);

	auto finish = std::chrono::high_resolution_clock::now();

    std::cout << "# Timing: "
              << std::chrono::duration_cast<std::chrono::milliseconds>(finish-start).count()
              << " ms\n";

    vec reference = {0.0117494,0.0103206,0.00886947}; // EXP_ORDER 3

	uword nvals = dipole.t.n_elem;
   	REQUIRE( dipole.value(nvals-1) == Approx( reference(0) ).epsilon(0.001) );
   	REQUIRE( dipole.value(nvals-2) == Approx( reference(1) ).epsilon(0.001) );
   	REQUIRE( dipole.value(nvals-3) == Approx( reference(2) ).epsilon(0.001) );
}
