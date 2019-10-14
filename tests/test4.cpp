// -------------------------------------------------------
// Copyright (C) 2019 by Dmitry R. Gulevich
// You may use or modify this code but not distribute it.
// -------------------------------------------------------
#include <iostream>
#include <chrono>  // timing
#include "catch2/catch.hpp" // testing framework
#include "tdcdft-omxc/xc.hpp"
#include "tdcdft-omxc/qwmodels.hpp" // supplementary QuantumWell structs
#include "tdcdft-omxc/dft.hpp"
#include "tdcdft-omxc/tddft.hpp"
#include "tdcdft-omxc/tools.hpp"
#include "tdcdft-omxc/fxcmodels.hpp"

using namespace std; 
using namespace arma; 

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
	recommend_num_threads(1);

	InfQW qwell;
	Mesh<QuantumWell> mesh = qwell.mesh(100);
	dft::KsGs ks;
	dft::KsArgs ksargs = {10,100,0.1}; // { number of bands, iterations, smoothness }

	qwell.Efield = qwell.effau.to_au(0.5,"mV/nm");
	ks = dft::Ks(qwell, mesh, ksargs); // Solve KS equation 

	qwell.Efield = qwell.effau.to_au(0.0,"mV/nm");
	tddft::Args args = {0.,500.,0.05,3};

	Fxc_ALDA fxc;

	auto start = std::chrono::high_resolution_clock::now();

	tddft::Result tdks = tddft::Tdks(qwell, mesh, fxc, ks, args);

	auto finish = std::chrono::high_resolution_clock::now();

/*	for(uword i=0;i<tdks.t.n_elem;++i)
		cout << tdks.t(i) << " " << tdks.dipole(i) << endl;*/

    std::cout << "# Timing: "
              << std::chrono::duration_cast<std::chrono::milliseconds>(finish-start).count()
              << " ms\n";

	vec result(3);
	uword npoints = tdks.t.n_elem;
    for(uword i=0;i<result.n_elem;++i) {
		result(i) = tdks.dipole(npoints-1-i);
	}

    vec reference = {0.0117494,0.0103206,0.00886947}; // EXP_ORDER 3

   	REQUIRE( result(0) == Approx( reference(0) ).epsilon(0.001) );
   	REQUIRE( result(1) == Approx( reference(1) ).epsilon(0.001) );
   	REQUIRE( result(2) == Approx( reference(2) ).epsilon(0.001) );
}
