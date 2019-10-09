// -------------------------------------------------------
// Copyright (C) 2019 by Dmitry R. Gulevich
// You may use or modify this code but not distribute it.
// -------------------------------------------------------
#include <iostream>
#include <chrono>  // timing
#include "catch2/catch.hpp" // testing framework
#include "../qwmodels.hpp" // supplementary QuantumWell structs
#include "../dft.hpp"
#include "../tddft.hpp"

// Important note: for small problem like this on it is best to execute 
// in serial by switching off OpenMP:
// 		$ export OMP_NUM_THREADS=1
// Check:
// 		$ echo $OMP_NUM_THREADS

using namespace std; 
using namespace arma; 

TEST_CASE( "Comparison to qwell.GaAs_mesh(100), ksargs = {10,100,0.1}, args = {0.,100.,0.2,3}", "[tddft]" ) {

	cout << "Test 4: " << endl;

	QWell_WU_PRL_2005 qwell;
	Mesh<QuantumWell> mesh = qwell.GaAs_mesh(100);
	dft::KsGs ks;
	dft::KsArgs ksargs = {10,100,0.1}; // { number of bands, iterations, smoothness }

	qwell.Efield = qwell.effau.to_au(0.5,"mV/nm");
	ks = dft::Ks(qwell, mesh, ksargs); // Solve KS equation 

	qwell.Efield = qwell.effau.to_au(0.0,"mV/nm");
	tddft::Args args = {0.,100.,0.2,3};

	auto start = std::chrono::high_resolution_clock::now();

	tddft::Result tdks = tddft::Tdks(qwell, mesh, ks, args);

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

	// Comment: these values (obtained with timestep 0.2) are not converged, use for test only
    //vec reference = {-0.0150219,-0.0148807,-0.0141548}; // EXP_ORDER not defined (i.e. full exp is used)
    vec reference = {-0.0150072,-0.0148585,-0.014131}; // EXP_ORDER 3

   	REQUIRE( result(0) == Approx( reference(0) ).epsilon(0.001) );
   	REQUIRE( result(1) == Approx( reference(1) ).epsilon(0.001) );
   	REQUIRE( result(2) == Approx( reference(2) ).epsilon(0.001) );
}
