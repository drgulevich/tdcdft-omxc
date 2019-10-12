// -------------------------------------------------------
// Copyright (C) 2019 by Dmitry R. Gulevich
// You may use or modify this code but not distribute it.
// -------------------------------------------------------
#include <iostream>
#include <chrono>  // timing
#include "catch2/catch.hpp" // testing framework
#include "tdcdft-omxc/qwmodels.hpp" // supplementary QuantumWell structs
#include "tdcdft-omxc/dft.hpp"

// Important note: for small problem like this on it is best to execute 
// in serial by switching off OpenMP:
// 		$ export OMP_NUM_THREADS=1
// Check:
// 		$ echo $OMP_NUM_THREADS

using namespace std; 
using namespace arma; 

TEST_CASE( "Comparison to qwell.GaAs_mesh(100), ksargs = {10,100,0.1} at field 0.5 mV/nm", "[dft]" ) {

	cout << "Test 1: " << endl;

	QWell_WU_PRL_2005 qwell;
	Mesh<QuantumWell> mesh = qwell.GaAs_mesh(100);
	dft::KsGs ks;
	dft::KsArgs ksargs = {10,100,0.1}; // { number of bands, iterations, smoothness }
	qwell.Efield = qwell.effau.to_au(0.5,"mV/nm");
	ks = dft::Ks(qwell, mesh, ksargs); // Solve KS equation 

	vec result(3);
    for(uword i=0;i<result.n_elem;++i) {
		result(i) = qwell.effau.Eh*(ks.ens(1+i)-ks.ens(0));
	}
	cout << "# Calculated subband spacings at Efield = " << qwell.effau.from_au( qwell.Efield,"mV/nm" ) << " mV/nm: "
		<< "E_01 = " << result(0) << " meV, "
		<< "E_02 = " << result(1) << " meV, "
		<< "E_03 = " << result(2) << " meV" << endl;

	vec reference = {8.95635,23.8527,44.8116};
	cout << "# Reference: "
		<< "E_01 = " << reference(0) << " meV, "
		<< "E_02 = " << reference(1) << " meV, "
		<< "E_03 = " << reference(2) << " meV" << endl;

   	REQUIRE( result(0) == Approx( reference(0) ).epsilon(0.0).margin(0.0001) );
   	REQUIRE( result(1) == Approx( reference(1) ).epsilon(0.0).margin(0.0001) );
   	REQUIRE( result(2) == Approx( reference(2) ).epsilon(0.0).margin(0.0001) );
}
