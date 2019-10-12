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

TEST_CASE( "Comparison to the Ullrich-Vignale-PRB-1998 result E01=8.18 meV for single quantum well", "[dft]" ) {

	cout << "Test 2: " << endl;

	QWell_UV_PRB_1998 qwell;
	Mesh<QuantumWell> mesh = qwell.GaAs_mesh(100);
	dft::KsGs ks;
	dft::KsArgs ksargs = {5,100,0.1}; // { number of bands, iterations, smoothness }
	qwell.Efield = qwell.effau.to_au(0.0,"mV/nm");
	ks = dft::Ks(qwell, mesh, ksargs); // Solve KS equation 

	double E01 = qwell.effau.Eh*(ks.ens(1)-ks.ens(0));
	cout << "# Calculated subband spacings at Efield = " << qwell.effau.from_au( qwell.Efield,"mV/nm" ) << " mV/nm: "
		<< "E01 = " << E01 << " meV" << endl;

	double E01_reference=8.18;
	cout << "# Reference (from Ullrich-Vignale-PRB-1998): " << E01_reference << " meV" << endl;

   	REQUIRE( E01 == Approx( E01_reference ).epsilon(0.0).margin(0.02) );
}
