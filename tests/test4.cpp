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

using namespace std; 
using namespace arma; 

struct Fxc_M1_test : xc::FXC {

	Fxc_M1_test() {
		Mosc=1;
    	cout << "# Number of oscillators: " << Mosc << endl;
	}

	cx_mat get_p(vec rho) {
		cx_double i(0.,1.);
		cx_mat p(rho.n_elem, Mosc);
		p.col(0) = (2.-2.*i)*omega_pl(rho);
		return p;
	}

	// Returns: n^(2/3) * Cm
	cx_mat get_n23Coeffs(cx_mat p, vec n13) {
		mat n23f0finf = xc::get_n23f0finf(n13);
		cx_mat Coeffs(p.n_rows,p.n_cols);
		Coeffs.col(0) = conj(p.col(0)) % (n23f0finf.col(1)-n23f0finf.col(0)) / real(p.col(0));
		return Coeffs;
	}
};

TEST_CASE( "Comparison to qwell.GaAs_mesh(100), ksargs = {10,100,0.1}, args = {0.,100.,0.2,3}", "[tddft]" ) {

	cout << "Test 4: " << endl;
	recommend_num_threads(1);

	QWell_WU_PRL_2005 qwell;
	Mesh<QuantumWell> mesh = qwell.GaAs_mesh(100);
	dft::KsGs ks;
	dft::KsArgs ksargs = {10,100,0.1}; // { number of bands, iterations, smoothness }

	qwell.Efield = qwell.effau.to_au(0.5,"mV/nm");
	ks = dft::Ks(qwell, mesh, ksargs); // Solve KS equation 

	qwell.Efield = qwell.effau.to_au(0.0,"mV/nm");
	tddft::Args args = {0.,100.,0.2,3};

	Fxc_M1_test fxc;

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

	// Comment: these values (obtained with timestep 0.2) are not converged, use for test only
    //vec reference = {-0.0150219,-0.0148807,-0.0141548}; // EXP_ORDER not defined (i.e. full exp is used)
    vec reference = {-0.0150072,-0.0148585,-0.014131}; // EXP_ORDER 3

   	REQUIRE( result(0) == Approx( reference(0) ).epsilon(0.001) );
   	REQUIRE( result(1) == Approx( reference(1) ).epsilon(0.001) );
   	REQUIRE( result(2) == Approx( reference(2) ).epsilon(0.001) );
}
