// -------------------------------------------------------
// Copyright (C) 2019 by Dmitry R. Gulevich
// You may use or modify this code but not distribute it.
// -------------------------------------------------------
#include <iostream>
#include <chrono>  // timing
#include "tdcdft-omxc/qwmodels.hpp" // supplementary QuantumWell structs
#include "tdcdft-omxc/dft.hpp"
#include "tdcdft-omxc/tddft.hpp"
#include "tdcdft-omxc/tools.hpp"

// Important note: for small problem like this on it is best to execute 
// in serial by switching off OpenMP:
// 		$ export OMP_NUM_THREADS=1
// Check:
// 		$ echo $OMP_NUM_THREADS

using namespace std; 
using namespace arma; 

#define OUTPUT

int main() {

	recommend_num_threads(1);

	QWell_WU_PRL_2008 qwell;
	Mesh<QuantumWell> mesh = qwell.GaAs_mesh(100);
	dft::KsGs ks;
	dft::KsArgs ksargs = {10,100,0.1}; // Set number of bands, iterations and smoothness (-std=c++17)

	qwell.Efield = qwell.effau.to_au(0.01,"mV/nm");
	ks = dft::Ks(qwell, mesh, ksargs); // Solve KS equation

	cout << "# Energies at Efield = " << qwell.Efield << ": "
		<< qwell.effau.Eh*ks.ens(0) << ", "
		<< qwell.effau.Eh*ks.ens(1) << ", "
		<< qwell.effau.Eh*ks.ens(2) << endl;


/*	cout << "# Lowest subband spacing at Efield = " << qwell.Efield << ": "
		<< qwell.effau.Eh*(ks.ens(1)-ks.ens(0)) << ", "
		<< qwell.effau.Eh*(ks.ens(2)-ks.ens(0)) << ", "
		<< qwell.effau.Eh*(ks.ens(3)-ks.ens(0)) << endl;
*/

	cout << "# qwell.Nocc: " << qwell.Nocc << endl;
	cout << "# qwell.EF: " << qwell.effau.Eh*qwell.EF << " meV" << endl;
	cout << "# qwell.ns: " << qwell.ns << endl;

	qwell.Efield = qwell.effau.to_au(0.0,"mV/nm");
	tddft::Args args = {0., 50., 0.05, 3}; // {0., 100., 0.05, 3} gives visually converged results when plotted

	auto start = std::chrono::high_resolution_clock::now();

	tddft::Result tdks = tddft::Tdks(qwell, mesh, ks, args);

	auto finish = std::chrono::high_resolution_clock::now();

#ifdef OUTPUT
	for(uword i=0;i<tdks.t.n_elem;++i)
		cout << tdks.t(i) << " " << tdks.dipole(i) << endl;
#endif

    std::cout << "# Timing: "
              << std::chrono::duration_cast<std::chrono::milliseconds>(finish-start).count()
              << " ms\n";

	return 0;
}
