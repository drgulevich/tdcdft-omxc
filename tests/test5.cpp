// -------------------------------------------------------
// Copyright (C) 2019 by Dmitry R. Gulevich
// You may use or modify this code but not distribute it.
// -------------------------------------------------------
#include <chrono>  // timing
#include <iostream>

#include "catch2/catch.hpp"  // testing framework
#include "tdcdft-omxc/fxcmodels.hpp"
#include "tdcdft-omxc/gs.hpp"
#include "tdcdft-omxc/qwmodels.hpp"  // supplementary QuantumWell structs
#include "tdcdft-omxc/td.hpp"
#include "tdcdft-omxc/tools.hpp"
#include "tdcdft-omxc/xc.hpp"

using namespace std;
using namespace arma;
using namespace tdcdft;

TEST_CASE("Comparison to qwell.GaAs_mesh(100), ksargs = {10,100,0.1}, args = {0.,100.,0.2,3}", "[tddft]") {
    cout << "Test 5: " << endl;
    tools::recommend_num_threads(1);

    QWell_WU_PRL_2005 qwell;
    Mesh<QuantumWell> mesh = qwell.GaAs_mesh(100);
    KsGs ks;
    KsArgs ksargs = {10, 100, 0.1};  // { number of bands, iterations, smoothness }

    qwell.Efield = qwell.effau.to_au(0.5, "mV/nm");
    ks = Ks(qwell, mesh, ksargs);  // Solve KS equation

    qwell.Efield = qwell.effau.to_au(0.0, "mV/nm");
    TdArgs args = {0., 100., 0.2, 3};

    xc::Fxc_M1_test fxc;

    auto start = std::chrono::high_resolution_clock::now();

    TdDipole dipole = Tdks(qwell, mesh, fxc, ks, args);

    auto finish = std::chrono::high_resolution_clock::now();

    std::cout << "# Timing: " << std::chrono::duration_cast<std::chrono::milliseconds>(finish - start).count()
              << " ms\n";

    // Comment: these values (obtained with timestep 0.2) are not converged, use for test only
    // vec reference = {-0.0150219,-0.0148807,-0.0141548}; // EXP_ORDER not defined (i.e. full exp is used)
    vec reference = {-0.0150072, -0.0148585, -0.014131};  // EXP_ORDER 3

    uword nvals = dipole.t.n_elem;
    REQUIRE(dipole.value(nvals - 1) == Approx(reference(0)).epsilon(0.001));
    REQUIRE(dipole.value(nvals - 2) == Approx(reference(1)).epsilon(0.001));
    REQUIRE(dipole.value(nvals - 3) == Approx(reference(2)).epsilon(0.001));
}
