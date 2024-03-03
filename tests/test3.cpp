// -------------------------------------------------------
// Copyright (C) 2019 by Dmitry R. Gulevich
// You may use or modify this code but not distribute it.
// -------------------------------------------------------
#include <chrono>  // timing
#include <iostream>

#include "catch2/catch.hpp"  // testing framework
#include "tdcdft-omxc/gs.hpp"
#include "tdcdft-omxc/qwmodels.hpp"  // supplementary QuantumWell structs
#include "tdcdft-omxc/tools.hpp"

using namespace std;
using namespace arma;
using namespace tdcdft;

TEST_CASE("Comparison to energy spacing in a double quantum well from the Ullrich-Vignale-PRB-1998", "[dft]") {
    cout << "Test 3: " << endl;
    tools::recommend_num_threads(1);

    DQWell_UV_PRB_1998 dqwell;
    Mesh<QuantumWell> mesh = dqwell.DQWmesh(300);
    KsGs ks;
    KsArgs ksargs = {10, 100, 0.1};  // { number of bands, iterations, smoothness }
    dqwell.Efield = dqwell.effau.to_au(0.0, "mV/nm");
    ks = Ks(dqwell, mesh, ksargs);  // Solve KS equation

    vec result(3);
    for (uword i = 0; i < result.n_elem; ++i) {
        result(i) = dqwell.effau.Eh * (ks.ens(1 + i) - ks.ens(0));
    }
    cout << "# Calculated subband spacings at Efield = " << dqwell.effau.from_au(dqwell.Efield, "mV/nm") << " mV/nm: "
         << "E_01 = " << result(0) << " meV, "
         << "E_02 = " << result(1) << " meV, "
         << "E_03 = " << result(2) << " meV" << endl;

    vec reference = {11.7, 109.2, 154.5};
    cout << "# Reference: "
         << "E_01 = " << reference(0) << " meV, "
         << "E_02 = " << reference(1) << " meV, "
         << "E_03 = " << reference(2) << " meV" << endl;

    REQUIRE(result(0) == Approx(reference(0)).epsilon(0.05));
    REQUIRE(result(1) == Approx(reference(1)).epsilon(0.05));
    REQUIRE(result(2) == Approx(reference(2)).epsilon(0.05));
}
