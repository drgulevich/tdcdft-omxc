// -------------------------------------------------------
// Copyright (C) 2019 by Dmitry R. Gulevich
// You may use or modify this code but not distribute it.
// -------------------------------------------------------
#include <chrono>  // timing
#include <iostream>

#include "tdcdft-omxc/fxcmodels.hpp"  // supplementary fxc structs
#include "tdcdft-omxc/gs.hpp"
#include "tdcdft-omxc/qwmodels.hpp"  // supplementary QuantumWell structs
#include "tdcdft-omxc/td.hpp"
#include "tdcdft-omxc/tools.hpp"
#include "tdcdft-omxc/xc.hpp"

using namespace std;
using namespace arma;
using namespace tdcdft;

#define OUTPUT

struct InfQW : QuantumWell {
    const double eps = 13.0;                         // dielectric constant, from WU-PRL-2005
    const double meff = 0.067;                       // effective mass in units m0, from WU-PRL-2005
    au::Effau effau = au::Effau(eps, meff);          // set up effective atomic units
    const double W = effau.to_au(400., "Angstrom");  // from WU-PRL-2005
    const double Ewell = 257.6 / effau.Eh;           // Quantum well depth, from WU-PRL-2005
    double Efield = effau.to_au(0., "mV/nm");
    InfQW() {
        cout << "# GaAs quantum well with infinite walls" << endl;
        ns = 1.0e11 * effau.a0cm * effau.a0cm;  // Sheet density, from WU-PRL-2005
    }

    // Quantum well potential, including electric field
    double Vext(double z) { return -Efield * z; }

    // Seed electron density for DFT calculations (initial guess)
    double rho_seed(double z) {
        double cosz = cos(M_PI * z / W);
        return 2. * (ns / W) * cosz * cosz;
    }

    // Mesh generator
    Mesh<QuantumWell> mesh(int M) {
        double dz = W / (M + 1);
        Mesh<QuantumWell> mesh(M, dz);
        return mesh;
    }
};

struct Fxc_custom : xc::Omxc {
    Fxc_custom() {
        Mosc = 1;
        cout << "# Number of oscillators: " << Mosc << endl;
    }

    cx_mat get_p(vec rho) {
        cx_double i(0., 1.);
        cx_mat p(rho.n_elem, Mosc);
        p.col(0) = (2. - 2. * i) * xc::omega_pl(rho);
        return p;
    }

    // Returns: n^(2/3) * Cm
    cx_mat get_n23C(cx_mat p, vec n13) {
        mat n23f0finf = xc::get_n23f0finf(n13);
        cx_mat Coeffs(p.n_rows, p.n_cols);
        Coeffs.col(0) = conj(p.col(0)) % (n23f0finf.col(1) - n23f0finf.col(0)) / real(p.col(0));
        return Coeffs;
    }
};

int main() {
    tools::recommend_num_threads(1);

    InfQW qwell;
    Mesh<QuantumWell> mesh = qwell.mesh(100);
    KsGs ks;
    KsArgs ksargs = {10, 100, 0.1};  // Set number of bands, iterations and smoothness (-std=c++17)

    qwell.Efield = qwell.effau.to_au(0.0, "mV/nm");
    ks = Ks(qwell, mesh, ksargs);  // Solve KS equation

    cout << "# Lowest subband spacing at Efield = " << qwell.Efield << ": " << qwell.effau.Eh * (ks.ens(1) - ks.ens(0))
         << ", " << qwell.effau.Eh * (ks.ens(2) - ks.ens(0)) << ", " << qwell.effau.Eh * (ks.ens(3) - ks.ens(0))
         << endl;

    qwell.Efield = qwell.effau.to_au(0.5, "mV/nm");
    ks = Ks(qwell, mesh, ksargs);  // Solve KS equation

    cout << "# Lowest subband spacing at Efield = " << qwell.Efield << ": " << qwell.effau.Eh * (ks.ens(1) - ks.ens(0))
         << ", " << qwell.effau.Eh * (ks.ens(2) - ks.ens(0)) << ", " << qwell.effau.Eh * (ks.ens(3) - ks.ens(0))
         << endl;

    cout << "# Eh: " << qwell.effau.Eh << endl;
    cout << "# ns: " << qwell.ns << endl;

    qwell.Efield = qwell.effau.to_au(0.0, "mV/nm");
    TdArgs args = {0., 1000., 0.05, 3};  // {0., 100., 0.05, 3} gives visually converged results when plotted

    Fxc_custom fxc;
    //	xc::Fxc_ALDA fxc;
    //	xc::Fxc_M1_test fxc;
    //	xc::Fxc_M1_D0 fxc;
    //	xc::Fxc_M1_DQV fxc;

    auto start = std::chrono::high_resolution_clock::now();

    TdDipole dipole = Tdks(qwell, mesh, fxc, ks, args);

    auto finish = std::chrono::high_resolution_clock::now();

#ifdef OUTPUT
    for (uword i = 0; i < dipole.t.n_elem; ++i) cout << dipole.t(i) << " " << dipole.value(i) << endl;
#endif

    std::cout << "# Timing: " << std::chrono::duration_cast<std::chrono::milliseconds>(finish - start).count()
              << " ms\n";

    return 0;
}
