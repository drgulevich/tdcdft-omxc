#ifndef AU_HEADER
#define AU_HEADER
// -------------------------------------------------------
// Copyright (C) 2019 by Dmitry R. Gulevich
// You may use or modify this code but not distribute it.
// -------------------------------------------------------
#include <cmath>

// Physical constants from NIST (https://physics.nist.gov)
namespace au {

	const double BOHR=0.52917721067; // Bohr radius in Angstrom
	const double HARTREE=27.21138602; // eV
	const double AUOT=2.418884326509e-17; // atomic unit of time, seconds
//	const double eps0=8.854187817e-12; // electric constant, F/m

	// Effective atomic units 
	struct Effau {
		Effau(double eps, double meff) : eps(eps), meff(meff) {}
		const double eps;
		const double meff;
		const double eeff=1./std::sqrt(eps); // effective charge in units of e
		const double a0=BOHR/(meff*eeff*eeff); // effective Bohr radius in Angstrom
		const double a0cm=a0*1.e-8; // effective Bohr radius in cm
		const double Eh=1.e3*HARTREE*meff*std::pow(eeff,4); // effective Hartree energy in meV
		const double Eau=1.e-3*Eh/(a0*1.e-10); // V/m, effective atomic unit of electric field by Eh/e0a0 (note e not eeff!)
		const double EmVnm=1.e-6*Eau; // mV/nm
		const double utfs=1.e15*AUOT*HARTREE/(Eh*1.e-3); // unit time in fs

		// convert from physical units to effective au
		double to_au( double value, std::string units ) {
			if(units=="meV")
				return value/Eh;
			else if(units=="Angstrom")
				return value/a0;
			else if(units=="cm")
				return value/a0cm;
			else if(units=="mV/nm")
				return value/EmVnm;
			else {
				// Through an exception here <--- !!!
				std::cout << "Error. Supplied units are not implemented." << std::endl;
				return 0.;
			}
		}

		// convert from effective au to physical units
		double from_au( double value, std::string units ) {
			if(units=="meV")
				return value*Eh;
			else if(units=="Angstrom")
				return value*a0;
			else if(units=="cm")
				return value*a0cm;
			else if(units=="mV/nm")
				return value*EmVnm;
			else {
				// Through an exception here <--- !!!
				std::cout << "Error. Supplied units are not implemented." << std::endl;
				return 0.;
			}
		}

	};

} // end of au namespace

#endif // AU_HEADER
