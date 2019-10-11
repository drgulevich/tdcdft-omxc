// -------------------------------------------------------
// Copyright (C) 2019 by Dmitry R. Gulevich
// You may use or modify this code but not distribute it.
// -------------------------------------------------------
#include "dft.hpp"
#include "xc.hpp"

namespace dft {

	// M+1 intervals of length dz=W/(M+1); infinite boundaries at m=-1 and m=M; M-1 internal nodes
	vec get_VHartree(Mesh<QuantumWell> &mesh, vec &rho, double ns) {
		double dzdzfactor = -4.*M_PI*mesh.dz*mesh.dz/12.;
		vec VH(mesh.M);
		VH(0) = 0.;
		VH(1) = 2.*VH(0) + dzdzfactor*(rho(1) + 10.*rho(0)); // assuming zero derivative, i.e. VH(0)=VH(-1)=0
		for(uword m=1; m<mesh.M-1; ++m)
			VH(m+1) = 2.*VH(m) - VH(m-1) + dzdzfactor*(rho(m+1) + 10.*rho(m) + rho(m-1));
		VH += 2.*M_PI*ns*mesh.z;
		return VH;
	}
	
	KsGs Ks(QuantumWell &qwell, Mesh<QuantumWell> &mesh, KsArgs args) {

		vec Vext(mesh.M);
		for(uword m=0; m<mesh.M; ++m)
			Vext(m) = qwell.Vext(mesh.z(m));

		vec rho0(mesh.M);
		for(uword m=0; m<mesh.M; ++m)
			rho0(m) = qwell.rho_seed(mesh.z(m));

		vec VHartree = get_VHartree(mesh,rho0,qwell.ns);

		vec VxcLDA = xc::get_VxcLDA(rho0);
		vec Veff = Vext + VHartree + VxcLDA;
	
		KsGs ks;
		ks.rhoarr.push_back(rho0);
	
		vec eigval;
		mat eigvec;
	
		double invdz2=1./(mesh.dz*mesh.dz);
		mat KSmatrix0 = diagmat( invdz2*ones<vec>(mesh.M) );
		KSmatrix0 += diagmat( -0.5*invdz2*ones<vec>(mesh.M-1), 1 );
		KSmatrix0 += diagmat( -0.5*invdz2*ones<vec>(mesh.M-1), -1 );
	
		for(int n=0; n<args.niters; ++n) {
			mat KSmatrix = KSmatrix0 + diagmat(Veff);
			eig_sym( eigval, eigvec, KSmatrix ); // can be improved using sparse matrices

			// Fill the orbitals, starting from the lowest one
			qwell.Nocc=1;
			qwell.EF = M_PI*qwell.ns + eigval(0);
			while(eigval(qwell.Nocc)<qwell.EF) {
				qwell.Nocc+=1;
				qwell.EF = ( M_PI*qwell.ns + sum(eigval.head(qwell.Nocc)) )/qwell.Nocc;
			}

			// Calculate electron density
			vec rho = zeros<vec>(mesh.M);
			for(uword m=0; m<qwell.Nocc; ++m)
				rho += (eigvec.col(m) % eigvec.col(m)) * (qwell.EF-eigval(m))/(M_PI*mesh.dz);  // use conj() if complex

			// Checks:
			// cout << "qwell.Nocc: " << qwell.Nocc << "  qwell.EF (in units Eh): " << qwell.EF  << endl;

			vec VHartree = get_VHartree(mesh,rho,qwell.ns);
			vec cr = pow(rho,1/3.);
			vec VxcLDA = xc::get_VxcLDA(cr);
			vec Veff_new = Vext + VHartree + VxcLDA;
			Veff = (1.-args.theta)*Veff + args.theta*Veff_new;
			ks.rhoarr.push_back(rho);
		}
	
		ks.rho = ks.rhoarr.back();
		ks.Veff = Veff;
		ks.ens = eigval.head(args.nbands);
		ks.orbs = eigvec.head_cols(args.nbands)/sqrt(mesh.dz); // orthogonality: \int \psi_n \psi_m dz = \delta_nm        
		return ks;
	}
	
} // end of namespace dft
