// -------------------------------------------------------
// Copyright (C) 2019 by Dmitry R. Gulevich
// You may use or modify this code but not distribute it.
// -------------------------------------------------------
#include <armadillo>
#include "tdcdft-omxc/xc.hpp"
#include "tdcdft-omxc/dft.hpp"
#include "tdcdft-omxc/tddft.hpp"

// Kernel versions:
//#include "tdcdft-omxc/kernels/M0.hpp"
//#include "tdcdft-omxc/kernels/M1_test1.hpp"
//#include "tdcdft-omxc/kernels/M1.hpp"
//#include "tdcdft-omxc/kernels/M1_new.hpp"
//#include "tdcdft-omxc/kernels/M5_new.hpp"

#define ALDAM 1 //  DEBUG
#define EXP_ORDER 3

namespace tddft {

// ksorbs assumed to be real
mat to_ksbasis(Mesh<QuantumWell> &mesh, dft::KsGs &ks, vec &v) {
	uword dim = ks.ens.n_elem;
	mat vmatrix(dim,dim);
	for(uword n1=0;n1<dim;++n1)
		for(uword n2=0;n2<dim;++n2)
			vmatrix(n1,n2) = mesh.dz*accu( ks.orbs.col(n1) % v % ks.orbs.col(n2) );
	return vmatrix;
}

cx_mat grad(Mesh<QuantumWell> &mesh, cx_mat &psi) {
	uword dim = psi.n_rows;
	cx_mat result(dim,psi.n_cols);
	double factor = 1./(2.*mesh.dz);
	for(uword m=0;m<psi.n_cols;++m) {
		result(span(1,dim-2),m) = factor*(psi.col(m).tail(dim-2)-psi.col(m).head(dim-2)); // central differences
		result(0,m) = factor*(-3.*psi(0,m)+4.*psi(1,m)-psi(2,m)); // 2nd order edge
		result(dim-1,m) = -factor*(-3.*psi(dim-1,m)+4.*psi(dim-2,m)-psi(dim-3,m)); // 2nd order edge
	}
	return result;
}

cx_vec grad(Mesh<QuantumWell> &mesh, cx_vec &psi) {
	uword dim = psi.n_elem;
	cx_vec result(dim);
	double factor = 1./(2.*mesh.dz);
	result(span(1,dim-2)) = factor*(psi.tail(dim-2)-psi.head(dim-2)); // central differences
	result(0) = factor*(-3.*psi(0)+4.*psi(1)-psi(2)); // 2nd order edge
	result(dim-1) = -factor*(-3.*psi(dim-1)+4.*psi(dim-2)-psi(dim-3)); // 2nd order edge
	return result;
}

vec grad(Mesh<QuantumWell> &mesh, vec &v) {
	uword dim = v.n_elem;
	vec result(dim);
	double factor = 1./(2.*mesh.dz);
	result(span(1,dim-2)) = factor*(v.tail(dim-2)-v.head(dim-2)); // central differences
	result(0) = factor*(-3.*v(0)+4.*v(1)-v(2)); // 2nd order edge
	result(dim-1) = -factor*(-3.*v(dim-1)+4.*v(dim-2)-v(dim-3)); // 2nd order edge
	return result;
}


cx_mat get_memory(cx_mat &p, cx_mat &previous_memory, vec &previous_gradv, vec &gradv, double dt) {

	uword nrows = previous_memory.n_rows;
	uword ncols = previous_memory.n_cols;

	cx_double i(0.,1.);
	double invdt = 1./dt;

	cx_mat memory(nrows,ncols);
	for(uword ncol=0;ncol<ncols;++ncol)
		for(uword m=0;m<nrows;++m) {

			cx_double invp = 1./p(m,ncol);
			#if EXP_ORDER == 1
				cx_double dep = -i*p(m,ncol)*dt; // 1st order in dt (not recommended)
				cx_double epm = invdt*invp*dep;
				cx_double ep = 1.+dep;
			#elif EXP_ORDER == 2
				cx_double x = -i*p(m,ncol)*dt;
				cx_double dep = x + 0.5*x*x; // 2nd order in dt
				cx_double epm = invdt*invp*dep;
				cx_double ep = 1.+dep;
			#elif EXP_ORDER == 3
				cx_double x = -i*p(m,ncol)*dt;
				cx_double dep = x + 0.5*x*x + (1./9.)*x*x*x; // 3rd order in dt
				cx_double epm = invdt*invp*dep;
				cx_double ep = 1.+dep;
			#else
				cx_double ep = exp(-i*p(m,ncol)*dt);
				cx_double epm = invdt*invp*(ep - 1.);
			#endif
				cx_double alpha = invp * ( i*ep + epm );
				cx_double beta = invp * ( -i - epm );
				memory(m,ncol) = ep * previous_memory(m,ncol) + alpha * previous_gradv(m) + beta * gradv(m);

//			memory(m,ncol) = previous_memory(m,ncol)-i*p(m,ncol)*dt*previous_memory(m,ncol) + dt*previous_gradv(m); // same results as EXP_ORDER=1
//			memory(m,ncol) = previous_memory(m,ncol)-i*p(m,ncol)*dt*previous_memory(m,ncol) + dt*gradv(m); // more stable but overdamped

		}

	return memory;
}


//mat get_VxcM(Mesh<QuantumWell> &mesh, cx_mat &memory, vec &n13, cx_mat &p) {
mat get_VxcM(Mesh<QuantumWell> &mesh, cx_mat &memory, vec &n13, cx_mat &n23Coeffs) {

	vec n23sxc = real(n23Coeffs.col(0) % memory.col(0));
	for(uword m=1;m<n23Coeffs.n_cols;++m)
		n23sxc += real(n23Coeffs.col(m) % memory.col(m));

	vec gradn13 = grad(mesh,n13);
	vec VxcM = -n13 % n23sxc - 3.*mesh.dz*cumsum(n23sxc % gradn13);
	return VxcM;
}


double get_dipole(Mesh<QuantumWell> &mesh, vec &rho) {
	return mesh.dz*accu( mesh.z % rho );
}


Result Tdks(QuantumWell &qwell, Mesh<QuantumWell> &mesh, xc::FXC &fxc, dft::KsGs &ks, Args args) {

#if defined(ALDAM) && ALDAM>0
    cout << "# ALDA+M" << endl;
    cout << "# Number of poles: " << ALDAM << endl;
#elif defined(ALDAM) && ALDAM==0
    cout << "# ALDA" << endl;
#else
    cout << "# Hartree" << endl; 
    cout << "# Warning! Ensure that Hartree is used in static Veff. Continuing." << endl;
#endif

	vec Vext(mesh.M);
	for(uword m=0; m<mesh.M; ++m)
		Vext(m) = qwell.Vext(mesh.z(m));

#if defined(ALDAM) && ALDAM>0
    struct State {
	    State(QuantumWell &qwell, dft::KsGs &ks, uword ncols) {
			gradv = zeros<vec>(ks.orbs.n_rows);
			memory = zeros<cx_mat>(ks.orbs.n_rows,ncols);
	    	Hmatrix = diagmat(ks.ens);
	    	Cpsi = eye<cx_mat>(ks.ens.n_elem,qwell.Nocc);
	    }
		vec gradv;
	    cx_mat memory;
	    cx_mat Cpsi;
	    mat Hmatrix;
    };
	State previous(qwell,ks,ALDAM);
#else
    struct State {
	    State(QuantumWell &qwell, dft::KsGs &ks) {
	    	Hmatrix = diagmat(ks.ens);
	    	Cpsi = eye<cx_mat>(ks.ens.n_elem,qwell.Nocc);
	    }
	    cx_mat Cpsi;
	    mat Hmatrix;
    };
    State previous(qwell,ks);
#endif

	State current = previous;

	std::vector<double> t_array;
	std::vector<double> dipole_array;

	double t=args.tmin;
	vec rho = ks.rho;
	Result result;

	do {

		double d = get_dipole(mesh,rho);
		t_array.push_back(t);
		dipole_array.push_back(d);


		// weights
		vec weights(qwell.Nocc); 
		for(uword m=0; m<qwell.Nocc; ++m)
			weights(m) = (qwell.EF-ks.ens(m))/M_PI;

		for(int nc=0; nc<args.ncorrections; ++nc) {

			mat Hmatrix_midpoint = 0.5*(previous.Hmatrix + current.Hmatrix); // Eq.(4.25) Ullrich book
			cx_mat H_left = eye<mat>(size(current.Hmatrix)) + args.dt*cx_double(0.,0.5)*Hmatrix_midpoint;
			cx_mat H_right = eye<mat>(size(current.Hmatrix)) - args.dt*cx_double(0.,0.5)*Hmatrix_midpoint;
			cx_mat rhs = H_right*previous.Cpsi;
			current.Cpsi = solve(H_left, rhs, solve_opts::fast); 
			cx_mat psi  = ks.orbs * current.Cpsi;

			//rho = qwell.ns*real(conj(psi)%psi);

			// Calculate electron density
			rho = zeros<vec>(mesh.M);
			for(uword m=0; m<qwell.Nocc; ++m)
				rho += weights(m) * real( conj(psi.col(m)) % psi.col(m) );

			vec VHartree = dft::get_VHartree(mesh,rho,qwell.ns);
			vec dV = Vext + VHartree - ks.Veff;
#ifdef ALDAM
			vec n13 = pow(rho,1/3.);
			vec VxcLDA = xc::get_VxcLDA(n13);
			dV += VxcLDA;
#endif
#if defined(ALDAM) && ALDAM>0
			cx_mat gradpsi = grad(mesh,psi);

			// Calculate current
			vec J = zeros(mesh.M);
			for(uword m=0;m<qwell.Nocc;++m)
				J += weights(m)*imag(conj(psi.col(m))%gradpsi.col(m));
			vec vfield = J/rho;

			current.gradv = grad(mesh, vfield);
			cx_mat p = fxc.get_p(rho);
			cx_mat n23Coeffs = fxc.get_n23Coeffs(p, n13);
			current.memory = get_memory(p, previous.memory, previous.gradv, current.gradv, args.dt);
			vec VxcM = get_VxcM(mesh, current.memory, n13, n23Coeffs);

			dV += VxcM;
#endif
			current.Hmatrix = to_ksbasis(mesh, ks, dV);
			current.Hmatrix.diag() += ks.ens;
		}
	
		previous = current;

		t += args.dt;

	} while(t<args.tmax);

	result.dipole = conv_to<vec>::from(dipole_array); 
	result.t = conv_to<vec>::from(t_array); 
	return result;
}

	
} // end of namespace tddft
