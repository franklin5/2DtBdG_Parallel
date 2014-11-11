/*
 * tBdG.cpp
 *
 *  Created on: Feb 5, 2014
 *      Author: ld7
 */

#include "tBdG.h"

void ctBdG :: input(){
	for (int ig = 0; ig < _size; ++ig) {
		if (ig ==_rank){
		  char dummyname[100];
		  double dummyvalue;
		  int intdummyvalue;
		  FILE *input;
		  input = fopen("input.txt","r");
		  assert(input != NULL);
		  if (ig == _root)  cout << "Starting to read in parameters from file input.txt" << endl;
		  fscanf(input,"%s %lf", dummyname, &dummyvalue);
		  _hi = dummyvalue;    if (ig == _root) cout << dummyname << "=" << _hi << endl;

		  FILE *sf_input;
		  	sf_input = fopen ("superfluid.dat","r"); // Zeeman, Delta, Mu, Eg
		  	double dummyH, dummyD, dummyM, dummyE;
		  	assert (sf_input != NULL);
		  	for (int nsf_input = 0; nsf_input <= int(_hi/0.001); ++nsf_input) {
		  		fscanf(sf_input, "%lf %lf %lf %lf", &dummyH, &dummyD, &dummyM, &dummyE);
		  	}
		  	fclose (sf_input);
		  	_delta.real() = dummyD; _delta.imag() = 0.0;

		  fscanf(input,"%s %lf", dummyname, &dummyvalue);
		_hf = dummyvalue;    if (ig == _root) cout << dummyname << "=" << _hf << endl;
		fscanf(input,"%s %lf", dummyname, &dummyvalue);
		_Eb = dummyvalue;    if (ig == _root) cout << dummyname << "=" << _Eb << endl;
		  fscanf(input,"%s %lf", dummyname, &dummyvalue);
		  _v = dummyvalue;    if (ig == _root) cout << dummyname << "=" << _v << endl;
		  fscanf(input,"%s %d", dummyname, &intdummyvalue);
		  _NK = intdummyvalue;    if (ig == _root) cout << dummyname << "=" << _NK << endl;
		  _NK2 = _NK*_NK;
		fscanf(input,"%s %lf", dummyname, &dummyvalue);
		_kc = dummyvalue;    if (ig == _root) cout << dummyname << "=" << _kc << endl;
		_Ueff = (-8.0*M_PI)/log(1.0+2.0*_kc*_kc/_Eb);

		_bdg_E.resize(4);
		_bdg.resize(4,4);
		_bdg_u.resize(_NK2,4);
		_bdg_a.resize(_NK2,4);
		_bdg_b.resize(_NK2,4);
		_bdg_v.resize(_NK2,4);
		}
	}

}

void ctBdG :: update_Delta(double dt, complex<double>& Delta, double& N0, double& N1){
	Delta.real() = 0.0; Delta.imag() = 0.0;
	N0 = 0.0; N1 =0.0;
    RK_Propagator(dt); // UPDATE wave functions for the next iteration.

    for (int nk = 0; nk < _NK; ++nk) {
//    		result += complex<double> (_gauss_w_k[nk],0.0) * integrand(nk);
    	Delta += _gauss_w_k[nk] * integrand(nk);
	N0    += _gauss_w_k[nk] * N0_CF(nk);
	N1    += _gauss_w_k[nk] * N1_CF(nk);
    }
    _delta = Delta; // UPDATE delta for the next iteration!!!
}

void ctBdG:: RK_Propagator(double dt){
	double normwv;
	complex<double> ncI(0.0,-1.0);
	VectorXcd wvVec(4), k1(4), k2(4), k3(4), k4(4);
	for (int nk = 0; nk < _NK; ++nk) {
		for (int eta = 0; eta < 4; ++eta) {
			wvVec(0) = _bdg_u(nk,eta);
			wvVec(1) = _bdg_a(nk,eta);
			wvVec(2) = _bdg_b(nk,eta);
			wvVec(3) = _bdg_v(nk,eta);

			MatrixXcd bdg(4,4);
			construct_BdG(bdg, _gauss_k[nk], 0.0); // TODO: set mu = 0 temperarily.
			k1 = ncI * (bdg *  wvVec);
			k2 = ncI * (bdg * (wvVec+dt*k1*0.5));
			k3 = ncI * (bdg * (wvVec+dt*k2*0.5));
			k4 = ncI * (bdg * (wvVec+dt*k3));
			wvVec += dt*(k1+2.0*k2+2.0*k3+k4)/6.0;
			normwv = wvVec.norm();
			wvVec /= normwv;
			//		cout << normwv << endl;
			_bdg_u(nk,eta) = wvVec(0);
			_bdg_a(nk,eta) = wvVec(1);
			_bdg_b(nk,eta) = wvVec(2);
			_bdg_v(nk,eta) = wvVec(3);
// _bdg_E(nk,eta) = wvVec.adjoint() * bdg * wvVec; // Don't update Ek(t). Use Ek(0) for all the time!

		}
	}
}


complex<double> ctBdG :: integrand(int nk){

    complex<double> result (0.0,0.0);
    complex<double> u, a, b, v, E;
//    complex<double> prefactor (-_Ueff*0.5/(4.0*M_PI*M_PI)*_gauss_k[nk]*(2*M_PI),0.0);
    double prefactor = -_Ueff/(4.0*M_PI)*_gauss_k[nk];
    for (int eta=0; eta<4; ++eta) {

        E = _bdg_E(nk,eta);
        u = _bdg_u(nk,eta);
        a = _bdg_a(nk,eta);
        b = _bdg_b(nk,eta);
        v = _bdg_v(nk,eta);

        switch (sgn(E.real())) {
            case -1:
                result += prefactor * u * conj(v);
                break;
            case 1:
                result += prefactor * a * conj(b);
                break;
            case 0:
                result += prefactor * (u * conj(v) + a * conj(b))/2.0;
                break;
            default:
                break;
        }
    }
    return result;
}

double ctBdG :: N0_CF(int nk){

    complex<double> result (0.0,0.0);
    complex<double> u, a, b, v, E;
    double prefactor = 1.0/(8.0*M_PI)*_gauss_k[nk];
    for (int eta=0; eta<4; ++eta) {

        E = _bdg_E(nk,eta);
        u = _bdg_u(nk,eta);
        a = _bdg_a(nk,eta);
        b = _bdg_b(nk,eta);
        v = _bdg_v(nk,eta);

        switch (sgn(E.real())) {
            case -1:
                result += u * conj(v);
                break;
            case 1:
                result += a * conj(b);
                break;
            case 0:
                result += (u * conj(v) + a * conj(b))/2.0;
                break;
            default:
                break;
        }
    }
    return prefactor * pow(abs(result),2.0);
}

double ctBdG :: N1_CF(int nk){

    complex<double> result (0.0,0.0);
    complex<double> u, a, b, v, E;
    double prefactor = 1.0/(8.0*M_PI)*_gauss_k[nk];
    for (int eta=0; eta<4; ++eta) {

        E = _bdg_E(nk,eta);
        u = _bdg_u(nk,eta);
        a = _bdg_a(nk,eta);
        b = _bdg_b(nk,eta);
        v = _bdg_v(nk,eta);

        switch (sgn(E.real())) {
            case -1:
                result += u * conj(b);
                break;
            case 1:
                result += a * conj(v);
                break;
            case 0:
                result += (u * conj(b) + a * conj(v))/2.0;
                break;
            default:
                break;
        }
    }
    return prefactor * pow(abs(result),2.0);
}

void ctBdG:: Initialize_Euabv(){

    MatrixXcd bdg(4,4);
    ComplexEigenSolver<MatrixXcd> ces;
    // could be possibly set to SelfAdjointEigenSolver<MatrixXcd> ces;
    for (int nk=0; nk<_NK; ++nk) {
    	construct_BdG(bdg, _gauss_k[nk], _mu);
    	ces.compute(bdg);
    	_bdg_E.row(nk) = ces.eigenvalues();
    	_bdg_u.row(nk) = ces.eigenvectors().row(0);
    	_bdg_a.row(nk) = ces.eigenvectors().row(1);
    	_bdg_b.row(nk) = ces.eigenvectors().row(2);
    	_bdg_v.row(nk) = ces.eigenvectors().row(3);
    }
}


void ctBdG :: construct_BdG(MatrixXcd& bdg, double k, double mu){
	double xi = k*k - mu;

	bdg(0,0) = complex<double> (xi+_h,0.0);
	bdg(0,1) = complex<double> (_v*k,0.0);
	bdg(0,2) = complex<double> (0.0,0.0);
	bdg(0,3) = -_delta;
	bdg(1,0) = complex<double> (_v*k,0.0);
	bdg(1,1) = complex<double> (xi-_h,0.0);
	bdg(1,2) = _delta;
	bdg(1,3) = complex<double> (0.0,0.0);
	bdg(2,0) = complex<double> (0.0,0.0);
	bdg(2,1) = conj(_delta);
	bdg(2,2) = complex<double> (-(xi+_h),0.0);
	bdg(2,3) = complex<double> (_v*k,0.0);
	bdg(3,0) = -conj(_delta);
	bdg(3,1) = complex<double> (0.0,0.0);
	bdg(3,2) = complex<double> (_v*k,0.0);
	bdg(3,3) = complex<double> (-xi+_h,0.0);
}



