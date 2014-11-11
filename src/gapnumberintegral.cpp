/*
 * gapnumberintegral.cpp
 *
 *  Created on: Feb 5, 2014
 *      Author: ld7
 */
#include "gapnumberintegral.h"
//#include "stdcpp.h"
void cSeek_Gap_Number::gapnumber(){
 double dD = 0.01, dm = 0.01;
 int g = 0; // gap
 int n = 1; // number
 double temp1 = compute(_Delta, _mu, g), temp2 = compute(_Delta, _mu, n);
 double temp = sqrt(pow(temp1,2)+pow(temp2,2));
 double alpha = 0.4;
 double dgd, dgm, dnd, dnm, det;
 	while ( temp > _myeps){
 		// central finite difference
 		// /*
 		dgd = (compute(_Delta+dD, _mu,g) - compute(_Delta-dD, _mu,g) )/(2*dD);
 		dgm = (compute(_Delta, _mu+dm, g) - compute(_Delta, _mu-dm,g) )/(2*dm);
 		dnd = (compute(_Delta+dD, _mu,n) - compute(_Delta-dD, _mu,n) )/(2*dD);
 		dnm = (compute(_Delta, _mu+dm,n) - compute(_Delta, _mu-dm,n) )/(2*dm);
 		det = (dgd*dnm-dgm*dnd);
 		_Delta -=  alpha * ( dnm*temp1-dgm*temp2);
 		_mu 	   -=  alpha * (-dnd*temp1+dgd*temp2);
 		temp1 = compute(_Delta, _mu,g), temp2 = compute(_Delta, _mu,n);
 		temp = sqrt(pow(temp1,2)+pow(temp2,2));
// 		cout.precision(16);
// 		cout << "Delta = " << _Delta << endl;
// 		cout << "mu = " << _mu << endl;
 	}
// 	cout.precision(16);
// 	cout << "Delta = " << _Delta << endl;
// 	cout << "mu = " << _mu << endl;
}

double cSeek_Gap_Number::compute(double delta, double mu, int i){
	double result = 0.0;
	double k;
	double deno ;
	double k1, k2, k3;
	double y1, y2, y3;
	double b1, b2, b3;
	int j;
	switch (i) {
			case 0: //gap
				j = 9;
				break;
			case 1: //number
				j = 7;
				break;
			default:
				break;
	}
	for (int nk = 0; nk < _NK; ++nk) {
		k = _gauss_k[nk];
		result += _gauss_w_k[nk] * integrand(delta,mu,k,i);
	}
	k1=_kc;k2=1.5*_kc;k3=2*_kc;
	deno=(k1*k1-k2*k2)*(k2*k2-k3*k3)*(k3*k3-k1*k1);
	y1=pow(k1,j)*integrand(delta,mu, k1,i);
	y2=pow(k2,j)*integrand(delta,mu, k2,i);
	y3=pow(k3,j)*integrand(delta,mu, k3,i);
	b1=(k3*k3-k2*k2)*y1+(k1*k1-k3*k3)*y2+(k2*k2-k1*k1)*y3;
	b2=(pow(k2,4)-pow(k3,4))*y1+(pow(k3,4)-pow(k1,4))*y2+(pow(k1,4)-pow(k2,4))*y3;
	b3=(k2*k2*pow(k3,4)-pow(k2,4)*k3*k3)*y1+(pow(k1,4)*k3*k3-k1*k1*pow(k3,4))*y2+(k1*k1*pow(k2,4)-pow(k1,4)*k2*k2)*y3;
	switch (i) {
		case 0: //gap
			result += (b1/(4*pow(_kc,4))+b2/(6*pow(_kc,6))+b3/(8*pow(_kc,8)))/deno;
			result = result/ (2*M_PI);
			break;
		case 1: //number
			result += (b1/(2*pow(_kc,2))+b2/(4*pow(_kc,4))+b3/(6*pow(_kc,6)))/deno;
			result = result/ (2*M_PI)-1/(2*M_PI);
			break;
		default:
			break;
	}
	return result;
}


double cSeek_Gap_Number::integrand(double delta, double mu,
		double k, int i){
	double result;
	double xi = k*k - mu;
	double temp = _h*_h+_v*_v*k*k,
		  temp1 = sqrt(temp*xi*xi+_h*_h*delta*delta);
	double ekp = sqrt(xi*xi+delta*delta+ temp + 2* temp1),
		   ekm = sqrt(xi*xi+delta*delta+ temp - 2* temp1);
	switch (i) {
		case 0: //gap
			result = k*(1.0/(2*k*k+_Eb)-1.0/4*( (1+_h*_h/temp1)/ekp + (1-_h*_h/temp1)/ekm ));
			break;
		case 1: //number
			result = k*(1-xi/2*( (1+temp/temp1)/ekp + (1-temp/temp1)/ekm ));
			break;
		default:
			break;
	}
	return result;
}



