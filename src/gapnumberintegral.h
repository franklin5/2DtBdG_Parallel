//
//  gapnumberintegral.h
//  tBdG
//
//  Created by Dong Lin on 1/29/14.
//  Copyright (c) 2014 Dong Lin. All rights reserved.
//

#ifndef tBdG_gapnumberintegral_h
#define tBdG_gapnumberintegral_h
#include "stdcpp.h"

class cSeek_Gap_Number {
private:
	double _myeps;
    double _Eb, _h, _v;
    int _NK;
    double _kc;
    double *_gauss_k, *_gauss_w_k;
    double _Delta, _mu;
public:
    cSeek_Gap_Number(const double Delta, const double Mu, const sPara& para, const sGauss& gauss)
    : _myeps(1e-7), _Eb(para.t), _h(para.h), _v(para.v), _NK(gauss.N), _kc(gauss.kc),
      _gauss_k(gauss.gauss_x), _gauss_w_k(gauss.gauss_w),
      _Delta(Delta), _mu(Mu){}
    ~cSeek_Gap_Number(){}

    void gapnumber();
    double compute(double, double, int);
    double integrand(double, double, double, double, int);
    double integrand_k(double, double, double, int);
    void getresult(double& delta, double& mu, double& Eg)
    {delta = _Delta; mu = _mu; Eg = abs(_h-sqrt(_mu*_mu+_Delta*_Delta));}
    void printEg(){ cout << "Eg = " << abs(_h-sqrt(_mu*_mu+_Delta*_Delta)) <<endl;}
};


#endif
