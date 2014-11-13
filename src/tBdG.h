//
//  tBdG.h
//  tBdG
//
//  Created by Dong Lin on 1/30/14.
//  Copyright (c) 2014 Dong Lin. All rights reserved.
//

#ifndef tBdG_tBdG_h
#define tBdG_tBdG_h

#include "stdcpp.h"
#include "dist.h"
class ctBdG: public cDistribute{
private:
    double _Eb, _hi, _hf, _v, _mu, _dt, _total_t;
    complex<double> _delta;
    int _NK, _NK2;
    double _kc, _Ueff;
    double *_gauss_k, *_gauss_w_k;
    complex<double> myI;
    MatrixXd _bdg_E;
    MatrixXcd _bdg, _bdg_u, _bdg_a, _bdg_b, _bdg_v;
public:
    ctBdG (const int rank, const int size, const int root) : cDistribute(rank,size,root){}
    ~ctBdG(){
    	delete []_gauss_k;
    	delete []_gauss_w_k;
    }
    void input();
    complex<double> DELTA_K(int,int);
    void RK_Propagator(int, int);
    void compute_DeltaK(complex<double>&);
    void quench();
    void Initialize_Euabv();
    void construct_BdG(double, double, double, double);

};





#endif

