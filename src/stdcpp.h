/*! handy standard header files
  stdcpp.h
//
//
//  Created by Dong Lin on 9/2/14.
//  Copyright (c) 2014 Dong Lin. All rights reserved.
*/

#ifndef STDCPP_H_
#define STDCPP_H_
#include <mpi.h>
#include <iostream>
#include <cstdio>
#include <new>
#include <math.h>
#include <iomanip>
#include <fstream>
#include <assert.h>
#include <Eigen/Eigenvalues>
#include <time.h>
using namespace std;
using namespace Eigen;
using namespace MPI;
template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}
struct sPhys {
    double mu;
    complex<double> delta;
};

struct sPara {
    double t; // inverse scattering length, or bound state energy Eb in 2D
    double h; // zeeman field
    double v; // rashba SOC strength
};

struct sGauss {
    int N;
    double kc;
    double *gauss_x;
    double *gauss_w;
};

#endif
