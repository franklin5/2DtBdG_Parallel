/*
 * lgwt.cpp
 *
 *  Created on: Feb 6, 2014
 *      Author: ld7
 */
#include "lgwt.h"
double maxfunabs(const int N, const double *y, const double *y0){
    double z[N];
    for (int i=0; i<N; i++) z[i]=fabs(y[i]-y0[i]);
    return *std::max_element(z, z+N);
}


void gauss_lgwt(int N, const double a, const double b, double x[], double w[]){
    N = N-1;
    int N1(N+1), N2(N+2);
    double xu[N1],y[N1];
    for (int i=0; i<N1-1; i++) {
        xu[i]=-1+2/(floor(N1)-1)*i;
        y[i]=cos((2*i+1)*M_PI/(2*N+2))+(0.27/N1)*sin(M_PI*xu[i]*N/N2);

    }
    xu[N1-1]=1;
    y[N1-1] =cos((2*N+1)*M_PI/(2*N+2))+(0.27/N1)*sin(M_PI*xu[N]*N/N2);
//    double L[N1][N2]; ====>>>> leads to segmentation fault.
    double *L = new double[N1*N2];
    double y0[N1],Lp[N1];
    double flag=1;
    while (flag>1e-15) {
        for (int i=0; i<N1; i++) {
//            L[i][0]=1;
//            L[i][1]=y[i];
        	L[N2*i] = 1;
        	L[N2*i+1] = y[i];
        }
        int ik;
        for (int k=2; k<N1+1; k++) {
            ik=k-1;
            for (int i=0;i<N1;i++){
//                L[i][ik+1]=((2*k-1)*y[i]*L[i][ik]-(k-1)*L[i][ik-1])/k;
            	L[i*N2+ik+1]=((2*k-1)*y[i]*L[i*N2+ik]-(k-1)*L[i*N2+ik-1])/k;
            }
        }

        for (int i=0; i<N1; i++) {
//            Lp[i]=N2*(L[i][N1-1]-y[i]*L[i][N2-1])/(1-y[i]*y[i]);
        	Lp[i]=N2*(L[i*N2+N1-1]-y[i]*L[i*N2+N2-1])/(1-y[i]*y[i]);
        }
        for (int i=0;i<N1;i++) {
            y0[i]=y[i];
//            y[i]=y0[i]-L[i][N2-1]/Lp[i];
            y[i]=y0[i]-L[i*N2+N2-1]/Lp[i];
        }

        flag=maxfunabs(N,y,y0);
    }
    for (int i=0;i<N1;i++) {
        x[i]=(a*(1-y[i])+b*(1+y[i]))/2;
        //cout << setprecision(16) << x[i]<< '\t';
        w[i]=(b-a)/((1-y[i]*y[i])*Lp[i]*Lp[i])*N2*N2/N1/N1;
        //cout << setprecision(16) << w[i] << '\n';
    }

}

