/*
 * grst.cpp
 *
 *  Created on: Feb 6, 2014
 *      Author: ld7
 */
#include "grst.h"
void cGRST :: compute(){
	ofstream superfluid_output;
	superfluid_output.open("superfluid.dat");
	assert(superfluid_output.is_open());

	int    NK = 500;    double kc = 10;
	double *gauss_k = new double [NK];double *gauss_w_k = new double [NK];
	gauss_lgwt(NK,-kc,kc,gauss_k,gauss_w_k);
	sGauss gauss;
	gauss.N=NK; gauss.kc=kc;
	gauss.gauss_x = gauss_k;gauss.gauss_w = gauss_w_k;
	double mu0= 0.05;    double delta0 = 0.05; double Eg;
	double Zeeman = 0.0;
	do {

		// *********************************** //
		// --> input parameters of the system.
		sPara para = {0.2, Zeeman, 1.2}; // E_b, Zeeman, SOC
		// *********************************** //
		cSeek_Gap_Number GapNumber(delta0, mu0 ,para,gauss);
		GapNumber.gapnumber();
		GapNumber.getresult(delta0, mu0, Eg);
		cout.precision(16);
		cout << "Zeeman = " << Zeeman << endl;
		cout << "Delta = " << delta0 << endl;
		cout << "Mu = "    << mu0    << endl;
		GapNumber.printEg();
		superfluid_output.precision(16);
		superfluid_output << Zeeman << '\t' << delta0 << '\t'
				<< mu0 << '\t' << Eg << endl;
		Zeeman += 0.001;

	} while(Zeeman < 2.5);
	delete []gauss_k;
	delete []gauss_w_k;
	superfluid_output.close();
}


