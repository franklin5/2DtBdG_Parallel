/*
 * main.cpp
 *
 *  Created on: Nov 11, 2014
 *      Author: ld7
 */
#include "stdcpp.h"
#include "grst.h"
#include "lgwt.h"
#include "tBdG.h"
#define root 0
#define char_length	100
int main(int argc, char** argv){
  Init(argc,argv);
  int rank, size;
  rank = COMM_WORLD.Get_rank();
  size = COMM_WORLD.Get_size();
  /******************************I/O for main.cpp************************************/
  if (rank == root) {
    cout << "======================================================================\n"
	 << "The purpose of this program is to study\n"
	 << "the quench dynamics of 2D Fermi gas with spin-orbit coupling.\n"
	 << "Motivated, proposed, designed, implemented and researched \n"
	 << "by Lin Dong at Rice University. \n"
	 << "at " << __TIME__ << ", on " << __DATE__ << endl
	 << "MPI is initialized and program starts from " << __FILE__ << endl
	 << "======================================================================\n";
  }

  int grst_flag = 0;
  switch (grst_flag) {
	case 1:
		if (rank==root) {
			// *********************************** //
			// --> grst solution from gap and number equation
			// *********************************** //
			cGRST grst; // zeeman field is up to 2.416
			grst.compute();
		}
		break;
	case 0:
		ctBdG quench(rank,size,root);
		quench.input();
		quench.Initialize_Euabv();
		quench.quench();
		quench.output();
		break;
	default:
		break;
  }
  Finalize();
  return 0;
}
