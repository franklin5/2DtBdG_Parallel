//! \file
#include "dist.h"
int cDistribute::compute_count(int rank, int size, int NJOBS){
  int result;
  if (rank != size-1) {
    result = int(NJOBS/size);
  } else {
    result = int(NJOBS/size) + NJOBS % size;
  }
  return result;
}

void cDistribute::distribution(int NJOBS){
  if (_rank == _root){ // send process is only root significant
    sendbuf = new int[NJOBS];
    for(int i = 0; i< NJOBS; ++i){
      sendbuf[i] = i;
    }
    sendcounts = new int[_size];
    displs = new int[_size];
    for(int i=0; i<_size; i++){
      sendcounts[i] = compute_count(i,_size,NJOBS);
      displs[i] = i*int(NJOBS/_size);
    }
  }
  recvcount = compute_count(_rank,_size,NJOBS); // This is a rank dependent variable.
  recvbuf = new int[recvcount]; // So is this array: rank dependent size
  MPI_Scatterv(sendbuf,sendcounts,displs,MPI_INT,recvbuf,recvcount,MPI_INT,_root,COMM_WORLD);
}

void cDistribute::print_rank(){
  int testing_jobs_number = 100;
  distribution(testing_jobs_number);
  for(int ig = 0; ig<_size; ++ig) {
    if (ig ==_rank){
      cout << "Hello! This is rank " << ig << " and I have " << recvcount << " tasks in total." << endl;
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }
}
