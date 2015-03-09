#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <mpi.h>
#include <time.h>
//#include "util.h"

int main( int argc, char *argv[])
{
  int rank, tag, origin, destination, np, N, i;
//  timestamp_type time1, time2;
  clock_t tstart, tend;

  if (argc != 2) {
      fprintf(stderr, "need an argument (times the message go around). \n");
      abort();
    }


  MPI_Status status;

  char hostname[1024];
  gethostname(hostname, 1024);

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  MPI_Comm_size(MPI_COMM_WORLD,&np);
  N = atoi(argv[1]);

  int *message_out;
  int *message_in;
  unsigned int L = 2e6 / sizeof(int);
  message_out = (int *) malloc(2e6);
  message_in = (int *) malloc(2e6);
  tag = 99;
  
//  get_timestamp(&time1);

  tstart = clock();

  for (i = 0; i < N; i++) {

//  Compute the origin and destination

    if (rank != np - 1) {
      destination = rank + 1;
    }
    else {
      destination = 0;
    }
    if (rank != 0) {
      origin = rank - 1;
    }
    else {
      origin = np - 1;
    }


//  Passing message


    if ( (i != 0) || (rank != 0) ) {
      MPI_Recv(message_in, L, MPI_INT, origin, tag, MPI_COMM_WORLD, &status);
      message_out[0] = message_in[0] + rank;
      printf("rank %d hosted on %s received from %d the message %d\n", rank, hostname, origin, message_in[0]);
    }
    else {
      message_out[0] = 0;
      printf("Start with initial value %i \n", message_out[0]);

    }

    if ( (i != N - 1) || (rank != np - 1) ) {

      MPI_Send(message_out, L, MPI_INT, destination, tag, MPI_COMM_WORLD);
    }

  }

  if (rank == np - 1) {
//    get_timestamp(&time2);
//    double tepc = timestamp_diff_in_seconds(time1,time2)/(N*np-1);
//
    tend = clock();
    double elapsed = ((double) (tend - tstart)) / CLOCKS_PER_SEC;
    double tepc = elapsed / (N*np-1);
    printf("Time elapsed per communication is %e seconds. \n",tepc);
    }

  MPI_Finalize();
  return 0;
}
