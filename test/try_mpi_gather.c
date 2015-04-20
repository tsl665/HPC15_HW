#include <stdio.h>
#include <unistd.h>
#include <mpi.h>
#include <stdlib.h>


int main( int argc, char *argv[])
{
  int rank;
  int i, N, s, size;
  int *vec,*rec;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  N = 100;
  s = 5;
  MPI_Comm_size(MPI_COMM_WORLD, &size);


  vec = calloc(N, sizeof(int));
  rec = calloc(N, sizeof(int));
  for (i = 0; i<N; i++) {
    vec[i] = i+1;
    rec[i] = 0;
  }

  for (i = 0; i<s; i++) {
    MPI_Gather(&vec[i],1,MPI_INT,&rec[i],s,MPI_INT,0,MPI_COMM_WORLD);
  }

//    MPI_Gather(&vec[i],1,MPI_INT,&rec[size*i],size,MPI_INT,0,MPI_COMM_WORLD);

  if (rank == 0){
    printf("rec = ");
    for (i = 0; i<s*size; i++) {
      printf("%i ",rec[i]);
    }
    printf("; \n");
  }

  free(vec);
  free(rec);
  MPI_Finalize();

  return 0;

}


