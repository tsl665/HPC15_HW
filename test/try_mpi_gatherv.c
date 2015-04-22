#include <stdio.h>
#include <unistd.h>
#include <mpi.h>
#include <stdlib.h>


int main( int argc, char *argv[])
{
  int rank;
  int i, N, s, P;
  int *vec,*rec,*len,*displs;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  N = 100;
  s = rank+1;
  MPI_Comm_size(MPI_COMM_WORLD, &P);


  vec = calloc(N, sizeof(int));
  rec = calloc(N, sizeof(int));
  len = calloc(P, sizeof(int));
  displs = calloc(P, sizeof(int));
  for (i = 0; i<N; i++) {
    vec[i] = i+1;
    rec[i] = -1;
  }
/*
  for (i = 0; i<1; i++) {
    MPI_Gather(&vec[i],2,MPI_INT,&rec[i],size,MPI_INT,0,MPI_COMM_WORLD);
  }
*/

//    MPI_Gather(&vec[0],s,MPI_INT,&rec[0],s,MPI_INT,0,MPI_COMM_WORLD);
  if (rank == 0){
    printf("rec = ");
    for (i = 0; i<(P+1)*P/2; i++) {
      printf("%i ",rec[i]);
    }
    printf("; \n");
  }

    MPI_Gather(&s,1,MPI_INT,len,1,MPI_INT,0,MPI_COMM_WORLD);
    
    displs[0] = 0;
    for (i = 1; i<P; i++) {
      displs[i] = displs[i-1] + len[i-1];
    }

    MPI_Gatherv(&vec[0],s,MPI_INT,rec,len,displs,MPI_INT,0,MPI_COMM_WORLD);
    printf("rank %i is finished. \n",rank);


  if (rank == 0){
    printf("len = ");
    for (i = 0; i<P; i++) {
      printf("%i ", len[i]);
    }
    printf("; \n");

    printf("displs = ");
    for (i = 0; i<P; i++) {
      printf("%i ", displs[i]);
    }
    printf("; \n");


    printf("rec = ");
    for (i = 0; i<(P+1)*(P+2)/2; i++) {
      printf("%i ",rec[i]);
    }
    printf("; \n");
  }


//  s = (5+4)/2;

//  printf("%i \n",s);


  free(vec);
  free(rec);
  MPI_Finalize();

  return 0;

}


