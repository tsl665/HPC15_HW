// Iterative Solver for Laplace Equation using Jacobi Method.

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <mpi.h>
#include <time.h>
// #include "util.h"

int main (int argc, char **argv)
{
    long  i,j, nMatrix, iterMax;
    int np, rank;
    double *u, *f,*u_new;

    if (argc != 3) {
      fprintf(stderr, "need two arguments (size of matrix, number of iterations). \n");
      return 1;
    }

    nMatrix = atol(argv[1]);
    iterMax = atol(argv[2]);

    MPI_Init(&argc, &argv);
    MPI_Status status[2];
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD,&np);
    int tag = 99;

    if (nMatrix % np != 0) {
      fprintf(stderr, "The number of processors cannot divide the size of matrix! \n");
      return 2;
    }

    long n = nMatrix / np + 2;
//    printf("length of submatrix = %ld. \n",n - 2);

    u = (double *) malloc(sizeof(double)*n);
    f = (double *) malloc(sizeof(double)*n);
    u_new = (double *) malloc(sizeof(double)*n);
    double h = 1.0/(nMatrix+1.0);
//    printf("h = %f, Rank = %i \n",h,rank);

    for (i = 0; i < n; ++i) {
        f[i] = 1;
        u[i] = 1e-3;
    }
    
    if (rank == 0) {
      u[0] = 0;
    }
    if (rank == np - 1) {
      u[n-1] = 0;
    }
    
    double message_out[2];

    for (i = 0; i < iterMax; i++) {



      for (j = 1; j < n - 1 ; j++) {
        u_new[j] = (f[j]*h*h + u[j-1] + u[j+1])/2;
      }
      message_out[0] = u_new[1];
      message_out[1] = u_new[n-2];


      if (rank != 0) {
        MPI_Send(&message_out[0], 1, MPI_DOUBLE, rank - 1, tag, MPI_COMM_WORLD);
      }
      if (rank != np - 1) {
        MPI_Recv(&u[n-1], 1, MPI_DOUBLE, rank + 1, tag, MPI_COMM_WORLD, &status[0]);
      }




      if (rank != np - 1) {
        MPI_Send(&message_out[1], 1, MPI_DOUBLE, rank + 1, tag, MPI_COMM_WORLD);
      }
      if (rank != 0) {
        MPI_Recv(&u[0], 1, MPI_DOUBLE, rank - 1, tag, MPI_COMM_WORLD, &status[1]);
      }
      
      for (j = 1; j < n - 1; j++) {
        u[j] = u_new[j];
      }
      

    }
    
//    printf("End of loop. Rank = %i \n", rank);

    for (i = 0; i < n; i++) {
      int pt = (rank * (n - 2) + i);
      double x = ((double) (rank * (n - 2) + i))*h;
      double true_sol = 0.5 * x * x - 0.5 * x;

      printf("pt = %i, x = %f, u(x) = %f; p.w.error = %e \n", pt, x, u[i],u[i]-true_sol);
    }


    free(u);
    free(f);
    free(u_new);

//    printf("End of free. Rank = %i \n",rank);

    MPI_Finalize();

//    printf("End of MPI_Finalize(). Rank = %i \n",rank);

    return 0;
}




























