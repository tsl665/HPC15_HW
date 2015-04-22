// Iterative Solver for Laplace Equation using Red-Black Gauss-Seidel Method.


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "util.h"

int main (int argc, char **argv)
{
    long  i, n, iter;
    int nthreads;
    double *u, *f;
    double h, res,resInit, tol;
    
    if (argc != 2) {
      fprintf(stderr, "Function needs matrix size as the input argument!\n");
      abort();
    }
    
//    n = 100;
    n = atoi(argv[1]);
    u = (double *) malloc(sizeof(double)*n);
    f = (double *) malloc(sizeof(double)*n);
    h = 1/(n+1);

    for (i = 0; i < n; ++i) {
        f[i] = 1;
        u[i] = 1E-3;
    }
    u[0] = 0;
    u[n-1] = 0;
    

    for (i = 1; i < n - 1; ++i) {
        res = res + fabs(f[i]*h*h - (-u[i-1] + 2*u[i] - u[i+1]));
    }
    resInit = res;
    tol = resInit*1E-6;

    #pragma omp parallel
      #pragma omp master 
      {
        nthreads = omp_get_num_threads();
        printf("Number of threads = %d\n ",nthreads);
      }

    timestamp_type time1, time2;
    get_timestamp(&time1);
    
    iter = 0;
    while (res > tol) {
        // Update Odd Points
        #pragma omp parallel for shared(u,h,f,n) private(i)
        for (i = 1; i < n - 1; i+=2) {
            u[i] = (u[i-1]+u[i+1]+h*h*f[i])/2;
        }
        
        // Update Even Points
        #pragma omp parallel for shared(u,h,f,n) private(i)
        for (i = 2; i < n - 1; i+=2) {
            u[i] = (u[i-1]+u[i+1]+h*h*f[i])/2;
        }
        // Calculate Residual
        res = 0;
        #pragma omp parallel for shared(u,h,f,n) private(i) reduction(+:res)
        for (i = 1; i < n - 1; ++i) {
            res = res + fabs(f[i]*h*h - (-u[i-1] + 2*u[i] - u[i+1]));
        }
        iter = iter + 1;
    }
    
    get_timestamp(&time2);
    double elapsed = timestamp_diff_in_seconds(time1,time2);
    
  {
    FILE* fd = NULL;
    char filename[256];
    snprintf(filename, 256, "RBGS_output%02d.txt", nthreads);
    fd = fopen(filename,"w+");

    if(NULL == fd)
    {
      printf("Error opening file \n");
      return 1;
    }


    fprintf(fd, "========= Red-Black Gauss-Seidel Method ========= \n");
    fprintf(fd, "N =  %ld \n", n);
    fprintf(fd, "Number of threads = %d\n ",nthreads);
    fprintf(fd, "Initial Residual is %e \n", resInit);
    fprintf(fd, "Final Residual is %e \n", res);
    fprintf(fd, "Converge in %ld Iterations. \n", iter);
    fprintf(fd, "Time elapsed is %f seconds.\n", elapsed);

    fclose(fd);
  }
    
    printf("========= Red-Black Gauss-Seidel Method ======== \n");
    printf("N =  %ld \n", n);    
    printf("Initial Residual is %e \n", resInit);
    printf("Final Residual is %e \n", res);
    printf("Converge in %ld Iterations. \n", iter);
    printf("Time elapsed is %f seconds.\n", elapsed);

    free(u);
    free(f);
    
    return 0;
    
}
