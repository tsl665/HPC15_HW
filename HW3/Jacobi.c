// Iterative Solver for Laplace Equation using Jacobi Method.

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include "util.h"

int main (int argc, char **argv)
{
    long  i, n, iter;
    int nthreads;
    double *u, *f,*u_new;
    double h, res,resInit, tol;

    n = 100;
    u = (double *) malloc(sizeof(double)*n);
    f = (double *) malloc(sizeof(double)*n);
    u_new = (double *) malloc(sizeof(double)*n);
    h = 1/(n+1);

    for (i = 0; i < n; ++i) {
        f[i] = 1;
        u[i] = 1e-3;
    }
    u[0] = 0;
    u[n-1] = 0;
    

    for (i = 1; i < n - 1; ++i) {
        res = res + fabs(f[i]*h*h - (-u[i-1] + 2*u[i] - u[i+1]));
    }
    resInit = res;
    tol = resInit*1E-6;
    
    timestamp_type time1, time2;
    get_timestamp(&time1);

    #pragma omp parallel
    {
      #pragma omp master
      {
      nthreads = omp_get_num_threads();
      printf("Number of threads = %d\n ",nthreads);
      }
    }
    
    /*
      nthreads = omp_get_num_threads();
      printf("Number of threads = %d\n ",nthreads);
    */
    
    
/*    printf("========= Jacobi Method ========= \n");
    printf("N =  %ld \n", n); 
    printf("================================= \n");
    printf("Iter 0: res = %e \n", res);
    */
    iter = 0;
    while (res > tol) {
        
        #pragma omp parallel for default(none) shared(u_new,u,f,n,h) private(i) 
          for (i = 1; i < n - 1; ++i) {
            u_new[i] = (f[i]*h*h + u[i-1] + u[i+1])/2;
          }
        

        #pragma omp parallel for shared(u_new,u,n) private(i) 
          for (i = 1; i < n - 1; ++i) {
            u[i] = u_new[i];
          }
        
        // Calculate Residual
        res = 0;
        #pragma omp parallel for shared(u,f,n) private(i) \
        reduction(+:res) 
          for (i = 1; i < n - 1; ++i) {
            res += fabs(f[i]*h*h - (-u[i-1] + 2*u[i] - u[i+1]));
          }
        
        iter = iter + 1;
        //printf("Iter %ld: res = %e \n", iter,res);
    }
    get_timestamp(&time2);
    double elapsed = timestamp_diff_in_seconds(time1,time2);
    
    
    printf("========= Jacobi Method ========= \n");
    printf("N =  %ld \n", n);    
    printf("Initial Residual is %e \n", resInit);
    printf("Final Residual is %e \n", res);
    printf("Converge in %ld Iterations. \n", iter);
    printf("Time elapsed is %f seconds.\n", elapsed);
    
    
    free(u);
    free(f);
    free(u_new);
    
    return 0;
    
}
