/* Parallel sample sort
 */
#include <stdio.h>
#include <unistd.h>
#include <mpi.h>
#include <stdlib.h>


static int compare(const void *a, const void *b)
{
  int *da = (int *)a;
  int *db = (int *)b;

  if (*da > *db)
    return 1;
  else if (*da < *db)
    return -1;
  else
    return 0;
}

static int binsearch(int *a, int N, int b, int *ind)
{
  int iL, iR, iM;
  iL = 0;
  iR = N-1;
  *ind = -1;

  while (iL <= iR) {
    iM = iL + (iR - iL)/2;
    if (a[iM] < b) {
      iL = iM + 1;
    }
    else if (a[iM] > b) {
      iR = iM - 1;
    }
    else {
      *ind = iM;
      break;
    }
  }

  if (*ind == -1) {
    *ind = iL;
  }

  return 0;
}




int main( int argc, char *argv[])
{
  int rank;
  int i, N, S, P,*iSpl,iSold;
  int *vec,*rec,*spl;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  /* Number of random numbers per processor (this should be increased
   * for actual tests or could be passed in through the command line */
  N = 100;
  S = 10;
  MPI_Comm_size(MPI_COMM_WORLD, &P);


  vec = calloc(N, sizeof(int));
  rec = calloc(N, sizeof(int));
  spl = calloc(P-1,sizeof(int));
  iSpl = calloc(P-1,sizeof(int));
  /* seed random number generator differently on every core */
  srand((unsigned int) (rank + 393919));

  /* fill vector with random integers */
  for (i = 0; i < N; ++i) {
    vec[i] = rand();
    rec[i] = 0;
  }
  printf("rank: %d, first entry: %d\n", rank, vec[0]);

  MPI_Gather(&vec[0],S,MPI_INT,&rec[0],S,MPI_INT,0,MPI_COMM_WORLD);

  /* sort locally */
  qsort(vec, N, sizeof(int), compare);

/*
  printf("vec = ");
  for (i = 0; i<N; i++){
    printf("%i, ",vec[i]);
  }
  printf("\n");

  printf("S = %i, ind(S) = %i \n",S,P);

*/


  if (rank == 0) {
    qsort(rec, S*P, sizeof(int), compare);
    for (i = 0; i < P-1; i++) {
      spl[i] = rec[i*S];
    }
    MPI_Bcast(&spl[0],P-1,MPI_INT,0,MPI_COMM_WORLD);
  }

  iSold = 0;

/*

  k = 0;
  for (i = 0; i < P-1; i++) {
    // find the i-th splitter and send the length to i-th core.
    binsearch(vec, N, spl[i], &iSnew);
    len = iSnew - iSold + 1;
    MPI_Gather(&len,1,MPI_INT,&vecLen[0],1,MPI_INT,i,MPI_COMM_WORLD);
    if (rank != i) {
      // Not at i-th core: just send the piece of data need by i-th
      // core.
      MPI_Isend(&vec[iSold],len,MPI_INT,i,rank,MPI_COMM_WORLD,&req);
      iSold = iSnew;
    } else {
      //At i-th core: recieve data from each other core; copy data
      //from itself.
      for (j = 0; j < i; j++) {
        MPI_Irecv(&rec[k],vecLen[j],MPI_INT,j,j,MPI_COMM_WORLD,&req);
        k = k + vecLen[j];
      }

      for (j = 0; j < vecLen[i]; j++) {
        rec[k+j] = 


*/











  /* randomly sample s entries from vector or select local splitters,
   * i.e., every N/P-th entry of the sorted vector */

  /* every processor communicates the selected entries
   * to the root processor; use for instance an MPI_Gather */

  /* root processor does a sort, determinates splitters that
   * split the data into P buckets of approximately the same size */

  /* root process broadcasts splitters */

  /* every processor uses the obtained splitters to decide
   * which integers need to be sent to which other processor (local bins) */

  /* send and receive: either you use MPI_AlltoallV, or
   * (and that might be easier), use an MPI_Alltoall to share
   * with every processor how many integers it should expect,
   * and then use MPI_Send and MPI_Recv to exchange the data */

  /* do a local sort */

  /* every processor writes its result to a file */

  free(vec);
  MPI_Finalize();
  return 0;
}
