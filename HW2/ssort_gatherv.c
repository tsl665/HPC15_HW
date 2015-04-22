/* Parallel sample sort
 */
#include <stdio.h>
#include <unistd.h>
#include <mpi.h>
#include <stdlib.h>
#include <math.h>


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
  int i,j, N,M, S, P,iSold,iSnew,Nrec,len;
  int *vec,*rec,*spl,*vecLen,*displs;
//  int rangeMax = 100;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  /* Number of random numbers per processor (this should be increased
   * for actual tests or could be passed in through the command line */

  if (argc != 2) {
    fprintf(stderr, "need one argument (size of vector on each core). \n");
    abort();
  }

  N = atoi(argv[1]);
  M = N*2;
  S = sqrt(N);
  MPI_Comm_size(MPI_COMM_WORLD, &P);


  vec = calloc(N, sizeof(int));
  rec = calloc(M, sizeof(int));
  spl = calloc(P-1,sizeof(int));
  vecLen = calloc(P,sizeof(int));
  displs = calloc(P,sizeof(int));
  /* seed random number generator differently on every core */
  srand((unsigned int) (rank + 393919));

  /* fill vector with random integers */
  for (i = 0; i < N; ++i) {
//    vec[i] = rangeMax*cos(rand())/2+rangeMax/2;
    vec[i] = rand();
//    rec[i] = -1;
  }

// Take the first S elements as random samples
  MPI_Gather(vec,S,MPI_INT,&rec[0],S,MPI_INT,0,MPI_COMM_WORLD);


  /* sort locally */
  qsort(vec, N, sizeof(int), compare);


/*
  for (i = 0; i < P; i++) {
    if (rank == i) {
      printf("Rank: %d: vec = ", rank);
      for (j = 0; j<N; j++) {
        printf("%d ",vec[j]);
      }
      printf("\n");
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }

  printf("\n");
  */

//  Compute splitter in root node and broadcast to everyone
  if (rank == 0) {
    qsort(rec, S*P, sizeof(int), compare);
    for (i = 0; i < P-1; i++) {
      spl[i] = rec[(i+1)*S];
    }
  }

  MPI_Bcast(spl,P-1,MPI_INT,0,MPI_COMM_WORLD);
/*
  if (rank == 0) {
//    printf("rec = ");
//    for (j = 0; j < S*P; j++) {
//      printf("%d ",rec[j]);
//    }
    printf("\nspl = ");
    for (j = 0; j < P-1; j++) {
      printf("%d ",spl[j]);
    }
    printf("\n");
  }
*/
/*
  for (i = 0; i < P; i++) {
    if (rank == i) {
      printf("Rank: %d: spl = ", rank);
      for (j = 0; j<P-1; j++) {
        printf("%d ",spl[j]);
      }
      printf("\n");
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }
*/


  //  Use MPI_Gatherv to collect different pieces of vectors.
  iSold = 0;
  for (i = 0; i < P-1; i++) {
    // find the i-th splitter and send the length to i-th core.
    binsearch(vec, N, spl[i], &iSnew);
    len = iSnew - iSold + 1;
    MPI_Gather(&len,1,MPI_INT,&vecLen[0],1,MPI_INT,i,MPI_COMM_WORLD);

    // Compute the displacement for MPI_Gatherv
    if (rank == i) {
      displs[0] = 0;
      for (j = 1; j < P; j++) {
        displs[j] = displs[j-1] + vecLen[j-1];
      }
      Nrec = displs[P-1] + vecLen[P-1];
    }

    MPI_Gatherv(&vec[iSold],len,MPI_INT,&rec[0],vecLen,displs,MPI_INT,i,MPI_COMM_WORLD);

    iSold = iSnew + 1;
//    printf("rank = %d, loop = %d: iSold = %d / %d \n",rank,i,iSold,N);

  }

// Put the rest in the last core
  len = N - iSold;
//  printf("rank = %d: len = %d \n",rank,len);
  MPI_Gather(&len,1,MPI_INT,&vecLen[0],1,MPI_INT,P-1,MPI_COMM_WORLD);
  if (rank == P-1) {
    displs[0] = 0;
    for (j = 1; j < P; j++) {
      displs[j] = displs[j-1] + vecLen[j-1];
    }
    Nrec = displs[P-1] + vecLen[P-1];
  }

  MPI_Gatherv(&vec[iSold],len,MPI_INT,&rec[0],vecLen,displs,MPI_INT,P-1,MPI_COMM_WORLD);


// Local sort.
  qsort(rec,Nrec,sizeof(int),compare);

// Print result to screen
  printf("Rank: %d: Length = %d; \n", rank, Nrec);

//
//
/*
  for (i = 0; i < P; i++) {
    if (rank == i) {
      printf("Rank: %d: Length = %d; \n", rank, Nrec);
      
      //printf("rec = ")
      //for (j = 0; j<Nrec; j++) {
      //  printf("%d ",rec[j]);
      //}
      //printf("\n\n");
      
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }
*/
//  Output to a file
  {
    FILE* fd = NULL;
    char filename[256];
    snprintf(filename, 256, "output%02d.txt", rank);
    fd = fopen(filename,"w+");

    if(NULL == fd)
    {
      printf("Error opening file \n");
      return 1;
    }

    fprintf(fd, "rank %d stores vector of length %d. \n", rank, Nrec);
    for(j = 0; j < Nrec; j++)
      fprintf(fd, "%d  ", rec[j]);

    fclose(fd);
  }

  free(vec);
  free(rec);
  free(spl);
  free(vecLen);
  free(displs);
  MPI_Finalize();
  return 0;
}
