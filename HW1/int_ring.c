/* Communication ping-pong:
 * Exchange between messages between mpirank
 * 0 <-> 1, 2 <-> 3, ....
 */
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <mpi.h>

int main( int argc, char *argv[])
{
  int rank, tag, origin, destination, np, N, i;

  if (argc != 2) {
      fprintf(stderr, "need an argument (number of summands)\n");
      abort();
    }


  MPI_Status status;

  char hostname[1024];
  gethostname(hostname, 1024);

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  MPI_Comm_size(MPI_COMM_WORLD,&np);
  N = atoi(argv[1]);

  int message_out = rank;
  int message_in = -1;
  tag = 99;

  for (i = 0; i < N; i++) {
    if (rank != np - 1) {
      destination = rank + 1;
    }
    else {
      destination = 0;
    }

    MPI_Send(&message_out, 1, MPI_INT, destination, tag, MPI_COMM_WORLD);

    if (rank != 0) {
      origin = rank - 1;
    }
    else {
      origin = np - 1;
    }
    
    MPI_Recv(&message_in, 1, MPI_INT, origin, tag, MPI_COMM_WORLD, &status);

//  printf("rank %d hosted on %s received from %d the message %d\n", rank, hostname, origin, message_in);    
  printf("rank %d received the message %d in loop %i \n", rank, message_in, i); 
  }

  MPI_Finalize();
  return 0;
}



//  int np;
//  np = atoi(argv[1]);
//  printf("The number of processors is %i. \n",np);

//  printf("%i \n",argc);
//  printf("%s \n",argv[0]);
//  
//
//  printf("The number of processors is %i. \n",np);
/*
  if(rank % 2 == 0)
  {
    destination = rank + 1;
    origin = rank + 1;

    MPI_Send(&message_out, 1, MPI_INT, destination, tag, MPI_COMM_WORLD);
    MPI_Recv(&message_in,  1, MPI_INT, origin,      tag, MPI_COMM_WORLD, &status);
  }
  else
  {
    destination = rank - 1;
    origin = rank - 1;

    MPI_Recv(&message_in,  1, MPI_INT, origin,      tag, MPI_COMM_WORLD, &status);
    MPI_Send(&message_out, 1, MPI_INT, destination, tag, MPI_COMM_WORLD);
  }

  printf("rank %d hosted on %s received from %d the message %d\n", rank, hostname, origin, message_in);

  MPI_Finalize();
  return 0;
}
*/
