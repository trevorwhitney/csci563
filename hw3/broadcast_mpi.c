#include <mpi.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define TAG 1

int main(int argc, char *argv[])
{
	//define variables
	double elapsed_time;
	int p;
	int id;
	int i;
  int j;
	int source;
	int destination;
  int *package;
  int package_size;
  MPI_Status *status;

	//Initialize MPI
  MPI_Init (&argc, &argv);
  MPI_Barrier(MPI_COMM_WORLD);
  elapsed_time = -MPI_Wtime();
  MPI_Comm_rank (MPI_COMM_WORLD, &id);
  MPI_Comm_size (MPI_COMM_WORLD, &p);

  //create 10MB array to send and recieve
  package_size = 2621440;
  package = malloc(sizeof(int[package_size]));
  if (!id) {
    for (i = 0; i < package_size; i++) {
      package[i] = i % 64; //use mod to fill array with different ints that won't overflow
    }
  }
  
  MPI_Bcast(package, package_size, MPI_INT, 0, MPI_COMM_WORLD);

  //verify arrays
  for (i = 0; i < p; i++) {
    if (id == i) {
      for (j = 0; j < package_size; j++) {
        if (package[j] != j % 64) {
          printf("Invalid array element at address %d on id %d\n", j, i);
          exit(1);
        }
      }
    }
  }

  //get execution time
  elapsed_time += MPI_Wtime();
  
  //print results on main processor
  if (!id)
  {
    printf ("Total elapsed time: %10.6f\n", elapsed_time);
  }
  

  //Finalize MPI and exit
  MPI_Finalize();
  return 0;
}