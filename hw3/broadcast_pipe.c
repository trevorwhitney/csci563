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
  double j;
  double k;
	int source;
	int destination;
  int *package;
  int package_size;
  int count;
  int send_id;
  MPI_Status *status;

	//Initialize MPI
  MPI_Init (&argc, &argv);
  MPI_Barrier(MPI_COMM_WORLD);
  elapsed_time = -MPI_Wtime();
  MPI_Comm_rank (MPI_COMM_WORLD, &id);
  MPI_Comm_size (MPI_COMM_WORLD, &p);

  //create 10MB array to send and recieve
  package_size = 2621440;
  //package_size = 10;
  package = malloc(sizeof(int[package_size]));
  if (!id) {
    for (i = 0; i < package_size; i++) {
      package[i] = i % 64; //use mod to fill array with different ints that won't overflow
    }
  }
  
  for (i = 0; i < package_size + p - 1; i++) {
    if (i + 1 < p) {
      //only send from procs with data
      //do i + 1 sends, starting at proc 0
      count = i;
      for (j = 0; j < i + 1; j++) {
        if (id == (int)j && (int)j+1 < p) {
          //send package[count] to j+1
          //printf("Iteration: %d. Sending package[%d] from %d to %d in loop 1\n", i, count, (int)j, (int)j+1);
          MPI_Send(&package[count], 1, MPI_INT, (int)j + 1, TAG, MPI_COMM_WORLD);
        }
        if (id == (int)j+1) {
          //recv package[count] from j
          //printf("Iteration: %d. Recieving package[%d] on %d from %d in loop 1\n", i, count, (int)j+1, (int)j);
          MPI_Recv(&package[count], 1, MPI_INT, (int)j, TAG, MPI_COMM_WORLD, status);
        }
        count--;
      }
    }
    else if (i >= package_size) {
      //proc 0 is done, so not all procs will be sending on each iteration
      send_id = i - package_size + 1;
      count = package_size - 1;
      for (j = 0; j < p - 1; j++) {
        /*if (!id) {
          printf("I: %d. Send id is : %d\n", i, send_id);
          printf("I: %d. Count is: %d\n", i, count);
        }*/
        if (id == send_id && send_id + 1 < p) {
          //send package[count] to start_id + 1
          //printf("Iteration: %d. Sending package[%d] from %d to %d in loop 2\n", i, count, send_id, send_id + 1);
          MPI_Send(&package[count], 1, MPI_INT, send_id + 1, TAG, MPI_COMM_WORLD);
        }
        if (id == send_id + 1) {
          //recv package[count] from start_id
          //printf("Iteration: %d. Recieving package[%d] on %d from %d in loop 2\n", i, count, send_id + 1, send_id);
          MPI_Recv(&package[count], 1, MPI_INT, send_id, TAG, MPI_COMM_WORLD, status);
        }
        count--;
        send_id++;
      }
    }
    else {
      //do p - 1 sends
      //0 to 1 is i, each sucessive is i - 1
      count = i;
      for (j = 0; j < p - 1; j++) {
        if (id == (int)j && (int)j + 1 < p) {
          //send package[count] to j+1
          //printf("Iteration: %d. Sending package[%d] from %d to %d in loop 3\n", i, count, (int)j, (int)j+1);
          MPI_Send(&package[count], 1, MPI_INT, (int)j + 1, TAG, MPI_COMM_WORLD);
        }
        if (id == (int)j + 1) {
          //recv pacakge[count] from j
          //printf("Iteration: %d. Recieving package[%d] on %d from %d in loop 3\n", i, count, (int)j+1, (int)j);
          MPI_Recv(package + count, 1, MPI_INT, (int)j, TAG, MPI_COMM_WORLD, status);
        }
        count--;
      }
    }
  }

  //verify arrays
  for (i = 0; i < p; i++) {
    if (id == i) {
      for (k = 0; k < package_size; k++) {
        if (package[(int)k] != (int)k % 64) {
          printf("Invalid array element at address %d on id %d\n", (int)k, i);
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