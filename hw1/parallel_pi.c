#include <stdio.h>
#include <mpi.h>
#include <math.h>

int main (int argc, char *argv[])
{
  long double i;
  int id; /*Processor rank */
  int p; /*Number of Processors*/
  long double n;
  long double sum;
  long double global_sum;
  long double pi;
  long elapsed_time;
  long double calculate_pi (long double, long double);

  MPI_Init (&argc, &argv);
  MPI_Barrier (MPI_COMM_WORLD);
  elapsed_time = - MPI_Wtime();

  MPI_Comm_rank (MPI_COMM_WORLD, &id);
  MPI_Comm_size (MPI_COMM_WORLD, &p);

  n = pow(10,10);
  if (!id) printf ("n is %li\n", (long)n);

  sum = 0;
  for (i = id; i < n; i += p)
    sum += calculate_pi(i, n);

  printf("Process %d is done\n", id);
  fflush(stdout);

  MPI_Reduce (&sum, &global_sum, 1, MPI_LONG_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  if (!id) {
    pi = 3.1415926535897932384626;
    printf("The value of Pi is: %1.22Le\n", pi);
    printf("My result is %1.22Le\n", global_sum);
    printf("Difference between Pi and my calculation is: %1.22Le\n", pi - global_sum);
  }
  
  MPI_Barrier (MPI_COMM_WORLD);
  elapsed_time += MPI_Wtime();
  MPI_Finalize();

  if (!id) printf("Executed in %li seconds\n", elapsed_time);
  return 0;
}

long double calculate_pi(long double i, long double n)
{
  long double sum;
  //i += 0.5;
  //base = pow(i, 2);
  //base += 1;
  //sum = (1/n) * (4/base);
  sum = (1/n) * (4 /(1 + pow((i + 0.5), 2)) );

  return sum;
}
