#include <mpi.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#define MIN(a,b)  ((a)<(b)?(a):(b))
#define MAX(a,b)  ((a)<(b)?(b):(a))
#define BLOCK_LOW(id,p,n)  ((id)*(n)/(p))

#define BLOCK_HIGH(id,p,n) \
        ( BLOCK_LOW((id)+1,p,n)-1 ) 

#define BLOCK_SIZE(id,p,n) \
        (BLOCK_LOW( (id)+1, p, n) - BLOCK_LOW( (id), p, n) )

#define BLOCK_OWNER(index,p,n) \
        ( ( ((p)*(index)+1)-1 ) / (n) )
#define ARRAY_INDEX(i, prime, iteration) (i - iteration*prime)

int main (int argc, char *argv[])
{
  //define variables
  int n;
  double elapsed_time;
  int p;
  int id;
  int low_value;
  int high_value;
  int size;
  int array_size;
  int proc0_size;
  char* marked;
  int index;
  int prime;
  int first;
  int count;
  int global_count;
  int i;
  int iteration;

  //Initialize MPI
  MPI_Init (&argc, &argv);
  MPI_Barrier(MPI_COMM_WORLD);
  elapsed_time = -MPI_Wtime();
  MPI_Comm_rank (MPI_COMM_WORLD, &id);
  MPI_Comm_size (MPI_COMM_WORLD, &p);
  
  //Check for proper command line parameters, must include N
  if (argc != 2) {
      if (!id) printf ("Command line: %s <m>\n", argv[0]);
      MPI_Finalize();
      exit(1);
  }

  //Convert parameter string to integer
  //N represents the number up to which we need to calculate primes
  n = atoi(argv[1]);

  //Low and high values for each processor
  low_value = 2 + BLOCK_LOW(id,p,n-1);
  high_value = 2 + BLOCK_HIGH(id,p,n-1);
  size = BLOCK_SIZE(id,p,n-1);
  array_size = size/2;
  
  //largest prime is sqrt(n), so first processor has all primes if
  //p is less than sqrt(n). We need to check we don't have more processors
  //than we need.
  proc0_size = (n-1)/p;
  if ((2 + proc0_size) < (int) sqrt((double) n)) {
    if (!id) printf ("Too many processes\n");
    MPI_Finalize();
    exit (1);
  }

  //allocate memory for block, error if unable to
  marked = (char *) malloc (array_size);
  if (marked == NULL) {
    printf ("Cannot allocate enough memory\n");
    MPI_Finalize();
    exit (1);
  }

  /* Begin Sieve of Eratosthenes Algorithm */
  //First fill marked[] with zero/false for all items in block
  for (i = 0; i < array_size; i++) marked[i] = 0;
  
  if (!id) index = 0;
  
  //first prime is 2
  prime = 3;
  do {
    if (prime * prime > low_value)
       first = prime * prime - low_value;
    else {
       if (!(low_value % prime)) first = 3;
       else first = prime - (low_value % prime);
    }
    printf("Id: %d, Low value: %d, First: %d, Size: %d, Array Size: %d\n", low_value, id, first, size, array_size);
    //increment by prime, marking the non-primes with 1, or true
    iteration = 1;
    for (i = first; i < size; i += 2*prime) {
      //translate i into array index
      marked[ARRAY_INDEX(i, prime, iteration)] = 1;
      iteration++;
    }
    if (!id) {
       while (marked[++index]);
       prime = index + 2;
    }
    MPI_Bcast (&prime,  1, MPI_INT, 0, MPI_COMM_WORLD);
  } while (prime * prime <= n);
  /* End Sieve of Eratosthenes Algorithm */

  /*Begin count of primes */
  count = 0;

  //for all elements in block, if prime is 1/true, increment count
  for (i = 0; i < array_size; i++)
    if (!marked[i]) {
      count++;
      printf("Location %d is prime\n", i);
      fflush(stdout);
    }
  
  //Sum count of primes from each process
  MPI_Reduce (&count, &global_count, 1, MPI_INT, MPI_SUM,
    0, MPI_COMM_WORLD);
  elapsed_time += MPI_Wtime();
  
  //print results on main processor
  if (!id) {
    printf ("%d primes are less than or equal to %d\n",
       global_count + 1, n);
    printf ("Total elapsed time: %10.6f\n", elapsed_time);
  }
  MPI_Finalize();
  return 0;
}