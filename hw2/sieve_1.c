#include <mpi.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#define MIN(a,b)  ((a)<(b)?(a):(b))

#define BLOCK_LOW(id,p,n)  ((id)*(n)/(p))

#define BLOCK_HIGH(id,p,n) \
        ( BLOCK_LOW((id)+1,p,n)-2 ) 

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
  int local_index;

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

  //Low and high values for each processor, from 3 to n
  low_value = 3 + BLOCK_LOW(id,p,n-1);
  if (low_value % 2 == 0) low_value -= 1;
  high_value = 3 + BLOCK_HIGH(id,p,n-1);
  if (high_value % 2 == 0) high_value -= 1;
  size = BLOCK_SIZE(id,p,n-1);

  //remove even integers
  array_size = ceil((float)size/2);
  
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
  
  //first prime is 3
  prime = 3;
  do {
    printf("New prime is %d\n", prime);
    if (prime * prime > low_value)
       first = ceil((float)(prime*prime)/2) - 2;
    else {
       if (!(low_value % prime)) first = 0;
       else first = (prime - (low_value % prime))/2;
    }

    //printf("ID: %d, size: %d, array_size: %d, low value: %d, high value: %d\n", id, size, array_size, low_value, high_value);

    //increment by prime, marking the non-primes with 1, or true
    //printf("ID: %d, Starting at %d\n", id, first);
    for (i = first; i < array_size; i += prime) {
      marked[i] = 1;
      //printf("ID: %d, Marking local array at %d as not prime\n", id, i);
    }
    if (!id) {
       while (marked[++index]);
       prime = index*2 + 3;
    }
    MPI_Bcast (&prime,  1, MPI_INT, 0, MPI_COMM_WORLD);
  } while (prime * prime <= n);
  /* End Sieve of Eratosthenes Algorithm */

  /*Begin count of primes */
  count = 0;

  //for all elements in block, if prime is 1/true, increment count
  for (i = 0; i < array_size; i++) {
    if (!marked[i]) count++;
    //printf("ID: %d, Location %d is %d\n", id, i, marked[i]);
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