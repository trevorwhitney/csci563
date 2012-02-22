#include <mpi.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

//Macros
#define MIN(a,b)  ((a)<(b)?(a):(b))

#define BLOCK_LOW(id,p,n)  (( (id)*(n)/(p) % 2 == 0)?( ((id)*(n)/(p)) ):( (id)*(n)/(p) - 1 ))

#define BLOCK_HIGH(id,p,n) \
        ( ((BLOCK_LOW((id)+1,p,n)-2) % 2 == 0)?((BLOCK_LOW((id)+1,p,n)-2)):((BLOCK_LOW((id)+1,p,n)-2) - 1) ) 

#define BLOCK_SIZE(id,p,n) \
        (BLOCK_LOW( (id)+1, p, n) - BLOCK_LOW( (id), p, n) )

#define BLOCK_OWNER(index,p,n) \
        ( ( ((p)*(index)+1)-1 ) / (n) )
#define ARRAY_INDEX(i, prime, iteration) (i - iteration*prime)
#define GLOBAL_TO_LOCAL(value, low_value) \
        ( (ceil((float)value/2) - 2) - (ceil((float)low_value/2) - 2) )
#define INDEX_TO_VALUE(index) \
        ( index*2 + 3 )

//function prototypes
void decompose_data(int, int, int, int*, int*, int*);
void processor_count_check(int, int, int);
char* allocate_memory(int);
void find_first_index(int*, int, int);
int count_local_primes(char*, int);
char* discover_primes(int);


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
  char* marked;
  int index;
  int prime;
  int first;
  int count;
  int global_count;
  int i;
  char* primes;

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

  //Add a function call to discover primes
  primes = discover_primes(n);

  decompose_data(id, p, n, &low_value, &high_value, &size);
  processor_count_check(p, n, id);

  /* Begin Sieve of Eratosthenes Algorithm */
  marked = allocate_memory(size);
  if (!id) index = 0;
  
  //first prime is 3, since we're skipping 2
  prime = 3;
  do {
    find_first_index(&first, prime, low_value);

    //increment by prime, marking the non-primes with 1, or 'marked'
    for (i = first; i < size; i += prime) {
      marked[i] = 1;
    }

    //calculate new prime on id 0
    if (!id) {
       while (marked[++index]);
       prime = INDEX_TO_VALUE(index);
    }
    //broadcast new prime to other processors
    MPI_Bcast (&prime,  1, MPI_INT, 0, MPI_COMM_WORLD);
  } while (prime * prime <= n);
  /* End Sieve of Eratosthenes Algorithm */

  
  count = count_local_primes(marked, size);
  
  //Sum count of primes from each process
  MPI_Reduce (&count, &global_count, 1, MPI_INT, MPI_SUM,
    0, MPI_COMM_WORLD);
  elapsed_time += MPI_Wtime();

  //add one for 2, which is the only even prime we skipped
  global_count += 1;
  
  //print results on main processor
  if (!id) {
    printf ("%d primes are less than or equal to %d\n",
       global_count, n);
    printf ("Total elapsed time: %10.6f\n", elapsed_time);
  }

  //free memory and finish
  free(marked);
  MPI_Finalize();
  return 0;
}


void decompose_data(int id, int p, int n, int *low_value, int *high_value, int *size) {
  //Low and high values for each processor, from 3 to n
  int true_size;

  *low_value = 3 + BLOCK_LOW(id,p,n-1);
  *high_value = 3 + BLOCK_HIGH(id,p,n-1);
  true_size = BLOCK_SIZE(id,p,n-1);

  //remove even integers
  *size = true_size/2;
}


void processor_count_check(int p, int n, int id) {
  //largest prime is sqrt(n), so first processor has all primes if
  //p is less than sqrt(n). We need to check we don't have more processors
  //than we need.
  int proc0_size;

  proc0_size = (n-1)/p;
  if ((2 + proc0_size) < (int) sqrt((double) n)) {
    if (!id) printf ("Too many processes\n");
    MPI_Finalize();
    exit (1);
  }
}


char* allocate_memory(int size) {
  //allocate memory for block, error if unable to
  char* marked;
  int i;
  
  marked = (char *) malloc (size);
  if (marked == NULL) {
    printf ("Cannot allocate enough memory\n");
    MPI_Finalize();
    exit (1);
  }

  //Fill marked[] with zero/false for all items in block
  for (i = 0; i < size; i++) marked[i] = 0;
  return marked;
}


void find_first_index(int *first, int prime, int low_value) {
  int mod_prime;
  int value;

  if (prime * prime > low_value) {
      //If low_value is less than prime*prime, then we need to start at prime*prime
      *first = GLOBAL_TO_LOCAL(prime*prime, low_value);
  }
  else {
    if (!(low_value % prime)) *first = 0;
    else {
      mod_prime = low_value % prime;
      if (mod_prime % 2 == 0) value = low_value - mod_prime + 2*prime;
      else value = low_value - mod_prime + prime;
      *first = GLOBAL_TO_LOCAL(value, low_value);
    }
  }
}


int count_local_primes(char* marked, int size) {
  int count = 0;
  int i;

  //for all elements in block, if prime is 1/true, increment count
  for (i = 0; i < size; i++) {
    if (!marked[i]) count++;
  }

  return count;
}

char* discover_primes(int n) {
  char* primes;
  char* marked;
  int prime;
  int index;
  int i;
  int first;
  int size;
  int array_size;
  int num_primes;
  int counter;

  //use sequential algorithm to find primes from 3 to sqrt(n)
  size = sqrt((double) n);
  array_size = ceil((float)size/2) - 1;
  printf("Size is %d, with array size of %d\n", size, array_size);

  //1. create a list of natural, odd numbers, up to sqrt(n), all unmarked
  marked = allocate_memory(array_size);

  //2. Set K to 3, the first prime in the list
  prime = 3;

  //3. Repeat
  //  a) Mark all multiples of k between k^2 and sqrt(n)
  //  b) Find the smallest number greater than k that is unmarked, set k to this new value
  //Until k^2 > sqrt(n)
  index = 0;
  num_primes = 0;
  
  do {
    //keep track of the number of primes
    num_primes++;

    //low value is 3
    find_first_index(&first, prime, 3);

    //increment by prime, marking the non-primes with 1, or 'marked'
    for (i = first; i < array_size; i += prime) {
      printf("First: %d\n", first);
      marked[i] = 1;
      printf("Marking %d at location %d as not prime\n", INDEX_TO_VALUE(i), i);
    }

    while (marked[++index]);
    prime = INDEX_TO_VALUE(index);

  } while (prime * prime <= n);

  printf("Found %d primes under %d\n", num_primes, size);

  //4. The unmarked numbers are primes, store in an array

  primes = (char *) malloc(num_primes);
  counter = 0;

  for (index = 0; index < array_size; index++) {
    if(!marked[index]) {
      prime = INDEX_TO_VALUE(index);
      printf("Prime #%d: %d\n", index, prime);
      primes[counter++] = prime;
    }
  }

  for (i = 0; i < num_primes; i++) {
    printf("primes[%d] = %d\n", i, primes[i]);
  }

  return primes;
  
}