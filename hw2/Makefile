EXECS = sieve_0 sieve_1 sieve_2

all: $(EXECS)

sieve_0: sieve_0.c
	mpicc -o sieve_0 sieve_0.c

sieve_1: sieve_1.c
	mpicc -o sieve_1 sieve_1.c

sieve_2: sieve_2.c
	mpicc -o sieve_2 sieve_2.c

clean:
	rm $(EXECS)