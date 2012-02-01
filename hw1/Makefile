EXECS = parallel_pi circuits

all: $(EXECS)

parallel_pi: parallel_pi.c
	mpicc -o parallel_pi parallel_pi.c
	
circuits: circuits.c
	mpicc -o circuits circuits.c

clean:
	rm $(EXECS)