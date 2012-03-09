EXECS = broadcast broadcast_mpi broadcast_pipe

all: $(EXECS)

broadcast: broadcast.c
	mpicc -o broadcast broadcast.c

broadcast_mpi: broadcast_mpi.c
	mpicc -o broadcast_mpi broadcast_mpi.c

broadcast_pipe: broadcast_pipe.c
	mpicc -o broadcast_pipe broadcast_pipe.c

clean:
	rm $(EXECS)