8.times do |num|
  system("mpiexec -np #{8 - num} --hostfile hostfile broadcast >> results/broadcast_#{8 - num}.txt")
  system("mpiexec -np #{8 - num} --hostfile hostfile broadcast_mpi >> results/broadcast_mpi_#{8 -num}.txt")
  system("mpiexec -np #{8 - num} --hostfile hostfile broadcast_pipe >> results/broadcast_pipe_#{8 -num}.txt")
end