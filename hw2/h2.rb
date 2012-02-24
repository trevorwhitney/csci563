8.times do |num|
  system("mpiexec -np #{8 - num} --hostfile hostfile sieve_0 1000 >> results/0_#{8 - num}.txt")
  system("mpiexec -np #{8 - num} --hostfile hostfile sieve_1 1000 >> results/1_#{8 -num}.txt")
  system("mpiexec -np #{8 - num} --hostfile hostfile sieve_2 1000 >> results/2_#{8 -num}.txt")
end