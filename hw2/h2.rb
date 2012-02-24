16.times do |num|
  system("mpiexec -np #{16 - num} --hostfile hostfile ./parallel_pi >> results/pi_#{16 - num}.txt")
  system("mpiexec -np #{16 - num} --hostfile hostfile ./circuits >> results/circuits_#{16 -num}.txt")
end
