16.times do |num|
  system("mpiexec -np #{16 - num} --hostfile hostfile parallel_pi >> pi_#{num + 1}.txt")
  system("mpiexec -np #{16 - num} --hostfile hostfile circuits >> circuits_#{num + 1}.txt")
end