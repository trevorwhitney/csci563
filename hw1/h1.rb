2.times do |num|
  system("mpiexec -np #{num + 1} --hostfile hostfile parallel_pi >> pi_#{num + 1}.txt")
  system("mpiexec -np #{num + 1} --hostfile hostfile circuits >> circuits_#{num + 1}.txt")
end