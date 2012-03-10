require 'csv'

broadcast = Array.new(8)
broadcast_mpi = Array.new(8)
broadcast_pipe = Array.new(8)

8.times do |num|
  broadcast_results = File.open("broadcast_#{num + 1}.txt")
  broadcast_results.each do |line|
    seconds = line[/^Total elapsed time:\s*\d+[\.\d]*$/]
    broadcast[num] = seconds[/\d[\.\d]*/] unless seconds.nil?
  end
  
  broadcast_mpi_results = File.open("broadcast_mpi_#{num + 1}.txt")
  broadcast_mpi_results.each do |line|
    seconds = line[/^Total elapsed time:\s*\d+[\.\d]*$/]
    broadcast_mpi[num] = seconds[/\d[\.\d]*/] unless seconds.nil?
  end

  broadcast_pipe_results = File.open("broadcast_pipe_#{num + 1}.txt")
  broadcast_pipe_results.each do |line|
    seconds = line[/^Total elapsed time:\s*\d+[\.\d]*$/]
    broadcast_pipe[num] = seconds[/\d[\.\d]*/] unless seconds.nil?
  end

end

CSV.open("broadcast_results.csv", "wb") do |csv|
  csv << %w(num_procs time)
  broadcast.each_with_index do |result, index|
    csv << [index + 1, result]
  end
end

CSV.open("broadcast_mpi_results.csv", "wb") do |csv|
  csv << %w(num_procs time)
  broadcast_mpi.each_with_index do |result, index|
    csv << [index + 1, result]
  end
end

CSV.open("broadcast_pipe_results.csv", "wb") do |csv|
  csv << %w(num_procs time)
  broadcast_pipe.each_with_index do |result, index|
    csv << [index + 1, result]
  end
end
