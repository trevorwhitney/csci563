require 'csv'

pi = Array.new(16)
circuits = Array.new(16)

16.times do |num|
  pi_results = File.open("pi_#{num + 1}.txt")
  pi_results.each do |line|
    seconds = line[/^Executed\ in\ \d+ seconds$/]
    pi[num] = seconds[/\d+/] unless seconds.nil?
  end

  c_results = File.open("circuits_#{num + 1}.txt")
  c_results.each do |line|
    seconds = line[/^Executed\ in\ \d[\.\d]*\ seconds$/]
    circuits[num] = seconds[/\d[\.\d]*/] unless seconds.nil?
  end
end

CSV.open("circuit_results.csv", "wb") do |csv|
  csv << %w(num_procs time)
  circuits.each_with_index do |result, index|
    csv << [index + 1, result]
  end
end

CSV.open("pi_results.csv", "wb") do |csv|
  csv << %w(num_procs time)
  pi.each_with_index do |result, index|
    csv << [index + 1, result]
  end
end