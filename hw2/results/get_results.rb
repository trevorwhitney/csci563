require 'csv'

script_0 = Array.new(8)
script_1 = Array.new(8)
script_2 = Array.new(8)

8.times do |num|
  script_0_results = File.open("0_#{num + 1}.txt")
  script_0_results.each do |line|
    seconds = line[/^Total elapsed time:\s*\d+[\.\d]*$/]
    script_0[num] = seconds[/\d[\.\d]*/] unless seconds.nil?
  end
  
  script_1_results = File.open("1_#{num + 1}.txt")
  script_1_results.each do |line|
    seconds = line[/^Total elapsed time:\s*\d+[\.\d]*$/]
    script_1[num] = seconds[/\d[\.\d]*/] unless seconds.nil?
  end

  script_2_results = File.open("2_#{num + 1}.txt")
  script_2_results.each do |line|
    seconds = line[/^Total elapsed time:\s*\d+[\.\d]*$/]
    script_2[num] = seconds[/\d[\.\d]*/] unless seconds.nil?
  end

end

CSV.open("0_results.csv", "wb") do |csv|
  csv << %w(num_procs time)
  script_0.each_with_index do |result, index|
    csv << [index + 1, result]
  end
end

CSV.open("1_results.csv", "wb") do |csv|
  csv << %w(num_procs time)
  script_1.each_with_index do |result, index|
    csv << [index + 1, result]
  end
end

CSV.open("2_results.csv", "wb") do |csv|
  csv << %w(num_procs time)
  script_2.each_with_index do |result, index|
    csv << [index + 1, result]
  end
end
