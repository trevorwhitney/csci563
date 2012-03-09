procs = ARGV[0].to_i

procs.times do |i|
	index = 2**i
	index.times do |j|
		destination = 2**i + j
		puts "Sending from #{j} to #{destination}" if destination < procs
		puts "Recieving on #{destination} from #{j}" if destination < procs
	end
end
