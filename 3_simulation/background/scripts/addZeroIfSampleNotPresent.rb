
samples = []
samples_vcf = []



def first(line)
	return line[0..8].join("\t")
end

File.open(ARGV[0]) do |f|
	samples = f.read.chomp.split("\t")
end

while line = STDIN.gets
	if line.start_with? "##"
		puts line
		next
	end
	line_split = line.chomp.split("\t")
	if line.start_with? "#"
		samples_vcf = line_split[9..line_split.size]
		print first(line_split)+"\t"+samples.join("\t")
		next
	end
	puts
	print first(line_split)
	genotypes = line_split[9..line_split.size]
	i = 0
	while i < samples.size
		if (!samples_vcf.find_index(samples[i]).nil?)
			pos =  samples_vcf.find_index(samples[i])
			print "\t"+genotypes[pos]
		else
			print "\t0"
		end
		i+=1
	end
end
