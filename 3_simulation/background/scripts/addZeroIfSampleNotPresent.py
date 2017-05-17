import sys

samples = []
samples_vcf = []

def first(line):
	return "\t".join(line[0:9])


f = open(sys.argv[1])
samples = f.read().strip().split("\t")

for line in sys.stdin:

	if line.startswith("##"):
		sys.stdout.write(line)
		continue

	line_split = line.strip().split("\t")

	if line.startswith("#"):
		samples_vcf = line_split[9:len(line_split)]
		sys.stdout.write(first(line_split)+"\t" + "\t".join(samples))
		continue
	
	sys.stdout.write("\n")
	sys.stdout.write(first(line_split))
	genotypes = line_split[9:len(line_split)]
	i = 0

	while i < len(samples):
		if samples[i] in samples_vcf:
			pos =  samples_vcf.index(samples[i])
			sys.stdout.write("\t"+genotypes[pos])
		else:
			sys.stdout.write("\t0")

		i+=1
