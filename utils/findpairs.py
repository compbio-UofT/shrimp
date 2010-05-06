import sys

if len(sys.argv) < 6:
	print >> sys.stderr, 'usage %s <suffix 1> <suffix 2> <paired output file> <unpaired output file> <readfile1> <readfile2> [...]' % (sys.argv[0])
	sys.exit(1)

suffix1 = sys.argv[1] + '\n'
suffix2 = sys.argv[2] + '\n'

paired_output_file = open(sys.argv[3],'w')
unpaired_output_file = open(sys.argv[4],'w')

def pout(x):
	paired_output_file.write(x)

def upout(x):
	unpaired_output_file.write(x)

reads = []

for f in sys.argv[5:]:
	for line in open(f,'r').readlines():
		if line[0] == '>':
			reads.append([line])
		elif line[0] == '#':
			continue
		elif reads:
			reads[-1].append(line)

def compare(x,y):
	return cmp(x[0],y[0])

reads.sort(cmp=compare)
i = 0

while i < len(reads):
	if reads[i][0].endswith(suffix1):
		if i+1<len(reads) and reads[i+1][0].endswith(suffix2) and reads[i][0][:-1*len(suffix1)] == reads[i + 1][0][:-1*len(suffix2)]:
			pout(''.join(reads[i]))
			pout(''.join(reads[i + 1]))
			i += 2
			continue
	elif reads[i][0].endswith(suffix2):
		if i+1<len(reads) and reads[i+1][0].endswith(suffix1) and reads[i][0][:-1*len(suffix2)] == reads[i + 1][0][:-1*len(suffix1)]:
			pout(''.join(reads[i + 1]))
			pout(''.join(reads[i]))
			i += 2
			continue
	upout(''.join(reads[i]))
	i += 1
