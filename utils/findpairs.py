import sys

if len(sys.argv) < 6:
	print >> sys.stderr, 'usage %s <suffix 1> <suffix 2> <paired output file> <unpaired output file> <readfile1> <readfile2> [...]' % (sys.argv[0])
	exit(1)

suffix1 = sys.argv[1]
suffix2 = sys.argv[2]

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
		if i+1<len(reads) and reads[i+1][0].endswith(suffix2):
			pout('\n'.join(reads[i]))
			pout('\n'.join(reads[i + 1]))
			i += 2
			continue
	elif reads[i][0].endswith(suffix2):
		if i+1<len(reads) and reads[i+1][0].endswith(suffix1):
			pout('\n'.join(reads[i]))
			pout('\n'.join(reads[i + 1]))
			i += 2
			continue
	upout('\n'.join(reads[i]))
	i += 1
