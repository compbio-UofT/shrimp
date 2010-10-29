#!/usr/bin/python
import sys

if len(sys.argv) < 3:
	print >> sys.stderr, 'usage %s [-Q] <file1> <file2>\nnote: -Q (fastq mode) requires 4 lines per read' % (sys.argv[0])
	sys.exit(1)

if sys.argv[1] == '-Q':
	fastq = True
	mod = 1
else:
	fastq = False
	mod = 0



file1 = open(sys.argv[1 + mod],'r')
file2 = open(sys.argv[2 + mod],'r')

def out(x):
	sys.stdout.write(x)

files = [file1,file2]
savelines = ['','']
cfile = 0
onedone = False
done = False

if fastq:
	while not done:
		for i in range(4):
			line = files[cfile].readline()
			if line:
				out(line)
			elif onedone:
				done = True
			else:
				cfile = (cfile + 1) % 2
				onedone = True
				break
		cfile = (cfile + 1) % 2
				
else:
	while 1:
		line = files[cfile].readline()
		if not line:
			if not onedone:
				cfile = (cfile + 1) % 2
				out(savelines[cfile])
				onedone = True
			else:
				break
		elif line[0] == '>':
			savelines[cfile] = line
			cfile = (cfile + 1) % 2
			out(savelines[cfile])
		else:
			out(line)
