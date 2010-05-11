#!/usr/bin/python
import sys

if len(sys.argv) < 3:
	print >> sys.stderr, 'usage %s <file1> <file2>' % (sys.argv[0])
	sys.exit(1)
	
file1 = open(sys.argv[1],'r')
file2 = open(sys.argv[2],'r')

def out(x):
	sys.stdout.write(x)

files = [file1,file2]
savelines = ['','']
cfile = 0
onedone = False


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
