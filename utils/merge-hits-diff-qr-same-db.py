#!/usr/bin/python
import sys

if len(sys.argv) < 3:
	print >> sys.stderr, 'usage %s <readfile1> <readfile2> [...]' % (sys.argv[0])
	sys.exit(1)

mode = 'none'
header = []
firstfile = True
hits = []
for f in sys.argv[1:]:
	for line in open(f,'r').readlines():
		if line[0] == '#':
			if firstfile:
				mode = 'shrimp'
				header = [line]
				sys.stdout.write(line)
			elif mode == 'sam':
				print >> sys.stderr, 'Error. Are you mixing sam and shrimp output?'
				print >> sys.stderr, 'Offending file: %s' % (f)
				exit(1)
			elif header[0] != line:
				print >> sys.stderr, 'Error. Are you mixing shrimp and probcalc output?'
				print >> sys.stderr, 'Offending file: %s' % (f)
				exit(1)
		elif line[0] == '@':
			if firstfile:
				mode = 'sam'
				header.append(line)
				sys.stdout.write(line)
			elif mode == 'shrimp':
				print >> sys.stderr, 'Error. Are you mixing sam and shrimp output?'
				print >> sys.stderr, 'Offending file: %s' % (f)
				exit(1)
			elif not line in header and not line[0:3] == '@PG':
				print >> sys.stderr, 'Error. Are you mixing output from two different databases?'
				print >> sys.stderr, 'Offending file: %s' % (f)
				exit(1)
		else:
			sys.stdout.write(line)
	firstfile = False
