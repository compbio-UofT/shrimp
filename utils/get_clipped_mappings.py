#!/usr/bin/python
import sys

def usage():
	print "%s" % (sys.rgv[0])

if __name__=='__main__':
	line=sys.stdin.readline()
	while line:
		line_list=line.split()
		if line[0]=='@':
			pass
		elif (line_list[5].find("H")>=0):
			print "@"+line_list[0]
			print ":".join(line_list[14].split(":")[2:])+"\n+"
			print ":".join(line_list[13].split(":")[2:])
		line=sys.stdin.readline()
