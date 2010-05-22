#!/usr/bin/python
import getopt

import sys


def usage():
	print >> sys.stderr, sys.argv[0], "[optional arguments] file.fa offset length"
	print >> sys.stderr, '''extract-subsequence
----------

This script is used to extract part of a fasta sequence from a file.


Parameters:

<file.fa>
        Input file in fasta format.

<offset>
        The offset form where to begin capturing the subsequence.

<length>
        The length of subsequence to caputre and return.

--contig <contig>
        The name of the contig from which to extract sequence from.
        Default is first contig in fasta file.

--fasta 
        If specified returns the subsequence in fasta format.
'''

if __name__=='__main__':
        #Parse options
	try:
		opts, args = getopt.getopt(sys.argv[1:], "c:f", ["contig=","fasta"])
	except getopt.GetoptError, err:
		print str(err)
		usage()
		sys.exit(1)
	#default values
	contig=""
	fasta=False
	for o,a in opts:
		if o in ("-c","--contig"):
			contig=a
		elif o in ("-f","--fasta"):
			fasta=True

	if len(args)!=3:
		usage()
		sys.exit(1)
	fasta_filename=args[0]
	try:
		offset=int(args[1])
		length=int(args[2])
	except ValueError, err:
		print >> sys.stderr, str(err)
		usage()
		sys.exit(1)

	#open the file and start reading the contigs
	try:
	        h=open(fasta_filename,'r')
		current_contig_name=None
		current_contig=[]
		line=h.readline()
		while line:
			if line[0]=='>':
				if contig.strip()==current_contig_name:
					break
				current_contig=[]
				current_contig_name=line[1:].strip()
				if contig=="":
					contig=current_contig_name
			if contig.strip()==current_contig_name:
        			current_contig.append(line.strip())
			line=h.readline()
		if contig.strip()!=current_contig_name:
			print >> sys.stderr, "Contig \"%s\" not found!" % contig
			sys.exit(1)
		h.close()
	except IOError, err:
			print >> sys.stderr, str(err)
			sys.exit(1)

	#get the sequence as a string
	contig_sequence="".join(current_contig[1:])
	#get the coorinates, 1 based
	start_coord=offset
	end_coord=offset+length-1
	s=""
	contig_length=len(contig_sequence)
	#its 1 based cannot be 0!
	sgn=0
	if (start_coord<0 and end_coord>=0):
		sgn=1
	if (abs(start_coord)>contig_length or abs(end_coord+sgn)>contig_length):
		print >> sys.stderr, "%d or %d is out of bounds, contig only %d base pairs long" % (start_coord,end_coord,contig_length)
		sys.exit(1)
	if start_coord>0 and end_coord>0:
		s=contig_sequence[start_coord-1:end_coord]
	elif start_coord<0 and end_coord<0:
		try:
			s=contig_sequence[contig_length+start_coord:contig_length+end_coord+1]
		except IndexError:
			print >> sys.stderr, "%d or %d is out of bounds, contig only %d long" % (start_coord,end_coord,contig_length)
	else:
		end_coord+=1
		s=contig_sequence[start_coord:]
		s+=contig_sequence[:end_coord]

	if not fasta:
		print s	
	else:
		print ">subsequence_%d_to_%d_LENGTH_%d" % (start_coord,end_coord,len(s))
		index=0
		perline=80
		while index<len(s):
			print s[index:index+perline]
			index+=perline
