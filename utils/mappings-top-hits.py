#!/usr/bin/python
import sys

paired=False
strata=False
def usage():
	print >> sys.stderr, sys.argv[0], "#alignments_per_read strata[yes/no] paired[yes/no] input_filename"

def get_next_entry(h):
	score=0
	sam_line=h.readline()
	while sam_line and sam_line[0]=='@':
		print sam_line,
		sam_line=h.readline()
	if not sam_line:
		return None
	sam_line_list=sam_line.split()
	for field in sam_line_list[11:]:
		if (field[0:3]=="AS:"):
			score=int(field[5:])
			break
	read_name=sam_line_list[0]
	return score,read_name,sam_line

def get_next_alignment(h):
	if not paired:
		return get_next_entry(h)
	else:
		entries=[get_next_entry(h),get_next_entry(h)]
		if entries[0]==None and entries[1]==None:
			return None
		if entries[0]==None or entries[1]==None:
			print >> sys.stderr, "An error has ocured in reading the next read pair"
			sys.exit(1)
		if (entries[0][1]!=entries[1][1]):
			print >> sys.stderr, "An error has occured in matching read pair names"
			print >> sys.stderr, entries
			sys.exit(1)
		return entries[0][0]+entries[1][0],entries[0][1],entries[0][2]+entries[1][2]		
			

if __name__=='__main__':
	if len(sys.argv)!=5:
		usage()
		sys.exit(1)
	alignments_per_read=int(sys.argv[1])
	if (sys.argv[2].lower()=="yes"):
		strata=True
	if (sys.argv[3].lower()=="yes"):
		paired=True

	file_name=sys.argv[4]
	file_handle=sys.stdin
	if (file_name!='-'):
		file_handle=open(file_name,'r')

	current_read=get_next_alignment(file_handle)
	while current_read:
		current_reads=[current_read]
		current_readname=current_read[1]
		current_read=get_next_alignment(file_handle)
		while current_read and current_read[1]==current_readname:
			current_reads.append(current_read)
			current_read=get_next_alignment(file_handle)
		current_reads.sort(reverse=True)
		best_alignments=0
		if strata:
			for x in range(len(current_reads)):
				if current_reads[0][0]==current_reads[x][0]:
					best_alignments+=1
				else:
					break
		if (len(current_reads)<=alignments_per_read or (strata and best_alignments<=alignments_per_read)):
			for x in range(0,min(len(current_reads),alignments_per_read)):
				print current_reads[x][2],
	
		
	
	
	
