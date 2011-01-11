#!/usr/bin/python
import os
import getopt
import sys

paired=False
strata=False
unaligned=False
half_paired=False
def usage(shrimp_folder):
        print >> sys.stderr, 'Usage: %s [--strata] [--sam-unaligned] [--paired] [--half-paired] [--max-alignments=100] [-o/--report 10] input_filename' % sys.argv[0]
        if os.path.exists(shrimp_folder):
                try:
                        handle=open(shrimp_folder+'/utils/MAPPINGS-TOP-HITS')
                        print >> sys.stderr, handle.read()
                        handle.close()
                except IOError, err:
                        print >> sys.stderr, str(err)
        else:   
                print >> sys.stderr, "Please set the SHRIMP_FOLDER environment variable, to your shrimp install folder"

def get_next_entry(h):
	score=0
	sam_line=h.readline()
	while sam_line and sam_line[0]=='@':
		print sam_line,
		sam_line=h.readline()
	if not sam_line:
		return None
	sam_line_list=sam_line.split()
	pos=int(sam_line_list[3])
	new_sam_line_list=sam_line_list[:11]
	for field in sam_line_list[11:]:
		if (field[0:3]=="AS:"):
			score=int(field[5:])
			new_sam_line_list.append(field)
		elif field[0:3] in ("NH:","IH:","H0:","H1:","H2:"):
			pass
		else:
			new_sam_line_list.append(field)
	read_name=sam_line_list[0]
	return [score,read_name,"\t".join(new_sam_line_list),pos]

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
		if not half_paired:
			if entries[0][0]==0:
				entries[1][0]=0
			if entries[1][0]==0:
				entries[0][0]=0
		return entries[0][0]+entries[1][0],entries[0][1],entries[0][2]+'\n'+entries[1][2],entries[0][3]*entries[1][3]	


def unaligned_read(x):
	new_line=[]
	read_list=x.split('\n')
	for read in read_list:
		line_list=read.split('\t')
		read_line=[line_list[0],"4","*","0","0","*","*","0","0","*","*"]
		for field in line_list[11:]:
			if field[0:3] in ("AS:","XX:","CM:","NM:","NH:","IH:","H0:","H1:","H2:"):
				pass
			else:
				read_line.append(field)
		new_line.append("\t".join(read_line))
	return "\n".join(new_line)	

if __name__=='__main__':
	shrimp_folder=os.environ.get("SHRIMP_FOLDER","")
	#Parse options
	try:
		opts, args = getopt.getopt(sys.argv[1:], "o:", ["unpaired","strata", "paired","sam-unaligned","half-paired","report=","max-alignments="])
	except getopt.GetoptError, err:
		print str(err)
		usage(shrimp_folder)
		sys.exit(1)
	#default parameters
	unaligned=False
	strata=False
	number_outputs=10
	max_alignments=0
	half_paired=False
	paired=False
	for o,a in opts:
		if o in ("--strata"):
			strata=True
		elif o in ("--sam-unaligned"):
			unaligned=True
		elif o in ("--half-paired"):
			half_paired=True
		elif o in ("--paired"):
			paired=True
		elif o in ("-o","--report"):
			number_outputs=int(a)
		elif o in ("--max-alignments"):
			max_alignments=int(a)
		elif o in ("--unpaired"):
			pass
		else:
			print "Invalid option : %s" % o
			usage(shrimp_folder)
			exit(1)

	if half_paired and not paired:
		print "Cannot half half_paired without paired"
		usage(shrimp_folder)
		sys.exit(1)

	if len(args)!=1:
		usage(shrimp_folder)
		sys.exit(1)

	file_name=args[0]
	file_handle=sys.stdin
	if (file_name!='-'):
		file_handle=open(file_name,'r')

	current_read=get_next_alignment(file_handle)
	while current_read:
		#find out what the read name is and read in all sequential alignments for this read name
		current_reads=[current_read]
		current_readname=current_read[1]
		current_read=get_next_alignment(file_handle)
		while current_read and current_read[1]==current_readname:
			current_reads.append(current_read)
			current_read=get_next_alignment(file_handle)
		#sort the alginments based on the score, largest first
		current_reads.sort(reverse=True)
		#find out how many best alginments there are
		alignments_to_print=len(current_reads)
		if (strata):
			alignments_to_print=0
			for x in range(len(current_reads)):
				if current_reads[0][0]==current_reads[x][0]:
					alignments_to_print+=1
				else:
					break
		if current_reads[0][0]==0:
			alignments_to_print=0 
			#falls through to unpaired printing
		if alignments_to_print!=0 and (max_alignments==0 or alignments_to_print<=max_alignments):
			alignments_to_print=min(alignments_to_print,number_outputs)
			#want to output the top , alignmnets_to_print	
			for x in range(0,min(len(current_reads),alignments_to_print)):
				print current_reads[x][2]
		elif unaligned:
			print unaligned_read(current_reads[0][2])
