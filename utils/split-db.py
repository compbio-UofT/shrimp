#!/usr/bin/python
import subprocess
import sys
import os
import getopt
import re
from get_contigs import get_contigs


def usage(shrimp_folder):
	print >> sys.stderr, 'Usage: %s --ram-size ramsize file1.fa file2.fa ...' % sys.argv[0]
	if os.path.exists(shrimp_folder):
		try:
			handle=open(shrimp_folder+'/utils/SPLIT-DB')
			print >> sys.stderr, handle.read()
			handle.close()
		except IOError, err:
			print >> sys.stderr, str(err)
	else:
		print >> sys.stderr, "Please set the SHRIMP_FOLDER environment variable, to your shrimp install folder"


def remove_trailing(s):
	while len(s)>0 and s[-1]==os.sep:
		s=s[:-1]
	return s

def append_fasta_file(source_handle,destination_handle):
	#make sure that there are not extra new lines between
	#fasta files, when writting last chunk, strip whitespace
	#from the right
	buffer_size=1024*1024*60 #60 megs
	data_current=source_handle.read(buffer_size)
	data_next=source_handle.read(buffer_size)
	while data_next:
		destination_handle.write(data_current)
		data_current=data_next
		data_next=source_handle.read(buffer_size)
	destination_handle.write(data_current.rstrip()+'\n')
	source_handle.close()

def main(argv):
	shrimp_folder=os.environ.get("SHRIMP_FOLDER","")
	#Parse options
	try:
		opts, args = getopt.getopt(argv, "r:d:p:t:s:h", ["ram-size=", "dest-dir=","prefix="\
,"tmp-dir=","seed=","h-flag","shrimp-folder="])
	except getopt.GetoptError, err:
		print str(err)
		usage(shrimp_folder)
		sys.exit(1)
	#default parameters
	ram_size=-1
	dest_dir="."
	prefix="db"
	tmp_dir="/tmp/"+str(os.getpid())
	seed=""
	h_flag=False
	for o,a in opts:
		if o in ("-r","--ram-size"):
			ram_size=a
		elif o in ("-d","--dest-dir"):
			dest_dir=remove_trailing(a)
			if not dest_dir:
				usage(shrimp_folder)
				sys.exit(1)
		elif o in ("-p","--prefix"):
			prefix=a
		elif o in ("-t","--tmp-dir"):
			tmp_dir=remove_trailing(a)
			if not tmp_dir:
				usage(shrimp_folder)
				sys.exit(1)
		elif o in ("-s","--seed"):
			seed=a
		elif o in ("-h","--h-flag"):
			h_flag=True
		elif o in ("--shrimp-folder"):
			shrimp_folder=a

	#get the shrimp folder
	if not os.path.exists(shrimp_folder):
		usage(shrimp_folder)
		sys.exit(1)
	shrimp_folder=remove_trailing(shrimp_folder)

	#check the ram size
	if float(ram_size)<0:
		usage(shrimp_folder)
		sys.exit(1)

	#check if temp directory exists
	if os.path.exists(tmp_dir):
		print >> sys.stderr, "Temp directory %s exists!" % tmp_dir
		usage(shrimp_folder)
		sys.exit(1)

        #check h flag
        valid_seed=re.compile('[01]+$')
        seeds=seed.split(',')
        if len(seeds[0])>0:
        	max_seed_length=0
		for s in seeds:
                	if not valid_seed.match(s):
                       		print >> sys.stderr, "Invalid seed %s" % s
                        	sys.exit(1)
                	max_seed_length=max(max_seed_length,len(s.replace('0','')))
        	if not h_flag and max_seed_length>14:
                	print >> sys.stderr, "For seeds of weight greater then 14, h-flag is required"
                	usage(shrimp_folder)
                	sys.exit(1)

	#check that can run split-contigs
	split_contigs_executable=shrimp_folder+'/utils/split-contigs'
	if not (os.path.exists(split_contigs_executable) and os.access(split_contigs_executable,os.X_OK)):
		print >> sys.stderr, "Cannot find split-contigs in %s" % split_contigs_executable
		sys.exit(1)

	#check for files already in the destination
	matching=[]
	try:
		listing=os.listdir(dest_dir+"/")
	except OSError, err:
		print >> sys.stderr, str(err)
		sys.exit(1)
	
	r=re.compile('%s-%sgb-.+[.]fa$' % (prefix,ram_size))
	for filename in listing:
		if r.match(filename):
			matching.append(filename)
	if matching:
		matching=", ".join(map(lambda x : dest_dir+"/"+x,matching))
		print >> sys.stderr, 'splitting already done in files %s' % matching
		sys.exit(0)

	#check genome_files
	genome_files=args
	if not genome_files:
		print >> sys.stderr, "No genome files given..."
		usage(shrimp_folder)
		sys.exit(1)

	#open genome files
	#could use os.access but this could cause a security hole? by python docs
	#just going to open since need to 'cat' together anyway
	handles=[]
	for filename in genome_files:
		try:
			handles.append(open(filename,'rU'))
		except IOError, err:
			print >> sys.stderr, str(err)
			for handle in handles:
				handle.close()
			usage(shrimp_folder)
			sys.exit(1)
	#all files are open for reading

	#make temp directory
	os.makedirs(tmp_dir)

	#make all.fa
	tmp_filename=tmp_dir+'/all.fa'
	tmp_file=open(tmp_filename,'w')
	#cat all the fasta files
	append_fasta_file(handles[0],tmp_file)
	for handle in handles[1:]:
		append_fasta_file(handle,tmp_file)	
	tmp_file.close()

	#do the seeds - seeds are gauranteed to be valid by h_flag check
	default_number_seeds=3
	default_seed_weights="12,12,12"
	hash_table_weight=12
	seed_weights=""
	if not h_flag:
		if seed=="":
			#this means no H flag and no seeds given
			#SEED_WEIGHTS := DEF_SEED_WEIGHTS
			seed_weights=default_seed_weights
		else:
			seeds=seed.split(',')
			seed_weights=[]
			for seed in seeds:
				seed_weights.append(len(seed.replace('0','')))
			seed_weights=",".join(map(str, seed_weights))
	else:
		if seed=="":
			seed_weights=",".join([str(hash_table_weight)]*default_number_seeds)
		else:
			seeds=seed.split(',')
			seed_weights=[]
			for seed in seeds:
				seed_weights.append(str(hash_table_weight))
			seed_weights=",".join(seed_weights)
	print seed_weights

	#run split-contigs
	split_contigs_output_filename=tmp_dir+'/split_contigs_out'
	split_contigs_output_handle=open(split_contigs_output_filename,'w')
	command=('%s %s %s %s' % (split_contigs_executable,tmp_filename,ram_size,seed_weights)).split()
	print " ".join(command)
	split_contigs_process=subprocess.Popen(command,stdout=split_contigs_output_handle)
	if split_contigs_process.wait()!=0:
		print >> sys.stderr, "An error has occured."
		sys.exit(1)
	split_contigs_output_handle.close()

	#extract the contigs	
	output_prefix='%s/%s-%sgb-%sseeds-' % (dest_dir,prefix,ram_size,seed_weights.replace(",","_"))
	ret=get_contigs(tmp_filename,split_contigs_output_filename,output_prefix)
	
	#clean up 
	#remove the tmp_dir
	os.remove(tmp_filename)
	os.remove(split_contigs_output_filename)
	os.rmdir(tmp_dir)
	return ret

if __name__=='__main__':
	main(sys.argv[1:])

