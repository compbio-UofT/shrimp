#!/usr/bin/python
import subprocess
import sys
import os
import getopt
import re
from get_contigs import get_contigs

def usage():
	print >> sys.stderr, 'Usage: %s --ram-size ramsize file1.fa file2.fa ...' % sys.argv[0]
	print >> sys.stderr, '''
This script is used to split a given reference genome into a set of fasta files
that fit into a target RAM size. gmapper can then be run independently on each
of these files.

Parameters:

<file1.fa> <file2.fa> ...
        List of reference genome files in fasta format.

--ram-size <ram-size>
        Target RAM size, in GB. This parameter is required.

--dest-dir <dest-dir>
        Destination directory where to place the database files. If not given,
        files are placed in the current working directory.

--prefix <prefix>
        Prefix for database files. Default is "db".

--tmp-dir <tmp-dir>
        Directory to store temporary files into. This defaults to
        /tmp/<PID>. Note, the script requires 1x(genome size) temporary space.

--seed <seed0,seed1,...>
        Comma-separated list of seeds that gmapper will use. This list is passed
        on directly to gmapper as argument of parameter -s. See README for more
        details. If absent, gmapper will not be given explicitly any seeds, so
        it will run with its default set of seeds.

--h-flag
        This corresponds to giving gmapper the flag -H, telling it to use
        hashing to index spaced kmers. For seeds of weight greater than 14, this
        is required. See README for more details.

Output:

<dest-dir>/<prefix>-<ram-size>gb-<i>of<n>.fa
        This is the <i>-th of <n> pieces in fasta format.
'''

#need to fix this if going to run on windows
rt=re.compile("/*.*[^/]+")
def remove_trailing(s):
	m=rt.match(s)
	if m:
		return m.group(0)
	return None

#copied from    
#http://stackoverflow.com/questions/377017/test-if-executable-exists-in-python
#but hard to condense?
def which(program):
	def is_exe(fpath):
		return os.path.exists(fpath) and os.access(fpath, os.X_OK)
        
	fpath, fname = os.path.split(program)
	if is_exe("./"+program):
		return "./"+program
	if is_exe("../bin/"+program):
		return "../bin/"+program
	for path in os.environ["PATH"].split(os.pathsep):
		exe_file = os.path.join(path, program)
		if is_exe(exe_file):
			return exe_file
	return None

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
	#Parse options
	try:
		opts, args = getopt.getopt(argv, "r:d:p:t:s:h", ["ram-size=", "dest-dir=","prefix="\
,"tmp-dir=","seed=","h-flag"])
	except getopt.GetoptError, err:
		print str(err)
		usage()
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
			ram_size=int(a)
		elif o in ("-d","--dest-dir"):
			dest_dir=remove_trailing(a)
			if not dest_dir:
				usage()
				sys.exit(1)
		elif o in ("-p","--prefix"):
			prefix=a
		elif o in ("-t","--tmp-dir"):
			tmp_dir=remove_tailing(a)
			if not tmp_dir:
				usage()
				sys.exit(1)
		elif o in ("-s","--seed"):
			seed=a
		elif o in ("-h","--h-flag"):
			h_flag=True

	#check the ram size
	if ram_size<0:
		usage()
		sys.exit(1)

	#check if temp directory exists
	if os.path.exists(tmp_dir):
		for handle in handles:
			handle.close()
		print >> sys.stderr, "Temp directory %s exists!" % tmp_dir
		usage()
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
                	usage()
                	sys.exit(1)

	#check that can run split-contigs
	split_contigs_executable=which('split-contigs')
	if not split_contigs_executable:	
		print >> sys.stderr, "Cannot find split-contigs in PATH"
		sys.exit(1)

	#check for files already in the destination
	matching=[]
	try:
		listing=os.listdir(dest_dir+"/")
	except OSError, err:
		print >> sys.stderr, str(err)
		sys.exit(1)
	
	r=re.compile('%s-%dgb-.+[.]fa$' % (prefix,ram_size))
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
		usage()
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
			usage()
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
	default_number_seeds=4
	default_seed_weights="12,12,12,12"
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
	command=('%s %s %d %s' % (split_contigs_executable,tmp_filename,ram_size,seed_weights)).split()
	split_contigs_process=subprocess.Popen(command,stdout=split_contigs_output_handle)
	if split_contigs_process.wait()!=0:
		print >> sys.stderr, "An error has occured."
		sys.exit(1)
	split_contigs_output_handle.close()

	#extract the contigs	
	output_prefix='%s/%s-%dgb-%sseeds-' % (dest_dir,prefix,ram_size,seed_weights.replace(",","_"))
	ret=get_contigs(tmp_filename,split_contigs_output_filename,output_prefix)
	
	#clean up 
	#remove the tmp_dir
	os.remove(tmp_filename)
	os.remove(split_contigs_output_filename)
	os.rmdir(tmp_dir)
	return ret

if __name__=='__main__':
	main(sys.argv[1:])

