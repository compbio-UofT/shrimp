#!/usr/bin/python
import subprocess
import sys
import os
import getopt
import re

def usage():
	print '%s --shrimp-mode [ls|cs] file1.fa file2.fa ...' % sys.argv[0]
	print '''
This script is used to split a given reference genome into a set of database
files that fit into a target RAM size. gmapper can then be run independently on
each of the database files.

Parameters:

<file1.fa> <file2.fa> ...
        Input files in fasta format.

--dest-dir <dest-dir>
        Destination directory where to place the database files. If not given,
        files are placed in the current working directory.

--shrimp-mode <mode>
        This is "ls" or "cs", for letter space or color space,
        respectively. This is a required parameter.

--seed <seed0,seed1,...>
        Comma-separated list of seeds that gmapper will use. This list is passed
        on directly to gmapper as argument of parameter -s. See README for more
        details. If absent, gmapper will not be given explicitly any seeds, so
        it will run with its default set of seeds.

--h-flag
        This corresponds to giving gmapper the flag -H, telling it to use
        hashing to index spaced kmers. For seeds of weight greater than 14, this
        is required. See README for more details.

Notes:

1) The gmapper calls issued by this script can be parallelized across machines.

2) A single machine needs about (1.5)*<ram-size> to create an index to be used
in a machine with <ram-size> RAM. Thus, to create an index for a 16GB machine,
it is highly recommended to use a machine with 32GB of RAM in order to avoid the
use of the swap.


Output:

<file>-<mode>.genome
<file>-<mode>.seed.*
        This is the projection of <file.fa>
'''

#need to fix this if going to run on windows
rt=re.compile("/*.*[^/]+")
def remove_tailing(s):
	m=rt.match(s)
	if m:
		return m.group(0)
	return None

#copied from 
#http://stackoverflow.com/questions/377017/test-if-executable-exists-in-python
#but hard to condense?
def which(program):
    import os
    def is_exe(fpath):
        return os.path.exists(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file

    return None


def main(argv):
	#Parse options
	try:
		opts, args = getopt.getopt(argv, "d:m:s:h", ["dest-dir=",\
"shrimp-mode=","seed=","h-flag"])
	except getopt.GetoptError, err:
		print str(err)
		usage()
		sys.exit(1)
	#default parameters
	dest_dir="."
	shrimp_mode=""
	seed=""
	h_flag=False
	for o,a in opts:
		if o in ("-m","--shrimp-mode"):
			shrimp_mode=a
		elif o in ("-d","--dest-dir"):
			dest_dir=remove_trailing(a)
			if not dest_dir:
				usage()
				sys.exit(1)
		elif o in ("-s","--seed"):
			seed=a
		elif o in ("-h","--h-flag"):
			h_flag=True
	#check that shrimp_mode is set
	if shrimp_mode not in ("ls","cs"):
		if len(shrimp_mode)>0:
			print "%s is not a valid shrimp mode" % shrimp_mode
		else:
			print "Please specify a shrimp_mode"
		usage()
		sys.exit(1)

	#check h flag
	valid_seed=re.compile('[01]+$')
	seeds=seed.split(',')
	mx=0
	for s in seeds:
		m=valid_seed.match(s)
		if not m:
			print "Invalid seed %s" % s
			sys.exit(1)
		mx=max(mx,len(s.replace('0','')))
	if not h_flag and mx>14:
		print "For seeds of weight greater then 14, h-flag is required"
		usage()
		sys.exit(1)

	#get non option parameters, aka genome files
	genome_files=args
	if not genome_files:
		print "No genome files given..."
		usage()
		sys.exit(1)

	#check gmapper in path
	gmapper_executable="gmapper-"+shrimp_mode
	if not which(gmapper_executable):	
		print "Cannot find %s executable" % gmapper_executable
		sys.exit(1)
	command_template=[gmapper_executable]
	if seed!="":
		command_template.append('-s')
		command_template.append(seed)
	if h_flag:
		command_template.append('-H')

	#check for files already in the destination
	matching=[]
	listing=os.listdir(dest_dir+"/")
	for fasta_filename in genome_files:
		target_filename=".".join(fasta_filename.split('.')[:-1])+'-'+shrimp_mode+".genome"
		if listing.__contains__(target_filename):
			print "database already exists in file %s" % target_filename
		elif not os.path.exists(fasta_filename):
			print '%s does not exist!' % fasta_filename
			usage()
			sys.exit(1)
		else:
			command=command_template[:]
			command.append('-S')
			command.append('%s/%s' % (dest_dir,".".join(fasta_filename.split('.')[:-1])+'-'+shrimp_mode))
			command.append(fasta_filename)
			#run the thing
			gmapper_process=subprocess.Popen(command)
			if gmapper_process.wait()!=0:
				print "An error has occured"
				sys.exit(1)
			

if __name__=='__main__':
	main(sys.argv[1:])


