#!/usr/bin/python
import subprocess
import sys
import os
import getopt
import re
split_db=__import__('split-db')
project_db=__import__('project-db')

def usage():
	print 'Usage: %s --shrimp-mode [ls|cs] file1.fa file2.fa ...' % sys.argv[0]
	print '''
This script is used to split a given reference genome into a set of database
files that fit into a target RAM size, then project the files into indexes that
can be reused by gmapper.

Parameters:

<file1.fa> <file2.fa> ..
        List of reference genome files in fasta format.

--ram-size <ram-size>
        Target RAM size, in GB. This parameter is required.

--dest-dir <dest-dir>
        Destination directory where to place the database files. If not given,
        files are placed in the current working directory. Note: both fasta
        files and projections are placed here.

--prefix <prefix>
        Prefix for database files. Default is "db".

--tmp-dir <tmp-dir>
        Directory to store temporary files into. This defaults to
        /tmp/<PID>. Note, the script requires 1x(genome size) temporary space.

--shrimp-mode <mode>
        This is "ls" or "cs", for letter space or color space,
        respectively. This is a required parameter.

--seed <seed0,seed1,..>
        Comma-separated list of seeds that gmapper will use. This list is passed
        on directly to gmapper as argument of parameter -s. See README for more
        details. If absent, gmapper will not be given explicitly any seeds, so
        it will run with its default set of seeds.

--h-flag
        This corresponds to giving gmapper the flag -H, telling it to use
        hashing to index spaced kmers. For seeds of weight greater than 14, this
        is required. See README for more details.


Output:

<dest-dir>/<prefix>-<ram-size>gb-<weights>seeds-<i>of<n>.fa
        This is the <i>-th of <n> pieces in fasta format

<dest-dir>/<prefix>-<ram-size>gb-<weights>seeds-<i>of<n>-<mode>.genome
<dest-dir>/<prefix>-<ram-size>gb-<weights>seeds-<i>of<n>-<mode>.seed.*
        This is the projection of the <i>-th piece.
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
	for path in os.environ["PATH"].split(os.pathsep):
		exe_file = os.path.join(path, program)
		if is_exe(exe_file):
			return exe_file

	return None

def main(argv):
	#Parse options
	try:
		opts, args = getopt.getopt(argv, "r:d:p:t:s:hm:", ["ram-size=", "dest-dir=","prefix="\
,"tmp-dir=","seed=","h-flag","shrimp-mode="])
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
	shrimp_mode=""

	split_db_args=[]
	project_db_args=[]

	for o,a in opts:
		if o in ("-r","--ram-size"):
			ram_size=a
			split_db_args+=[o,a]
		elif o in ("-d","--dest-dir"):
			dest_dir=remove_trailing(a)
			if not dest_dir:
				usage()
				sys.exit(1)
			split_db_args+=[o,a]
			project_db_args+=[o,a]
		elif o in ("-p","--prefix"):
			prefix=a
			split_db_args+=[o,a]
		elif o in ("-t","--tmp-dir"):
			tmp_dir=remove_trailing(a)
			if not tmp_dir:
				usage()
				sys.exit(1)
			split_db_args+=[o,a]
		elif o in ("-s","--seed"):
			seed=a
			split_db_args+=[o,a]
			project_db_args+=[o,a]
		elif o in ("-h","--h-flag"):
			h_flag=True
			split_db_args+=[o,a]
			project_db_args+=[o,a]
		elif o in ("-m","--shrimp-mode"):
			shrimp_mode=a
			project_db_args+=[o,a]
	#check ram size
	if flost(ram_size)<0:
		usage()
		sys.exit(1)


        #check that shrimp_mode is set
        if shrimp_mode not in ("ls","cs"):
                if len(shrimp_mode)>0:
                        print >> sys.stderr, "%s is not a valid shrimp mode" % shrimp_mode
                else:
                        print >> sys.stderr, "Please specify a shrimp_mode"
                usage()
                sys.exit(1)
            
        #check h flag
        valid_seed=re.compile('[01]+$')
        seeds=seed.split(',')
        mx=0
	if len(seeds[0])>0:
       		for s in seeds:
                	m=valid_seed.match(s)
                	if not m:
                        	print >> sys.stderr, "Invalid seed %s" % s
                        	sys.exit(1)
                	mx=max(mx,len(s.replace('0','')))
        	if not h_flag and mx>14:
                	print >> sys.stderr, "For seeds of weight greater then 14, h-flag is required"
                	usage()
                	sys.exit(1)

	if which('gmapper-'+shrimp_mode)==None:
		print >> sys.stderr, "Cannot find gmapper-" + shrimp_mode +" executable"
		sys.exit(1)
	if which('split-contigs')==None:
		print >> sys.stderr, "Cannot find split-contigs executable"
		sys.exit(1)

	genome_files=args
	output_filenames=split_db.main(split_db_args+genome_files)
	project_db.main(project_db_args+map(lambda x : x.split(os.sep)[-1],output_filenames))

if __name__=='__main__':
	main(sys.argv[1:])
	
