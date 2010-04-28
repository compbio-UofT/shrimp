#!/usr/bin/python
import subprocess
import sys
import os
import getopt
split_db=__import__('split-db')
project_db=__import__('project-db')

def usage():
	print '%s rest of usage here' % sys.argv[0]

def main(argv):
	#Parse options
	try:
		opts, args = getopt.getopt(argv, "r:d:p:t:s:hm:", ["ram-size=", "dest-dir=","prefix="\
,"tmp-dir=","seed=","h-flag","shrimp-mode"])
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
	shrimp_mode="ls"

	split_db_args=[]
	project_db_args=[]

	for o,a in opts:
		if o in ("-r","--ram-size"):
			ram_size=int(a)
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
			tmp_dir=remove_tailing(a)
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
	if ram_size<0:
		usage()
		sys.exit(1)

	genome_files=args
	for fasta_filename in genome_files:
		print split_db_args+[fasta_filename]
		output_filenames=split_db.main(split_db_args+[fasta_filename])
		print output_filenames

if __name__=='__main__':
	main(sys.argv[1:])
	
