#!/usr/bin/python
import subprocess
import sys
import os
import getopt
import re

def usage(shrimp_folder):
	print >> sys.stderr, 'Usage: %s --shrimp-mode [ls|cs] --ram-size <size> file1.fa file2.fa ...' % sys.argv[0]
        if os.path.exists(shrimp_folder):
		try:
			handle=open(shrimp_folder+'/utils/SPLIT-PROJECT-DB')
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

def main(argv):
	shrimp_folder=os.environ.get("SHRIMP_FOLDER","")
	#Parse options
	try:
		opts, args = getopt.getopt(argv, "r:d:p:t:s:hm:", ["ram-size=", "dest-dir=","prefix="\
,"tmp-dir=","seed=","h-flag","shrimp-mode=","shrimp-folder=","print-script"])
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
				usage(shrimp_folder)
				sys.exit(1)
			split_db_args+=[o,a]
			project_db_args+=[o,a]
		elif o in ("-p","--prefix"):
			prefix=a
			split_db_args+=[o,a]
		elif o in ("-t","--tmp-dir"):
			tmp_dir=remove_trailing(a)
			if not tmp_dir:
				usage(shrimp_folder)
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
		elif o in ("--shrimp-folder"):
			shrimp_folder=a
			project_db_args+=[o,a]
			split_db_args+=[o,a]
		elif o in ("--print-script"):
			project_db_args+=[o,a]

        #get the shrimp folder
        if not os.path.exists(shrimp_folder):
                usage(shrimp_folder)
                sys.exit(1)
        shrimp_folder=remove_trailing(shrimp_folder)


	#check ram size
	if float(ram_size)<0:
		usage(shrimp_folder)
		sys.exit(1)


        #check that shrimp_mode is set
        if shrimp_mode not in ("ls","cs"):
                if len(shrimp_mode)>0:
                        print >> sys.stderr, "%s is not a valid shrimp mode" % shrimp_mode
                else:
                        print >> sys.stderr, "Please specify a shrimp_mode"
                usage(shrimp_folder)
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
                	usage(shrimp_folder)
                	sys.exit(1)

	#import split-db and project-db
	if not os.path.exists(shrimp_folder+'/utils/split-db.py'):
		print >> sys.stderr, "Cannot find %s" % (shrimp_folder+'/utils/split-db.py')
		usage(shrimp_folder)
		sys.exit(1)
	if not os.path.exists(shrimp_folder+'/utils/project-db.py'):
		print >> sys.stderr, "Cannot find %s" % (shrimp_folder+'/utils/project-db.py')
		usage(shrimp_folder)
		sys.exit(1)

	sys.path.insert(0, shrimp_folder+'/utils/')
	split_db=__import__('split-db')
	project_db=__import__('project-db')

	gmapper_executable=shrimp_folder+'/bin/gmapper-'+shrimp_mode
	if not (os.path.exists(gmapper_executable) and os.access(gmapper_executable,os.X_OK)):
		print >> sys.stderr, "Cannot find gmapper-"+shrimp_mode+" in %s" % gmapper_executable
		sys.exit(1)
	split_contigs_executable=shrimp_folder+'/utils/split-contigs'
	if not (os.path.exists(split_contigs_executable) and os.access(split_contigs_executable,os.X_OK)):
		print >> sys.stderr, "Cannot find split-contigs in %s" % split_contigs_executable
		sys.exit(1)
	genome_files=args
	output_filenames=split_db.main(split_db_args+genome_files)
	project_db.main(project_db_args+map(lambda x : x.split(os.sep)[-1],output_filenames))

if __name__=='__main__':
	main(sys.argv[1:])
	
