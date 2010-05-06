#!/usr/bin/python
import subprocess
import sys
import os
import getopt
import re

def usage(shrimp_folder):
	print >> sys.stderr, 'Usage: %s --shrimp-mode [ls|cs] file1.fa file2.fa ...' % sys.argv[0]
	if os.path.exists(shrimp_folder):
		try:
			handle=open(shrimp_folder+'/utils/PROJECT-DB')
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

def basename(s):
	s=remove_trailing(s)
	return s.split(os.sep)[-1]

def main(argv):
	shrimp_folder=os.environ.get("SHRIMP_FOLDER","")
	#Parse options
	try:
		opts, args = getopt.getopt(argv, "d:m:s:hp",\
                  ["shrimp-folder=","dest-dir=","shrimp-mode=","seed=","h-flag","print-script"])
	except getopt.GetoptError, err:
		print >> sys.stderr, str(err)
		usage(shrimp_folder)
		sys.exit(1)
	#default parameters
	dest_dir="."
	shrimp_mode=""
	seed=""
	script=False
	h_flag=False
	for o,a in opts:
		if o in ("-m","--shrimp-mode"):
			shrimp_mode=a
		elif o in ("-d","--dest-dir"):
			dest_dir=remove_trailing(a)
			if not dest_dir:
				usage(shrimp_folder)
				sys.exit(1)
		elif o in ("-s","--seed"):
			seed=a
		elif o in ("-p","--print-script"):
			script=True
		elif o in ("-h","--h-flag"):
			h_flag=True
		elif o in ("--shrimp-folder"):
			shrimp_folder=a

	#get the shrimp folder
	if not os.path.exists(shrimp_folder):
		usage(shrimp_folder)
		sys.exit(1)
	shrimp_folder=remove_trailing(shrimp_folder)

	#check that shrimp_mode is set
	if shrimp_mode not in ("ls","cs"):
		if len(shrimp_mode)>0:
			print >> sys.stderr, "%s is not a valid shrimp mode" % shrimp_mode
		else:
			print >> sys.stderr, "Please specify a shrimp-mode"
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

	#check genome files
	genome_files=args
	if not genome_files:
		print >> sys.stderr, "No genome files given..."
		usage(shrimp_folder)
		sys.exit(1)

	#check gmapper in path
	gmapper_executable=shrimp_folder+"/bin/gmapper-"+shrimp_mode
	if not (os.path.exists(gmapper_executable) and os.access(gmapper_executable,os.X_OK)):
		print >> sys.stderr, "Cannot find gmapper-"+shrimp_mode+" in %s" % gmapper_executable
		sys.exit(1)
	command_template=[gmapper_executable]
	if seed!="":
		command_template.append('-s')
		command_template.append(seed)
	if h_flag:
		command_template.append('-H')

	#check for files already in the destination
	matching=[]
	try:
		listing=os.listdir(dest_dir+"/")
	except OSError, err:
		print >> sys.stderr, str(err)
		sys.exit(1)
	for fasta_filename in genome_files:
		if fasta_filename[-3:]!=".fa":
			print >> sys.stderr, "Skipping file %s, does not have '.fa' extension" % fasta_filename
			continue
		target_filename=".".join(basename(fasta_filename).split('.')[:-1])+'-'+shrimp_mode
		if listing.__contains__(target_filename+".genome"):
			print >> sys.stderr, "Database already exists in file %s" % (target_filename+".genome")
		elif not os.path.exists(fasta_filename):
			print >> sys.stderr, '%s does not exist!' % fasta_filename
			usage(shrimp_folder)
			sys.exit(1)
		else:
			command=command_template[:]
			command.append('-S')
			command.append('%s/%s' % (dest_dir,target_filename))
			command.append(fasta_filename)
			#run the thing
			if not script:
				print >> sys.stderr, " ".join(command)
				gmapper_process=subprocess.Popen(command)
				if gmapper_process.wait()!=0:
					print >> sys.stderr, "An error has occured"
					sys.exit(1)
			else:
				print " ".join(command)
			

if __name__=='__main__':
	main(sys.argv[1:])


