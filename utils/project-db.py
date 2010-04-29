#!/usr/bin/python
import subprocess
import sys
import os
import getopt
import re

def usage():
	print '%s rest of usage here' % sys.argv[0]

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
	shrimp_mode="ls"
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


