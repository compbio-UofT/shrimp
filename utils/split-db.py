#!/usr/bin/python
import subprocess
import sys
import os
import getopt
import re
from get_contigs import get_contigs

def usage():
	print '%s rest of usage here' % sys.argv[0]

#need to fix this if going to run on windows
rt=re.compile("/*.*[^/]+")
def remove_tailing(s):
	m=rt.match(s)
	if m:
		return m.group(0)
	return None


def append_fasta_file(source_handle,destination_handle):
	#make sure that there are not extra new lines between
	#fasta files, when writting last chunk, strip whitespace
	#from the right
	buffer_size=1024*1024*10 #10 megs
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
	print 'ram-size %dGB, dest-dir %s, prefix %s' % (ram_size,dest_dir,prefix)
	print 'tmp-dir %s, seed %s, h flag %s' % (tmp_dir, seed, str(h_flag))
	if ram_size<0:
		usage()
		sys.exit(1)

	#check for files already in the destination
	#<dest-dir>/<prefix>-<ram-size>gb-*.fa
	matching=[]
	listing=os.listdir(dest_dir+"/")
	r=re.compile('%s-%dgb-.+[.]fa$' % (prefix,ram_size))
	for filename in listing:
		m=r.match(filename)
		if m:
			matching.append(filename)
	if matching:
		matching=", ".join(map(lambda x : dest_dir+"/"+x,matching))
		print 'splitting already done in files %s' % matching
		sys.exit(0)

	#get non option parameters, aka genome files
	genome_files=args
	if not genome_files:
		print "No genome files given..."
		usage()
		sys.exit(1)
	#could use os.access but this could cause a security hole? by python docs
	#just going to open since need to 'cat' together anyway
	handles=[]
	for filename in genome_files:
		try:
			handles.append(open(filename,'rU'))
		except IOError, err:
			print str(err)
			for handle in handles:
				handle.close()
			usage()
			sys.exit(1)
	#all files are open for reading

	#test if temp directory exists
	if os.path.exists(tmp_dir):
		for handle in handles:
			handle.close()
		print "Temp directory %s exists!" % tmp_dir
		usage()
		sys.exit(1)
	#make the directory
	os.makedirs(tmp_dir)
	#make all.fa
	tmp_filename=tmp_dir+'/all.fa'
	tmp_file=open(tmp_filename,'w')
	#cat all the fasta files
	append_fasta_file(handles[0],tmp_file)
	for handle in handles[1:]:
		append_fasta_file(handle,tmp_file)	
	tmp_file.close()

	#do the seeds
	default_number_seeds=4
	default_seed_weights="12,12,12,12"
	hash_table_weight=12
	seed_weights=""
	valid_seed=re.compile('[01]+$')
	if not h_flag:
		if seed=="":
			#this means no H flag and no seeds given
			#SEED_WEIGHTS := DEF_SEED_WEIGHTS
			seed_weights=default_seed_weights
		else:
			seeds=seed.split(',')
			seed_weights=[]
			for seed in seeds:
				m=valid_seed.match(seed)
				if m:
					seed_weights.append(len(seed.replace('0','')))
				else:
					print "Invalid seed, %s" % seed
					#remove the tmp_dir
					os.rempve(tmp_filename)
					os.rmdir(tmp_dir)
					usage()
					sys.exit(1)		
			seed_weights=",".join(map( lambda x : str(x), seed_weights))
	else:
		if seed=="":
			seed_weights=",".join([str(hash_table_weight)]*default_number_seeds)
		else:
			seeds=seed.split(',')
			seed_weights=[]
			for seed in seeds:
				m=valid_seed.match(seed)
				if m:
					seed_weights.append(str(hash_table_weight))
				else:
					print "Invalid seed, %s" % seed
					#remove the tmp_dir
					os.remove(tmp_filename)
					os.rmdir(tmp_dir)
					usage()
					sys.exit(1)		
			seed_weights=",".join(seed_weights)
	print seed_weights

	#want to run split-contigs
	if not os.access('./split-contigs',os.X_OK):	
		print "Cannot find split-contigs in current directory"
		#clean up and exit
		os.remove(tmp_filename)
		os.rmdir(tmp_dir)
		sys.exit(1)
	split_contigs_output_filename=tmp_dir+'/split_contigs_out'
	split_contigs_output_handle=open(split_contigs_output_filename,'w')
	command=('./split-contigs %s %d %s' % (tmp_filename,ram_size,seed_weights)).split()
	split_contigs_process=subprocess.Popen(command,stdout=split_contigs_output_handle)
	split_contigs_process.wait()
	split_contigs_output_handle.close()	

	#extract the contigs	
	output_prefix='%s/%s-%dgb-%sseeds-' % (dest_dir,prefix,ram_size,seed_weights)
	ret=get_contigs(tmp_filename,split_contigs_output_filename,output_prefix)
	
	#clean up 
	#remove the tmp_dir
	os.remove(tmp_filename)
	os.remove(split_contigs_output_filename)
	os.rmdir(tmp_dir)
	return ret

if __name__=='__main__':
	main(sys.argv[1:])

