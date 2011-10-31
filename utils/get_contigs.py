import sys, re




def get_contigs(reference_filename,splits_filename,prefix=""):
	#Read in the reference
	ref_contigs={}
	h=open(reference_filename,'r')
	current_contig=[]
	line=h.readline()
	while line:
		if line[0]=='>' and len(current_contig)>0:
			ref_contigs[current_contig[0][1:].split()[0]]="".join(current_contig)
			current_contig=[]
		current_contig.append(line)
		line=h.readline()
	if len(current_contig)>0:
		ref_contigs[current_contig[0][1:].split()[0]]="".join(current_contig)
		current_contig=[]
	h.close()

	#Read in the splits
	h=open(splits_filename,'r')
	chunks={}
	current_chunk=[]
	z=re.compile('^chunk \d+:$')  
	f=lambda x: ref_contigs[" ".join(x.split()[:-1])]
	for line in h.readlines():
		m=z.match(line)
		if m and len(current_chunk)>0:
			#parse the whole chunk
			chunks[current_chunk[0]]="".join(map(f,current_chunk[1:]))
			current_chunk=[]
		current_chunk.append(line)
	if len(current_chunk)>0:
		#parse the whole chunk
		chunks[current_chunk[0]]="".join(map(f,current_chunk[1:]))
	h.close()
	#write it out
	output_filenames=[]
	for k in chunks.keys():
		filename=prefix+k[:-2].split()[1]+"of"+str(len(chunks.keys()))+".fa"
		output_filenames.append(filename)
		h=open(filename,'w')
		h.write(chunks[k])
		h.close()
	return output_filenames


if __name__=='__main__':
	if len(sys.argv)!=3:
		print "%s reference_filename splits_filename" % (sys.argv[0])
		sys.exit(1)

	reference_filename=sys.argv[1]
	splits_filename=sys.argv[2]
	get_contigs(reference_filename,splits_filename,"")


