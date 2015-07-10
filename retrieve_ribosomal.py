#!/usr/bin/env python2.7
import sys, re

#First file input must be list of organisms you want to retrieve ribosomal proteins from.
orgs= [line.strip() for line in open(sys.argv[1])]
print orgs
#Next file is ribosomal protein file from ggkbase.
seq_dict= {}
for entry in open(sys.argv[2]).read().split(">")[1:]:
	seq_dict[entry.split('\n')[0]]= ''.join(entry.split('\n')[1:])

#Generate fasta output for each ribosomal protein
rp_list= ["L2", "L3", "L4", "L5", "L6", "L14", "L15", "L16", "L18", "L22", "L24", "S3", "S8", "S10", "S17", "S19"]
for rp in rp_list:
	print rp
	file_handle= open('pulled_ribosomal-{0}.fasta'.format(rp), 'w')
	for header in seq_dict:
		Bin= header.split("bin=")[1].split()[0]
		annotation= header.split("bin=")[0]
		if re.search(rp, annotation) and Bin in orgs:
			if re.search("{0}[0-9]".format(rp), annotation):
				continue
			file_handle.write(">{0}\n{1}\n".format(Bin, seq_dict[header]))
	file_handle.close()
