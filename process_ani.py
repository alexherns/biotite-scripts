#!/usr/bin/env python
import sys, os, getopt, glob

opts, args = getopt.getopt(sys.argv[1:], 'f:i:d:o:g:ht', ['help'])

def usage():
        print """
Runs ANI analysis on table of organisms.

Usage: process_ani.py [OPTIONS]

OPTIONS:
		-f <string>	Directory to FASTA files
		-i <string>	Input table of organisms by group (required if -f is passed)
		-g <glob> glob to list of FASTA files
		-o <string>	Alternative base name for output directory
		-d <string>	Delimiter (Default: \\t)
		-t 		Flag if input table contains header row (Default: off)
"""
        exit()

fasta_directory= ''
input_table= ''
output_directory= 'ani_output'
delim= '\t'
header_row= False
glb= ''

for o, a in opts:
    if o in ('-h', '--help'):
        usage()
    elif o in ('-f'):
        fasta_directory= a
    elif o in ('-i'):
        input_table= a
    elif o in ('-d'):
        delim= a
    elif o in ('-o'):
        output_directory= a
    elif o in ('-t'):
        header_row= True
    elif o in ('-g'):
        glb= a
        

if '' in [fasta_directory, input_table] and glb == '':
        print """
Please specify appropriate parameters. See -h (--help) for help
"""
        exit()

#Simple circuiting for glob path
if glob != '':
	output_directory= output_directory.strip("/")
	temp_dir= "ani_temporary_storage"
	if not os.path.exists(temp_dir):
		os.makedirs(temp_dir)
	else:
		os.system("rm {0}/*".format(temp_dir))
	os.system("cp {0} {1}".format(glb, temp_dir))
	os.system("average_nucleotide_identity.py -i {0} -o {1} -g --gformat png".format(temp_dir, output_directory))
	exit()
	
#Format input and output directories to strip trailing /
fasta_directory= fasta_directory.strip("/")
output_directory= output_directory.strip("/")


fasta_files= os.listdir(fasta_directory)
seqs= [file.split('.fasta')[0] for file in fasta_files]
seqset= set(seqs)

orgset= set()
data= []
#Confirm the presence of all organisms in the folder before moving forward with analysis
table_handle= open(input_table)
if header_row:
	table_handle.readline()
for line in table_handle:
	line= [str.replace(" ", "_").replace("-", "_").replace(".", "_") for str in line.strip().split(delim)]
	data.append(line)
	for org in [str for str in line[1:] if str != '']:
		orgset.add(org)
		if org not in seqset:
			print org, '\tabsent'
		else:
			print org, '\tpresent'

while True:
	answer= input("Continue with ANI analysis (y/n)?")
	if answer == 'y':
		break
	elif answer == 'n':
		exit()
	print 'Incorrect answer'

temp_dir= "ani_temporary_storage"

for line in data:
	ID= int(line[0])
	files= [str + ".fasta" for str in line[1:] if str != '']
	if len(files) <= 1:
		continue

	if not os.path.exists(temp_dir):
		os.makedirs(temp_dir)
	else:
		os.system("rm {0}/*".format(temp_dir))

	for file in files:
		os.system("cp {0}/{1} {2}".format(fasta_directory, file, temp_dir))
		#print "cp {0}/{1} {2}".format(fasta_directory, file, temp_dir)
	output_dir= "{1}_{0}".format(ID, output_directory)
	
	print "Processing bin #{0}".format(ID)
	os.system("average_nucleotide_identity.py -i {0} -o {1} -g --gformat png".format(temp_dir, output_dir))
	#print "average_nucleotide_identity.py -i {0} -o {1} -g --gformat png".format(temp_dir, output_dir)
