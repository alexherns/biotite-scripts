#!/usr/bin/env python2.7
import sys, os, getopt, glob, argparse

opts, args = getopt.getopt(sys.argv[1:], 'f:i:d:o:g:ht', ['help'])

parser = argparse.ArgumentParser(description='Will create a directory to store contigs binned by ESOM mapping from glob or directory', formatter_class=argparse.ArgumentDefaultsHelpFormatter, add_help=False)

#Required arguments
required = parser.add_argument_group('REQUIRED')

exclusive= parser.add_mutually_exclusive_group(required=True)
exclusive.add_argument('-d', help= 'FASTA directory', type=str)
exclusive.add_argument('-g', help= 'glob path to list of FASTA files', type=str)
required.add_argument('-i', help= 'input table of organisms by group', required= True, type=argparse.FileType('r'), default= '')

#Optional arguments
optional = parser.add_argument_group('OPTIONAL')
optional.add_argument('-h', action="help", help="show this help message and exit")
optional.add_argument('-o', help= 'Alternative base name for output directory', default= 'ani_output', type=str)
optional.add_argument('--delimiter', help= 'Delimiter', default= '\\t', type=str)
optional.add_argument('--header', help= 'Flag if input table contains header row', action='store_true')

args = parser.parse_args()

fasta_directory= args.d
input_table= args.i
output_directory= args.o
delim= args.delimiter
header_row= args.header
glb= args.g        

#Simple circuiting for glob path
if glob != None:
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
table_handle= input_table
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
