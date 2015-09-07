#!/usr/bin/env python2.7
import sys, argparse, glob, os, subprocess

parser = argparse.ArgumentParser(description='Extracts single copy genes using HMM profiles from phylosift', formatter_class=argparse.ArgumentDefaultsHelpFormatter, add_help=False)

#Required arguments
required = parser.add_argument_group('REQUIRED')
required.add_argument('-g', help= 'input glob path to individual genomes', required=True, type=str)

#Optional arguments
optional = parser.add_argument_group('OPTIONAL')
optional.add_argument('-h', action="help", help="show this help message and exit")
optional.add_argument('-o', help= "alternative output directory for ribosomal protein files", type= str)

args = parser.parse_args()

hmm_dict= {
'DNGNGWU00014': 'rpL14',
'DNGNGWU00021': 'rpL15',
'DNGNGWU00033': 'rpL18',
'DNGNGWU00007': 'rpL22',
'DNGNGWU00040': 'rpL24',
'DNGNGWU00010': 'rpL2',
'DNGNGWU00012': 'rpL3',
'DNGNGWU00009': 'rpL4',
'DNGNGWU00025': 'rpL5',
'DNGNGWU00023': 'rpL6',
'DNGNGWU00002': 'rpS10',
'DNGNGWU00036': 'rpS17',
'DNGNGWU00016': 'rpS19',
'DNGNGWU00028': 'rpS3',
'DNGNGWU00031': 'rpS8'
}

def run_phylosift(fasta):
	command= 'phylosift search -f --threads 10 --unique --output {0}-phylosift {0} >log'.format(fasta)
	try:
		print "Running phylosift on: {0}".format(fasta)
		subprocess.call(command, shell=True)
	except:
		print "Phylosift failed!!!"
		exit()

def move_files(fasta):
	files= os.listdir('{0}-phylosift/blastDir/'.format(fasta))
	files= [file for file in files if 'DNGNGWU000' in file]
	
	#Read the lookup table into a dictionary
	lookup_handle= open('{0}-phylosift/blastDir/lookup_ID.1.tbl'.format(fasta))
	lookup_dict= {}
	for line in lookup_handle:
		line= line.strip().split('\t')
		lookup_dict[line[1].split('/')[0]]= line[0]
	lookup_handle.close()
	
	
	print "Moving files..."
	for file in files:
		#Retrieve the full length sequence from each ribosomal protein in the HMM dictionary
		lookup= file.split('.')[0]
		if lookup not in hmm_dict:
			continue
		rp= hmm_dict[lookup]
		
		sequence= lookup_dict[open('{0}-phylosift/blastDir/{1}'.format(fasta, file)).readline().strip()[1:].split('.')[0]]
		try:
			subprocess.call('echo "{0}" | pullseq -i {1} -N > {1}-{2}.faa'.format(sequence, fasta, rp), shell=True)
		except:
			exit()
		output_handle= open('{0}-{1}.fixed.faa'.format(fasta, rp), 'w')
		for line in open('{0}-{1}.faa'.format(fasta, rp)):
			if line[0] == '>':
				output_handle.write('>{0}\n'.format(fasta.split('.')[0]))
			else:
				output_handle.write(line)
		output_handle.close()

for fasta in glob.glob(args.g):
	run_phylosift(fasta)
	move_files(fasta)
	
for rp in hmm_dict.values():
	os.system('cat *{0}.fixed.faa > {0}.query.faa'.format(rp))
os.system('rm *-rp*.faa')
os.system('rm -r *-phylosift')




	
