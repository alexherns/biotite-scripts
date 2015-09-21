#!/usr/bin/env python2.7
import sys, argparse

parser = argparse.ArgumentParser(description='Reformats fasta headers for tree-building/alignment programs.', formatter_class=argparse.ArgumentDefaultsHelpFormatter, add_help=False)

#Required arguments
required = parser.add_argument_group('REQUIRED')
required.add_argument('-f', help= 'input.fasta', type=argparse.FileType('r'), default= sys.stdin)

#Optional arguments
optional = parser.add_argument_group('OPTIONAL')
optional.add_argument('-h', action="help", help="show this help message and exit")
optional.add_argument('--ggkbase', action= 'store_true', help= 'automatic fixes if file downloaded from ggkbase')
optional.add_argument('--prodigal', action= 'store_true', help= 'automatic fixes if file generateed from prodigal')
optional.add_argument('--uniprot', action= 'store_true', help= 'automatic fixes if file downloaded from uniprot')
optional.add_argument('--silva', action= 'store_true', help= 'automatic fixes if file downloaded from silva')
optional.add_argument('--seed', action= 'store_true', help= 'labels sequences as SEED')

args = parser.parse_args()

replacement_tuples= [('\t', '_'), ('[', '_'), (']', '_'), ('(', '_'), (')', '_'), ('.', '_'), (' ', '_'), (';', '_')]

if args.ggkbase:
	for line in args.f:
		if ">" in line:
			if "bin=" in line:
				pre_bin, post_bin= line.strip().split("bin=")[0:2]
				bin= post_bin.split(" ")[0]
				gene= pre_bin.split(" ")[0]
				print gene + "-" + bin
			else:
				print line.strip()
		else:
			if "-" in line:
				line.replace('*', '-')
			print line.strip()
	
elif args.prodigal:
	for line in args.f:
		if ">" in line:
			print line.strip().split(' ')[0]
		else:
			print line.strip()
			
elif args.uniprot:
	for line in args.f:
		if ">" in line:
			start, species= line.strip().split("OS=")
			accession= start.split(' ')[0]
			description= '_'.join(start.split(' ')[1:-1])
			species= "_".join(species.split("=")[0].split(" ")[:-1])
			line= "---".join([accession, description, species])
			for tuple in replacement_tuples:
				line= line.replace(tuple[0], tuple[1])
			if args.seed:
				line+= '---SEED'
		print line.strip()

elif args.silva:
	for line in args.f:
		if ">" in line:
			replacement_tuples= [('\t', '_'), ('[', '_'), (']', '_'), ('(', '_'), (')', '_'), ('.', '_'), (' ', '_'), (';', '|')]
			for tuple in replacement_tuples:
				line= line.replace(tuple[0], tuple[1])
			if args.seed:
				line+= '---SEED'
		print line.strip()
		
else:
	for line in args.f:
		if ">" in line:
			for tuple in replacement_tuples:
				line= line.replace(tuple[0], tuple[1])
		elif "-" in line:
			line.replace('*', '-')
		if args.seed:
				line+= '---SEED'
		print line.strip()

