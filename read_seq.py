#!/usr/bin/env python2.7
import argparse, os

parser = argparse.ArgumentParser(description='Faster implementation of pullseq ideas for reads files when large RAM available.', formatter_class=argparse.ArgumentDefaultsHelpFormatter, add_help=False)

#Required arguments
required = parser.add_argument_group('REQUIRED')
required.add_argument('-i', help= 'Read list', required=True, type=str)
required.add_argument('-f', help= 'FASTQ file', required=True, type=str)

#Optional arguments
optional = parser.add_argument_group('OPTIONAL')
optional.add_argument('-h', action="help", help="show this help message and exit")
optional.add_argument('--exclude', action="store_true", help="exclude, rather than select, the reads from file")

args = parser.parse_args()

select_set= set([])
with open(args.i) as input_handle:
	for line in input_handle:
		select_set.add(line.strip())
		
switch= False
with open(args.f) as reads_handle:
	for i, line in enumerate(reads_handle):
		line=line.strip()
		if i%4 == 0:
			read= line[1:].split(' ')[0]
			if read not in select_set:
				switch= False
			else:
				switch= True
			if not args.exclude and switch:
				print line
			elif args.exclude and not switch:
				print line
		elif i%4 == 1:
			if not args.exclude and switch:
				print line
			elif args.exclude and not switch:
				print line
		