#!/usr/bin/env python2.7
import sys, argparse

parser = argparse.ArgumentParser(description='Reformats fasta headers for tree-building/alignment programs.', formatter_class=argparse.ArgumentDefaultsHelpFormatter, add_help=False)

#Required arguments
required = parser.add_argument_group('REQUIRED')
required.add_argument('-f', help= 'input.fasta', type=argparse.FileType('r'), default= sys.stdin)

#Optional arguments
optional = parser.add_argument_group('OPTIONAL')
optional.add_argument('-h', action="help", help="show this help message and exit")

args = parser.parse_args()

replacement_tuples= [('\t', '_'), ('[', '_'), (']', '_'), ('(', '_'), (')', '_'), ('.', '_'), (' ', '_'), (';', '_')]
for line in args.f:
	if ">" in line:
		for tuple in replacement_tuples:
			line= line.replace(tuple[0], tuple[1])
	elif "-" in line:
		line.replace('*', '-')
	print line.strip()
