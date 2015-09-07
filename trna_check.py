#!/usr/bin/env python2.7
import sys, os, getopt, glob, argparse

parser = argparse.ArgumentParser(description='Runs a check of the tRNAs predicted by tRNAscan-SE to determine which are missing', formatter_class=argparse.ArgumentDefaultsHelpFormatter, add_help=False)

#Required arguments
required = parser.add_argument_group('REQUIRED')

required.add_argument('-t', help= 'tRNAscane-SE output', type=argparse.FileType('r'), default= sys.stdin)

#Optional arguments
optional = parser.add_argument_group('OPTIONAL')
optional.add_argument('-h', action="help", help="show this help message and exit")

args = parser.parse_args()

amino_acids= ['Gly', 'Ala', 'Val', 'Leu', 'Ile', 'Pro', 
'Phe', 'Tyr', 'Trp',
'Ser', 'Thr', 'Cys', 'Met', 'Asn', 'Gln',
'Lys', 'Arg', 'His',
'Asp', 'Glu']

for line in args.t:
	line= line.strip().split('\t')
	if line[4] in amino_acids:
		amino_acids.pop(amino_acids.index(line[4]))
		
print ('None missing!' if amino_acids == [] else 'Missing:\n'+'\n'.join(amino_acids))