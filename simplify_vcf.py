#!/usr/bin/env python2.7
import argparse, sys

parser = argparse.ArgumentParser(description='Rapid scan of protein primary structure motifs', formatter_class=argparse.ArgumentDefaultsHelpFormatter, add_help=False)

#Required arguments
required = parser.add_argument_group('REQUIRED')
required.add_argument('-i', help= 'input VCF', default= sys.stdin, type=str)

#Optional arguments
optional = parser.add_argument_group('OPTIONAL')
optional.add_argument('-h', action="help", help="show this help message and exit")
optional.add_argument('--fields', help= 'VCF annotation fields to retain', nargs= '+', type=str)

args = parser.parse_args()

def ann_dict(s):
	return {annotation.split('=')[0]:annotation.split('=')[1] for annotation in s.split(';')}

if args.i == sys.stdin:
	pass
elif args.i == None:
	parser.error('Input VCF must either be passed via STDIN or -i')
else:
	args.i= open(args.i)
for line in args.i:
	if line[0] == '#':
		print line.strip()
		continue
	line= line.strip().split('\t')
	annotations= ann_dict(line[7])
	select_annotations= ';'.join([field+'='+annotations[field] for field in args.fields])
	print '{0}\t{1}'.format('\t'.join(line[:7]), select_annotations)