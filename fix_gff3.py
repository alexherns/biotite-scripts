#!/usr/bin/env python2.7
import argparse, sys

parser = argparse.ArgumentParser(description='Reformats GFF3 files produced from Prodigal to include ID and Name', formatter_class=argparse.ArgumentDefaultsHelpFormatter, add_help=False)

#Required arguments
required = parser.add_argument_group('REQUIRED')
required.add_argument('-i', help= 'input gff file', required=True, type=argparse.FileType('r'), default= sys.stdin)

#Optional arguments
optional = parser.add_argument_group('OPTIONAL')
optional.add_argument('-h', action="help", help="show this help message and exit")


args = parser.parse_args()
id= ''

for line in args.i:
	if line[0] == '#':
		print line.strip()
		continue
	line= line.strip().split('\t')
	new_id= line[0]
	if id != new_id:
		id= new_id
		i= 1
	line[-1]= 'ID={0}_{2};NAME={0}_{2};{1}'.format(id, line[-1], str(i))
	print '\t'.join(line)
	i+=1