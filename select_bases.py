#!/usr/bin/env python2.7
import pysam, argparse, sys

parser = argparse.ArgumentParser(description='Retrieves regions of fasta file as specified', formatter_class=argparse.ArgumentDefaultsHelpFormatter, add_help=False)

#Required arguments
required = parser.add_argument_group('REQUIRED')
required.add_argument('-f', help= 'fasta file for retrieval', type=str, required=True)
required.add_argument('-c', help= 'contig/scaffold for read detection', type=str, required=True)

#Optional arguments
optional = parser.add_argument_group('OPTIONAL')
optional.add_argument('-h', '--help', action="help", help="show this help message and exit")
optional.add_argument('--min', help= 'minimum on scaffold (inclusive)', type=int, default=0)
optional.add_argument('--max', help= 'maximum on scaffold (exclusive)', type=int, default=sys.maxint)
optional.add_argument('--base', help= 'base for coordinate access', type=int, choices= [0,1], default=0)

args = parser.parse_args()

args.min, args.max= args.min-args.base, args.max-args.base

if args.min < 0:
	parser.error('Minimum specified comes to less than 0')

fastafile= open(args.f, 'r').read()[1:]
for id_read in fastafile.split('\n>'):
    if id_read.split('\n')[0] == args.c:
    #   Found the sequence of interest
			args.max= min(args.max, len(''.join(id_read.split('\n')[1:])))
			if args.min > args.max:
				parser.error('Minimum specified is higher than length of scaffold')
			print '>{0}:{1}-{2}'.format(args.c, args.min, args.max)
			print ''.join(id_read.split('\n')[1:])[int(args.min):int(args.max)]
			exit()

print 'No sequences were found with the designated identifier'
