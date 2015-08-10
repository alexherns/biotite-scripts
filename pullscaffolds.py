#!/usr/bin/env python2.7
import argparse, sys

parser = argparse.ArgumentParser(description='''Can retrieve scaffolds according
to gc and coverage specifications. Currently only works with *fa.genes format as
provided in assembly.d/sample/idba_ud_full folder.''', 
	formatter_class=argparse.ArgumentDefaultsHelpFormatter, add_help=False)

#Required arguments
required = parser.add_argument_group('REQUIRED')
required.add_argument('-i', help= 'input fasta', required=True, 
	type=argparse.FileType('r'))

#Optional arguments
optional = parser.add_argument_group('OPTIONAL')
optional.add_argument('-h', action="help", 
	help="show this help message and exit")
optional.add_argument('--max_cov', help= 'Maximum coverage of scaffold', 
	metavar= 'X', type=float, default= float('inf'))
optional.add_argument('--min_cov', help= 'Minimum coverage of scaffold', 
	metavar= 'X', type=float, default= 0.)
optional.add_argument('--max_gc', help= 'Maximum GC of scaffold', metavar= 'X',
	type=float, default= 100.)
optional.add_argument('--min_gc', help= 'Minimum GC of scaffold', metavar= 'X',
	type=float, default= 0.)

args = parser.parse_args()

input_file= args.i
max_cov= args.max_cov
min_cov= args.min_cov
max_gc= args.max_gc
min_gc= args.min_gc

if args.min_cov < 0.:
	parser.error("Coverage must be greater than or equal to 0")
if args.min_gc < 0.:
	parser.error("GC must be greater than or equal to 0")
if args.max_gc > 100.:
	parser.error("GC must be less than or equal to 100")


for o, a in opts:
    if o in ('-h', '--help'):
        usage()
    elif o in ('-i', '--input'):
        input_file= a
    elif o in ('--max_cov'):
        max_cov= float(a)
    elif o in ('--min_cov'):
        min_cov= float(a)
    elif o in ('--max_gc'):
        max_gc= float(a)
    elif o in ('--min_gc'):
        min_gc= float(a)

for line in input_file:
    if line.strip().split()[0] == "DEFINITION":
        seq_header= line.split('seqhdr="')[1].split('"')[0]
        sepped= seq_header.split()
        read_length= float(sepped[1].split("read_length_")[1])
        contig_length= float(sepped[2].split("length_")[1])
        read_count= float(sepped[3].split("read_count_")[1])
        cov= read_length*read_count/contig_length
        if cov < min_cov or cov > max_cov:
            continue
        gc= float(line.split('gc_cont=')[1].split(';')[0])
        if gc < min_gc or gc > max_gc:
            continue
        print "{0} cov={1} gc={2}".format(seq_header, cov, gc)

