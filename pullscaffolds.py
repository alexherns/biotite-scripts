#!/usr/bin/env python2.7
import getopt, sys

opts, args = getopt.getopt(sys.argv[1:], 'i:h', ['help', 'input', 'max_cov=', 'min_cov=', 'max_gc=', 'min_gc='])

def usage():
        print """
Can retrieve scaffolds according to gc and coverage specifications.
Currently only works with *fa.genes format as provided in assembly.d/sample/idba_ud_full folder

Usage: pullscaffolds.py -i <INPUT>  [OPTIONS]

OPTIONS:

--max_cov=
--min_cov=
--max_gc=
--min_gc=
"""
        exit()

input_file= ''
max_cov= float('inf')
min_cov= 0.
max_gc= 100.
min_gc= 0.

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

if input_file == '':
        print """
Please specify input file. See -h (--help) for help
"""
        exit()

for line in open(input_file):
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

