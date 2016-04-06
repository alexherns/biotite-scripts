#!/usr/bin/env python
import common_objects
import sys

def simplify_fasta(fhandle, fout = sys.stdout, save_pickle = False):
    d = common_objects.uniquify_fasta_names(fhandle, output_handle = fout)
    if save_pickle:
        import cPickle
        cPickle.dump(d, save_pickle)

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description = \
            'Uniquify names of sequences in FASTA, saving pickled map upon request')
    parser.add_argument(\
            '-f', type = str, required = False, default = '', \
            help = 'specify fasta file')
    parser.add_argument(\
            '-p', type = str, required = False, default = '', \
            help = 'save pickle file of sequence name maps')
    parser.add_argument(\
            '-o', type = str, default = '', \
            help = 'write to output file (default:stdout)')
    args = parser.parse_args()
    if args.f == '':
        fin = sys.stdin
    else:
        fin = open(args.f, 'r')
    if args.p:
        args.p = open(args.p, 'wb')
    if args.o:
        fout = open(args.o, 'wb')
    else:
        fout = sys.stdout
    simplify_fasta(fin, save_pickle = args.p, fout = fout)
