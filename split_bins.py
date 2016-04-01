#!/usr/bin/env python2.7
import argparse, os, sys
import re
import common_objects

def split_bins(fhandle):
    bin_name= ''
    for fasta in common_objects.parse_fasta(fhandle):
        tmp_bin = re.search('bin=(\S+)$', fasta.header).groups()[0]
        if tmp_bin == bin_name

def main(fname):
    with open(fname) as fhandle:
        split_bins(fhandle)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='SCRIPT DESCRIPTION', formatter_class=argparse.ArgumentDefaultsHelpFormatter, add_help=False)

    #Required arguments
    required = parser.add_argument_group('REQUIRED')
    required.add_argument('-f', help= 'fasta file', required=True, type=str)

    #Optional arguments
    optional = parser.add_argument_group('OPTIONAL')
    optional.add_argument('-h', action="help", help="show this help message and exit")

    args = parser.parse_args()
    
    main(args.f)
