#!/usr/bin/env python2.7
import argparse 
from Bio import SeqIO

parser = argparse.ArgumentParser(description='Reformats lines properly for use with samtools or other picky programs.', formatter_class=argparse.ArgumentDefaultsHelpFormatter, add_help=False)

#Required arguments
required = parser.add_argument_group('REQUIRED')
required.add_argument('-f', help= 'input.fasta', required=True, type=str)

#Optional arguments
optional = parser.add_argument_group('OPTIONAL')
optional.add_argument('-h', action="help", help="show this help message and exit")

args = parser.parse_args()

SeqIO.convert(args.f, "fasta", args.f+"-reformatted", "fasta")
