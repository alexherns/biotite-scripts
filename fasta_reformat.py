#!/usr/bin/env python2.7
import sys
from Bio import SeqIO

if sys.argv[1] == "-h" or len(sys.argv) < 3:
    print """
Reformats lines properly for use with samtools or other picky programs.

Usage:  fasta_reformat.py   input.fasta output.fasta
"""
    exit()

SeqIO.convert(sys.argv[1], "fasta", sys.argv[2], "fasta")
