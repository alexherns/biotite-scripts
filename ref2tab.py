#!/usr/bin/env python2.7
from Bio import SeqIO
import sys
'''
Prints length of each sequence in fasta file
Usage: ref2tab.py <file1>
'''

handle = open(sys.argv[1], "r")
for record in SeqIO.parse(handle, "fasta") :
    print "{0}\t{1}".format(record.id, len(record.seq))
handle.close()
