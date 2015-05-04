#!/usr/bin/env python2.7
import sys

if '-h' in sys.argv:
	print """
grep_scaffolds is a simple wrapper to match lines from a scafs2bins file.
Usage: grep_scaffolds.py scaffolds.fasta scafs2bins.tsv

Outputs lines from scafs2bins file that match sequence headers from fasta file
"""
	exit()

scaffolds= set([line.strip().split()[0][1:] for line in open(sys.argv[1]) if line[0]=='>'])
for line in open(sys.argv[2]):
    if line.strip().split("\t")[0] in scaffolds:
        print line.strip()
