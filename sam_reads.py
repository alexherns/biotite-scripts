#!/usr/bin/env python2.7
"""
Simple script to pull all reads matching specified contigs from SAM file

Usage: sam_reads.py <input.sam> <contigs.txt>
Output: input.sliced.sam

Notes: contigs are separated by line and are formatted identically to contig definitions in SAM header
"""

import pysam, sys

samfile= pysam.AlignmentFile(sys.argv[1], "r")
scaffoldfile= sys.argv[2]

samoutput= pysam.AlignmentFile(sys.argv[1].split(".sam")[0]+".sliced.sam", "wh", header= samfile.header)

scaffolds= set([line.strip() for line in open(scaffoldfile)])

id_dict= {}
for SQ in samfile.header["SQ"]:
    id_dict[samfile.gettid(SQ["SN"])]= SQ["SN"]

for i, read in enumerate(samfile.fetch()):
    if id_dict[read.reference_id] in scaffolds:
        samoutput.write(read)

samfile.close()
samoutput.close()


