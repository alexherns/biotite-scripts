#!/usr/bin/env python2.7
import sys

if len(sys.argv) < 4:
    print """
Usage: bin_coverage.py scaffold_hits.tsv header.sam scaffolds2bins.tsv

Output: Coverage per bin.
"""
    exit()

read_size= 150 #read size in bp

#Read the coverage file first
scafCov= {} #an entry maintains information for a scaffold. scafCov[contig]= [hits, length, coverage]
for line in open(sys.argv[1]):
    line= line.strip().split("\t")
    scafCov[line[1]]= [float(line[0])]

#Read the SAM file for line lengths next (or simply SAM header)
for line in open(sys.argv[2]):
    if "@SQ" != line[:3]:
        continue
    line= line.strip().split()
    scaffold= line[1].split("SN:")[1]
    length= float(line[2].split("LN:")[1])
    if scaffold in scafCov:
        scafCov[scaffold].append(length)
        scafCov[scaffold].append(scafCov[scaffold][0]*read_size/length)

#Read the scaffolds2bins file and compute bin length and hits
binCov= {} #an entry maintains information for a bin. binCov[bin]= [hits, length, coverage]
for line in open(sys.argv[3]):
    if "scaffold_name" in line:
        continue
    scaffold, Bin= line.strip().split()[:2]
    if scaffold in scafCov:
        binCov[Bin]= [0., 0., 0.]
        binCov[Bin][0]+= scafCov[scaffold][0]
        binCov[Bin][1]+= scafCov[scaffold][1]

#Compute and print bin coverage
for Bin in binCov:
    coverage= binCov[Bin][0]*read_size/binCov[Bin][1]
    binCov[Bin][2]= coverage
    print "{0}\t{1}".format(Bin, coverage)
