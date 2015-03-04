#!/usr/bin/env python2.7
import sys, operator
from Bio import SeqIO

if len(sys.argv)<3 or sys.argv[1]=="-h" or sys.argv[1]=="--help":
    print """
Usage: abundance_sort.py <features.fasta> <features.tsv>

TSV of features as downloaded from ggkbase.
Scaffold_gene is in column 2.
Coverage value is in column 5.
"""
    exit()

#Create a dictionary of feature:coverage values
#Read in the tsv of features
handle= open(sys.argv[2], "r")
feat2cov= {}
for line in handle:
    contig_features= line.strip().split("\t")
    feature, coverage= contig_features[1], contig_features[4]
    feat2cov[feature]= float(coverage)
handle.close()

#Read in the fasta file of features
feat2seq= {}
handle= open(sys.argv[1], "rU")
record_dict= SeqIO.to_dict(SeqIO.parse(handle, "fasta"))

#Print the features sorted by coverage value
for feature_tuple in sorted(feat2cov.items(), key=operator.itemgetter(1))[::-1]:
    print record_dict[feature_tuple[0]].format("fasta").strip()
