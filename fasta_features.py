#!/usr/bin/env python2.7
import sys

"""
Simple script to print all unique features of desired type from FASTA file.

Example usage: fasta_features.py <input.fasta> genus

Prints out list of all genera predicted in fasta file formatted as per usual from ggkbase
Can be used to select any identifier as output in labeling scheme used by ggkbase
"""

first_pass= [line.strip().split(sys.argv[2])[1][1:] for line in open(sys.argv[1]) if line[0]==">" and sys.argv[2] in line]
second_pass= []
for line in first_pass:
    if len(line.split("="))==1:
        second_pass.append(line)
    else:
        second_pass.append(" ".join(line.split("=")[0].split(" ")[:-1]))
for unique in set(second_pass):
    print unique
