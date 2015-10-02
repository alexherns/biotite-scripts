#!/usr/bin/env python2.7
import sys, pickle
from Bio import SeqIO


def get_bin(line):
    genome= line.split('bin=')[1].split()[0]
    return genome

def get_annotation(line):
    annotation= ' '.join(line.split(' id=')[0].split(' ')[1:])
    if annotation == ' ' or not annotation:
        return 'unk'
    return annotation

def get_scaffold(line):
    return '_'.join(line.split(' ')[0].split('_')[:-1])

gene_d= {}
for seq_record in SeqIO.parse(sys.argv[1], 'fasta'):
    gene= seq_record.id
    annotation= get_annotation(seq_record.description)
    genome= get_bin(seq_record.description)
    print "{0}\t{1}".format(gene, annotation)
