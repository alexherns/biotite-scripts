#!/usr/bin/env python2.7
import sys, pickle
from Bio import SeqIO


def get_bin(line):
    if len(line.split('bin=')) == 1:
        raise Exception('No bin found in line: {0}'.format(line))
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
    try:
        genome= get_bin(seq_record.description)
    except Exception as inst:
        sys.stderr.write(str(inst))
        sys.stderr.write(str(inst.args))
        exit()
    print "{0}\t{1}".format(gene, annotation)
