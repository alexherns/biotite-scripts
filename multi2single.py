#!/usr/bin/env python2.7
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
import sys

'''
Concatenates multi-fasta file into single-fasta
Recommended usage with ref2tab.py in order to determine ordering of sequences and their lengths

Usage: multi2single.py <input_fasta> <output_fasta> <label_suffix>
'''


out_seq = ''
for seq_record in SeqIO.parse(open(sys.argv[1], 'r'), 'fasta'):
	out_seq += str(seq_record.seq)
out_record = SeqRecord(Seq(out_seq, generic_dna), id=sys.argv[1]+sys.argv[3], description='')
SeqIO.write(out_record, open(sys.argv[2], 'w'), 'fasta')
