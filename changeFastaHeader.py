#!/usr/bin/env python
from Bio import SeqIO
import sys

'''
Appends all sequence identifiers in a fasta file with a desired
string.

Example usage: rewriteFastaHeaders.py <input.fasta> <string> <delimiter>
If no delimiter is provided, default to single pipe "|"
Output is to stdout
'''


records = SeqIO.parse(sys.argv[1], 'fasta')
identifier = sys.argv[2]
if len(sys.argv) == 3:
	delimiter = '|'
else:
	delimiter = sys.argv[3]

for seq_record in records:
	seq_record.description = identifier + delimiter + seq_record.description
	print '>{0}\n{1}'.format(seq_record.description, str(seq_record.seq))


