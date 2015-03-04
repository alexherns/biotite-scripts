#!/usr/bin/env python
from Bio import SeqIO
import sys

if len(sys.argv) == 1:
	print '''
Please see rewrite_contigs.py --help for information
on correct formatting of command.
'''
	exit()
if sys.argv[1] == '-h' or sys.argv[1] == '--help':
	print '''
This program will rewrite the contigs of a .fa file 
containing multiple genomes so that each contig can 
be correctly reassigned to its original genome.

The output contigs will be labeled prefix_[original], 
where prefix is an input string and [original] 
is the sequence identifier as included in the 
original file.

Example usage:
	rewrite_contigs.py <input_file> 'prefix1|M' 
		'prefix2|N' ... <output_file>

<input_file>	Contains contigs to be rewritten
'prefix1|M'...	List of space separated strings 
		where:
	Prefix1 = prefix for each fasta file **(must not have "_" in it)**
	M = Number of contigs in X's genome file
<output_file>	Name of file to write to

'''
	exit()

count = 0
records = SeqIO.parse(sys.argv[1], 'fasta')
add_inputs = sys.argv[2:-1]

genome_list = []
for input in add_inputs:
	id, reps = input.split('|')
	reps = int(reps)
	genome_list.extend([id]*reps)

out_records = []
for i, seq_record in enumerate(records):
	prefix = genome_list[i]
	seq_record.id = prefix + '_' + seq_record.id
	seq_record.description = ''
	count += 1
	out_records.append(seq_record)
	print seq_record

SeqIO.write(out_records, sys.argv[-1], 'fasta')
