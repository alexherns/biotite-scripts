#!/usr/bin/env python
import sys, getopt

def usage():
	print '''
Will take in a sam file and bin reads according to reference genome to which they best map.
Outputs read files in FASTQ format.
Requires SAM file to be built from reference contigs with sequence headers in format:
	>pseudoX_N
	ACTG...
where X is the number of the isolate, and N is a unique string which identifies the contig

Usage: assign2bin.py <sam_file>
'''
	exit()

opts, args = getopt.getopt(sys.argv[1:], 'h', ['help'])
if opts == []:
	print '''
This script requires more arguments. Please pass -h or --help for assistance
'''

for o, a in opts:
	if o in ('-h', '--help'):
		usage()

from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

time_ints = [time.time()]
samfile = pysam.AlignmentFile(sys.argv[1], 'r')

genomes = set([line.split()[1].split('pseudo')[1].split('_')[0] for line in samfile.text.strip().split('\n') if line.split()[0] == '@SQ'])
read_by_genome = {}
for genome_num in genomes:
	read_by_genome[genome_num] = []

print 'Processing reads...'
for read in samfile.fetch():
	if not read.is_unmapped:
		contig = samfile.getrname(read.reference_id)
		genome_num = contig.split('pseudo')[1].split('_')[0]
		if not read.is_reverse:
			seq_obj = Seq(read.query_sequence, generic_dna)
			seq_record = SeqRecord(seq_obj, id = read.query_name, description = '', letter_annotations = {'phred_quality': read.query_qualities.tolist()})
#			print seq_record.letter_annotations
		else:
			seq_obj = Seq(read.query_sequence, generic_dna).reverse_complement()
			seq_record = SeqRecord(seq_obj, id = read.query_name, description = '', letter_annotations = {'phred_quality': read.query_qualities.tolist()[::-1]})
		read_by_genome[genome_num].append(seq_record)

time_ints.append(time.time())
print 'Processing of {1} reads finished in {0}s'.format(str(time_ints[-1]-time_ints[-2]), str(sum([len(read_by_genome[key]) for key in read_by_genome]))) 
samfile.close()

for key in read_by_genome:
	out_file = open(sys.argv[1][:-3]+'reads_'+key+'.fastq', 'w')
	print 'Writing files to {0}...'.format(sys.argv[1][:-3]+'reads_'+key+'.fastq')
	SeqIO.write(read_by_genome[key], out_file, 'fastq')
#	out_file = open(sys.argv[1][:-3]+'reads_'+key+'.txt', 'w')
#	out_file.write('\n'.join(read_by_genome[key]))
	out_file.close()
	time_ints.append(time.time())
	print '\t{0}s'.format(str(time_ints[-1]-time_ints[-2]))
