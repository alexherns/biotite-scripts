#!/usr/bin/env python
import argparse

parser = argparse.ArgumentParser(description='Will take in a sam file and bin reads according to reference genome to which they best map. NOTE: CURRENTLY BROKEN, ONLY WORKS FOR PSEUDO POP', formatter_class=argparse.ArgumentDefaultsHelpFormatter, add_help=False)

#Required arguments
required = parser.add_argument_group('REQUIRED')
required.add_argument('-s', help= 'mapping.sam', required=True, type=string)

#Optional arguments
optional = parser.add_argument_group('OPTIONAL')
optional.add_argument('-h', action="help", help="show this help message and exit")

args = parser.parse_args()

from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

time_ints = [time.time()]
samfile = pysam.AlignmentFile(args.s, 'r')

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
	out_file = open(args.c[:-3]+'reads_'+key+'.fastq', 'w')
	print 'Writing files to {0}...'.format(args.c[:-3]+'reads_'+key+'.fastq')
	SeqIO.write(read_by_genome[key], out_file, 'fastq')
	out_file.close()
	time_ints.append(time.time())
	print '\t{0}s'.format(str(time_ints[-1]-time_ints[-2]))
