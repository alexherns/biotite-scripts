#!/usr/bin/env python2.7
import pysam, argparse, sys

parser = argparse.ArgumentParser(description='Helps with genome curation using linked paired-end read information from subsetted regions of a genome.', formatter_class=argparse.ArgumentDefaultsHelpFormatter, add_help=False)

#Required arguments
required = parser.add_argument_group('REQUIRED')
required.add_argument('-b', help= 'indexed bam file required for random read access', type=str, required=True)
#required.add_argument('-s', help= 'QNAME sorted SAM file required for mate pair lookups', type=str, required=True)
required.add_argument('-f', help= 'indexed fasta file required for random sequence access', type=str, required=True)

#Optional arguments
optional = parser.add_argument_group('OPTIONAL')
optional.add_argument('-h', '--help', action="help", help="show this help message and exit")
optional.add_argument('-c', help= 'contig/scaffold for read detection', type=str)
optional.add_argument('--min', help= 'minimum on scaffold', type=int, default=0)
optional.add_argument('--max', help= 'maximum on scaffold', type=int, default=sys.maxint)
optional.add_argument('--base', help= 'N-base for coordinate access', type=int, choices= [0,1], default=0)
optional.add_argument('--connectivity', help= 'outputs reads that should span multiple scaffolds and the regions they occur', action="store_false")
optional.add_argument('--threshold', help= 'minimum number of reads bridging scaffolds required to print', type=int, default=0)

args = parser.parse_args()

args.min, args.max= args.min-args.base, args.max-args.base

bamfile = pysam.AlignmentFile(args.b, "rb")
fastafile= pysam.FastaFile(args.f)

sys.stderr.write("Reading BAM file into mate-based dictionary...\t")
read_dict= {}
for read in bamfile:
	read_dict[read.query_name + ("-1" if read.is_read1 else "-2")]= read
sys.stderr.write("finished!\n")

sys.stderr.write("Piling BAM file...\t")
piled_columns = bamfile.pileup(args.c, args.min, (args.max if args.max != sys.maxint else fastafile.get_reference_length(args.c)))
sys.stderr.write("finished!\n")

if args.connectivity:
	for piled_column in piled_columns:
		if piled_column.reference_pos < args.min or piled_column.reference_pos > args.max:
			continue
		if piled_column.reference_pos % 1000 == 0:
			sys.stderr.write("Processing column {0}\n".format(str(piled_column.reference_pos)))
		if args.threshold == 0:
			print "Scaffold: {0}\tPosition: {1}\tReference base: {2}\tDepth: {3}".format(
				bamfile.getrname(piled_column.reference_id),
				piled_column.reference_pos,
				fastafile.fetch(bamfile.getrname(piled_column.reference_id), piled_column.reference_pos, piled_column.reference_pos+1),
				piled_column.nsegments)
			for piled_read in piled_column.pileups:
				seg= piled_read.alignment
				mate= read_dict[seg.query_name + ("-2" if read.is_read1 else "-1")]
				print "\t".join([
					seg.query_name, 
					seg.query_sequence[piled_read.query_position], 
					seg.cigarstring,
					mate.query_name,
					("unmapped" if mate.is_unmapped else bamfile.getrname(mate.reference_id)),
					("NA" if mate.is_unmapped else str(mate.reference_start)),
					("NA" if mate.is_unmapped else mate.cigarstring)
					])
			continue

		
		bridge_count= 0
		segs= []
		mates= []
		for piled_read in piled_column.pileups:
			seg= piled_read.alignment
			if seg.is_proper_pair:
				continue
			mate= read_dict[seg.query_name + ("-2" if read.is_read1 else "-1")]
			if bamfile.getrname(mate.reference_id) != bamfile.getrname(seg.reference_id):
				bridge_count+= 1
			segs.append(seg), mates.append(mate)
	
		if bridge_count >= args.threshold:
			print "Scaffold: {0}\tPosition: {1}\tReference base: {2}\tDepth: {3}".format(
				bamfile.getrname(piled_column.reference_id),
				piled_column.reference_pos,
				fastafile.fetch(bamfile.getrname(piled_column.reference_id), piled_column.reference_pos, piled_column.reference_pos+1),
				piled_column.nsegments)

			print "READ\tBASE_CALL\tMAPPED?"
			for seg, mate in zip(segs, mates):
				print "\t".join([
					seg.query_name, 
					seg.query_sequence[piled_read.query_position], 
					seg.cigarstring,
					mate.query_name,
					("unmapped" if mate.is_unmapped else bamfile.getrname(mate.reference_id)),
					("NA" if mate.is_unmapped else str(mate.reference_start)),
					("NA" if mate.is_unmapped else mate.cigarstring)
					])
			
	
"""
for piled_column in piled_columns:
	if piled_column.reference_pos < args.min or piled_column.reference_pos > args.max:
		continue
	print "Scaffold: {0}\tPosition: {1}\tReference base: {2}\tDepth: {3}".format(
		bamfile.getrname(piled_column.reference_id),
		piled_column.reference_pos,
		fastafile.fetch(bamfile.getrname(piled_column.reference_id), piled_column.reference_pos, piled_column.reference_pos+1),
		piled_column.nsegments)
	print "READ\tBASE_CALL\tMAPPED?"
	for piled_read in piled_column.pileups:
		seg= piled_read.alignment
		pos = bamfile.tell()
		try:
			mate= bamfile.mate(seg)
		except ValueError:
			1
		finally:
			bamfile.seek(pos)
		print "{0}\t{1}\t{2}".format(
			seg.query_name, 
			seg.query_sequence[piled_read.query_position], 
			("unmapped" if seg.mate_is_unmapped else "{0} is on {1}".format(mate.query_name, bamfile.getrname(mate.reference_id))))
	print "\n\n\n"
	"""
