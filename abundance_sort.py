#!/usr/bin/env python2.7
import sys, operator, argparse
from Bio import SeqIO

parser = argparse.ArgumentParser(description='Sorts fasta features according to coverage information.', formatter_class=argparse.ArgumentDefaultsHelpFormatter, add_help=False)

#Required arguments
required = parser.add_argument_group('REQUIRED')
required.add_argument('-f', help= 'FASTA file containing features of interest', required=True, type=argparse.FileType('r'))
required.add_argument('-t', help= 'TSV containing features as downloaded from ggkbase', required=True, type=argparse.FileType('rU'))

#Optional arguments
optional = parser.add_argument_group('OPTIONAL')
optional.add_argument('-h', action="help", help="show this help message and exit")
optional.add_argument('--scaf_col', help= '0-based column with scaffold_gene', type=int, default=1)
optional.add_argument('--cov_col', help= '0-based column with coverage values', type=int, default=4)

args = parser.parse_args()

#Create a dictionary of feature:coverage values
#Read in the tsv of features
handle= args.t
feat2cov= {}
for line in handle:
    contig_features= line.strip().split("\t")
    feature, coverage= contig_features[args.scaf_col], contig_features[args.cov_col]
    feat2cov[feature]= float(coverage)
handle.close()

#Read in the fasta file of features
feat2seq= {}
record_dict= SeqIO.to_dict(SeqIO.parse(args.f, "fasta"))

#Print the features sorted by coverage value
for feature_tuple in sorted(feat2cov.items(), key=operator.itemgetter(1))[::-1]:
    print record_dict[feature_tuple[0]].format("fasta").strip()
