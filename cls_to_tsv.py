#!/usr/bin/env python2.7
import argparse, os

parser = argparse.ArgumentParser(description='Writes tsv assignments from an ESOM class file', formatter_class=argparse.ArgumentDefaultsHelpFormatter, add_help=False)

#Required arguments
required = parser.add_argument_group('REQUIRED')
required.add_argument('-c', help= 'esom.class', required=True, type=argparse.FileType('r'))
required.add_argument('-n', help= 'esom.names', required=True, type=argparse.FileType('r'))
required.add_argument('-f', help= 'input.fasta', required=True, type=str)

#Optional arguments
optional = parser.add_argument_group('OPTIONAL')
optional.add_argument('-h', action="help", help="show this help message and exit")

args = parser.parse_args()

class_file= args.c
names_file= args.n
fasta_file= args.f

bins= {}
class_assignments= {}
scaffold_bins= {}
for line in class_file:
    if line[0] == "%" and len(line.strip().split()) > 1:
        binID, classBin= line.strip().split()[:2] 
        binID= binID[1:]
        bins[binID]= {}
        bins[binID]["bin"]= classBin
    elif line[0] != "%":
        seqName, binID= line.strip().split()
        class_assignments[seqName]= bins[binID]["bin"]
        if class_assignments[seqName] not in scaffold_bins:
            scaffold_bins[class_assignments[seqName]]= set()

for line in names_file:
    if line[0] == "%":
        continue
    seqName, seqHeader= line.strip().split()[:2]
    seqHeader= "_".join(seqHeader.split("_")[:-1])
    scaffold_bins[class_assignments[seqName]].add(seqHeader)

for bin_name in scaffold_bins:
    print "\n".join(["\t".join((scaffold, bin_name)) for scaffold in scaffold_bins[bin_name]])
