#!/usr/bin/env python2.7
import argparse, os

parser = argparse.ArgumentParser(description='Will create a directory to store contigs binned by ESOM mapping.', formatter_class=argparse.ArgumentDefaultsHelpFormatter, add_help=False)

#Required arguments
required = parser.add_argument_group('REQUIRED')
required.add_argument('-c', help= 'esom.class', required=True, type=argparse.FileType('r'))
required.add_argument('-n', help= 'esom.names', required=True, type=argparse.FileType('r'))
required.add_argument('-f', help= 'input.fasta', required=True, type=argparse.FileType('r'))
required.add_argument('-d', help= 'output directory', required=True)

#Optional arguments
optional = parser.add_argument_group('OPTIONAL')
optional.add_argument('-h', action="help", help="show this help message and exit")
optional.add_argument('-s', help= 'Name of class for selection')

args = parser.parse_args()

class_file= args.c
names_file= args.n
fasta_file= args.f
output_dir= args.d
selection= args.s

bins= {}
class_assignments= {}
scaffold_bins= {}
for line in open(class_file):
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

for line in open(names_file):
    if line[0] == "%":
        continue
    seqName, seqHeader= line.strip().split()[:2]
    seqHeader= "_".join(seqHeader.split("_")[:-1])
    scaffold_bins[class_assignments[seqName]].add(seqHeader)

if selection != '':
    temp_scaffold_bins= {}
    temp_scaffold_bins[selection]= scaffold_bins[selection]
    scaffold_bins= temp_scaffold_bins

os.system("mkdir {0}".format(output_dir))
for Bin in scaffold_bins:
    temp_handle= open('temp', 'w')
    temp_handle.write("\n".join(scaffold_bins[Bin]))
    temp_handle.close()
    os.system("cat temp")
    pullseq_command= "pullseq -i {0} -n {1} > {2}/{3}.fasta".format(fasta_file, "temp", output_dir, Bin)
    print pullseq_command
    os.system(pullseq_command)
