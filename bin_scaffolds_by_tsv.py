#!/usr/bin/env python2.7
import sys, os, getopt, argparse

parser = argparse.ArgumentParser(description='Will create FASTA files for bins as provided from a scaffolds2bins.tsv file.', formatter_class=argparse.ArgumentDefaultsHelpFormatter, add_help=False)

#Required arguments
required = parser.add_argument_group('REQUIRED')
required.add_argument('-t', help= 'scaf2bin file', type=str, required=True)
required.add_argument('-f', help= 'input fasta file', type=str, required=True)

#Optional arguments
optional = parser.add_argument_group('OPTIONAL')
optional.add_argument('-h', action="help", help="show this help message and exit")
optional.add_argument('-d', help= 'Output directory', type=str, default= '')
optional.add_argument('-i', help= 'Read input bins from file', type=str, default= '')

args = parser.parse_args()

tsv_file= args.t
fasta_file= args.f
output_dir= args.d

#Format output directory properly
if output_dir != '':
        if output_dir[-1] != '/':
                output_dir+= "/"

#scaffolds2bins file
tsv_handle= open(tsv_file)

#skip header row
tsv_handle.readline()

bin_dict= {}
for line in tsv_handle:
    scaffold, Bin= line.strip().split()[:2]
    if Bin not in bin_dict:
        bin_dict[Bin]= []
    bin_dict[Bin].append(scaffold)

if args.i != '':
    selected_bins = set()
    for line in open(args.i):
        line = line.strip()
        selected_bins.add(line)

for Bin in bin_dict:
    if args.i != '':
        if Bin not in selected_bins:
            continue
    scaffolds= bin_dict[Bin]
    open("temp.list", "w").write("\n".join(scaffolds))
    print "pullseq -i {0} -n temp.list > {2}{1}.fasta".format(fasta_file, Bin, output_dir)
    os.system("pullseq -i {0} -n temp.list > {2}{1}.fasta".format(fasta_file, Bin, output_dir))
