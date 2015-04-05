#!/usr/bin/env python2.7
import sys, os, getopt

opts, args = getopt.getopt(sys.argv[1:], 't:d:f:h', ['help', 'tsv', 'dir', 'fasta'])

def usage():
        print """
Will create FASTA files for bins as provided from a scaffolds2bins.tsv file

Usage: -t file.tsv -f input.fasta [OPTIONS]

OPTIONS:
-d, --dir=      Output directory

"""
        exit()

tsv_file= ''
fasta_file= ''
output_dir= ''

for o, a in opts:
    if o in ('-h', '--help'):
        usage()
    if o in ('-f', '--fasta='):
        fasta_file= a
    if o in ('-t', '--tsv='):
        tsv_file= a
    if o in ('-d', '--dir='):
        output_dir= a

if '' in [tsv_file, fasta_file]:
        print """
Please specify appropriate parameters. See -h (--help) for help
"""
        exit()

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

for Bin in bin_dict:
    scaffolds= bin_dict[Bin]
    open("temp.list", "w").write("\n".join(scaffolds))
    os.system("pullseq -i {0} -n temp.list > {2}{1}.fasta".format(fasta_file, Bin, output_dir))
