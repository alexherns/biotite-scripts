#!/usr/bin/env python2.7
import sys, getopt, os

opts, args = getopt.getopt(sys.argv[1:], 'c:n:d:f:h', ['help', 'class', 'names', 'dir', 'fasta'])

def usage():
        print """
Can retrieve scaffolds according to gc and coverage specifications.
Currently only works with *fa.genes format as provided in assembly.d/sample/idba_ud_full folder

Usage: pullscaffolds.py -n esom.names -c esom.class -f input.fasta -d output_directory  [OPTIONS]

OPTIONS:

"""
        exit()

class_file= ''
names_file= ''
fasta_file= ''
output_dir= ''

for o, a in opts:
    if o in ('-h', '--help'):
        usage()
    elif o in ('-n', '--names'):
        names_file= a
    elif o in ('-c', '--class'):
        class_file= a
    elif o in ('-d', '--dir'):
        output_dir= a
    elif o in ('-f', '--fasta'):
        fasta_file= a

if '' in [class_file, names_file, fasta_file, output_dir]:
        print """
Please specify appropriate parameters. See -h (--help) for help
"""
        exit()

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

os.system("mkdir {0}".format(output_dir))
for Bin in scaffold_bins:
    temp_handle= open('temp', 'w')
    temp_handle.write("\n".join(scaffold_bins[Bin]))
    temp_handle.close()
    os.system("cat temp")
    pullseq_command= "pullseq -i {0} -n {1} > {2}/{3}.fasta".format(fasta_file, "temp", output_dir, Bin)
    print pullseq_command
    os.system(pullseq_command)
