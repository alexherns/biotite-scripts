#!/usr/bin/env python2.7
import argparse, os

parser = argparse.ArgumentParser(description='Second step (muscle alignment) in three-step tree-building pipeline.', formatter_class=argparse.ArgumentDefaultsHelpFormatter, add_help=False)

#Required arguments
required = parser.add_argument_group('REQUIRED')
required.add_argument('-a', help= 'protein alignment', required=True, type=str)

#Optional arguments
optional = parser.add_argument_group('OPTIONAL')
optional.add_argument('-h', action="help", help="show this help message and exit")
optional.add_argument('-f', help= 'minimum fraction of aligned length to retain a sequence', default=0.5, type=float)

args = parser.parse_args()

fraction= args.f
basename= ".".join(args.a.split(".")[:-1])
os.system("java -jar /home/alexh/bin/BMGE-1.12/BMGE.jar -i {0}.afa -m BLOSUM62 -of {0}.BLOSUM62.afa -t AA".format(basename))

count= 0
length= 0
for line in open("{0}.BLOSUM62.afa".format(basename)):
	if line[0]==">":
		count+=1
	else:
		length+= len(line.strip())
	if count>1:
		break
	
min_length= int(length*fraction)
print """
The minimum aligned length for processing is {0}
""".format(min_length)

os.command("pullseq -i {0}.BLOSUM62.afa -l {1} > {0}.BLOSUM62.trimmed.afa".format(basename, min_length))
os.command("FastTree {0}.BLOSUM62.trimmed.afa > {0}.BLOSUM62.trimmed.tree".format(basename))

print """
FastTree finished.

For RAxML, run ProteinModelSelection.pl as:
perl ~/bin/ProteinModelSelection.pl {0}.BLOSUM62.trimmed.afa
...

and then run RAxML as:
raxmlHPC-PTHREADS-SSE3 -f a -m <BEST_MODEL> -p 12345 -x 12345 -# 100 -s {0}.BLOSUM62.trimmed.afa -n T20 -T 10
"""