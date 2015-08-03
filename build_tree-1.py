#!/usr/bin/env python2.7
import argparse, os

parser = argparse.ArgumentParser(description='Step 1 (initial alignment against database) in three-step tree-building pipeline.', formatter_class=argparse.ArgumentDefaultsHelpFormatter, add_help=False)

#Required arguments
required = parser.add_argument_group('REQUIRED')
required.add_argument('-f', help= 'input.fasta', required=True, type=str)
required.add_argument('-d', help= 'protein database for tree', required=True, type=str)

#Optional arguments
optional = parser.add_argument_group('OPTIONAL')
optional.add_argument('-h', action="help", help="show this help message and exit")

args = parser.parse_args()

	
database_locs= {
"mcrA": "/data3/borehole/genes/mcrA/mcrAdb2/gis.cleaned.fasta",
"rubisco": "sample",
"dsrA": "sample"
}

database= args.d
if database not in database_locs:
	print "Database not available. Please try:\n\t{0}".format("\n\t".join(database_locs.keys()))
	exit()
database_location= database_locs[database]

base_name= ".".join(args.f.split(".")[:-1])
os.system("cat {0} {1} > {2}.wr.fasta".format(database_location, args.f, base_name))
os.system("fix_header.py {0} > {1}.wr.clean.fasta".format(args.f, base_name))
os.system("mafft --thread 6 {0}.wr.clean.fasta > {0}.wr.clean.afa".format(base_name))

print """
Initial steps finished... please manually select sequences or proceed directly to next step
"""