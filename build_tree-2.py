#!/usr/bin/env python2.7
import argparse, os

parser = argparse.ArgumentParser(description='Second step (muscle alignment) in three-step tree-building pipeline.', formatter_class=argparse.ArgumentDefaultsHelpFormatter, add_help=False)

#Required arguments
required = parser.add_argument_group('REQUIRED')
required.add_argument('-f', help= 'input.fasta', required=True, type=str)

#Optional arguments
optional = parser.add_argument_group('OPTIONAL')
optional.add_argument('-h', action="help", help="show this help message and exit")

args = parser.parse_args()
	
os.system("muscle -in {0} -out {1}.afa".format(args.f, ".".join(args.f.split(".")[:-1])))

print """
Muscle alignment finished... please manually select sequences and trim alignment or proceed directly to next step
"""