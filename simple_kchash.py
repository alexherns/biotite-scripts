#!/usr/bin/env python
import kyotocabinet as kc
import sys
import argparse, os

parser = argparse.ArgumentParser(description='Will create a directory to store contigs binned by ESOM mapping.', formatter_class=argparse.ArgumentDefaultsHelpFormatter, add_help=False)

#Required arguments
required = parser.add_argument_group('REQUIRED')
required.add_argument('-t', help= 'input tsv', required=True, type=str)
required.add_argument('-k', help= 'output kyotocabinet', required=True, type=str)

#Optional arguments
optional = parser.add_argument_group('OPTIONAL')
optional.add_argument('-h', action="help", help="show this help message and exit")

args = parser.parse_args()

db= kc.DB()
if not db.open(args.k, db.OWRITER | db.OCREATE):
    sys.stderr.write("opening error: " + str(db.error()))

for line in open(args.t):
    key, value= line.split('\t')[:2]
    value= value.strip()
    db.set(key, value)

if not db.close():
    sys.stderr.write("close error: " + str(db.error()))
