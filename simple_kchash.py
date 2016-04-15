#!/usr/bin/env python
import kyotocabinet as kc
import sys
import argparse, os

parser = argparse.ArgumentParser(description='Simple generation of kyotocabinet hash table from tab-separated table', 
        formatter_class=argparse.ArgumentDefaultsHelpFormatter, add_help=False)

#Required arguments
required = parser.add_argument_group('REQUIRED')
required.add_argument('-t', help= 'input tsv', required=True, type=str)
required.add_argument('-k', help= 'output kyotocabinet', required=True, type=str)

#Optional arguments
optional = parser.add_argument_group('OPTIONAL')
optional.add_argument('-h', action="help", help="show this help message and exit")
optional.add_argument('-b', type = int, default = 1000, help = 'bulk size')

args = parser.parse_args()

db= kc.DB()
if not db.open(args.k, db.OWRITER | db.OCREATE):
    sys.stderr.write("opening error: " + str(db.error()))

count = 0
d = {}
bulk = 0
for line in open(args.t):
    key, value= line.split('\t')[:2]
    value= value.strip()
    d[key] = value
    bulk += 1
    count += 1
    if bulk % args.b == 0:
        db.set_bulk(d)
        d = {}
        bulk = 0
    if not count % 1000000:
        sys.stderr.write('Processed {0} lines\n'.format(count))

if bulk:
    db.set_bulk(d)

if not db.close():
    sys.stderr.write("close error: " + str(db.error()))
