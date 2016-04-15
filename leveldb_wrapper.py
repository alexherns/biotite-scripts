#!/usr/bin/env python
import leveldb as leveldb
import sys
import argparse, os

parser = argparse.ArgumentParser(description='Simple generation of leveldb hash table from tab-separated table', 
        formatter_class=argparse.ArgumentDefaultsHelpFormatter, add_help=False)

#Required arguments
required = parser.add_argument_group('REQUIRED')
required.add_argument('-t', help= 'input tsv', required=True, type=str)
required.add_argument('-d', help= 'output leveldb', required=True, type=str)

#Optional arguments
optional = parser.add_argument_group('OPTIONAL')
optional.add_argument('-h', action="help", help="show this help message and exit")
optional.add_argument('-b', type = int, default = 10000, help = 'batch size')

args = parser.parse_args()

db = leveldb.LevelDB(args.d)

count = 0
batch = leveldb.WriteBatch()
for line in open(args.t):
    key, value= line.split('\t')[:2]
    value= value.strip()
    batch.Put(key, value)
    count += 1
    if count % args.b == 0:
        db.Write(batch, sync = False)
        batch = leveldb.WriteBatch()
    if count % 1000000 == 0:
        sys.stderr.write('Processed {0} lines\n'.format(count))

if count % args.b:
    db.Write(batch, sync = False)
