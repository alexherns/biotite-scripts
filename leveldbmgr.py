#!/usr/bin/env python
import leveldb
import sys
import argparse
import os

LDB_UTILITIES = [
        'import',
        'get',
        'getbulk'
        ]

def ldb_writer(db, tsv_handle, batch_size = 10000, delim = '\t'):
    count = 0
    batch = leveldb.WriteBatch()
    for line in tsv_handle:
        key = line.split(delim)[0]
        value = delim.join(line.split(delim)[1:])
        batch.Put(key, value)
        count += 1
        if count % batch_size == 0:
            db.Write(batch, sync = False)
            batch = leveldb.WriteBatch()
        if count % 1000000 == 0:
            sys.stderr.write('Processed {0} lines\n'.format(count))
    if count % args.b:
        db.Write(batch, sync = False)

def ldb_get(db, key):
    try:
        return db.Get(key)
    except KeyError:
        raise KeyError, '{0} not found in db!'.format(key)

def ldb_get_iter(db, keys):
    for key in keys:
        yield key, ldb_get(db, key)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Simple generation of leveldb hash table from tab-separated table', 
            formatter_class=argparse.ArgumentDefaultsHelpFormatter, add_help=False)

    #Required arguments
    required = parser.add_argument_group('REQUIRED')
    required.add_argument(\
            '-d', type = str, required = True, \
            help = 'leveldb for read/write')
    required.add_argument(\
            '-p', type = str, required = True, \
            help = 'program')

    #Optional arguments
    optional = parser.add_argument_group('OPTIONAL')
    optional.add_argument(\
            '-h', action = 'help', \
            help = 'show this help message and exit')
    optional.add_argument(\
            '-t', type = str, \
            help= 'input file')
    optional.add_argument(\
            '-k', type = str, \
            help = 'key for database query')
    optional.add_argument(\
            '-b', type = int, default = 10000, \
            help = 'batch size for writes')
    optional.add_argument(\
            '--delim', type = str, default = '\t', \
            help = 'delimiter for input file')

    args = parser.parse_args()
    db = leveldb.LevelDB(args.d)
    if args.p not in LDB_UTILITIES:
        sys.exit('Please choose one of the available utilities')
    if args.p == 'import':
        tsv_handle = open(args.t)
        ldb_writer(db, tsv_handle, batch_size = args.b, delim = args.delim)
    if args.p == 'get':
        if not args.k:
            sys.exit('Please supply a key for querying the database')
        print ldb_get(db, args.k)
    if args.p == 'getbulk':
        pass
