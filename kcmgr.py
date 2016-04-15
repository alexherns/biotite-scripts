#!/usr/bin/env python
import kyotocabinet as kc
import sys
import argparse
import os
import glob

kchashmgr_help = \
'''
usage:
  kchashmgr create [-otr] [-onl|-otl|-onr] [-apow num] [-fpow num] [-ts] [-tl] [-tc] [-bnum num] path
  kchashmgr inform [-onl|-otl|-onr] [-st] path
  kchashmgr set [-onl|-otl|-onr] [-add|-rep|-app|-inci|-incd] [-sx] path key value
  kchashmgr remove [-onl|-otl|-onr] [-sx] path key
  kchashmgr get [-onl|-otl|-onr] [-rm] [-sx] [-px] [-pz] path key
  kchashmgr list [-onl|-otl|-onr] [-max num] [-rm] [-sx] [-pv] [-px] path [key]
  kchashmgr clear [-onl|-otl|-onr] path
  kchashmgr import [-onl|-otl|-onr] [-sx] path [file]
  kchashmgr copy [-onl|-otl|-onr] path file
  kchashmgr dump [-onl|-otl|-onr] path [file]
  kchashmgr load [-otr] [-onl|-otl|-onr] path [file]
  kchashmgr defrag [-onl|-otl|-onr] path
  kchashmgr setbulk [-onl|-otl|-onr] [-sx] path key value ...
  kchashmgr removebulk [-onl|-otl|-onr] [-sx] path key ...
  kchashmgr getbulk [-onl|-otl|-onr] [-sx] [-px] path key ...
  kchashmgr check [-onl|-otl|-onr] path
'''

KC_UTILITIES = [
    'create',
    'import',
    'get',
    'dump',
    'getbulk'
]

def open_multi(db_names, mode = kc.DB.OREADER):
    dbs = []
    for db_name in db_names:
        db= kc.DB()
        if not db.open(db_name, mode):
            sys.stderr.write("opening error: " + str(db.error()))
        dbs.append(db)
    return dbs

def close_multi(dbs):
    for db in dbs:
        if not db.close():
            sys.stderr.write('close error: ' + str(db.error()))

def distribute(dbs, fread, batch_size = 1000, delim = '\t'):
    num_dbs = len(dbs)
    line_batch = [{} for _ in xrange(num_dbs)]
    lines_read = 0
    for line in fread:
        lines_read += 1
        split = line.strip().split(delim)
        key = split[0]
        val = delim.join(split[1:])
        db_num = hash(key) % num_dbs
        this_map = line_batch[db_num]
        this_map[key] = val
        if len(this_map) >= batch_size:
            dbs[db_num].set_bulk(this_map)
            line_batch[db_num] = {}
        if lines_read % 100000 == 0:
            sys.stderr.write('Processed: {0}\n'.format(lines_read))
    for db_num in xrange(num_dbs):
        this_map = line_batch[db_num]
        if len(this_map) > 0:
            dbs[db_num].set_bulk(this_map)

def iter_batch(dbs, fread):
    num_dbs = len(dbs)
    for line in fread:
        query = line.strip()
        db_num = hash(query) % num_dbs
        try:
            val = dbs[db_num].get(query)
        except KeyError:
            val = None
        yield query, val

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Tool for managing distributed key-value stores using kyotocabinet', 
            formatter_class=argparse.ArgumentDefaultsHelpFormatter, add_help=False)

    #Required arguments
    required = parser.add_argument_group('REQUIRED')
    required.add_argument(\
            '-d', type = str, required = True, \
            help = 'database base-name for read/write')
    required.add_argument(\
            '-p', type = str, required = True, choices = KC_UTILITIES, \
            help = 'program')

    #Optional arguments
    optional = parser.add_argument_group('OPTIONAL')
    optional.add_argument(\
            '-h', action = 'help', \
            help = 'show this help message and exit')
    optional.add_argument(\
            '--tsv', type = str, \
            help= 'tab-delimited file')
    optional.add_argument(\
            '-k', type = str, \
            help = 'single key for database query')
    optional.add_argument(\
            '-n', type = str, \
            help = 'keys file for database query')
    optional.add_argument(\
            '-N', action = 'store_true', \
            help = 'read keys from stdin')
    optional.add_argument(\
            '-b', type = int, default = 10000, \
            help = 'batch size for reads and writes')
    optional.add_argument(\
            '--count', type = int, default = 100, \
            help = 'number of hash-tables to split over')
    optional.add_argument(\
            '--delim', type = str, default = '\t', \
            help = 'delimiter for input file')

    #Argument parsing
    args = parser.parse_args()
    db_basename = args.d
    program = args.p
    delim = args.delim
    tsv = args.tsv
    query = args.k
    batch_size = args.b

    #Program selection logic
    if program not in KC_UTILITIES:
        parser.error('Please choose one of the available utilities')

    elif program in ['create', 'import']:

        db_names = ['{0}.{1}.kch'.format(db_basename, num) for num in xrange(args.count)]
        dbs = open_multi(db_names, kc.DB.OWRITER | kc.DB.OCREATE)

        try:
            if program == 'import':
                fread = open(tsv)
                distribute(dbs, fread, batch_size = batch_size, delim = delim)

        finally:
            close_multi(dbs)

    elif program in ['get', 'dump', 'getbulk']:

        num_dbs = len(glob.glob('{0}.*.kch'.format(db_basename)))
        db_names = ['{0}.{1}.kch'.format(db_basename, num) for num in xrange(num_dbs)]
        dbs = open_multi(db_names, kc.DB.OREADER)

        try:

            if program == 'get':
                if not query:
                    parser.error('Please supply a key for querying the database')
                value = dbs[hash(query)%len(dbs)].get(query)
                if value != None:
                    print value
                else:
                    raise KeyError, '{0} not found in {1}'.format(query, db_basename)

            elif program == 'getbulk':
                if args.n:
                    fread = open(args.n)
                elif args.N:
                    fread = sys.stdin
                for key, val in iter_batch(dbs, fread):
                    if val != None:
                        print '{0}\t{1}'.format(key, val)

            elif program == 'dump':
                raise Exception, 'This seems dangerous... and you may not like the result. Not implemented yet.'

        finally:
            close_multi(dbs)
