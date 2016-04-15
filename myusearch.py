#!/usr/bin/env python
import subprocess
import sys
import os

DB_SHORTCUTS = {
        'nr': '/work/blastdb/nr.fa.udb',
        'refseq': '/work/blastdb/refseq_protein.udb'
        }
USEARCH_FORMATS = [
        'alnout',
        'blast6out',
        'dbmatched'
        ]


def gen_usearch_cmd(args):
    base = ['usearch7.0.1001_i86linux64', '-ublast']
    base.append(args.f)
    base.extend(['-db', args.db])
    base.extend(['-evalue', args.evalue])
    if args.maxhits:
        base.extend(['-maxhits', str(args.maxhits)])
    base.extend(['-threads', str(args.t)])
    base.extend(['-{0}'.format(args.format), args.o])
    return base


def run_usearch(args):
    usearch_cmd = gen_usearch_cmd(args)
    print ' '.join(usearch_cmd)
    p = subprocess.Popen(usearch_cmd)
    p.wait()


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description = \
            'simple generation of usearch commands')
    parser.add_argument(\
            '-f', type = str, required = True, \
                help = 'fasta (required)')
    parser.add_argument(\
            '-o', type = str, required = True, \
            help = 'output results file')
    parser.add_argument(\
            '-t', type = int, default = 10, \
            help = '# threads for usearch')
    parser.add_argument(\
            '--db', type = str, default = 'nr', \
            help = 'path to udb file')
    parser.add_argument(\
            '--evalue', type = str, default = '1e-8', \
            help = 'e-value threshold for reporting hits')
    parser.add_argument(\
            '--maxhits', type = int, default = 0, \
            help = 'maximum number hits to report')
    parser.add_argument(\
            '--format', type = str, default = 'blast6out', \
            help = 'output format type')
    parser.add_argument(\
            '--verbose', action = 'store_true', \
            help = 'print progress metrics to stderr')
    parser.add_argument(\
            '--log', type = str, default = 'log.txt', \
            help = 'log file for stderr logging when not run with verbose')
    parser.add_argument(\
            '--clean', action = 'store_true', \
            help = 'clean up temporary files after hits found')
    args = parser.parse_args()
    if args.db in DB_SHORTCUTS:
        args.db = DB_SHORTCUTS[args.db]
    if args.format not in USEARCH_FORMATS:
        sys.exit('Output format {0} not one of [{1}]'.format(args.format, ', '.join(USEARCH_FORMATS)))
    run_usearch(args)





