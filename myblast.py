#!/usr/bin/env python
import subprocess
import sys
import os
import multiprocessing
from common_objects import parse_fasta, write_fastas

REMOTE_DBS = [
        '16SMicrobial',
        'Representative_Genomes',
        'env_nr',
        'env_nt',
        'nr',
        'nt',
        'pdbaa',
        'refseq_protein',
        'refseq_genomic',
        'swissprot',
        'wgs'
        ]

LOCAL_DBS = {
        'nr':'/work/blastdb/nr',
        'nt':'/work/blastdb/nt',
        'refseq':'/work/blastdb/refseq_protein',
        'env_nr':'/work/blastdb/env_nr',
        'env_nt':'/work/blastdb/env_nt'
        }

def create_batches(fasta_handle, batch_size):
    count = 0
    fastas = []
    batch_num = 0
    batch_names = []
    for fasta in parse_fasta(fasta_handle):
        count += 1
        fastas.append(fasta)
        if count >= batch_size:
            batch_name = make_batch_name(fasta_handle.name, batch_num)
            batch_names.append(batch_name)
            batch_handle = open(batch_name, 'wb')
            write_fastas(batch_handle, fastas)
            count = 0
            fastas = []
            batch_num += 1
    if count > 0:
        batch_name = make_batch_name(fasta_handle.name, batch_num)
        batch_names.append(batch_name)
        batch_handle = open(batch_name, 'wb')
        write_fastas(batch_handle, fastas)
        batch_num += 1
    return batch_names

def make_batch_name(file_name, batch_num):
    return '{0}.{1}'.format(file_name, batch_num)

def gen_blast_cmd(fasta, args):
    base = ['blastp']
    base.extend(['-query', fasta])
    base.extend(['-db', args.db])
    base.extend(['-evalue', args.evalue])
    base.extend(['-max_target_seqs', args.maxhits])
    base.extend(['-outfmt', args.format])
    base.extend(['-out', args.o])
    if args.remote:
        base.append('-remote')
    else:
        base.extend(['-num_threads', args.threads])
    base = map(str, base)
    return base


def run_blast(fasta, args):
    blast_cmd = gen_blast_cmd(fasta, args)
    sys.stderr.write(' '.join(map(str, blast_cmd)) + '\n')
    p = subprocess.Popen(blast_cmd)
    return p

def merge_results(batch_names):
    pass


def clean_batches(batch_names):
    for batch in batch_names:
        if os.path.isfile(batch):
            os.remove(batch)

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description = \
            'simple generation of blast searches')
    parser.add_argument(\
            '-f', type = str, required = True, \
                help = 'fasta (required)')
    parser.add_argument(\
            '-o', type = str, required = True, \
            help = 'output results file (required)')
    parser.add_argument(\
            '--db', type = str, default = 'nr', \
            help = 'search database')
    parser.add_argument(\
            '--remote', action = 'store_true', \
            help = 'run blast on remote server')
    parser.add_argument(\
            '--threads', type = int, default = 10, \
            help = 'threads for search (local only)')
    parser.add_argument(\
            '--evalue', type = str, default = '1e-8', \
            help = 'e-value threshold for reporting hits')
    parser.add_argument(\
            '--maxhits', type = int, default = 100, \
            help = 'maximum number hits to report')
    parser.add_argument(\
            '--format', type = int, default = 6, \
            help = 'output format type')
    parser.add_argument(\
            '--batchsize', type = int, default = 5, \
            help = 'number of proteins to search per submission (remote only)')
    parser.add_argument(\
            '--concurrent', type = int, default = 5, \
            help = 'number of concurrent batches to submit (remote only)')
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
    if args.remote:
        sys.exit('Remote BLAST not working yet! Please check back soon :)')
        if args.db not in BLAST_DBS:
            sys.exit('Please select from these databases:\n\t' + '\n\t'.join(REMOTE_DBS))
        batch_names = create_batches(open(args.f, 'r'), args.batchsize)
        for batch in batch_names:
            run_blast(batch, args).wait()
        merge_results(batch_names, args.o)
        if args.clean:
            clean_batches(batch_names)
    else:
        assert 0 < args.threads < 20, '--threads must be integer between 0 and 20'
        if args.db not in LOCAL_DBS:
            sys.exit('Please select from these databases:\n\t' + '\n\t'.join(LOCAL_DBS.keys()))
        args.db = LOCAL_DBS[args.db]
        run_blast(args.f, args).wait()
