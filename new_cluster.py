#!/usr/bin/env python

import networkx as nx
import os
import argparse
import subprocess

def ani(fastas, sim_threshold, mash_file, threads = 6):
    """
    Use mash to estimate ANI of genomes
    Input:
        fastas <list[str]>      -- paths to fasta files
        sim_threshold <float>   -- fractional cutoff % identity for cluster joining
        mash_file <str>         -- pasted sketch file of all fastas being compared
        threads <int>           -- number threads for distance estimation
    """
    ANI = nx.Graph()
    # use Mash to estimate ANI
    for fasta in fastas:
        indiv_mash = fasta + '.msh'
        if os.path.isfile(indiv_mash):
            cmp_file = indiv_mash
        else:
            cmp_file = fasta
        mash_cmd = ['mash', 'dist', '-p', str(threads), cmp_file,
                    mash_file]
        process = subprocess.Popen(mash_cmd, stdout = subprocess.PIPE)
        for pair in process.communicate()[0].splitlines():
            a, b, dist, p, shared = pair.strip().split()
            p = float(p)
            similarity = (1 - float(dist)) * 100
            if similarity >= sim_threshold:
                ANI.add_edge(a, b, si = similarity, pval = p, sharedK = shared)
        process.wait()
    return ANI

def make_mashes(fastas, mash_file, threads = 6, kmer = 21, force = False):
    """
    Create mash files for multiple fasta files
    Input:
        fastas <list[str]>  -- paths to fasta files
        mash_file <str>     -- path to output mash file 
        threads <int>       -- # threads for parallelization
        kmer <int>          -- kmer size for mash sketching
        force <boolean>     -- force overwrite of all mash files
    """
    mash_processes = set()
    sketches = [fasta + '.msh' for fasta in fastas]
    # Perform the sketching
    for fasta, sketch in zip(fastas, sketches):
        if os.path.isfile(sketch):
            continue
        mash_cmd = ['mash', 'sketch', '-o', fasta, '-k', str(kmer), fasta]
        mash_processes.add(subprocess.Popen(mash_cmd))
        if len(mash_processes) >= threads:
            os.wait()
            mash_processes.difference_update([mp for mp in mash_processes if mp.poll() is not None])
    # Collect stragglers
    for mp in mash_processes:
        if mp.poll() is None:
            mp.wait()
    # Paste sketches into single mash
    paste_mashes(sketches, mash_file, force = force)
    return


def paste_mashes(sketches, pasted_mash, force = False):
    """
    Combine mash files into single sketch
    Input:
        sketches <list[str]>  -- paths to sketch files
        pasted_mash <str>     -- path to output mash file 
        force <boolean>     -- force overwrite of all mash file
    """
    if os.path.isfile(pasted_mash):
        if force:
            subprocess.Popen(['rm', pasted_mash]).wait()
        else:
            return
    pasted_mash = pasted_mash.rsplit('.msh')[0]
    mash_cmd = ['mash', 'paste', pasted_mash]
    mash_cmd.extend(sketches)
    process = subprocess.Popen(mash_cmd)
    process.wait()
    return


def print_clusters(ANI, info = None):
    """
    Print line-separated representation of clusters. Will attempt to guess header
    attributes from fields provided in the info dictionary.
    Input:
        ANI <nx.Graph>      -- networkx graph of pairwise connections
        info <dict[dict]>   -- dictionary of dictionaries to maintain genome information
    """
    header = ['#cluster', 'num. genomes', 'genome']
    if info:
        sample_item = info.itervalues().next()
        fields = sample_item.keys()
        header.extend(fields)
    else:
        fields = None
    header.append('list')
    print '\t'.join(header)
    for cluster_num, cluster in enumerate(nx.connected_components(ANI)):
        cluster = list(cluster)
        size = len(cluster)
        for genome in cluster:
            attrs = []
            if info and genome in info:
                genome_info = info[genome]
                for field in fields:
                    if field in genome_info:
                        attrs.append(str(genome_info[field]))
                    else:
                        attrs.append('n/a')
            output_elements = [cluster_num, size, genome]
            output_elements.extend(attrs)
            output_elements.append(cluster)
            output_elements = map(str, output_elements)
            print '\t'.join(output_elements)

def get_genome_lengths(fastas):
    """
    Get genome lengths of all fasta files
    TODO: multi-thread
    """
    info = {}
    for genome in fastas:
        try:
            size = quick_len(genome)
        except RuntimeError:
            sys.exit("Genome length retrieval failed on {0}".format(genome))
        else:
            info[genome] = {'size':size}
    return info 

def quick_len(path):
    """
    Return <int> length of genome accessed at <str> path
    """
    len_cmd = ['tot-bp', path]
    proc = subprocess.Popen(len_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = proc.communicate()
    if stderr == '':
        return int(stdout)
    raise RuntimeError

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description = \
            '# cluster genomes based on average nucleotide identity (ani)')
    parser.add_argument(\
            '-f', nargs = '*', action = 'store', required = True, \
                help = 'fastas')
    parser.add_argument(\
            '-s', default = 98, type = float, required = False, \
                help = 'percent similarity (default = 98)')
    parser.add_argument(\
            '-m', action = 'store', required = True, type = str, \
            help = 'mash file to compare each fasta file against (will be created if doesn\'t exist)')
    parser.add_argument(\
            '-t', required = False, default = 6, type = int, \
            help = 'threads (default = 6)')
    parser.add_argument(\
            '--force', action = 'store_true', default = False, \
            help = 'force creation of mash sketches')
    parser.add_argument(\
            '--kmer', default = 21, type = int, required = False, \
            help = 'kmer used for mash sketch generation')
    args = parser.parse_args()
    genome_lengths = get_genome_lengths(args.f)
    make_mashes(args.f, args.m, threads = args.t, force = args.force, kmer = args.kmer)
    ANI = ani(args.f, args.s, args.m, threads = args.t)
    print_clusters(ANI, info = genome_lengths)
