import os
import subprocess

def run_fasta2phy(fasta_file, phylip_file):
    fasta2phy_cmd = gen_fasta2phy_cmd(fasta_file, phylip_file)
    return subprocess.Popen(fasta2phy_cmd)

def gen_fasta2phy_cmd(fasta_file, phylip_file):
    return ['fasta2phy', '-i', fasta_file, '-o', phylip_file]

def run_uniquify_fasta(fasta_file, output_name, pickle_name):
    uniquify_cmd = gen_uniquify_fasta_cmd(fasta_file, output_name, pickle_name = pickle_name)
    return subprocess.Popen(uniquify_cmd)

def gen_uniquify_fasta_cmd(fasta_file, output_name, pickle_name = None):
    if pickle_name == None:
        pickle_name = os.path.splitext(fasta_file)[0]
        pickle_name += '.pkl'
    return ['uniquify_fasta.py', '-f', fasta_file, '-o', output_name, '-p', pickle_name]

def run_seqboot(args, phylip_file, output_file):
    if os.path.isfile('outfile') is False:
        subprocess.Popen(['touch', 'outfile']).wait()
    if os.path.isfile(output_file) is True:
        subprocess.Popen(['rm', output_file]).wait()
    seq_boot_str = '\n'.join(gen_seqboot_str(args, phylip_file, output_file))
    print_cmd = ['printf', seq_boot_str]
    print_process = subprocess.Popen(print_cmd, stdout=subprocess.PIPE)
    seq_boot_process = subprocess.Popen(['seqboot'], stdin=print_process.stdout)
    return seq_boot_process

def gen_seqboot_args(args, phylip_file, output_file):
    params = {}
    params['J'] = ('Jackknife' if args.jackknife else 'Bootstrap')
    params['R'] = str(args.n)
    params['2'] = ('Yes' if args.verbose else 'No')
    seq_boot_args = [phylip_file]
    for key, val in params.iteritems():
        seq_boot_args.extend([key, val])
    seq_boot_args.extend(['Y', args.s, output_file])
    return seq_boot_args

def run_FastTree(args, aln_file, tree_file, bootstrap = True):
    fasttree_cmd = gen_fasttree_cmd(args, aln_file, tree_file)
    return subprocess.Popen(fasttree_cmd)

def gen_fasttree_cmd(args, aln_file, tree_file, bootstrap = True):
    params = {}
    params['-cat'] = args.cat
    if bootstrap:
        params['-n'] = str(args.n)
    if args.gamma: params['-gamma'] = ''
    if args.wag: params['-wag'] = ''
    if args.nt: params['-nt'] = ''
    if not args.verbose: 
        params['-quiet'] = ''
        params['-nopr'] = ''
    params['-out'] = tree_file
    fastree_args = [('FastTreeMP' if args.mp else 'FastTree')]
    for key, val in params.iteritems():
        fastree_args.extend([key, val])
    fastree_args.append(aln_file)
    return fastree_args

def run_compare_bootstraps(single_tree, bootstrapped_trees):
    compare_cmd = gen_compare_cmd(single_tree, bootstrapped_trees)
    return subprocess.Popen(compare_cmd)

def gen_compare_cmd(single_tree, bootstrapped_trees):
    cmp_prog_path = '/home/alexh/bin/MOTreeComparison/CompareToBootstrap.pl'
    compare_cmd = ['perl', cmp_prog_path, '-tree', single_tree, '-boot', bootstrapped_trees]
    return compare_cmd

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description = \
            'generate bootstrapped trees using FastTree')
    parser.add_argument(\
            '-f', type = str, required = True, \
                help = 'aligned fasta')
    parser.add_argument(\
            '-n', type = int, default = 100, \
            help = 'number of bootstrap replicates')
    parser.add_argument(\
            '-s', type = int, default = 12345, \
            help = 'seed for bootstrapping')
    parser.add_argument(\
            '--mp', action = 'store_true', \
            help = 'Use MP version of FastTree for parallel computation')
    parser.add_argument(\
            '--jackknife', action = 'store_true', \
            help = 'generate jackknife rather than bootstrap replicates')
    parser.add_argument(\
            '--cat', type = int, default = 20, \
            help = 'number of rate categories for sites')
    parser.add_argument(\
            '--gamma', action = 'store_true', \
            help = 'rescale branch lengths to optimize gamma likelihood')
    parser.add_argument(\
            '--wag', action = 'store_true', \
            help = 'Whelan-And-Goldman 2001 model instead of (default) Jones-Taylor-Thorton 1992 model (aa only)')
    parser.add_argument(\
            '--nt', action = 'store_true', \
            help = 'nucleotide alignment (default:False)')
    parser.add_argument(\
            '--gtr', action = 'store_true', \
            help = 'generalized time-reversible model (default:False) (nt only)')
    parser.add_argument(\
            '--verbose', action = 'store_true', \
            help = 'print progress metrics to stderr')
    args = parser.parse_args()
    if args.wag and any(args.nt, args.gtr):
        sys.exit('WAG model incompatible with nt alignments')
    if args.gtr and not args.nt:
        sys.exit('GTR model incompatible with aa alignments')










