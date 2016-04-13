#!/usr/bin/env python
import os
import sys
import subprocess
import cPickle

def run_uniquify_fasta(fasta_file, output_file, pickle_file):
    """
    Returns a subprocess from command-line script "uniquify_fasta.py" <subprocess.Popen>
    Input:
        fasta_file <str>    -- path to fasta file
        output_file <str>   -- path to output (simplified) fasta file
        pickle_file <str>   -- path to pickle file containing name-mapping of new->old
    """ 
    uniquify_cmd = gen_uniquify_fasta_cmd(fasta_file, output_file, pickle_file = pickle_file)
    return subprocess.Popen(uniquify_cmd)

def gen_uniquify_fasta_cmd(fasta_file, output_file, pickle_file = None):
    """
    Returns a "uniquify_fasta.py" command <list>
    Input:
        fasta_file <str>    -- path to fasta file
        output_file <str>   -- path to output (simplified) fasta file
        pickle_file <str>   -- path to pickle file containing name-mapping of new->old
    """ 
    if pickle_file == None:
        pickle_file = os.path.splitext(fasta_file)[0]
        pickle_file += '.pkl'
    return ['uniquify_fasta.py', '-f', fasta_file, '-o', output_file, '-p', pickle_file]

def run_fasta2phy(fasta_file, phylip_file):
    """
    Returns a subprocess from command-line script "fasta2phy" <subprocess.Popen>
    Input:
        fasta_file <str>    -- path to fasta file
        phylip_file <str>   -- path to output phylip file
    """
    fasta2phy_cmd = gen_fasta2phy_cmd(fasta_file, phylip_file)
    return subprocess.Popen(fasta2phy_cmd)

def gen_fasta2phy_cmd(fasta_file, phylip_file):
    """
    Returns a "fasta2phy" command <list>
    Input:
        fasta_file <str>    -- path to fasta file
        phylip_file <str>   -- path to output phylip file
    """ 
    return ['fasta2phy', '-i', fasta_file, '-o', phylip_file]

def run_seqboot(phylip_file, output_file, args):
    """
    Returns a subprocess from Joe Felsenstein's "seqboot" utility <subprocess.Popen>
    Input:
        phylip_file <str>   -- path to phylip alignment
        output_file <str>   -- path to bootstrapped alignments
        args <Namespace>    -- keyword arguments for utility
    TO-DO:
        -> Refactor the args parameter to make use of *args
    """ 
    if os.path.isfile('outfile') is False:
        subprocess.Popen(['touch', 'outfile']).wait()
    if os.path.isfile(output_file) is True:
        subprocess.Popen(['rm', output_file]).wait()
    seq_boot_str = '\n'.join(gen_seqboot_args(phylip_file, output_file, args))
    print_cmd = ['printf', seq_boot_str]
    print_process = subprocess.Popen(print_cmd, stdout=subprocess.PIPE)
    if args.verbose:
        stderr = sys.stderr
    else:
        stderr = open(args.log, 'a')
    seq_boot_process = subprocess.Popen(['seqboot'], stdin=print_process.stdout, stderr=stderr, stdout=stderr)
    return seq_boot_process

def gen_seqboot_args(phylip_file, output_file, args):
    """
    Returns a command list for Joe Felsenstein's "seqboot" utility <list>
    Input:
        phylip_file <str>   -- path to phylip alignment
        output_file <str>   -- path to bootstrapped alignments
        args <Namespace>    -- keyword arguments for utility
    TO-DO:
        -> Refactor the args parameter to make use of *args
    """ 
    params = {}
    params['J'] = ('Jackknife' if args.jackknife else 'Bootstrap')
    params['R'] = str(args.n)
    seq_boot_args = [phylip_file]
    for key, val in params.iteritems():
        seq_boot_args.extend([str(key), str(val)])
    seq_boot_args.extend(['Y', str(args.s), 'F', output_file])
    return seq_boot_args

def run_FastTree(args, aln_file, tree_file, bootstrap = True):
    """
    Returns a subprocess of the FastTree program <subprocess.Popen>
    Input:
        args <Namespace>    -- keyword arguments for utility
        aln_file <str>      -- path to alignment (FASTA or PHYLIP format)
        tree_file <str>     -- path to output tree (Newick format)
        bootstrap <boolean> -- tree should be run on a bootstrapped alignment
    TO-DO:
        -> Refactor the args parameter to make use of *args
    """ 
    fasttree_cmd = gen_fasttree_cmd(args, aln_file, tree_file)
    if args.verbose:
        stderr = sys.stderr
    else:
        stderr = open(args.log, 'a')
    return subprocess.Popen(fasttree_cmd, stderr=stderr)

def gen_fasttree_cmd(args, aln_file, tree_file, bootstrap = True):
    """
    Returns a command list for the FastTree program <list>
    Input:
        args <Namespace>    -- keyword arguments for utility
        aln_file <str>      -- path to alignment (FASTA or PHYLIP format)
        tree_file <str>     -- path to output tree (Newick format)
        bootstrap <boolean> -- tree should be run on a bootstrapped alignment
    TO-DO:
        -> Refactor the args parameter to make use of *args
    """ 
    params = {}
    params['-cat'] = str(args.cat)
    if bootstrap:
        params['-n'] = str(args.n)
    if args.gamma: params['-gamma'] = ''
    if args.wag: params['-wag'] = ''
    if args.nt: params['-nt'] = ''
    if not args.verbose: 
        params['-quiet'] = ''
        params['-nopr'] = ''
        params['-log'] = args.log
    fastree_args = [('FastTreeMP' if args.mp else 'FastTree')]
    for key, val in params.iteritems():
        if val:
            fastree_args.extend([key, val])
        else:
            fastree_args.append(key)
    fastree_args.extend(['-out', tree_file, aln_file])
    return fastree_args

def run_compare_bootstraps(single_tree, bootstrapped_trees, output_file):
    """
    Returns a subprocess of Morgan Price's "CompareToBootstrap" perl script <subprocess.Popen>
    Input:
        single_tree <str>           -- path to reference tree file (Newick)
        bootstrapped_trees <str>    -- path to bootstrapped trees file (Newick)
        output_file <str>           -- path to write tree annotated with bootstrap support
    """ 
    if not check_perl_module('MMOTree'):
        raise OSError, "Check to make sure your PERL5LIB path includes the MOTree.pm module"
    output_handle = open(output_file, 'wb')
    compare_cmd = gen_compare_cmd(single_tree, bootstrapped_trees)
    if args.verbose:
        stderr = sys.stderr
    else:
        stderr = open(args.log, 'a')
    return subprocess.Popen(compare_cmd, stdout=output_handle, stderr=stderr)

def gen_compare_cmd(single_tree, bootstrapped_trees):
    """
    Returns a command list for Morgan Price's "CompareToBootstrap" perl script <list>
    Input:
        single_tree <str>           -- path to reference tree file (Newick)
        bootstrapped_trees <str>    -- path to bootstrapped trees file (Newick)
    """ 
    cmp_prog_path = '/home/alexh/bin/MOTreeComparison/CompareToBootstrap.pl'
    compare_cmd = ['perl', cmp_prog_path, '-tree', single_tree, '-boot', bootstrapped_trees]
    return compare_cmd

def check_perl_module(perl_module):
    """
    Returns perl_module is accessible in PERL5LIB
    Input:
        perl_module <str>   -- name of perl module
    """
    cmd = ['perl', '-M{0}'.format(perl_module), '-e1']
    try:
        exit_status = subprocess.Popen(cmd).wait()
        return int(exit_status) == 0
    except OSError:
        sys.exit("Checking for perl module {0} failed!".format(perl_module))

def relabel_nodes(tree_file, pickle_file, output_file):
    """
    Relabel nodes on tree according to mapping saved in pickle file
    Input:
        tree_file <str>     -- path to input tree file (Newick)
        pickle_file <str>   -- path to pickle containing new->old mapping
        output_file <str>   -- path to write renamed tree (Newick)
    """
    import ete3
    pickle_handle = open(pickle_file, 'rb')
    name_map = cPickle.load(pickle_handle)
    pickle_handle.close()
    tree = ete3.Tree(tree_file)
    for node in tree.iter_leaves():
        if node.name in name_map:
            node.name = name_map[node.name]
        else:
            sys.stderr.write("{0} from tree ({1}) not found in pickle ({2})".format(node.name, tree_file, pickle_file))
    tree.write(outfile = output_file)

def cleanup(files):
    """
    Deletes files from input
    Input:
        files <list[str]>   -- paths to files that should be removed
    """
    for f_path in files:
        if os.path.isfile(f_path):
            try:
                os.remove(f_path)
            except OSError:
                continue

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description = \
            'generate bootstrapped trees using FastTree')
    parser.add_argument(\
            '-f', type = str, required = True, \
                help = 'aligned fasta (required)')
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
    parser.add_argument(\
            '--log', type = str, default = 'log.txt', \
            help = 'log file for stderr logging when not run with verbose')
    parser.add_argument(\
            '--clean', action = 'store_true', \
            help = 'clean up temporary files after generation of final tree')
    args = parser.parse_args()
    if args.wag and any([args.nt, args.gtr]):
        sys.exit('WAG model incompatible with nt alignments')
    if args.gtr and not args.nt:
        sys.exit('GTR model incompatible with aa alignments')
    if args.s % 2 == 0:
        sys.exit('Seed must be odd for seqboot')
    base_name, ext = os.path.splitext(args.f)
    simple_name = base_name + '.simple' + ext
    pickle_name = base_name + '.pkl'
    phylip_name = base_name + '.phylip'
    bootstrap_name = base_name + '.boot.phylip'
    tree_name = base_name + '.tree'
    boottree_name = base_name + '.boots.tree'
    compared_name = base_name + '.bootvals.tree'
    relabeled_name = base_name + '.final.tree'
    run_uniquify_fasta(args.f, simple_name, pickle_name).wait()
    run_fasta2phy(simple_name, phylip_name).wait()
    run_seqboot(phylip_name, bootstrap_name, args).wait()
    run_FastTree(args, simple_name, tree_name, bootstrap = False).wait()
    run_FastTree(args, bootstrap_name, boottree_name, bootstrap = True).wait()
    run_compare_bootstraps(tree_name, boottree_name, compared_name).wait()
    relabel_nodes(compared_name, pickle_name, relabeled_name)
    if args.clean:
        cleanup([simple_name, pickle_name, phylip_name, bootstrap_name, tree_name, boottree_name, compared_name])
