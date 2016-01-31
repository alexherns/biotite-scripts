import numpy as np
import kyotocabinet as kc
import biotite_config

def mutate_sequence(seq, percent_id=0.9):
    """Mutates the StringBuffer sequence seq, retaining input percent id with
    the original"""
    import random
    from copy import deepcopy
    aminos= list('ACDEFGHIKLMNPQRSTVWY')
    seq= deepcopy(seq)
    num_mutations= int((1-percent_id)*len(seq))
    for _ in xrange(num_mutations):
        mut_position= random.randint(0, len(seq)-1)
        new_amino= aminos[random.randint(0, len(aminos)-1)]
        seq[mut_position]= new_amino
    return seq

def scaf2bin_dictionary(file_handle, header=True, sep="\t"):
    """"Creates a scaffold2bin dictionary object from the respective file"""
    if header==True:
        file_handle.readline()
    scaf2bin= {}
    for line in file_handle:
        scaffold, Bin= line.split(sep)[:2]
        scaf2bin[scaffold]= Bin
    return scaf2bin

def lrn2np(file_handle):
    """Reads an esom.lrn file into a numpy array"""
    return np.asarray([[float(i) for i in line.strip().split()[1:]] for line in file_handle if "%" not in line])
	
def writelrn(X, file_handle):
    """Writes an esom.lrn file from a numpy array given an output file handle"""
    num_rows, num_cols= np.shape(X)
    print num_rows
    print num_cols
    file_handle.write('%\t{0}\n%\t{1}\n% key \t{2}\n'.format(str(num_rows), str(num_cols), '\t'.join(["D"+str(d) for d in range(num_cols)])))
    for i, line in enumerate(X.tolist()):
        file_handle.write('{0}\t{1}\n'.format(str(i+1), '\t'.join([str(j) for j in line])))

def motif2regex(motif):
    """Converts a protein motif query to a regex search"""
    return motif.replace('X', '.').replace('x', '.').replace('(', '{').replace(')', '}')

class UClustSeq():
    def __init__(self, line):
        self.cluster= int(line[1])
        self.query= line[8]

def uclust2dict(file_handle):
    from collections import defaultdict
    d= defaultdict(list)
    for line in file_handle:
        line= line.strip().split('\t')
        d[line[1]].append(line[8])
    return d

def retrieve_scaf2bin_hash(hash_path= biotite_config.FULLSCAF2BINPATH):
    """Connect to kyotocabinet represenation of scaffolds-to-bins hash"""
    db= kc.DB()
    if not db.open(hash_path, kc.DB.OREADER):
        return db.error()
    return db

def get_best_bin(other_bin, bin_to_bin_mapping=None):
    """Returns name of best bin matching to other bin"""
    if not bin_to_bin_mapping:
        bin_to_bin_mapping= get_bin_to_bin_mapping()
    return bin_to_bin_mapping[other_bin]

def get_bin_to_bin_mapping(file_path= biotite_config.BIN_TO_BIN_CSV):
    """Returns defaultdict representation of bin_to_bin_mapping for best bin
    choices"""
    from collections import defaultdict
    bin_to_bin_mapping= defaultdict(list)
    for line in open(file_path, 'rU'):
        line= line.strip().split(',')
        best_bin= line[0]
        other_bins= [other_bin for other_bin in line if other_bin != '']
        for other_bin in other_bins:
            bin_to_bin_mapping[other_bin]= best_bin
    return bin_to_bin_mapping

def get_publication_name(best_bin, bin_name_mapping=None):
    """Returns publication name of bin"""
    if not bin_name_mapping:
        bin_name_mapping= get_publication_name_mapping()
    return bin_name_mapping[best_bin]

def get_publication_name_mapping(file_path= biotite_config.PUBLICATION_NAMES):
    """Returns dict represenation of publication names mapped to best bins"""
    bin_name_mapping= {}
    for line in open(file_path, 'rU'):
        line= line.strip().split(',')
        bin_name_mapping[line[0]]= line[1]
    return bin_name_mapping

def lookup_scaf2bin(scaffold, scaf2bin_hash=None):
    """Retrieve bin name of scaffold. Will use kc.db object if provided, otherwise creates a new connection to the kc"""
    if not scaf2bin_hash:
        scaf2bin_hash= retrieve_scaf2bin_hash()
    if isinstance(scaffold, str):
        return db.get(scaffold)
    elif isinstance(scaffold, list):
        return db.getbulk(scaffold)

def seq_list(file_path, file_type='fasta'):
    """Returns a list of SeqRecord objects from file_path using Bio.SeqIO.parse"""
    from Bio.SeqIO import parse
    return list(parse(open(file_path, 'rU'), file_type))

def scaffold_from_gene(gene):
    """Returns the scaffold from a gene (ggkbase style!)"""
    return '_'.join(gene.split('_')[:-1])

def rename_scaffold_to_bin(scaffold, **kwargs):
    """Renames scaffold to bin name, otherwise returns scaffold
    
    Keyword Args:
        scaf2bin_hash (kc.db) -- kyotocabint implementation of a scaffolds-to-bins
            mapping
        header_type (str) -- choose whether headers are genes or scaffolds
            [Default: scaffolds]
        desired_name (str) -- choose whether output bins are [Default: default, 
            best_bin, or publication] in format
        header_type='scaffold', desired_name='default'):
        bin_name_mapping= 
        bin_to_bin_mapping= 
    """
    if 'scaf2bin_hash' not in kwargs:
        kwargs['scaf2bin_hash']= retrieve_scaf2bin_hash()
    scaf2bin_hash= kwargs['scaf2bin_hash']
    genome_bin= scaf2bin_hash.get(scaffold)
    if not genome_bin:
        return scaffold
    elif 'UNK' in genome_bin.upper():
        return scaffold + '_unbinned'
    elif 'desired_name' not in kwargs or kwargs['desired_name'] == 'default':
        pass
    else:
        if 'bin_to_bin_mapping' not in kwargs:
            kwargs['bin_to_bin_mapping']= get_bin_to_bin_mapping()
        bin_to_bin_mapping= kwargs['bin_to_bin_mapping']
        if not bin_to_bin_mapping[genome_bin]:
            return genome_bin
        genome_bin= bin_to_bin_mapping[genome_bin]
        if desired_name == 'publication':
            if 'bin_name_mapping' not in kwargs:
                bin_name_mapping= get_publication_name_mapping()
            else:
                bin_name_mapping= kwargs['bin_name_mapping']
            if genome_bin not in bin_name_mapping:
                return genome_bin
            genome_bin= bin_name_mapping[genome_bin]
    return genome_bin

def rename_fasta_by_bin(multifasta, **kwargs):
    """Renames header in each SeqRecord object in the multifasta list according
    to bin
    
    Keyword Arguments:
        scaf2bin_hash (kc.db) -- kyotocabint implementation of a scaffolds-to-bins
            mapping
        header_type (str) -- choose whether headers are genes or scaffolds
            [Default: scaffolds]
    """
    if 'scaf2bin_hash' not in kwargs:
        kwargs['scaf2bin_hash']= retrieve_scaf2bin_hash()
    for seq_record in multifasta:
        if header_type == 'gene':
            seq_record.id= scaffold_from_gene(seq_record.id)
        seq_record.id= rename_to_bin(seq_record.id, **kwargs)
    return multifasta

def rename_tree_nodes_by_bin(tree, **kwargs):
    """Renames header in each Node in Tree according to bin
    
    Keyword Arguments:
        scaf2bin_hash (kc.db) -- kyotocabint implementation of a scaffolds-to-bins
            mapping
        header_type (str) -- choose whether headers are genes or scaffolds
            [Default: scaffolds]
        desired_name (str) -- choose whether output bins are [Default: default, 
            best_bin, or publication] in format
    header_type='scaffold', desired_name='default'):
    """
    if 'scaf2bin_hash' not in kwargs:
        kwargs['scaf2bin_hash']= retrieve_scaf2bin_hash()
    bin_to_bin_mapping= get_bin_to_bin_mapping()
    bin_name_mapping= get_publication_name_mapping()
    for node in tree.leaf_node_iter():
        if 'header_type' in kwargs and kwargs['header_type'] == 'genes':
            node.taxon.label= scaffold_from_gene(node.taxon.label)
        node.taxon.label= rename_scaffold_to_bin(node.taxon.label, bin_name_mapping=bin_name_mapping,
                bin_to_bin_mapping= bin_to_bin_mapping, **kwargs)
        print node.taxon.label
    return tree
