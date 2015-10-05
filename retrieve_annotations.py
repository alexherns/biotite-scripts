#!/usr/bin/env python2.7
import kyotocabinet as kc
import scipy as sp
import numpy as np
import pandas as pd
import warnings, argparse, sys

warnings.filterwarnings('ignore')
from Bio import SearchIO

scaf2bin_default= '/Users/alexh/Documents/berkeley/banfield_lab/scaffolds2bins/horonobe_nr.scaf2bin.kch'
ggkb_default= '/Users/alexh/Documents/berkeley/banfield_lab/proteins/horonobe.proteins.kch'
metab_default= '/Users/alexh/Documents/berkeley/banfield_lab/hmm/processed_table-names.csv'
default_print_options= 'o'

parser = argparse.ArgumentParser(description='Helps with retrieva of protein annotations from databases', 
        formatter_class=argparse.ArgumentDefaultsHelpFormatter, add_help=False)

#Optional arguments
optional = parser.add_argument_group('OPTIONAL')
optional.add_argument('-h', action="help", help="show this help message and exit")
optional.add_argument('-g', help= 'gene name', type=str)
optional.add_argument('-s', help= 'scaffolds-to-bins kyotocabinet', type=str, default=scaf2bin_default)
optional.add_argument('-a', help= 'annotations kyotocabinet', type=str, default=ggkb_default)
optional.add_argument('-m', help= 'matchup table for genomes and genes', type=str, default=metab_default)
optional.add_argument('-b', help= 'read genes as batch from STDIN', action='store_true')
optional.add_argument('--print_options', help= 'list of printing options', type=str, default=default_print_options)

args = parser.parse_args()

scaf2bin_db= kc.DB()
scaf2bin_db.open(args.s, kc.DB.OREADER)
ggkb_genes= kc.DB()
ggkb_genes.open(args.a, kc.DB.OREADER)

metab_df= pd.read_csv(args.m, index_col=0)
def gene2scaf(gene):
    return '_'.join(gene.split('_')[:-1])

pos, neg= 0, 0
if args.g:
    for genes in metab_df[args.g][metab_df[args.g].notnull()]: 
        genes= genes.strip().split(' ')
        if genes != ['']:
            for gene in genes:
                line= ''
                if 'o' in args.print_options:
                    line+= gene+'\t'
                if 'g' in args.print_options:
                    line+= scaf2bin_db.get(gene2scaf(gene))+'\t'
                if 'a' in args.print_options:
                    line+= str(ggkb_genes.get(gene))+'\t'
                line= line.strip()
                sys.stdout.write(line+'\n')

elif args.b:
    for gene in sys.stdin.readlines():
        gene= gene.strip()
        line= ''
        if 'o' in args.print_options:
            line+= gene+'\t'
        if 'g' in args.print_options:
            line+= scaf2bin_db.get(gene2scaf(gene))+'\t'
        if 'a' in args.print_options:
            line+= str(ggkb_genes.get(gene))+'\t'
        line= line.strip()
        sys.stdout.write(line+'\n')
