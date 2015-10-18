#!/usr/bin/env python2.7
import argparse
import pandas as pd

parser = argparse.ArgumentParser(description='Extracts the taxonomy string from ggkbase organism_info files into a columnar format', formatter_class=argparse.ArgumentDefaultsHelpFormatter, add_help=False)

#Required arguments
required = parser.add_argument_group('REQUIRED')
required.add_argument('-i', help= 'input organism_info sheet', required=True, type=str)
required.add_argument('-o', help= 'output csv', required=True, type=str)

#Optional arguments
optional = parser.add_argument_group('OPTIONAL')
optional.add_argument('-h', action="help", help="show this help message and exit")

args = parser.parse_args()

input_file= args.i
output_file= args.o
df= pd.read_csv(input_file, sep='\t', index_col='name')

levels= ('domain', 'phylum', 'class', 'order', 'family', 'genus', 'species')
for level in levels:
    df[level]= ''

domains= set(['Bacteria', 'Phage', 'Archaea', 'Plasmid-like', 'Eukaryota'])
for org in df.index:
    tax_string= df.loc[org].taxonomy
    tax_list= tax_string.split(", ")[::-1]
    if tax_list[0] == 'Eukaryote':
        tax_list[0]= 'Eukaryota'
    if tax_list[0] not in domains:
        tax_list.insert(0, 'Bacteria')
    for i in range(len(tax_list)):
        df.set_value(org,levels[i], tax_list[i])

df.to_csv(output_file)
