#!/usr/bin/env python2.7
import sys, colorsys, getopt
import numpy as np

def usage():
        print '''
Automatically generates ESOM class file (complete with color-coded bin assignments)
from ESOM names file and tab-delimited scaffold-bin assignments

Usage: class_wrapper -n esom.names -s scaffolds2bins.tsv

        -h, --help      Print this help dialog

        <REQUIRED ARGUMENTS>
        -n      ESOM Names file
        -s      List of tab-delimited scaffold-bin assignments

'''
        exit()

opts, args = getopt.getopt(sys.argv[1:], 'hn:s:', ['help'])
scaffold_file= ""
names_file= ""

for o, a in opts:
        if o in ('-h', '--help'):
                usage()
        elif o == '-s':
                scaffold_file = a
        elif o == '-n':
                names_file= a

if len([o for o, a, in opts])<2:
        print '''
This script requires more arguments. Please pass -h or --help for assistance
'''
        exit()

#Function to create an RGB colorspace from a defined number of colors
def get_colors(num_colors):
    colors=[]
    for i in np.arange(0., 360., 360. / num_colors):
                hue = i/360.
                lightness = (50 + np.random.rand() * 10)/100.
                saturation = (90 + np.random.rand() * 10)/100.
                temp_col= colorsys.hls_to_rgb(hue, lightness, saturation)
                colors.append([int(col*255) for col in temp_col])
    return colors

"""Dictionary of bins accessed by scaffold
    scaf_dict= {scaffold: bin}"""
scaf_list= [line.strip().split() for line in open(scaffold_file) if "scaffold_name" not in line]
scaf_dict= {}
for scaf in scaf_list:
    scaf_dict[scaf[0]]= scaf[1]

"""Define list of bins and color space"""
bins= list(set([scaf[1] for scaf in scaf_list]))
color_space= get_colors(len(bins))

"""Dictionary of bin information
    bin_dict= {bin: {"rgb": [r, g, b], "index": bin index}}"""
bin_dict= {}
for org, color in zip(bins, color_space):
    bin_dict[org]= {}
    bin_dict[org]["rgb"]= color
    bin_dict[org]["index"]= bins.index(org)+1

#Read in the names file and generate a header for the class file
file_handle= open(names_file)
print file_handle.readline().strip()
for org in bin_dict:
    print "%\t{0}\t{1}\t{2}".format(str(bin_dict[org]["index"]), org, "\t".join([str(col) for col in bin_dict[org]["rgb"]]))


"""Flesh out the class file with all assignments"""
for line in file_handle:
    i, scaf= line.strip().split()[:2]
    scaf= "_".join(scaf.split("_")[:-1])
    print i, bin_dict[scaf_dict[scaf]]["index"]

file_handle.close()
