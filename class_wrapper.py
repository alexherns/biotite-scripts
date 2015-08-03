#!/usr/bin/env python2.7
import colorsys, argparse
import numpy as np

parser = argparse.ArgumentParser(description='generates ESOM class file (with color-coded bins) from ESOM names file and scaffold-bin assignments.', formatter_class=argparse.ArgumentDefaultsHelpFormatter, add_help=False)

#Required arguments
required = parser.add_argument_group('REQUIRED')
required.add_argument('-n', help= 'esom.names', required=True, type=argparse.FileType('r'))
required.add_argument('-s', help= 'scaffold-bin assignments', required=True, type=argparse.FileType('r'))

#Optional arguments
optional = parser.add_argument_group('OPTIONAL')
optional.add_argument('-h', action="help", help="show this help message and exit")
optional.add_argument('--header', help= 'Scaffold-bin file contains header', action='store_true')

args = parser.parse_args()

scaffold_file= args.s
names_file= args.n
header= args.header

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
file_handle= scaffold_file
if not header:
	file_handle.readline()
scaf_list= [line.strip().split() for line in file_handle]
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
file_handle= names_file
print file_handle.readline().strip()
for org in bin_dict:
    print "%\t{0}\t{1}\t{2}".format(str(bin_dict[org]["index"]), org, "\t".join([str(col) for col in bin_dict[org]["rgb"]]))


"""Flesh out the class file with all assignments"""
for line in file_handle:
    i, scaf= line.strip().split()[:2]
    scaf= "_".join(scaf.split("_")[:-1])
    if scaf in scaf_dict:
    	print i, bin_dict[scaf_dict[scaf]]["index"]

file_handle.close()
