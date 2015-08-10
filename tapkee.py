#!/usr/bin/env python2.7
import sys, os, argparse

parser = argparse.ArgumentParser(
	description='''Runs Non-linear Dimensionality Reduction of lrn files using 
		the command-line interface "tapkee"''', 
	formatter_class=argparse.ArgumentDefaultsHelpFormatter, add_help=False)

#Required arguments
required = parser.add_argument_group('REQUIRED')
required.add_argument('-i', help= 'input.lrn', required=True, type=str)
required.add_argument('-o', help= 'output.lrn', required=True, type=str)

#Optional arguments
optional = parser.add_argument_group('OPTIONAL')
optional.add_argument('-h', action="help", help="show this help message and exit")
optional.add_argument('-m', help= 'input.fasta', type=str,
	choices= ['lle', 'npe', 'ltsa', 'lltsa', 'hlle', 'lpp', 'dm', 'isomap',
		'l-isomap', 'mds', 'l-mds', 'spe', 'kpca', 'pca', 'ra', 'fa', 't-sne'],
	default= 'pca')
optional.add_argument('-d', help= 'delimiter', type=str, default= '\t')
optional.add_argument('-n', help= 'output dimensions', type= int, default= 2)
optional.add_argument('-x', help= 'additional tapkee arguments', type= str, default= '')

args = parser.parse_args()

input_lrn= args.i
output_lrn= args.o
method= args.m
delim= args.d
dims= args.n

#Modify input lrn file for usage with tapkee
tmp_handle= open(input_lrn+'_tmp_file', 'w')
for line in open(input_lrn):
	if "%" in line:
		continue
	tmp_handle.write(delim.join(line.strip().split(delim)[1:])+'\n')
tmp_handle.close()

#Run tapkee on input lrn
command= 'tapkee -i {0} -o {1} -m {2} -td {3} -d "{4}" {5} --verbose'.format(input_lrn+'_tmp_file', output_lrn+'_tmp_file', method, dims, delim, args.x)
print command
os.system(command)

#Reformat output
tmp_handle= open(output_lrn, 'w')
output_lines= [line for line in open(output_lrn+'_tmp_file')]
tmp_handle.write('% {0}\n% {1}\n% 9\t{2}\n% key\t{3}\n'.format(str(len(output_lines)), str(dims), '\t'.join(["1" for i in range(dims)]), '\t'.join(["D"+str(i) for i in range(dims)])))
for i, line in enumerate(output_lines):
	tmp_handle.write('{0}\t{1}'.format(str(i+1), line))
tmp_handle.close()

os.system('rm {0}; rm {1}'.format(input_lrn+'_tmp_file', output_lrn+'_tmp_file'))
