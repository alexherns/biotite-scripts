#!/usr/bin/env python2.7
import sys, os, getopt

opts, args = getopt.getopt(sys.argv[1:], 'i:o:m:d:n:h', ['help'])

def usage():
        print """
Runs Non-linear Dimensionality Reduction of lrn files using the command-line interface "tapkee".

Usage: tapkee.py [OPTIONS]

OPTIONS:
		-i, --input	<string>	ESOM lrn file
		-o, --output	<string>	Output lrn file
		-m, --method	<string>	Dimensional reduction method
		-d, --delimiter	<string>	Delimiter (Default: \\t)
		-n, --dimensions	<integer>	Number of dimensions to output (Default: 2)
"""
        exit()

input_lrn= ''
output_lrn= ''
method= ''
delim= '\t'
dims= 2

for o, a in opts:
    if o in ('-h', '--help'):
        usage()
    elif o in ('-i', '--input'):
    	input_lrn= a
    elif o in ('-o', '--output'):
    	output_lrn= a
    elif o in ('-m', '--method'):
    	method= a
    elif o in ('-d', '--delimiter'):
    	delim= a
    elif o in ('-n', '--dimensions'):
    	dims= int(a)
        
if '' in [input_lrn, output_lrn, method]:
        usage()

#Modify input lrn file for usage with tapkee
tmp_handle= open(input_lrn+'_tmp_file', 'w')
for line in open(input_lrn):
	if "%" in line:
		continue
	tmp_handle.write(delim.join(line.strip().split(delim)[1:])+'\n')
tmp_handle.close()

#Run tapkee on input lrn
command= 'tapkee -i {0} -o {1} -m {2} -td {3} -d "{4}" --verbose'.format(input_lrn+'_tmp_file', output_lrn+'_tmp_file', method, dims, delim)
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
