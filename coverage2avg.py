#!/usr/bin/env python
'''
Reports simple metrics from ouput of bedtools genomecov.
Prints the average coverage of the genome, and percent of the genome covered by 
each coverage level.

Note: requires use of "bedtools genomecov" using the -d option
Usage: coverage2avg.py <coverage_file>
'''

import sys

coverage_list = [float(line.strip().split()[2]) for line in open(sys.argv[1], 'r')]
cov_dict = {}
for val in coverage_list:
	if val not in cov_dict:
		cov_dict[val] = 0
	cov_dict[val]+=1

for val in cov_dict:
	print '{0}:\t{1}%'.format(val, float(cov_dict[val])/len(coverage_list)*100)
print 'Average coverage:\t'.format(sum(coverage_list)/len(coverage_list))
