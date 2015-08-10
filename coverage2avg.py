#!/usr/bin/env python2.7
import argparse

parser = argparse.ArgumentParser(description='Prints the average coverage of the genome, and percent of the genome covered by each coverage level.', formatter_class=argparse.ArgumentDefaultsHelpFormatter, add_help=False)

#Required arguments
required = parser.add_argument_group('REQUIRED')
required.add_argument('-c', help= 'coverage file from "bedtools genomecov -d"', required=True, type=argparse.FileType('r'))

#Optional arguments
optional = parser.add_argument_group('OPTIONAL')
optional.add_argument('-h', action="help", help="show this help message and exit")

args = parser.parse_args()

coverage_list = [float(line.strip().split()[2]) for line in args.c]
cov_dict = {}
for val in coverage_list:
	if val not in cov_dict:
		cov_dict[val] = 0
	cov_dict[val]+=1

for val in cov_dict:
	print '{0}:\t{1}%'.format(val, float(cov_dict[val])/len(coverage_list)*100)
print 'Average coverage:\t'.format(sum(coverage_list)/len(coverage_list))
