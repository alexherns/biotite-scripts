#!/usr/bin/env python

############################################################################
# Copyright (c) Alex W Hernsdorf
# ahernsdorf@berkeley.edu
# All Rights Reserved
############################################################################

import sys
'''Compares the results of a VCF file with a SNPs type file from Mauve.

Required: vcf_snps.py <vcf_file> <snps_file>
Output: multiple numbers, including total number SNPs in VCF,
	number of SNPs from same location,
	and number of identical SNPs
Notes: Will only compare the single bp SNPs
'''

file1 = open(sys.argv[1], 'r')
file2 = open(sys.argv[2], 'r')

mut_snps = 0
mut_snp_list = []
for line in file1:
	if line[0] == '#':
		continue
	node, loc, dummy, site1, site2 = line.strip().split()[:5]
	if len(site1)==1 and len(site2) == 1:
		mut_snp_list.append((node, loc, site1, site2))
		mut_snps+=1

file2.readline()
dif_snp_list = []
dif_snp_locs = []
for line in file2:
	snp_set, node, loc = line.strip().split()[:3]
	site1 = snp_set[0].upper()
	site2 = snp_set[1].upper()
	dif_snp_list.append((node, loc, site1, site2))
	dif_snp_locs.append((node, loc))

mut_snp_set = set(mut_snp_list)
dif_snp_set = set(dif_snp_list)
dif_snp_locs_set = set(dif_snp_locs)
same_location = 0
reasonable_snp_count = 0
for snp in mut_snp_set:
	if snp[:2] in dif_snp_locs_set:
		same_location+=1
		if snp in dif_snp_set:
			reasonable_snp_count+=1
	else:
		print snp

print mut_snps, same_location, reasonable_snp_count
file1.close()
file2.close()
