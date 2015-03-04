#!/usr/bin/env python 
"""
Appends a breseq-input file with annotations of whether predicted intergenic mutations lie upstream of coding regions
"""

import sys, getopt

def usage():
        print '''

Usage: breseq_promoters.py <options>

	-h, --help	Print this help dialog

	<REQUIRED ARGUMENTS>
	-b	Breseq input file
	-g	Tabbed-CDS locations from genbank (see gbk_features.py for details)

	<OPTIONAL ARGUMENTS>
	-s	Specify length upstream of start site to start search [default: 100]
	-e	Specify length upstream of start site to end search (required s>e) [default: 1]
'''
        exit()

opts, args = getopt.getopt(sys.argv[1:], 'he:s:b:g:', ['help'])

pro_start = 100
pro_end = 1

for o, a in opts:
        if o in ('-h', '--help'):
                usage()
	elif o == '-s':
		pro_start = int(a)
	elif o == '-e':
		pro_end = int(a)
	elif o == '-b':
		breseq_file = a
	elif o == '-g':
		genbank_file = a

if len([o for o, a, in opts if o == '-b' or o == '-g'])<2:
        print '''
This script requires more arguments. Please pass -h or --help for assistance
'''
        exit()

if pro_end > pro_start:
	print '\n\nUpstream start site smaller than end site, please see help for assistance\n\n'
	exit()

promoter_dict = {}
for line in open(genbank_file, 'r'):
	if len(line.strip().split('\t'))!=6:
		print '\n\nGENBANK TAB FILE HAS INCORRECT NUMBER OF COLUMNS\n\n'
		exit()
	contig, contig_length, f_start, f_end, strand, product = line.strip().split('\t')
	contig_length, f_start, f_end = int(contig_length), int(f_start), int(f_end)
	if contig not in promoter_dict:
		promoter_dict[contig] = []
	if strand == '1':
		#CDS on forward strand
		this_p_start, this_p_end = max(1,f_start-pro_start+1), max(1,f_start-pro_end+1)
	if strand == '-1':
		#CDS on reverse strand
		this_p_start, this_p_end = min(contig_length, f_end+pro_end), min(contig_length, f_end+pro_start)
	promoter_dict[contig].append((this_p_start, this_p_end, strand, product))

for line in open(breseq_file, 'r'):

	if line.strip().split('\t')[4] == '':
	#Mutation is not of a funky formatted MC type
		print line.strip()
		continue
	if line.strip().split('\t')[4].split()[0] != 'intergenic':
		print line.strip()
		continue

	#Mutation is in an intergenic region
	
	contig = line.strip().split('\t')[1]
	start_loc = int(line.strip().split('\t')[2].replace(',', ''))
	mutation = line.strip().split('\t')[3]
	
	if contig not in promoter_dict:
		#Mutation is on a contig without any coding regions
		print line.strip()
		continue

	if mutation[0] == '+':
		#Mutation is of "+N base pairs" type
		stop_loc = start_loc
	elif mutation[:2] == '\xce\x94':
		#Mutation is of "delta base pairs" type
		stop_loc = int(mutation.split()[0][2:].replace(',', '')) + start_loc
	elif len(mutation) == 5:# and mutation[1:4] == '\xe2\x86\x92':
		#Mutation is of simple "N arrow N" type
		stop_loc = start_loc
	elif len(mutation) >= 9 and '\xe2\x86\x92' in mutation:
		stop_loc = int(mutation.split()[0].replace(',','')) + start_loc
	else:
		print '\n\n\nFAILURE TO IDENTIFY STOP LOCATION!\n\n\n'
		exit()

	for locus in [start_loc, stop_loc]:
		#print 'This is my locus:', locus
		success = 0
		if locus == start_loc:
			which_loc = 'starting_edge'
		else:
			which_loc = 'trailing_edge'
		for this_p_start, this_p_end, strand, product in promoter_dict[contig]:
			if locus >= this_p_start and locus <= this_p_end:
				if strand == '1':
					print '{0}\t-{1}\t{2}\t{3}'.format(line.strip(), this_p_end+pro_end-locus, which_loc, product)
				else:
					print '{0}\t-{1}\t{2}\t{3}'.format(line.strip(), locus-this_p_start+pro_end, which_loc, product)
				success = 1
				break

		if success == 1:
			break
		if start_loc == stop_loc:
			break
	if success == 1:
		continue
	else:
		print line.strip()
