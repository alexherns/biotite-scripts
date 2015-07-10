#!/usr/bin/env python2.7
import sys, getopt, os

opts, args = getopt.getopt(sys.argv[1:], 's:c:b:h', ['help'])

def usage():
        print """
Will create table to assist in linking 16S genes to bins by scaffold.

Usage: 16S_bins.py -s scaf2bin.tsv -c clusters.sorted.uc -b bin_table.csv 
"""
        exit()

scaf2bin_file= ''
cluster_file= ''
bin_file= ''

for o, a in opts:
    if o in ('-h', '--help'):
        usage()
    elif o in ('-s'):
        scaf2bin_file= a
    elif o in ('-c'):
        cluster_file= a
    elif o in ('-b'):
        bin_file= a

if '' in [scaf2bin_file, cluster_file, bin_file]:
        print """
Please specify appropriate parameters. See -h (--help) for help
"""
        exit()
        
from common_objects import *

#Create scaffold2bin dictionary
scaf2bin= scaf2bin_dictionary(open(scaf2bin_file), header=True, sep="\t")

#Create bin2bestbin dictionary
bestBins= {}
for line in open(bin_file):
	bins= [item for item in line.strip().split(",")[1:] if item != ""]
	best_bin= bins.pop(0)
	for Bin in bins:
		bestBins[Bin]= best_bin
	
cluster_dict= {}	
for line in open(cluster_file):
	i= line.strip().split('\t')[1]
	seqID= line.strip().split()[8]
	if "|" in seqID:
		seqID= seqID.split("|")[1]
	if i not in cluster_dict:
		cluster_dict[i]= []
	cluster_dict[i].append(seqID)


for cluster in cluster_dict:
	lineOut= []
	
	#Check if there were no scaffolds (and so must have been only from EMIRGE
	if len([seqID for seqID in cluster_dict[cluster] if "scaffold" in seqID and "ig2599" not in seqID and "Ig3401" not in seqID and "ig3402" not in seqID and "Ig3397" not in seqID and "Ig3399" not in seqID]) == 0:
		 lineOut = "\t".join([cluster, "", "", "", "", seqID])		 
	else:
		for seqID in cluster_dict[cluster]:
			scaffold= ""
			if "scaffold" in seqID and "ig2599" not in seqID and "Ig3401" not in seqID and "ig3402" not in seqID and "Ig3397" not in seqID and "Ig3399" not in seqID:
				scaffold= seqID.split("_16S")[0]
				Bin= scaf2bin[scaffold]
				if Bin in bestBins:
					if Bin == bestBins[Bin]:
						lineOut= "\t".join([cluster, scaffold, Bin])
						break
					else:
						lineOut.append("\t".join([cluster, scaffold, "", Bin, bestBins[Bin]]))
				else:
					lineOut.append("\t".join([cluster, scaffold, "", Bin, ""]))
		
	if isinstance(lineOut, basestring):
		print lineOut
	else:
		print "\n".join(lineOut)
	continue
	
exit()
	