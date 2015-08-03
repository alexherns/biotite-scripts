#!/usr/bin/env python2.7
import os, argparse

parser = argparse.ArgumentParser(description='Assists with linking 16s genes from EMIRGE to "best" bin.', formatter_class=argparse.ArgumentDefaultsHelpFormatter, add_help=False)

#Required arguments
required = parser.add_argument_group('REQUIRED')
required.add_argument('-s', help= 'scaf2bin file', type=argparse.FileType('r'), required=True)
required.add_argument('-c', help= 'sorted clusters file output from UCLUST', type=argparse.FileType('r'), required=True)
required.add_argument('-b', help= 'best bin choices in format "best_bin"\tbin\tbin', type=argparse.FileType('r'), required=True)

#Optional arguments
optional = parser.add_argument_group('OPTIONAL')
optional.add_argument('-h', action="help", help="show this help message and exit")

args = parser.parse_args()
        
from common_objects import *

#Create scaffold2bin dictionary
scaf2bin= scaf2bin_dictionary(args.s, header=True, sep="\t")

#Create bin2bestbin dictionary
bestBins= {}
for line in args.b:
	bins= [item for item in line.strip().split(",")[1:] if item != ""]
	best_bin= bins.pop(0)
	for Bin in bins:
		bestBins[Bin]= best_bin
	
cluster_dict= {}	
for line in args.c:
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
	