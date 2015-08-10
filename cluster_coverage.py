#!/usr/bin/env python2.7
import sys, operator, argparse
from Bio import SeqIO

parser = argparse.ArgumentParser(description='''Prints out the coverage values for each cluster, by sample and total.
Also lists number of hits in each cluster.''', formatter_class=argparse.ArgumentDefaultsHelpFormatter, add_help=False,
epilog= '''TSV of features and as downloaded from ggkbase.
Scaffold_gene is in column 2.
Coverage value is in column 5.

Clusters file as generated from USEARCH
''')

#Required arguments
required = parser.add_argument_group('REQUIRED')
required.add_argument('-c', help= 'clusters.uc', required=True, type=str)
required.add_argument('-t', help= 'features.tsv', required=True, type=str)

#Optional arguments
optional = parser.add_argument_group('OPTIONAL')
optional.add_argument('-h', action="help", help="show this help message and exit")

args = parser.parse_args()

cluster_file= args.c
tsv_file= args.t


#Create a dictionary of feature:coverage values
#Read in the tsv of features
handle= open(tsv_file, "r")
feat2cov= {}
samples= []
for line in handle:
    contig_features= line.strip().split("\t")
    samples.append(contig_features[1].split("_scaffold")[0])
    feature, coverage= contig_features[1], contig_features[4]
    feat2cov[feature]= float(coverage)
samples= list(set(samples))
handle.close()

#Select all non-redundant cluster lines from file
clusters= [line.strip().split("\t") for line in open(cluster_file) if line[0] in ["H", "C"]]

#Extract unique list of all clusters
cluster_names= list(set([line[1]for line in clusters]))

#Dictionary of clusters:
#   clust_dict[cluster_name: [clust1, ..., clustN]]
clust_dict= {}
for cluster in clusters:
    if cluster[1] not in clust_dict:
        clust_dict[cluster[1]]= []
    clust_dict[cluster[1]].append(cluster)

#List to contain output lines
cov_list= []
for cluster in clust_dict:
    #Each line in output, formatted as list
    clustercov= [cluster]+[0]*(len(samples)+3)
    for line in clust_dict[cluster]:
        scaf= line[8]
        #Append centroids
        if line[0]=="C":
            clustercov.append(scaf)
        sample= scaf.split("_scaffold")[0]
        if sample not in samples:
            print "FAIL: SCAF", scaf
        else:
            clustercov[samples.index(sample)+1]+=feat2cov[scaf.split(" ")[0]]
    #Number of samples with positive hits
    clustercov[-2]= len([i for i in clustercov[1:-4] if i > 0])
    #Number of hits
    clustercov[-3]= len(clust_dict[cluster])
    #Total (raw and not normalized) cluster coverage value
    clustercov[-4]= sum(clustercov[1:-4])
    
    cov_list.append(clustercov)

#Print header line
print "TAX\t"+"\t".join(samples)+"\tTotal\t#Hits\t#Samples\tCentroid"

#Print each line in output
print "\n".join(["\t".join([str(i) for i in row]) for row in cov_list])
