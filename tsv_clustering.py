#!/usr/bin/env python2.7
import pandas as pd
import numpy as np
import sys
import matplotlib
from scipy.spatial.distance import pdist, squareform
from scipy.cluster.hierarchy import linkage, dendrogram
from matplotlib import pyplot as plt
from scipy.cluster import hierarchy
import sys, getopt

opts, args = getopt.getopt(sys.argv[1:], 'i:o:d:hl', ['help', 'input', 'output', 'delimiter', 'log-normalize'])

def usage():
        print """
Example script for plotting hierarchical clusters on top of heatmap
using various Python modules. Output to SVG. Will need to move axes
around to better suit incoming data.

Assumes row identifiers in first column.
Assumes column identifiers in first row.
Assumes all input data is numerical

Usage: tsv_clustering.py -i <INPUT> -o <fig.svg> [OPTIONS]

OPTIONS:

-d, --delimiter=        (string)    Delimiter - default "\t"
-l, --log-normalize     (flag)      Option to log-normalize data matrix (will remove all rows with zeros)
"""
        exit()

input_file= ''
output_file= ''
delim= '\t'
logNormal= False

for o, a in opts:
        if o in ('-h', '--help'):
                usage()
        elif o in ('-i', '--input'):
                input_file= a
        elif o in ('-o', '--output'):
                output_file= a
        elif o in ('-d', '--delimiter'):
                delim= a
        elif o in ('-l', '--log-normalize'):
                logNormal= True;

if input_file == '' or output_file == '':
        print """
Please specify input and output files. See -h (--help) for help
"""
        exit()


#Use AGG backend for output to SVG
matplotlib.use("Agg")

#Load in data from file
df= pd.read_csv(input_file, sep=delim, header=0, index_col=0)

if logNormal == True:
    #Remove rows with zeros, as they cannot be log-transformed
    df.replace(to_replace=0., value=np.nan, inplace=True)   #first needs to be set elements to nan
    df.dropna(axis=0, inplace=True) #drop any row with any nan values
    df= np.log(df)

#Generate pair-wise distance matrix for each row
row_dist= pd.DataFrame(squareform(pdist(df, metric='euclidean')), columns= df.index, index= df.index)

#Apply clustering using complete linkage method
#Compute pairwise distances for columns
col_dist= squareform(pdist(df.T, metric= "euclidean"))
col_clusters= linkage(col_dist, method= "complete")

row_clusters= linkage(row_dist, method= 'complete')

#Plot column dendrogram
fig= plt.figure(figsize= (8,8))
axd2= fig.add_axes([0.25, 0.9, .56, 0.10])
col_dendr= dendrogram(col_clusters, orientation= "top", color_threshold= np.inf)
axd2.set_xticks([])
axd2.set_yticks([])

#Plot the row-dendrogram
hierarchy.set_link_color_palette(['black']) #makes dendrogram black (1)
axd1= fig.add_axes([0.05, 0.1, 0.15, 0.6])
row_dendr= dendrogram(row_clusters, orientation= 'right', color_threshold=np.inf) #makes dendrogram black (2)
axd1.set_xticks([])
axd1.set_yticks([])

#Remove axes spines from dendrogram
for i, j in zip(axd1.spines.values(), axd2.spines.values()):
    i.set_visible(False)
    j.set_visible(False)

#reorder colums and rows with respect to the clustering
df_rowclust= df.ix[row_dendr['leaves']]
df_colrowclust= df_rowclust.ix[:][col_dendr["leaves"]]

#Plot the heatmap
axm= fig.add_axes([0.25, 0.1, 0.7, 0.6])
cax= axm.matshow(df_rowclust, interpolation= 'nearest', cmap= 'hot_r',  aspect='auto')
cax.set_clim(0., 10.)
fig.colorbar(cax)

axm.set_yticklabels(['']*len(list(df_rowclust.index)))
plt.xticks(np.arange(len(list(df_colrowclust.columns))), list(df_colrowclust.columns), rotation="vertical")

for tick in axm.xaxis.get_minor_ticks():
    tick.label1.set_verticalalignment('top')

print ['']+list(df_colrowclust.columns)


plt.savefig(output_file)


