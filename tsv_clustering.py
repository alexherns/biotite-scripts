#!/usr/bin/env python2.7
import pandas as pd
import numpy as np
import sys
import matplotlib
import sys, argparse

parser = argparse.ArgumentParser(
	description='''A less than perfect script to generate hierarchical clusters
	using Python heatmaps and outputs to SVG''', 
	formatter_class=argparse.ArgumentDefaultsHelpFormatter, add_help=False,
	epilog= '''NOTE: Assumes row identifiers in first column.
Assumes column identifiers in first row.
Assumes all input data is numerical''')

#Required arguments
required = parser.add_argument_group('REQUIRED')
required.add_argument('-i', help= 'input data', required=True, type=str)
required.add_argument('-o', help= 'output figure', required=True, type=str)

#Optional arguments
optional = parser.add_argument_group('OPTIONAL')
optional.add_argument('-h', action="help", help="show this help message and exit")
optional.add_argument('-b', help= 'Matplotlib backend', default= 'AGG', 
	type= str, choices= ['AGG', 'PS', 'PDS', 'SVG', 'Cairo', 'GDK'])
optional.add_argument('-d', help= 'delimiter', type=str, default= '\t')
optional.add_argument('-l', help= '''log-normalize data (and remove all non-zero
	rows''', action= 'store_true')

args = parser.parse_args()

#Defaults
input_file= args.i
output_file= args.o
delim= args.d
logNormal= args.l
backend= args.b

#Use AGG backend for output to SVG
matplotlib.use(backend)

#Load additional modules (needs to be after backend selection for matplotlib)
from scipy.spatial.distance import pdist, squareform
from scipy.cluster.hierarchy import linkage, dendrogram
from matplotlib import pyplot as plt
from scipy.cluster import hierarchy

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

#axm.set_yticklabels(list(df_rowclust.index), )
#plt.yticks(np.arange(len(list(df_rowclust.columns))), list(df_rowclust.columns))
plt.yticks(range(0,len(list(df_rowclust.index))), list(df_rowclust.index))
plt.xticks(np.arange(len(list(df_colrowclust.columns))), list(df_colrowclust.columns), rotation="vertical")

plt.tick_params(axis="y", labelsize=1)

for tick in axm.xaxis.get_minor_ticks():
    tick.label1.set_verticalalignment('top')

plt.savefig(output_file)


