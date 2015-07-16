#!/usr/bin/env python
import networkx as nx
import matplotlib.pyplot as plt
import networkx_viewer as nv

import sys, getopt, os, re

opts, args = getopt.getopt(sys.argv[1:], 'c:o:m:h', ['help', 'connections', 'output', 'minimum'])

def usage():
        print """
Will create a directory to store contigs binned by ESOM mapping.

Usage: class2fasta.py -c sam.connections -o output.png

OPTIONS:
    -c, --connections	Connections file
    -o, --output		Output image file
    -m, --minimum		Minimum number of connections for scaffolds to be drawn together
"""
        exit()

connections= ''
output= ''
min_connect= 0

for o, a in opts:
    if o in ('-h', '--help'):
        usage()
    elif o in ('-c', '--connections'):
        connections= a
    elif o in ('-o', '--output'):
        output= a
    elif o in ('-m', '--minimum'):
    	min_connect= int(a)

if '' in [connections]:
        usage()

#Build the graph        
G= nx.Graph()
nodes= []
edges= {}
for line in open(connections):
	line= line.strip().split('\t')
	if 'accept' not in line or 'flanking' in line:
		continue
	attr= {}
	line[0]= re.search('(NODE_\d+)', line[0]).group()
	line[2]= re.search('(NODE_\d+)', line[2]).group()
	if line[0]==line[2]:
		nodes.append([line[0], {'self':'True', 'direction':" ".join(line[:4]), 'count':line[4]}])
		print line[0]+"\tSelf-edge"
		continue
	if line[0] not in nodes:
		nodes.append(line[0])
	if line[2] not in nodes:
		nodes.append(line[2])
	edge= sorted([line[0], line[2]])
	lookup= "\t".join(edge)
	if lookup in edges:
		continue
	if 'mid' in [line[1], line[3]]:
		attr= {'fill':'red'}
	attr['direction']= "  ".join(line[:4])
	attr['count']= line[4]
	if int(attr['count'])<min_connect:
		continue
	edge.append(attr)
	edges[lookup]= edge

G.add_nodes_from(nodes)
G.add_edges_from(edges.values())

#Draw the graph
app= nv.Viewer(G)
app.mainloop()
