#!/usr/bin/env python2.7
import networkx as nx
import matplotlib.pyplot as plt

import sys, getopt, os

opts, args = getopt.getopt(sys.argv[1:], 'c:o:h', ['help', 'connections', 'output'])

def usage():
        print """
Will create a directory to store contigs binned by ESOM mapping.

Usage: class2fasta.py -c sam.connections -o output.png

OPTIONS:
    -c, --connections	Connections file
    -o, --output		Output image file
"""
        exit()

connections= ''
output= ''

for o, a in opts:
    if o in ('-h', '--help'):
        usage()
    elif o in ('-c', '--connections'):
        connections= a
    elif o in ('-o', '--output'):
        output= a

if '' in [connections, output]:
        usage()

#Build the graph        
G= nx.Graph()
nodes= list(set([line.strip().split('\t')[0] for line in open(connections)]))
G.add_nodes_from(nodes)
input_lines= [line.strip().split('\t') for line in open(connections) if 'accept' in line and 'flanking' not in line]
edges= set()
for line in input_lines:
	edges.add(tuple(sorted([line[0], line[2]])))
for edge in edges:
	print edge
G.add_edges_from(edges)

#Draw the graph
nx.draw(G)
plt.savefig("{0}".format(output))