#!/usr/bin/env python
import networkx as nx
import matplotlib.pyplot as plt
import networkx_viewer as nv

import sys, argparse, os, re

parser = argparse.ArgumentParser(description='Visualizes connections in assembly using networkx_viewer module.', formatter_class=argparse.ArgumentDefaultsHelpFormatter, add_help=False)

#Required arguments
required = parser.add_argument_group('REQUIRED')
required.add_argument('-c', help= 'connections file', required=True, type=argparse.FileType('r'))

#Optional arguments
optional = parser.add_argument_group('OPTIONAL')
optional.add_argument('-h', action="help", help="show this help message and exit")
optional.add_argument('-o', metavar='<*.png>', type=argparse.FileType('w'))
optional.add_argument('-m', metavar='<int>', type=int, default=0)

args = parser.parse_args()

#Build the graph        
G= nx.Graph()
nodes= []
edges= {}
for line in args.c:
	line= line.strip().split('\t')
	if 'accept' not in line or 'flanking' in line:
		continue
	attr= {}
	#line[0]= re.search('(NODE_\d+)', line[0]).group()
	#line[2]= re.search('(NODE_\d+)', line[2]).group()
	if line[0]==line[2]:
		nodes.append([line[0], {'self':'True', 'fill':'blue', 'direction':" ".join(line[:4]), 'count':line[4]}])
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
	if int(attr['count'])<args.m:
		continue
	edge.append(attr)
	edges[lookup]= edge

G.add_nodes_from(nodes)
G.add_edges_from(edges.values())

#Draw the graph
app= nv.Viewer(G)
app.mainloop()
