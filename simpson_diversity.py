#!/usr/bin/env python2.7
from math import log
import argparse, sys

parser = argparse.ArgumentParser(description='''Prints out Simpson diversity
index (D) and Simpson's equitability E_D.
H= -sum(pi*ln(pi) for i: i=1 -> S)
E_H= H/H_max= H/ln(S)

where:  pi= proportion of species i relative to the total number of species
        S= total number of species in the community (richness)''', 
	formatter_class=argparse.ArgumentDefaultsHelpFormatter, add_help=False)

#Required arguments
required = parser.add_argument_group('REQUIRED')
required.add_argument('-i', nargs='?', type=argparse.FileType('r'), 
	help= 'Input list of species proportions', default=sys.stdin)

#Optional arguments
optional = parser.add_argument_group('OPTIONAL')
optional.add_argument('-h', action="help", help="show this help message and exit")

args = parser.parse_args()

abundances= [float(line.strip()) for line in args.i] #Read in abundance data
abundances= [i/sum(abundances) for i in abundances] #Calculate pi

D= 1/sum([pi**2 for pi in abundances])
E_D= D/len(abundances)

print "{0}\t{1}".format(D, E_D)

