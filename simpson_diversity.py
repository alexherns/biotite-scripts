#!/usr/bin/env python2.7
from fileinput import input
from math import log

"""
Prints out Simpson diversity index (D) and Simpson's equitability E_D.
Accepts input from STDIN or from file names.

D= 1/sum(pi^2 for i: i=1 -> S)
E_D= D/D_max= D/S

where:  pi= proportion of species i relative to the total number of species
        S= total number of species in the community (richness)
"""

abundances= [float(line.strip()) for line in input()] #Read in abundance data
abundances= [i/sum(abundances) for i in abundances] #Calculate pi

D= 1/sum([pi**2 for pi in abundances])
E_D= D/len(abundances)

print "{0}\t{1}".format(D, E_D)

