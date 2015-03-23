#!/usr/bin/env python2.7
from fileinput import input
from math import log

"""
Prints out Shannon diversity index (H) and Shannon's equitability E_H in ln form.
Accepts arguments from STDIN or as file names.

H= -sum(pi*ln(pi) for i: i=1 -> S)
E_H= H/H_max= H/ln(S)

where:  pi= proportion of species i relative to the total number of species
        S= total number of species in the community (richness)
"""

abundances= [float(line.strip()) for line in input()] #Read in abundance data
abundances= [i/sum(abundances) for i in abundances] #Calculate pi
H= (-1)*sum([pi*log(pi) for pi in abundances])
EH= H/log(len(abundances))
print "{0}\t{1}".format(H, EH)
