#!/usr/bin/env python2.7
import sys
import numpy as np


frac= []
for line in open(sys.argv[1]):
	if line[0]=="#":
		continue
	line= line.strip().split("\t")
	format= line[9]
	depth, ref, alt= [int(format.strip().split(":")[i]) for i in [1,2,4]]
	frac.append(float(depth)/alt)
	
print np.mean(frac)
print np.std(frac)