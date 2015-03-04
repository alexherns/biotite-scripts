#!/usr/bin/env python
'''
Repairs improperly formatted GFF3 files.
A common formatting error is to leave the 9th column of a feature line blank,
which causes issues when trying to import the features into Tablet. This simple
script repairs these errors by adding a "." in any blank columnes of a feature.
'''

import sys

for line in open(sys.argv[1], 'r'):
	if line[0] == '#' or len(line.strip().split('\t')) == 9:
		print line.strip()
	else:
		print line.strip()+'\t.'
