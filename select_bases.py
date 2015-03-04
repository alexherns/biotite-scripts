#!/usr/bin/env python2.7
import sys

'''Usage: select_bases.py <input.fasta> seq_id region_start region_end
Output: fasta file of sequencing matching region requested
Notes:	1-based numbering
	start and end are inclusive
	script terminates after finding first identifier
	prints error message if no sequences are found'''

line = open(sys.argv[1], 'r').read()[1:]
for id_read in line.split('\n>'):
    if id_read.split('\n')[0] == sys.argv[2]:
    #   Found the sequence of interest

        if len(sys.argv[2:]) == 1:
    #   No length parameters provided, select entire region
            print ">{0}".format(sys.argv[2])
            print "".join(id_read.split("\n")[1:])

        elif len(sys.argv[2:]) == 2:
    #   End site provided
            print ">{0}:{1}-{2}".format(sys.argv[2], 1, sys.argv[3])
            print "".join(id_read.split("\n")[1:])[0:int(sys.argv[3])]
        else:
    #   Start and end sites provided
            print '>{0}:{1}-{2}'.format(sys.argv[2], sys.argv[3], sys.argv[4])
            print ''.join(id_read.split('\n')[1:])[int(sys.argv[3])-1:int(sys.argv[4])]
        exit()

print 'No sequences were found with the designated identifier'
