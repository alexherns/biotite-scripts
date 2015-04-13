#!/usr/bin/env python2.7
import sys, getopt

if len(sys.argv) < 3:
    print """
Usage: compare_blast.py <input1.b6> <input2.b6>
"""

blast1= sys.argv[1]
blast2= sys.argv[2]


dict1= {}
for line in open(blast1):
    line= line.strip().split("\t")
    query= line[0].split()[0]
    dict1[query]= line

dict2= {}
for line in open(blast2):
    line= line.strip().split("\t")
    query= line[0].split()[0]
    dict2[query]= line

blast1better= blast1.split(".b6")[0] + "_good.b6"
blast1worse= blast1.split(".b6")[0] + "_bad.b6"
blast2better= blast2.split(".b6")[0] + "_good.b6"
blast2worse= blast2.split(".b6")[0] + "_bad.b6"

b1b_handle= open(blast1better, 'w')
b1w_handle= open(blast1worse, 'w')
b2b_handle= open(blast2better, 'w')
b2w_handle= open(blast2worse, 'w')

for query in dict1:
    line1= dict1[query]
    if query not in dict2:
        b1b_handle.write("\t".join(line1)+"\n")
    else:
        line2= dict2[query]
        if float(line1[-1]) > float(line2[-1]):
            b1b_handle.write("\t".join(line1)+"\n")
        else:
            b1w_handle.write("\t".join(line1)+"\n")

for query in dict2:
    line2= dict2[query]
    if query not in dict1:
        b2b_handle.write("\t".join(line2)+"\n")
    else:
        line1= dict1[query]
        if float(line2[-1]) > float(line1[-1]):
            b2b_handle.write("\t".join(line2)+"\n")
        else:
            b2w_handle.write("\t".join(line2)+"\n")

b1b_handle.close()
b1w_handle.close()
b2b_handle.close()
b2w_handle.close()
