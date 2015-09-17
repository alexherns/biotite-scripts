#!/usr/bin/env python
import sys

for line in open(sys.argv[1]):
    if line[:5] == '# Id:':
        scaffold= line.strip().split(' ')[2]
    elif line[0] != '#' and line[0] != 'S':
        line= line.strip().split()
        if line == []:
            continue
        begin= int(line[0])
        length= int(line[1])
        end= begin+length
        sequence= line[2]
        print '>{0}_crispr_{1}-{2}\n{3}'.format(scaffold, str(begin), str(end), sequence)
