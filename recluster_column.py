#!/usr/bin/env python
import sys

def recluster_columns(fpath, column, prefix='bin_', delim='\t'):
    uniq_cols = {}
    counter = 1
    for line in open(fpath, 'r'):
        split_line = line.strip('\n').split(delim)
        if len(split_line) <= column:
            print line.strip('\n')
            continue
        orig_data = split_line[column]
        if orig_data not in uniq_cols:
            uniq_cols[orig_data] = prefix + str(counter)
            counter += 1
        split_line[column] = uniq_cols[orig_data]
        line = delim.join(split_line)
        print line

if __name__ == '__main__':
    if len(sys.argv) < 3:
        print """recluster_columns.py <file> <column> [<prefix>]
        column -- 0-indexed
        prefix -- optional, default = bin_"""
        sys.exit(0)
    file_name = sys.argv[1]
    file_column = int(sys.argv[2])
    if len(sys.argv) > 3:
        prefix = sys.argv[3]
        recluster_columns(file_name, file_column, prefix=prefix)
    else:
        recluster_columns(file_name, file_column)


