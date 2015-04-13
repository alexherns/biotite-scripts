#!/usr/bin/env python2.7
import sys

if len(sys.argv) < 2:
    print """
Usage: sam_header.py <input.sam>
Output: header of <input.sam>

Notes:  Script will terminate upon first occurrence of a
non-header formatted line, so make sure your SAM header
is formatted properly. This ensures rapid processing of file.
"""
    exit()

for line in open(sys.argv[1]):
    if line[0] == "@":
        print line.strip()
    else:
        exit()
