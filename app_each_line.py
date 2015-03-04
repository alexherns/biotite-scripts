#!/usr/bin/env python2.7
"""Appends each line (first argument) by the string provided (second argument)"""
import sys
print "\n".join([line.strip("\n")+sys.argv[2] for line in open(sys.argv[1])])
