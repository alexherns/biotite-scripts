#!/usr/bin/env python2.7
import sys
"""
Creates a fasta file with X numers of N's
"""

print ">Ns\n{0}\n".format(''.join(["N" for i in range(int(sys.argv[1]))]))
