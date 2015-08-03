#!/usr/bin/env python2.7
import argparse
parser = argparse.ArgumentParser(description='Appends each line by the string provided', formatter_class=argparse.ArgumentDefaultsHelpFormatter, add_help=False)

#Required arguments
required = parser.add_argument_group('REQUIRED')
required.add_argument('-f', help= 'input file', required=True, type=argparse.FileType('r'))
required.add_argument('-s', help= 'string to append', required=True, type=str)

#Optional arguments
optional = parser.add_argument_group('OPTIONAL')
optional.add_argument('-h', action="help", help="show this help message and exit")

args = parser.parse_args()

print "\n".join([line.strip("\n")+args.s for line in args.f])
