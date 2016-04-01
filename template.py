#!/usr/bin/env python2.7
import argparse, os, sys


def main():
    sys.exit()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='SCRIPT DESCRIPTION', formatter_class=argparse.ArgumentDefaultsHelpFormatter, add_help=False)

    #Required arguments
    required = parser.add_argument_group('REQUIRED')
    required.add_argument('-X', help= 'required argument', required=True, type=None)

    #Optional arguments
    optional = parser.add_argument_group('OPTIONAL')
    optional.add_argument('-h', action="help", help="show this help message and exit")
    optional.add_argument('-x', help='optional argument', type=None)

    args = parser.parse_args()
    
    main()
