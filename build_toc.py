#!/usr/bin/env python2.7
import os, re, argparse

parser = argparse.ArgumentParser(description='''Prints details for all current scripts in ~/scripts/''', 
 formatter_class=argparse.ArgumentDefaultsHelpFormatter, add_help=False)

#Optional arguments
optional = parser.add_argument_group('OPTIONAL')
optional.add_argument('-h', action="help", help="show this help message and exit")
optional.add_argument('--colored', action= 'store_true', help= 'format colors STDOUT')

args = parser.parse_args()


class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'


ignore= [line.strip() for line in open('/home/alexh/scripts/toc/.toc_ignore')]

for script in sorted(os.listdir('/home/alexh/scripts/')):	
	if script in ignore:
		continue
	elif script[-2:] == 'py':
		handle= ''.join(open('/home/alexh/scripts/{0}'.format(script)).readlines())
		line= handle.split('ArgumentParser')[1].split(')')[0]
		match= re.search("ArgumentParser\(.*(description=.*\'),.*formatter_class", handle, re.DOTALL)
		print (bcolors.BOLD + bcolors.OKGREEN + script + bcolors.ENDC if args.colored else script)
		print '\t' + match.group(1) + '\n'
	elif script[-2:] == 'sh':
		handle= ''.join(open('/home/alexh/scripts/{0}'.format(script)).readlines())
		match= re.search("Description.*", handle)
		print (bcolors.BOLD + bcolors.OKGREEN + script + bcolors.ENDC if args.colored else script)
		print '\t' + match.group(0) + '\n'
