#!/usr/bin/env python2.7
import fileinput, numpy, sys, argparse

parser = argparse.ArgumentParser(description='Generates histogram bins from numerical data.', formatter_class=argparse.ArgumentDefaultsHelpFormatter, add_help=False)

#Required arguments
required = parser.add_argument_group('REQUIRED')
required.add_argument('-i', nargs='?', type=argparse.FileType('r'), default=sys.stdin)

#Optional arguments
optional = parser.add_argument_group('OPTIONAL')
optional.add_argument('-h', action="help", help="show this help message and exit")
optional.add_argument('-b', help= 'Number of bins', type=int, default=10)
optional.add_argument('-s', help= 'Size of each bin', type=float)
optional.add_argument('--min', help= 'Minimum for range of bins', type=float)
optional.add_argument('--max', help= 'Maximum for range of bins', type=float)

args = parser.parse_args()

if len([x for x in (args.min,args.max) if x is not None]) == 1:
   parser.error('--min and --max must be given together')
if args.s is not None and args.b != 10 > 1:
   parser.error('--s and -b are mutually exclusive')


data= []
for line in args.i:
	data.append(float(line.strip()))
	
if args.min != None:
	if args.s != None:
		bins= numpy.linspace(args.min, args.max, retstep=args.s)
	else:
		bins= numpy.linspace(args.min, args.max, num=args.b)
	heights, bins= numpy.histogram(data, bins)
elif args.s != None:
	bins= range(int(min(data)-args.s), int(max(data)+args.s), int(args.s))
	print bins
	heights, bins= numpy.histogram(data, bins)
else:
	heights, bins= numpy.histogram(data)

for a, b in zip(heights, bins):
	print "{0}\t{1}".format(a, b)