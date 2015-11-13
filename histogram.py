#!/usr/bin/env python2.7
import fileinput, numpy, sys, argparse
from collections import Counter

parser = argparse.ArgumentParser(description='Generates histogram bins from numerical data.', 
	formatter_class=argparse.ArgumentDefaultsHelpFormatter, add_help=False, 
	epilog= '''Bin definitions are as generated from numpy.histogram -- 
	All but the last (righthand-most) bin is half-open. In other words, if bins is: 
	[1, 2, 3, 4] then the first bin is [1, 2) (including 1, but excluding 2) and the second [2, 3). 
	The last bin, however, is [3, 4], which includes 4.''')

#Required arguments
required = parser.add_argument_group('REQUIRED')
required.add_argument('-i', nargs='?', type=argparse.FileType('r'), default=sys.stdin)

#Optional arguments
optional = parser.add_argument_group('OPTIONAL')
optional.add_argument('-h', action="help", help="show this help message and exit")
optional.add_argument('-o', action="store_true", help="list is non-numerical\
data")
optional.add_argument('-b', help= 'Number of bins', type=int, default=10)
optional.add_argument('-s', help= 'Size of each bin', type=float)
optional.add_argument('--min', help= 'Minimum for range of bins', type=float)
optional.add_argument('--max', help= 'Maximum for range of bins', type=float)
optional.add_argument('--pretty', action= "store_true", help= 'Prints x\'s to visualize height of bins')
optional.add_argument('--sort_numeric', action="store_true", help="Sorts a\
        non-numerical histogram by frequency")
optional.add_argument('--sort_alpha', action='store_true', help='Sorts a\
        non-numerical histogram by key')

args = parser.parse_args()

if len([x for x in (args.min,args.max) if x is not None]) == 1:
   parser.error('--min and --max must be given together')
if args.s is not None and args.b != 10 > 1:
   parser.error('-s and -b are mutually exclusive')

data= []
for line in args.i:
    if args.o:
        data.append(line.strip())
    else:
        data.append(float(line.strip()))

def print_bins(heights, bins):
    for a, b in zip(heights, bins):
	if not args.pretty:
            print "{0}\t{1}".format(b, a)
	else:
            print "{0}\t{1}\t{2}".format(b, a, 
                    ''.join(["x" for i in 
                        range(int(float(a)/max(100, max(heights))*100))]
                        )
                    )

if args.o:
    c= Counter(data)
    counted= c.items()
    if args.sort_numeric:
        counted.sort(key=lambda x: x[1])
    elif args.sort_alpha:
        counted.sort()
    bins, heights= zip(*counted)
    print_bins(heights, bins)
    exit()

if args.min != None:
	if args.s != None:
		bins= numpy.arange(args.min, args.max+1, step=args.s)
	elif args.b != None:
		bins= numpy.linspace(args.min, args.max, num=args.b)
	heights, bins= numpy.histogram(data, bins)
elif args.s != None:
	bins= range(int(min(data)-args.s), int(max(data)+args.s), int(args.s))
	heights, bins= numpy.histogram(data, bins)
elif args.b != None:
	bins= numpy.linspace(min(data), max(data), num=args.b)
	heights, bins= numpy.histogram(data, bins)
else:
	heights, bins= numpy.histogram(data)

print_bins(heights, bins)
