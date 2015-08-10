#!/usr/bin/env python2.7
import pysam, os, argparse, json, subprocess

parser = argparse.ArgumentParser(description='Maps reads to contigs on biotite.', 
	formatter_class=argparse.ArgumentDefaultsHelpFormatter, add_help=False)

#Required arguments
required = parser.add_argument_group('REQUIRED')
required.add_argument('-f', help= 'input fasta', required=True, type=str)

#Optional arguments
optional = parser.add_argument_group('OPTIONAL')
optional.add_argument('-h', action="help", help="show this help message and exit")
optional.add_argument('-p1', help= 'forward paired-end reads', type= str)
optional.add_argument('-p2', help= 'reverse paired-end reads', type= str)
optional.add_argument('-j', help= 'json containing read locations', type= str)
optional.add_argument('-r', help= 'read names if from json', type= str)
optional.add_argument('-t', help= 'Number of threads for running Bowtie mapping',
	default= 6, type= int)
optional.add_argument('-x', help= 'Additional bowtie2 parameters for mapping',
	default= '', type= str)

args = parser.parse_args()

if (args.p1 != None and args.p2 == None) or (args.p1 == None and args.p2 != None):
	parser.error('Paired end reads must be supplied together')

if args.p1 != None and args.j != None:
	parser.error('JSON file and direct path to reads cannot be supplied together')
	
if (args.j != None and args.r == None) or (args.j == None and args.r != None):
	parser.error('JSON file and read code must be supplied together')

###Create bowtie2 index
while True:
	try:
		os.mkdir('bt2')
	except OSError:
		break
		
bt2db= '.'.join(args.f.split('.')[:-1])
print 'Building bt2db from source fasta...'
os.system('bowtie2-build -q {0} bt2/{1}'.format(args.f, bt2db))


###Map to reference sequences using reads, or from json
if args.j != None:
	print args.j
	with open(args.j) as handle:
		read_json= json.load(handle)
	if args.r in read_json['Japanese Borehole Project']['reads']:
		args.p1= read_json['Japanese Borehole Project']['reads'][args.r][0]
		args.p2= read_json['Japanese Borehole Project']['reads'][args.r][1]
	else:
		parser.error('Reads not found in JSON!')

sam_base= '{0}-{1}'.format(args.f, ('vs-' + args.r if args.r != None else 'mapped'))
bowtie2_command= '''/opt/bin/bio/bowtie2 {5} -p {0} -x bt2/{1} -1 {2} -2 {3} \
| /opt/bin/bio/shrinksam -v > {4}.sam'''.format(str(args.t), bt2db, args.p1, args.p2, 
	sam_base, args.x)

try:
	subprocess.call(bowtie2_command, shell=True)
except:
	exit()

###Convert sam to sorted bam file
print 'Indexing fasta file...'
try:
	subprocess.call('samtools faidx {0}'.format(args.f), shell=True)
except:
	print 'Fasta indexing FAILED!'
	exit()
print 'Converting sam to sorted bam...'
try:
	subprocess.call(
		'samtools view -bS {0}.sam | samtools sort - {0}.sorted'.format(sam_base), 
		shell=True)
except:
	print 'SAM to BAM conversion FAILED!'
	exit()
print 'Indexing bam file...'
try:
	subprocess.call('samtools index {0}.sorted.bam'.format(sam_base), shell=True)
except:
	print 'BAM indexing FAILED!'
	exit()


