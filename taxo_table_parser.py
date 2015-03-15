#!/usr/bin/env python2.7
import sys, getopt

def usage():
        print '''
Useful script for extracting lines of interest from taxonomy table file as downloaded from ggkbase

Usage: taxo_table_parser.py <tab_sep_taxo_table> <options>

	-h, --help	Print this help dialog

	<ARGUMENTS>
	--min_gc=<float>  Default 0
	--max_gc=<float>  Default 100
	--min_cov=<float>   Default 0
	--max_cov=<float>   Default "inf"
	--min_size=<float>  Default 0
	--tax_def=<string>  e.g. species, genus, or class
	--tax_query=<string>    e.g. Geobacter uraniireducens
	--tax_cutoff=<float>    Default 0
        --exact
'''
        exit()

opts, args = getopt.getopt(sys.argv[2:], 'h', ['help', 'min_gc=', 'max_gc=', 'min_cov=', 'max_cov=', 'min_size=', 'tax_def=', 'tax_query=', 'tax_cutoff=', 'exact'])

min_gc = 0.
max_gc = 100.
min_cov = 0.
max_cov = float("inf")
min_size = 0.
tax_def, tax_query, tax_cutoff = [None, None, 0.]
exact= False

for o, a in opts:
	if o in ('-h', '--help'):
		usage()
	elif o == '--min_gc':
		min_gc = float(a)
	elif o == '--max_gc':
		max_gc = float(a)
	elif o == '--min_cov':
		min_cov = float(a)
	elif o == '--max_cov':
		max_cov = float(a)
	elif o == '--min_size':
		min_size = float(a)
	elif o == '--min_gc':
		min_gc = float(a)
	elif o == '--tax_def':
		tax_def = a
	elif o == '--tax_query':
		tax_query = a
	elif o == '--tax_cutoff':
		tax_cutoff = a
        elif o == '--exact':
                exact= True
print exact

if (tax_def is not None and tax_query is None) or (tax_query is not None and tax_def is None):
	sys.stderr.write("ERROR: If selecting by taxonomy, you MUST provide both a taxonomic group to search and a string to identify")
	exit()

class contig:
	def __init__(self, contig_detail_list):
		name, size, coverage, gc, taxo_winner, taxo_winner_pc, \
		species_winner, species_winner_pc, genus_winner, genus_winner_pc, \
		order_winner, order_winner_pc, class_winner, class_winner_pc, \
		phylum_winner, phylum_winner_pc, domain_winner, domain_winner_pc = contig_detail_list
		self.name = name
		self.size = float(size)
		self.coverage = float(coverage)/100
		self.gc = float(gc)/100
		if taxo_winner == '':
			self.taxo_winner = None
			self.taxo_winner_pc = None
		else:
			self.taxo_winner = taxo_winner
			self.taxo_winner_pc = float(taxo_winner_pc)
		self.species_winner = species_winner
		self.species_winner_pc = float(species_winner_pc)
		self.genus_winner = genus_winner
		self.genus_winner_pc = float(genus_winner_pc)
		self.order_winner = order_winner
		self.order_winner_pc = float(order_winner_pc)
		self.class_winner = class_winner
		self.class_winner_pc = float(class_winner_pc)
		self.phylum_winner = phylum_winner
		self.phylum_winner_pc = float(phylum_winner_pc)
		self.domain_winner = domain_winner
		self.domain_winner_pc = float(domain_winner_pc)
	
	def toString(self):
		return '\t'.join([self.name, str(self.size), str(self.coverage), str(self.gc), str(self.taxo_winner), str(self.taxo_winner_pc), \
		self.species_winner, str(self.species_winner_pc), self.genus_winner, str(self.genus_winner_pc), \
		self.order_winner, str(self.order_winner_pc), self.class_winner, str(self.class_winner_pc), \
		self.phylum_winner, str(self.phylum_winner_pc), self.domain_winner, str(self.domain_winner_pc)])

contig_list = []
print "Contig name	Size (bp)	Coverage	GC %	Taxonomy winner	Winner %	Species winner	Species winner %	Genus winner	Genus winner %	Order winner	Order winner %	Class winner	Class winner %	Phylum winner	Phylum winner %	Domain winner	Domain winner %"
for i, line in enumerate(open(sys.argv[1],'r')):
	if line.strip().split()[0] == 'Contig':
		continue
	contig_detail_list = line.strip().split('\t')
	if len(contig_detail_list)<18:
		sys.stderr.write("Contig {0} contains too few parameters to initialize\n".format(contig_detail_list[0]))
		continue
        this_contig = contig(contig_detail_list)
	if this_contig.gc >= min_gc and this_contig.gc <= max_gc \
		and this_contig.coverage >= min_cov and this_contig.coverage <= max_cov \
		and this_contig.size >= min_size:
		if tax_def is not None:
                        if exact:
                            exec("if this_contig.{0}_winner == '{1}' and this_contig.{0}_winner_pc >= {2}:\n\tprint this_contig.toString()\n\tcontig_list.append(this_contig)".format(tax_def, tax_query, tax_cutoff))
                        else:
			    exec("if '{1}' in this_contig.{0}_winner and this_contig.{0}_winner_pc >= {2}:\n\tprint this_contig.toString()\n\tcontig_list.append(this_contig)".format(tax_def, tax_query, tax_cutoff))
		else:
			print this_contig.toString()
			contig_list.append(this_contig)
			
print "Total Size:\t{0}\tNumber of contigs:{1}".format(sum([this_contig.size for this_contig in contig_list]), len(contig_list))
	
