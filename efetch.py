#!/usr/bin/env python2.7
import sys, argparse, glob, os, subprocess

parser = argparse.ArgumentParser(description='Extracts records from Entrez or UniProt databases', 
	formatter_class=argparse.ArgumentDefaultsHelpFormatter, add_help=False,
	epilog= 'Note: Using unknown_uniprot_id may return undesired results. Use true accession ID whenever possible')

entrez_dbs= ['protein', 'nuccore', 'pubmed', 'gene', 'taxonomy', 'uniprot']

#Required arguments
required = parser.add_argument_group('REQUIRED')
required.add_argument('--db', help= 'Database service used', required=True, type=str, choices= entrez_dbs)

#Optional arguments
optional = parser.add_argument_group('OPTIONAL')
optional.add_argument('-h', '--help', action="help", help="show this help message and exit")
optional.add_argument('--rettype', help= 'EFetch &rettype setting', type= str, default= 'fasta')
optional.add_argument('--retmode', help= 'EFetch &retmode setting', type= str, default= 'text')
optional.add_argument('--id', help= 'Query id for retrieval from EFetch', type= str)
optional.add_argument('--unknown_uniprot_id', help= 'Required if uniprot accession id is unknown... slow!', action= 'store_true')

args = parser.parse_args()




def efetch(id):
	from Bio import Entrez
	
	Entrez.email= 'example@other.com'
	
	results= ''
	if len(id) > 500:
		batch= 500
		for start in range(0, len(id), batch):
			stop= (start+batch if start+batch <= len(id) else len(id))
			results+= Entrez.efetch(db=args.db, rettype=args.rettype, retmode=args.retmode, id= id[start:stop]).read()
			sys.stderr.write("batch complete: {0} - {1}\n".format(str(start), str(stop-1)))
	else:
		results= Entrez.efetch(db=args.db, rettype=args.rettype, retmode=args.retmode, id= id).read()
	return results
	
def unifetch(id):
	import urllib2
	results= ''
	if args.unknown_uniprot_id:
		sys.stderr.write('Uniprot accessions will be queried individually!\n')
		for i, protein in enumerate(id):
			sys.stderr.write('Processing ({1}/{2})\t{0}...\r'.format(protein, i, len(id)))
			response= urllib2.urlopen('http://www.uniprot.org/uniprot/{0}.fasta'.format(protein))
			results+= response.read()
	else:
		query= '+OR+'.join(['id:'+protein for protein in id])
		results= urllib2.urlopen('http://www.uniprot.org/uniprot/?query={0}&format=fasta'.format(query)).read()
	return results
	
if args.id == None:
	#Read from stdin
	id= [line.strip() for line in sys.stdin]
else:
	id= [args.id]

if args.db == 'uniprot':
	results= unifetch(id)
else:
	results= efetch(id)
	
print '\n'.join([line.strip() for line in results.strip().split('\n') if line.strip() != ''])
