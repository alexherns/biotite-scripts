#!/usr/bin/env python2.7
import argparse, sys

parser = argparse.ArgumentParser(description='Rapid scan of protein primary structure motifs', formatter_class=argparse.ArgumentDefaultsHelpFormatter, add_help=False)

#Required arguments
required = parser.add_argument_group('REQUIRED')
required.add_argument('-i', help= 'input fasta', default= sys.stdin, type=str)

#Optional arguments
optional = parser.add_argument_group('OPTIONAL')
optional.add_argument('-h', action="help", help="show this help message and exit")
optional.add_argument('-s', '--string', help= 'Motif search string', type=str)
optional.add_argument('-m', '--motif', help= 'Choose from list of pre-formatted motifs', type=str)
optional.add_argument('--min', help= 'Minimum number of matches to motif to print result', type= int, default= 1)
optional.add_argument('--matches', help= 'Print the number of motif matches in header', action= 'store_true')
optional.add_argument('--list', help= 'List pre-formatted motifs available', action= 'store_true')
optional.add_argument('--regex', help= 'Supply regex rather than protein motif', action= 'store_true')
optional.add_argument('--prosite', help= 'Supply PROSITE motif rather than regex', action= 'store_true')

args = parser.parse_args()

motif_db= {
'heme_binding': 'CXXCH',
'hydrogenase_group1_N': '[EGMQS]RXC[GR][IV]CXXX[HT]XXX[AGS]X(0,4)[VANQD]',
'hydrogenase_group1_C': '[AFGIKLMV][HMR]XX[HR][AS][AFLY][DN]PC[FILMV]XC[AGS]XH',
'hydrogenase_group2a_N': 'PR[AIV]CGICX(1,3)HX(0,2)LXX[AST]',
'hydrogenase_group2a_C': 'Vx[KR]S[FHY]DxCxVC[ST][TV][HK]',
'hydrogenase_group2b_N': 'PR[IV]CGICS[IV][AS]Q[GS]xA',
'hydrogenase_group2b_C': 'H[IV]VRSFDPCMVCT[AV]H',
'hydrogenase_group3a_N': 'R[FIV]CG[ILV]C[PQ]x[APT]H[ACGT]x[AS][AGS]',
'hydrogenase_group3a_C': 'R[ACS]YD[IP]C[AILV][AS]Cx(2,3)Hx[ILMV]',
'hydrogenase_group3b_N': 'R[IV]C[AGS][FIL]Cxxx[HY]xx[AST][ANS]xx[AS][AILV]',
'hydrogenase_group3b_C': 'R[ANS][FHY]DPCISC[AS][ATV]H',
'hydrogenase_group3c_N': 'Px[FILV][TV][ADPST]x[IV]CG[IV]CxxxHxx[AC][AS]xxA',
'hydrogenase_group3c_C': 'E[FMV][AGLV][FIV]Rx[FY]DPCx[AS]C[AS][ST]Hx[AILV]',
'hydrogenase_group3d_N': 'Ex[APV]xxxxRxCG[IL]Cxx[AS]Hx[IL][ACS][AGS][AGNSV][KR][ATV]xD',
'hydrogenase_group3d_C': 'DPC[IL]SC[AS][AST]H[ASTV]x[AG]xx[APV]',
'hydrogenase_group4_N': 'C[GS][ILV]C[AGNS]xxH',
'hydrogenase_group4_C': '[DE][PL]Cx[AGST]Cx[DE][RL]',
'FeFe_L1': 'TSC(2,3)PxW',
'FeFe_L2': 'MPCxxKxxE',
'FeFe_L3': 'ExMACxxGCxxGGGxP',
'binuclear_copper': 'CxxxCxxxHxxM'
}

if args.list:
	for motif in sorted(motif_db.keys()):
		print motif, motif_db[motif]
	exit()

if args.string == None and args.motif == None:
	parser.error('You must supply either a search string or motif identifier')

import re
from Bio import SeqIO

def motif2regex(motif):
    """Converts a protein motif query to a regex search"""
    return motif.replace('X', '.').replace('x', '.').replace('(', '{').replace(')', '}')

search= args.string

if args.regex:
	pass
elif args.prosite:
	search= motif2regex(search)
elif args.motif != None:
	search= motif2regex(motif_db[args.motif])
else:
	parser.error('You must supply either a regex string or prosite motif')

re_search= re.compile(search)

for seq_record in SeqIO.parse(args.i, "fasta"):
	num_matches= len(re_search.findall(str(seq_record.seq)))
	if num_matches >= args.min:
		if args.matches:
			print ">{0}-matches={2}\n{1}".format(seq_record.id, re.sub("(.{70})", "\\1\n", str(seq_record.seq), 0, re.DOTALL), str(num_matches))
		else:
			print ">{0}\n{1}".format(seq_record.id, re.sub("(.{70})", "\\1\n", str(seq_record.seq), 0, re.DOTALL))
