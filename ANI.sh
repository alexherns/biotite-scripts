#!/bin/bash
#	Description:	Runs ANI on two input genomes and prints results to STDOUT
#	Usage:	ANI.sh	1.fasta	2.fasta
mkdir ANI_tmp
cp $1 ANI_tmp/
file1=`basename $1`
cp $2 ANI_tmp/
file2=`basename $2`

average_nucleotide_identity.py -i ANI_tmp -o ANI_$file1_$file2 --maxmatch
cat ANI_$file1_$file2/ANIm_percentage_identity.tab | cut -f 3 | awk 'NR==2{print "Percentage_identity:\t"$1}'
cat ANI_$file1_$file2/ANIm_alignment_coverage.tab | cut -f 3 | awk 'NR==2{print "Alignment_coverage:\t"$1}'
cat ANI_$file1_$file2/ANIm_alignment_lengths.tab | cut -f 3 | awk 'NR==2{print "Alignment_length:\t"$1}'
cat ANI_$file1_$file2/ANIm_similarity_errors.tab | cut -f 3 | awk 'NR==2{print "Similarity_error:\t"$1}'

rm -r ANI_tmp ANI_$file1_$file2
