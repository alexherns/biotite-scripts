#!/bin/bash
for fasta in $(ls -1 $1/*fasta)
do
	scaf=$(echo $fasta | cut -d. -f2)
	grep ">" $fasta | cut -c2- | sed 's/$/\t'$scaf'/g'
done
