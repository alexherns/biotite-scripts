#!/bin/bash
for fasta in $(ls -1 $1/*fa)
do
	scaf=$(echo $fasta | cut -d. -f4)
	grep ">" $fasta | cut -c2- | sed 's/$/\t'$scaf'/g'
done
