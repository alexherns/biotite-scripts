#!/bin/bash
if [ -z "$1" ]
then
	echo "Usage is as follows:	sam_slicer.sh <desired sam> <scaffolds.list> <output.sam>"
	exit 1
fi
cat $2 > temp.sliced.txt
app_each_line.py temp.sliced.txt "\s" > temp2.sliced.txt
sed -i "1i^@HD\n^@PG\n^@CO" temp2.sliced.txt
grep -f temp2.sliced.txt $1 > $3
rm temp2.sliced.txt temp.sliced.txt
