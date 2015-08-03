#!/bin/bash
#	Description:	Repairs linebreaks in CSV file and prints to STDOUT
#	Usage:	csvRepair.sh	bad_format.csv

tr '\r' '\n' < $1 > tempCSVfile
cat tempCSVfile 
rm tempCSVfile
