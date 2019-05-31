#!/bin/bash
#Teresita M. Porter, Sept. 12, 2018
#Script to get read counts from fasta headers for all *.uniques in a directory
#USAGE sh get_read_counts_uniques_3.sh

#set number of CPUs to use here if processing many infiles
NR_CPUS=5
count=0

echo -e 'sample\treadcount'

for f in *.uniques
do

base=${f%%.uniques*}
#echo $base

readcount=$(grep ">" $f | awk 'BEGIN {FS=";"} {print $3}' | sed 's/size=//g' | awk '{sum+=$1} END {print sum}')
#echo $readcount

echo -e $base'\t'$readcount

let count+=1 
[[ $((count%NR_CPUS)) -eq 0 ]] && wait

done
	
wait

echo "All jobs are done"
