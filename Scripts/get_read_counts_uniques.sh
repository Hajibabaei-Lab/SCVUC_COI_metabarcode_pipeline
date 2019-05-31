#!/bin/bash
#June 28, 2018 by Teresita M. Porter
#Script to get read counts from fasta headers after vsearch dereplication
#USAGE sh get_read_counts_uniques.sh

NR_CPUS=10
count=0

echo -e 'sample\treadcount'

for f in *.uniques
do

base=${f%%.uniques*}

readcount=$(grep ">" $f | awk 'BEGIN {FS=";"} {print $3}' | sed 's/size=//g' | awk '{sum+=$1} END {print sum}')

echo -e $base'\t'$readcount

let count+=1 
[[ $((count%NR_CPUS)) -eq 0 ]] && wait

done
	
wait

echo "All jobs are done"
