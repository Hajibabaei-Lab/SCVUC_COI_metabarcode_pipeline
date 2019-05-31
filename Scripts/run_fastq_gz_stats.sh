#!/bin/bash
#June 28, 2018 by Teresita M. Porter
#This script gets read stats from fastq.gz files
#Script to run jobs in parallel using bash
#stats2 links to fastaq_gz_stats.plx
#USAGE sh run_fastq_gz_stats.sh

echo -e sample'\t'numseqs'\t'minlength'\t'maxlength'\t'meanlength'\t'median'\t'modelength

NR_CPUS=20
count=0

EXT="$1"

for f in *$EXT
do

stats2 $f &

let count+=1
[[ $((count%NR_CPUS)) -eq 0 ]] && wait

done
	
wait

echo "All jobs are done"
