#!/bin/bash
#Script to run fasta_gz_stats.plx on a directory of fasta files
#symbolic link gz_stats points to this script
#USAGE gz_stats gz

echo -e sample'\t'numseqs'\t'minlength'\t'maxlength'\t'meanlength'\t'medianlength'\t'modelength

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
