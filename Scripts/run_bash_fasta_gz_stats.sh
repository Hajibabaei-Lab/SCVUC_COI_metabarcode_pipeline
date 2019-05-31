#!/bin/bash
#June 28, 2018 by Teresita M. Porter
#Scripts gets sequence stats from a directory of fasta.gz files
#stats3 links to fasta_gz_stats.plx
#be sure to include file extension that targets the fasta.gz files as a command-line argument
#USAGE sh run_bash_fasta_gz_stats.sh gz

echo -e sample'\t'numseqs'\t'minlength'\t'maxlength'\t'meanlength'\t'medianlength'\t'modelength

NR_CPUS=20
count=0

EXT="$1"

for f in *$EXT
do

stats3 $f &

let count+=1
[[ $((count%NR_CPUS)) -eq 0 ]] && wait

done
	
wait

echo "All jobs are done"
