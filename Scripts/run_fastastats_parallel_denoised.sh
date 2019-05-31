#!/bin/bash
#June 28, 2018 by Teresita M. Porter
#Script calculates read stats on FASTA-formatted files after unoise
#stats links to fasta_stats_parallel.sh
#USAGE sh run_fastastats_parallel.sh

echo -e sample'\t'numseq'\t'minlength'\t'maxlength'\t'meanlength'\t'modelength'\t'

NR_CPUS=10
count=0

for f in *.denoised
do

stats $f &

let count+=1 
[[ $((count%NR_CPUS)) -eq 0 ]] && wait

done
	
wait

echo "All jobs are done"
