#!/bin/bash
#June 28, 2018 by Teresita M. Porter
#Get sequence stats from a directory of files named *.uniques
#stats links to 
#USAGE sh run_fastastats_parallel_uniques.sh

echo -e sample'\t'numseq'\t'minlength'\t'maxlength'\t'meanlength'\t'modelength

#ngs-workflow only has 8 cores so don't run more than 7 at a once
NR_CPUS=10
count=0

for f in *.uniques
do

stats $f &

let count+=1 
[[ $((count%NR_CPUS)) -eq 0 ]] && wait

done
	
wait

echo "All jobs are done"

