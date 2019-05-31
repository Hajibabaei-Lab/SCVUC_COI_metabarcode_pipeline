#!/bin/bash
#June 28, 2018 by Teresita M. Porter
#Script grabs fasta.gz files in a directory and adds the sample parsed from a GRDI-formatted file to the fasta headers
#rename_fasta links to rename_fasta_gzip.plx
#Be sure to include the file extension that targetsthe fasta.gz files as a command-line argument
#USAGE sh run_rename_fasta.sh fasta.gz

NR_CPUS=10
count=0

EXT="$1"

for f in *$EXT
do

rename_fasta $f &

let count+=1
[[ $((count%NR_CPUS)) -eq 0 ]] && wait

done

wait

#Create a concatenated outfile
for g in *.renamed.fasta.gz
do

zcat $g >> cat.fasta
rm $g

done

gzip cat.fasta
rm cat.fasta

echo 'Job is done.'
