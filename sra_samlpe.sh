#!/bin/bash

fileItemString=$(cat  "$PWD/sra_files/sra_file_names.txt" |tr "\n" " ")
sample=($fileItemString)
len=${#sample[@]}
len=$((len-1))
o1="$PWD/Read_files"
for i in $(seq 0 $len); do
    in1="$PWD/sra_files/${sample[i]}"
    j=$((i+1))
    echo "File $j processing"
    fastq-dump -I --split-files -O $o1 $in1 
done
echo ""
echo ""
echo ""
echo "This script was successfully run."
echo ""
echo ""

