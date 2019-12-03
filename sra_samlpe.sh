#!/bin/bash

fileItemString=$(cat  "sra_file_names.txt" |tr "\n" " ")
sample=($fileItemString)
#echo "${sample[0]}"
len=${#sample[@]}

#fileItemString=$(cat  "$PWD/files/catalog_align.txt" |tr "\n" " ")
#parm=($fileItemString)
#echo "${parm[2]}"
#len=${#sample[@]}
len=$((len-1))
#in1="$PWD/Read_files/"
#o1="$PWD/fastqc_output"
#if [ ${parm[0]} != "q" ]; then
#   echo "This step only execute for fastq files"
#fi
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

