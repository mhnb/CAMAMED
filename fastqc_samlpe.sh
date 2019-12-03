#!/bin/bash

fileItemString=$(cat  "sample_file_names.txt" |tr "\n" " ")
sample=($fileItemString)
#echo "${sample[0]}"
len=${#sample[@]}

fileItemString=$(cat  "$PWD/files/catalog_align.txt" |tr "\n" " ")
parm=($fileItemString)
#echo "${parm[2]}"
#len=${#sample[@]}
len=$((len-1))
in1="$PWD/Read_files/"
o1="$PWD/fastqc_output"
if [ ${parm[0]} != "q" ]; then
   echo "This step only execute for fastq files"
fi
cd FastQC
if [ ${parm[0]} = "q" ]; then
   for i in $(seq 0 $len); do
      j=$((i+1))
      echo "File $j processing"
      in2="$in1${sample[i]}"
      ./fastqc -o $o1 -f fastq $in2
   done
fi

