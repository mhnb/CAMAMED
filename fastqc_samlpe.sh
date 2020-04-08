#!/bin/bash

fileItemString=$(cat  "$PWD/Read_files/sample_file_names.txt" |tr "\n" " ")
sample=($fileItemString)
len=${#sample[@]}

fileItemString=$(cat  "$PWD/files/catalog_align.txt" |tr "\n" " ")
parm=($fileItemString)
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

