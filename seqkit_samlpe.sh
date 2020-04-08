#!/bin/bash

fileItemString=$(cat  "$PWD/Read_files/sample_file_names.txt" |tr "\n" " ")
sample=($fileItemString)
len=${#sample[@]}

len=$((len-1))
in1="$PWD/Read_files/"
o1="$PWD/seqkit_output/"
for i in $(seq 0 $len); do
      j=$((i+1))
      echo "File $j processing"
      in2="$in1${sample[i]}"
      o2="$o1${sample[i]}"
      ./seqkit stats $in2 -o $o2
done
