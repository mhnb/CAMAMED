#!/bin/bash


fileItemString=$(cat  "$PWD/Read_files/sample_file_names.txt" |tr "\n" " ")
sample=($fileItemString)
len=${#sample[@]}

fileItemString=$(cat  "$PWD/files/catalog_align.txt" |tr "\n" " ")
parm=($fileItemString)
if [ ${parm[0]} = "q" ]; then
   t="fastq"
fi
if [ ${parm[0]} = "a" ]; then
   t="fasta"
fi
len=$((len-1))

for i in $(seq 0 $len); do
   echo "File $((i+1)) processing"
   in1="$PWD/Read_files/${sample[i]}"
   o1="$PWD/metaphlan_output/metaphlan_${sample[i]}"
   metaphlan2 $in1 $o1 --input_type $t --tax_lev $1 --bt2_ps very-sensitive-local  --min_alignment_len 21 --nproc $2 --ignore_viruses --ignore_eukaryotes --ignore_archaea
   # --bt2_ps very-sensitive-local  --min_alignment_len 21

   s="$PWD/Read_files/${sample[i]}.bowtie2out.txt"
   rm $s
done

echo ""
