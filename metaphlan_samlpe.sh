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
in1="$PWD/Read_files/"
o1="$PWD/metaphlan_output/metaphlan_"
s="$PWD/Read_files/"
cd MetaPhlAn-2.7.8
export PATH='pwd':$PATH
export mpa_dir='pwd'
for i in $(seq 0 $len); do
   echo "File $((i+1)) processing"
   in2="$in1${sample[i]}"
   o2="$o1${sample[i]}"
   ./metaphlan2.py $in2 $o2 --input_type $t --tax_lev $1 --bt2_ps very-sensitive-local  --min_alignment_len 21 --nproc $2 --ignore_viruses --ignore_eukaryotes --ignore_archaea
   # --bt2_ps very-sensitive-local  --min_alignment_len 21

   s1="$s${sample[i]}.bowtie2out.txt"
   rm $s1
done

echo ""
