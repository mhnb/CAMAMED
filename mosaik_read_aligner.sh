#!/bin/bash

fileItemString=$(cat  "$PWD/files/catalog_name.txt" |tr "\n" " ")
fileItemArray=($fileItemString)
ref="${fileItemArray[0]}"
d1="main_ref_${ref}.dat"
d2="ref_${ref}"


fileItemString=$(cat  "sample_file_names.txt" |tr "\n" " ")
sample=($fileItemString)
#echo "${sample[0]}"
len=${#sample[@]}
len=$((len-1))

fileItemString=$(cat  "$PWD/files/catalog_align.txt" |tr "\n" " ")
parm=($fileItemString)

#aligne fastq files

if [ ${parm[0]} = "q" ]; then
   if [ ${parm[1]} = "s" ]; then
      for i in $(seq 0 $len); do
	    j=$((i+1))
            echo "File $j processing"
            in1="$PWD/Read_files/${sample[i]}"
            o1="$PWD/mosaik_outputs/s_${sample[i]}.mkb"
	    o2="$PWD/mosaik_outputs/ss_${sample[i]}.mka"
            ./MosaikBuild -q $in1 -st illumina -out $o1
            ./MosaikAligner -in $o1 -out $o2 -ia $d1 -j $d2 -p ${parm[2]} -annpe 2.1.26.pe.100.0065.ann -annse 2.1.26.se.100.005.ann
            o3="$PWD/mosaik_outputs/ss_${sample[i]}.mka.bam"
            o4="$PWD/mosaik_outputs/${sample[i]}.sam"
            samtools view -h -F 4 $o3 > $o4
            rm $o1
            rm $o2
            rm $o3
            o1="$PWD/mosaik_outputs/s_${sample[i]}.mkb.bam"
            rm $o1
            o1="$PWD/mosaik_outputs/ss_${sample[i]}.mka.stat" 
            rm $o1
      done
   fi
fi

if [ ${parm[0]} = "q" ]; then
   if [ ${parm[1]} = "p" ]; then
      for i in $(seq 0 2 $len); do
            j=$((i+1))
            #echo "File $i processing"
            in1="$PWD/Read_files/${sample[i]}"
            echo $in1
            in2="$PWD/Read_files/${sample[j]}"
            o1="$PWD/mosaik_outputs/s_${sample[i]}.mkb"
	    o2="$PWD/mosaik_outputs/ss_${sample[i]}.mka"
            ./MosaikBuild -q $in1 -q2 $in2 -st illumina -mfl ${parm[3]} -out $o1
            ./MosaikAligner -in $o1 -out $o2 -ia $d1 -j $d2 -p ${parm[2]} -annpe 2.1.26.pe.100.0065.ann -annse 2.1.26.se.100.005.ann
            o3="$PWD/mosaik_outputs/ss_${sample[i]}.mka.bam"
            o4="$PWD/mosaik_outputs/${sample[i]}.sam"
            samtools view -h -F 4 $o3 > $o4
            rm $o1
            rm $o2
            rm $o3
            o1="$PWD/mosaik_outputs/s_${sample[i]}.mkb.bam"
            rm $o1
            o1="$PWD/mosaik_outputs/ss_${sample[i]}.mka.stat"
            rm $o1
      done
   fi
fi

#aligne fasta files

if [ ${parm[0]} = "a" ]; then
   if [ ${parm[1]} = "s" ]; then
      for i in $(seq 0 $len); do
	    j=$((i+1))
            echo "File $j processing"
            in1="$PWD/Read_files/${sample[i]}"
            o1="$PWD/mosaik_outputs/s_${sample[i]}.mkb"
	    o2="$PWD/mosaik_outputs/ss_${sample[i]}.mka"
            ./MosaikBuild -fr $in1 -st illumina -out $o1 -assignQual 99
            ./MosaikAligner -in $o1 -out $o2 -ia $d1 -j $d2 -p ${parm[2]} -annpe 2.1.26.pe.100.0065.ann -annse 2.1.26.se.100.005.ann
            o3="$PWD/mosaik_outputs/ss_${sample[i]}.mka.bam"
            o4="$PWD/mosaik_outputs/${sample[i]}.sam"
            samtools view -h -F 4 $o3 > $o4
            rm $o1
            rm $o2
            rm $o3
            o1="$PWD/mosaik_outputs/s_${sample[i]}.mkb.bam"
            rm $o1
            o1="$PWD/mosaik_outputs/ss_${sample[i]}.mka.stat" 
            rm $o1
      done
   fi
fi

if [ ${parm[0]} = "a" ]; then
   if [ ${parm[1]} = "p" ]; then
      for i in $(seq 0 2 $len); do
            j=$((i+1))
            #echo "File $i processing"
            in1="$PWD/Read_files/${sample[i]}"
            echo $in1
            in2="$PWD/Read_files/${sample[j]}"
            o1="$PWD/mosaik_outputs/s_${sample[i]}.mkb"
	    o2="$PWD/mosaik_outputs/ss_${sample[i]}.mka"
            ./MosaikBuild -fr $in1 -fr2 $in2 -st illumina -mfl ${parm[3]} -out $o1 -assignQual 99
            ./MosaikAligner -in $o1 -out $o2 -ia $d1 -j $d2 -p ${parm[2]} -annpe 2.1.26.pe.100.0065.ann -annse 2.1.26.se.100.005.ann
            o3="$PWD/mosaik_outputs/ss_${sample[i]}.mka.bam"
            o4="$PWD/mosaik_outputs/${sample[i]}.sam"
            samtools view -h -F 4 $o3 > $o4
            rm $o1
            rm $o2
            rm $o3
            o1="$PWD/mosaik_outputs/s_${sample[i]}.mkb.bam"
            rm $o1
            o1="$PWD/mosaik_outputs/ss_${sample[i]}.mka.stat"
            rm $o1
      done
   fi
fi


echo ""
echo ""
echo ""
echo "This script was successfully run."
echo ""
echo ""

