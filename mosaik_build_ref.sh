#!/bin/bash

fileItemString=$(cat  "$PWD/files/catalog_name.txt" |tr "\n" " ")
fileItemArray=($fileItemString)
ref="${fileItemArray[0]}"
d1="main_ref_${ref}.dat"
d2="ref_${ref}"
fileItemString=$(cat  "$PWD/files/catalog_hash.txt" |tr "\n" " ")
fileItemArray=($fileItemString)
ha1=${fileItemArray[0]}
ha2=${fileItemArray[1]}


./MosaikBuild -fr $ref -oa $d1
./MosaikJump -ia $d1 $ha1 $ha2 -out $d2 -kd -mem 4
echo ""
echo ""
echo ""
echo "This script was successfully run."
echo ""
echo ""

