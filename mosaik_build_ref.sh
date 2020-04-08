#!/bin/bash

fileItemString=$(cat  "$PWD/files/catalog_name.txt" |tr "\n" " ")
fileItemArray=($fileItemString)
ref="${fileItemArray[0]}"
d1="main_ref_${ref}.dat"
d2="ref_${ref}"

./MosaikBuild -fr $ref -oa $d1
./MosaikJump -ia $d1 -hs $1 -out $d2 -kd -mem 4
echo ""
echo ""
echo ""
echo "This script was successfully run."
echo ""
echo ""

