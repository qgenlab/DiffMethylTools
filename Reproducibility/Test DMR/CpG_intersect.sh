#!/bin/bash

for folder in "$@"; do
files=`ls "$folder"/*.bed`


# Replace ./CpG_With_Strand.bed with CpG positions file

for file in $files; do
    echo $file
    bedtools intersect -a $file -b ./CpG_With_Strand.bed -wa -c > $"${file%.bed}.filter.CpG.bed"
    awk '$NF > 0' $"${file%.bed}.filter.CpG.bed" > $"${file%.bed}.filter.CpG.not.null.bed"
done
done
