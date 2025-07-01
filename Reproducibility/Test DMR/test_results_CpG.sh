
input="$1"
output="${input%.bed}.count.CpG.bed"

bedtools intersect -b /mnt/analysis/derbelh/CpG_Without_Strand.bed -a $input -c > $output
