#!/bin/bash

l1=("$HOME/methylSig/AD/new.methylSig..hg38.window.1000.filtered.csv" "$HOME/methylSig/bismark_liver/new.methylSig..hg38.window.1000.filtered.csv"  "$HOME/methylSig/MO_ND/new.methylSig..hg38.window.1000.filtered.csv" "$HOME/methylSig/NK_B/new.methylSig..hg38.window.1000.filtered.csv" "$HOME/methylSig/KOB_WTA/new.methylSig..hg38.window.1000.filtered.csv")



l2=("$HOME/bsseq/AD/new.dsseq..hg38.DMR.csv" "$HOME/bsseq/bismark_liver/new.dsseq..hg38.DMR.csv" "$HOME/bsseq/MO_ND/new.dsseq..hg38.DMR.csv" "$HOME/bsseq/NK_B/new.dsseq..hg38.DMR.csv" "$HOME/bsseq/KOB_WTA/new.dsseq..hg38.DMR.csv"  )



l3=("$HOME/dss/AD/new.dss.CpG.hg38DMR.csv" "$HOME/dss/bismark_liver/new.dss.CpG.hg38DMR.csv" "$HOME/dss/MO_ND/new.dss.CpG.hg38DMR.csv" "$HOME/dss/NK_B/new.dss.CpG.hg38DMR.csv" "$HOME/dss/KOB_WTA/new.dss.CpG.hg38DMR.csv")



l4=("$HOME/methylkit/AD/new.methylkit..destranded.CpG.hg38.window.1000.step.500.cov.10.filtered.csv" "$HOME/methylkit/bismark_liver/new.methylkit..destranded.CpG.hg38.window.1000.step.500.cov.5.filtered.csv"  "$HOME/methylkit/MO_ND/new.methylkit.CpG.hg38.window.1000.step.500.cov.10.filtered.csv" "$HOME/methylkit/NK_B/new.methylkit..destranded.CpG.hg38.window.1000.step.500.cov.5.filtered.csv" "$HOME/methylkit/KOB_WTA/new.methylkit..destranded.CpG.hg38.window.1000.step.500.cov.10.filtered.csv" )




l5=("$HOME/AD/data/generate_DMR.csv_0.csv" "$HOME/new_liver5/data/generate_DMR.csv_0.csv" "$HOME/new_MO_ND4/data/generate_DMR.csv_0.csv" "$HOME/new_NK_B4/data/generate_DMR.csv_0.csv" "$HOME/KOB_WTA/data/generate_DMR.csv_0.csv")




#for arr in l1[@] l2[@] l3[@] l4[@] l5[@]; do
#  for f in "${!arr}"; do
#	echo $f
#    python to_tab.py "$f"
#  done
#done



n=${#l1[@]}

for ((i = 0; i < n; i++)); do
  elements=("${l1[i]}" "${l2[i]}" "${l3[i]}" "${l4[i]}" "${l5[i]}")
  mkdir -p triple_$i
  for ((j = 0; j < 5; j++)); do
     for ((k = j + 1; k < 5; k++)); do
        for ((l = k + 1; l < 5; l++)); do
    	  echo $"${elements[j]%.csv}.bed" , $"${elements[k]%.csv}.bed", $"${elements[l]%.csv}.bed"
	  filename=$(basename "${elements[j]}")
	  name1="${filename%.*}"
	  filename=$(basename "${elements[k]}")
	  name2="${filename%.*}"
          filename=$(basename "${elements[l]}")
          name3="${filename%.*}"
	  echo "generating file -> $name1.$name2.$name3.bed"
	  bedtools intersect -a $"${elements[j]%.csv}.bed" -b $"${elements[k]%.csv}.bed" | bedtools intersect -a - -b $"${elements[l]%.csv}.bed"  > triple_$i/$name1.$name2.$name3.bed
	done
     done
  done
done
