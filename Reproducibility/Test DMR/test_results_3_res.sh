#!/bin/bash


tools=("methylSig." "dsseq." "dss." "methylkit." "DiffMethylTools.")


l1=("$HOME/methylSig/AD/new.methylSig..hg38.window.1000.filtered.bed" "$HOME/methylSig/bismark_liver/new.methylSig..hg38.window.1000.filtered.bed"  "$HOME/methylSig/MO_ND/new.methylSig..hg38.window.1000.filtered.bed" "$HOME/methylSig/NK_B/new.methylSig..hg38.window.1000.filtered.bed" "$HOME/methylSig/KOB_WTA/new.methylSig..hg38.window.1000.filtered.bed")


l2=("$HOME/bsseq/AD/new.dsseq..hg38.DMR.bed" "$HOME/bsseq/bismark_liver/new.dsseq..hg38.DMR.bed" "$HOME/bsseq/MO_ND/new.dsseq..hg38.DMR.bed" "$HOME/bsseq/NK_B/new.dsseq..hg38.DMR.bed" "$HOME/bsseq/KOB_WTA/new.dsseq..hg38.DMR.bed"  )


l3=("$HOME/dss/AD/new.dss.CpG.hg38DMR.bed" "$HOME/dss/bismark_liver/new.dss.CpG.hg38DMR.bed" "$HOME/dss/MO_ND/new.dss.CpG.hg38DMR.bed" "$HOME/dss/NK_B/new.dss.CpG.hg38DMR.bed" "$HOME/dss/KOB_WTA/new.dss.CpG.hg38DMR.bed")


l4=("$HOME/methylkit/AD/new.methylkit..destranded.CpG.hg38.window.1000.step.500.cov.10.filtered.bed" "$HOME/methylkit/bismark_liver/new.methylkit..destranded.CpG.hg38.window.1000.step.500.cov.5.filtered.bed"  "$HOME/methylkit/MO_ND/new.methylkit.CpG.hg38.window.1000.step.500.cov.10.filtered.bed" "$HOME/methylkit/NK_B/new.methylkit..destranded.CpG.hg38.window.1000.step.500.cov.5.filtered.bed" "$HOME/methylkit/KOB_WTA/new.methylkit..destranded.CpG.hg38.window.1000.step.500.cov.10.filtered.bed" )


l5=("$HOME/AD/data/generate_DMR.csv_0.bed" "$HOME/new_liver5/data/generate_DMR.csv_0.bed" "$HOME/new_MO_ND4/data/generate_DMR.csv_0.bed" "$HOME/new_NK_B4/data/generate_DMR.csv_0.bed" "$HOME/KOB_WTA/data/generate_DMR.csv_0.bed")



all_lists=(l1 l2 l3 l4 l5)



for t in "${!all_lists[@]}"; do
    tool="${tools[$t]}"
    current_list="${all_lists[$t]}"
    # files=("${!current_list}")
    eval "files=(\"\${${current_list}[@]}\")"
    for i in "${!files[@]}"; do
        file="${files[$i]}"
        dataset_folder="triple_${i}"
        reference="${dataset_folder}/${tool}.test.bed"
        output="${dataset_folder}/intersect_${tool}_${i}"
        mkdir -p "$dataset_folder"
        echo "Running bedtools intersect on $file with $reference -> $output"
	sed -i 's/"//g' "$file"
        # bedtools intersect -a "$file" -b "$reference" > "$output"
	bedtools intersect -a "$file" -b "$reference" | awk 'BEGIN{OFS="\t"} {print $0, $3 - $2}' > "$output.TP.bed"
	bedtools intersect -a "$file" -b "$reference" -v | awk 'BEGIN{OFS="\t"} {print $0, $3 - $2}'> "$output.FP.bed"
	bedtools intersect -b "$file" -a "$reference" -v | awk 'BEGIN{OFS="\t"} {print $0, $3 - $2}'> "$output.FN.bed"
    done
done
