files=`ls /*.bed`

for f in $files; do echo $f; python new_destrand_bed.py $f; done


files=`ls ./new_destranded_filtered*.bam_5mc.new.methyl1_filtered_1.CpG.bed`
for f in $files; do echo $f; python bed2dss.py $f; done


# for AD data
###############################################################################
p="."


Rscript run_dss.r --control $p/AD_ctrlSample2/filtered/filtered_CpG/new_destranded_filtered_Methy5mC_AD_ctrlSample2.bam_5mc.new.methyl1_filtered_1.CpG_dss_stranded.txt,$p/AD_ctrlSample3/filtered/filtered_CpG/new_destranded_filtered_Methy5mC_AD_ctrlSample3.bam_5mc.new.methyl1_filtered_1.CpG_dss_stranded.txt,$p/AD_ctrlSample4/filtered/filtered_CpG/new_destranded_filtered_Methy5mC_AD_ctrlSample4.bam_5mc.new.methyl1_filtered_1.CpG_dss_stranded.txt,$p/AD_11_0562/filtered/filtered_CpG/new_destranded_filtered_Methy5mC_AD_11_0562.bam_5mc.new.methyl1_filtered_1.CpG_dss_stranded.txt,$p/AD_13_0582/filtered/filtered_CpG/new_destranded_filtered_Methy5mC_AD_13_0582.bam_5mc.new.methyl1_filtered_1.CpG_dss_stranded.txt,$p/AD_12_0575/filtered/filtered_CpG/new_destranded_filtered_Methy5mC_AD_12_0575.bam_5mc.new.methyl1_filtered_1.CpG_dss_stranded.txt,$p/AD_17_0604/filtered/filtered_CpG/new_destranded_filtered_Methy5mC_AD_17_0604.bam_5mc.new.methyl1_filtered_1.CpG_dss_stranded.txt,$p/AD_19_0608/filtered/filtered_CpG/new_destranded_filtered_Methy5mC_AD_19_0608.bam_5mc.new.methyl1_filtered_1.CpG_dss_stranded.txt --case $p/AD_ad0405/filtered/filtered_CpG/new_destranded_filtered_Methy5mC_AD_ad0405.bam_5mc.new.methyl1_filtered_1.CpG_dss_stranded.txt,$p/AD_ad0487/filtered/filtered_CpG/new_destranded_filtered_Methy5mC_AD_ad0487.bam_5mc.new.methyl1_filtered_1.CpG_dss_stranded.txt,$p/AD_7_0515/filtered/filtered_CpG/new_destranded_filtered_Methy5mC_AD_7_0515.bam_5mc.new.methyl1_filtered_1.CpG_dss_stranded.txt,$p/AD_8_0519/filtered/filtered_CpG/new_destranded_filtered_Methy5mC_AD_8_0519.bam_5mc.new.methyl1_filtered_1.CpG_dss_stranded.txt,$p/AD_9_0223/filtered/filtered_CpG/new_destranded_filtered_Methy5mC_AD_9_0223.bam_5mc.new.methyl1_filtered_1.CpG_dss_stranded.txt,$p/AD_10_0525/filtered/filtered_CpG/new_destranded_filtered_Methy5mC_AD_10_0525.bam_5mc.new.methyl1_filtered_1.CpG_dss_stranded.txt,$p/AD_15_0572/filtered/filtered_CpG/new_destranded_filtered_Methy5mC_AD_15_0572.bam_5mc.new.methyl1_filtered_1.CpG_dss_stranded.txt,$p/AD_16_0594/filtered/filtered_CpG/new_destranded_filtered_Methy5mC_AD_16_0594.bam_5mc.new.methyl1_filtered_1.CpG_dss_stranded.txt,$p/AD_18_ad0595/filtered/filtered_CpG/new_destranded_filtered_Methy5mC_AD_18_ad0595.bam_5mc.new.methyl1_filtered_1.CpG_dss_stranded.txt,$p/AD_20_0637/filtered/filtered_CpG/new_destranded_filtered_Methy5mC_AD_20_0637.bam_5mc.new.methyl1_filtered_1.CpG_dss_stranded.txt,$p/AD_21_0379/filtered/filtered_CpG/new_destranded_filtered_Methy5mC_AD_21_0379.bam_5mc.new.methyl1_filtered_1.CpG_dss_stranded.txt,$p/AD_adSample1/filtered/filtered_CpG/new_destranded_filtered_Methy5mC_AD_adSample1.bam_5mc.new.methyl1_filtered_1.CpG_dss_stranded.txt --output $p/dss/AD





###############################################################################
p="."

Rscript run_dss.r --control "$p/AD_ctrlSample2/filtered/filtered_CpG/filtered_Methy5mC_AD_ctrlSample2.bam_5mc.new.methyl1_filtered_1.CpG_dss_stranded.txt,$p/AD_ctrlSample3/filtered/filtered_CpG/filtered_Methy5mC_AD_ctrlSample3.bam_5mc.new.methyl1_filtered_1.CpG_dss_stranded.txt,$p/AD_ctrlSample4/filtered/filtered_CpG/filtered_Methy5mC_AD_ctrlSample4.bam_5mc.new.methyl1_filtered_1.CpG_dss_stranded.txt,$p/AD_11_0562/filtered/filtered_CpG/filtered_Methy5mC_AD_11_0562.bam_5mc.new.methyl1_filtered_1.CpG_dss_stranded.txt,$p/AD_13_0582/filtered/filtered_CpG/filtered_Methy5mC_AD_13_0582.bam_5mc.new.methyl1_filtered_1.CpG_dss_stranded.txt,$p/AD_12_0575/filtered/filtered_CpG/filtered_Methy5mC_AD_12_0575.bam_5mc.new.methyl1_filtered_1.CpG_dss_stranded.txt,$p/AD_17_0604/filtered/filtered_CpG/filtered_Methy5mC_AD_17_0604.bam_5mc.new.methyl1_filtered_1.CpG_dss_stranded.txt,$p/AD_19_0608/filtered/filtered_CpG/filtered_Methy5mC_AD_19_0608.bam_5mc.new.methyl1_filtered_1.CpG_dss_stranded.txt" --case "$p/AD_ad0405/filtered/filtered_CpG/filtered_Methy5mC_AD_ad0405.bam_5mc.new.methyl1_filtered_1.CpG_dss_stranded.txt,$p/AD_ad0487/filtered/filtered_CpG/filtered_Methy5mC_AD_ad0487.bam_5mc.new.methyl1_filtered_1.CpG_dss_stranded.txt,$p/AD_7_0515/filtered/filtered_CpG/filtered_Methy5mC_AD_7_0515.bam_5mc.new.methyl1_filtered_1.CpG_dss_stranded.txt,$p/AD_8_0519/filtered/filtered_CpG/filtered_Methy5mC_AD_8_0519.bam_5mc.new.methyl1_filtered_1.CpG_dss_stranded.txt,$p/AD_9_0223/filtered/filtered_CpG/filtered_Methy5mC_AD_9_0223.bam_5mc.new.methyl1_filtered_1.CpG_dss_stranded.txt,$p/AD_10_0525/filtered/filtered_CpG/filtered_Methy5mC_AD_10_0525.bam_5mc.new.methyl1_filtered_1.CpG_dss_stranded.txt,$p/AD_15_0572/filtered/filtered_CpG/filtered_Methy5mC_AD_15_0572.bam_5mc.new.methyl1_filtered_1.CpG_dss_stranded.txt,$p/AD_16_0594/filtered/filtered_CpG/filtered_Methy5mC_AD_16_0594.bam_5mc.new.methyl1_filtered_1.CpG_dss_stranded.txt,$p/AD_18_ad0595/filtered/filtered_CpG/filtered_Methy5mC_AD_18_ad0595.bam_5mc.new.methyl1_filtered_1.CpG_dss_stranded.txt,$p/AD_20_0637/filtered/filtered_CpG/filtered_Methy5mC_AD_20_0637.bam_5mc.new.methyl1_filtered_1.CpG_dss_stranded.txt,$p/AD_21_0379/filtered/filtered_CpG/filtered_Methy5mC_AD_21_0379.bam_5mc.new.methyl1_filtered_1.CpG_dss_stranded.txt,$p/AD_adSample1/filtered/filtered_CpG/filtered_Methy5mC_AD_adSample1.bam_5mc.new.methyl1_filtered_1.CpG_dss_stranded.txt" --output $p/dss/AD -d FALSE


###############################################################################


# for bismark_liver data

p="."

Rscript run_dss.r --control "$p/trimmed_SRR2074687/new_destranded_filtered_marked_output_forward_paired_bismark_bt2_pe.CpG_report_dss_bismark.txt,$p/trimmed_SRR2074683/new_destranded_filtered_marked_output_forward_paired_bismark_bt2_pe.CpG_report_dss_bismark.txt,$p/trimmed_SRR2074679/new_destranded_filtered_marked_output_forward_paired_bismark_bt2_pe.CpG_report_dss_bismark.txt,$p/trimmed_SRR2074675/new_destranded_filtered_marked_output_forward_paired_bismark_bt2_pe.CpG_report_dss_bismark.txt"  --case "$p/trimmed_SRR2074689/new_destranded_filtered_marked_output_forward_paired_bismark_bt2_pe.CpG_report_dss_bismark.txt,$p/trimmed_SRR2074685/new_destranded_filtered_marked_output_forward_paired_bismark_bt2_pe.CpG_report_dss_bismark.txt,$p/trimmed_SRR2074681/new_destranded_filtered_marked_output_forward_paired_bismark_bt2_pe.CpG_report_dss_bismark.txt,$p/trimmed_SRR2074677/new_destranded_filtered_marked_output_forward_paired_bismark_bt2_pe.CpG_report_dss_bismark.txt"  --output $p/dss/bismark_liver/ -d TRUE


# for MO_ND data (new)

Rscript run_dss.r --control "$p/MO1/filtered/filtered_CpG/new_destranded_MO1_5mc.new.methyl1_filtered_1.CpG_dss_stranded.txt,$p/MO2/filtered/filtered_CpG/new_destranded_MO2_5mc.new.methyl1_filtered_1.CpG_dss_stranded.txt,$p/MO3/filtered/filtered_CpG/new_destranded_MO3_5mc.new.methyl1_filtered_1.CpG_dss_stranded.txt"  --case "$p/ND1/filtered/filtered_CpG/new_destranded_ND1_5mc.new.methyl1_filtered_1.CpG_dss_stranded.txt,$p/ND2/filtered/filtered_CpG/new_destranded_ND2_5mc.new.methyl1_filtered_1.CpG_dss_stranded.txt,$p/ND3/filtered/filtered_CpG/new_destranded_ND3_5mc.new.methyl1_filtered_1.CpG_dss_stranded.txt"  --output $p/dss/MO_ND/ -d TRUE

# for bismark data

Rscript run_dss.r  --control $p/trimmed_SRR2074687/new_destranded_filtered_marked_output_forward_paired_bismark_bt2_pe.CpG_report_dss_bismark.txt,$p/trimmed_SRR2074683/new_destranded_filtered_marked_output_forward_paired_bismark_bt2_pe.CpG_report_dss_bismark.txt,$p/trimmed_SRR2074679/new_destranded_filtered_marked_output_forward_paired_bismark_bt2_pe.CpG_report_dss_bismark.txt,$p/trimmed_SRR2074675/new_destranded_filtered_marked_output_forward_paired_bismark_bt2_pe.CpG_report_dss_bismark.txt --case $p/trimmed_SRR2074689/new_destranded_filtered_marked_output_forward_paired_bismark_bt2_pe.CpG_report_dss_bismark.txt,$p/trimmed_SRR2074685/new_destranded_filtered_marked_output_forward_paired_bismark_bt2_pe.CpG_report_dss_bismark.txt,$p/trimmed_SRR2074681/new_destranded_filtered_marked_output_forward_paired_bismark_bt2_pe.CpG_report_dss_bismark.txt,$p/trimmed_SRR2074677/new_destranded_filtered_marked_output_forward_paired_bismark_bt2_pe.CpG_report_dss_bismark.txt --output 0/ --destranded TRUE


p="."

Rscript run_dss.r --case "$p/Blood-B-Mem/GSM5652319/data_CpG_report.txt,$p/Blood-B-Mem/GSM5652320/data_CpG_report.txt" --control "$pBlood-B/GSM5652316/data_CpG_report.txt,$p/Blood-B/GSM5652317/data_CpG_report.txt,$p/Blood-B/GSM5652318/data_CpG_report.txt" --destrand TRUE --output $p/dss/blood_b/
