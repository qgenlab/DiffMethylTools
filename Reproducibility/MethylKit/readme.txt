### position based AD data

p="."


Rscript run_methylkit.r --control $p/AD_ctrlSample2/filtered/filtered_CpG/filtered_Methy5mC_AD_ctrlSample2.bam_5mc.new.methyl1_filtered_1.CpG_amp.txt,$p/AD_ctrlSample3/filtered/filtered_CpG/filtered_Methy5mC_AD_ctrlSample3.bam_5mc.new.methyl1_filtered_1.CpG_amp.txt,$p/AD_ctrlSample4/filtered/filtered_CpG/filtered_Methy5mC_AD_ctrlSample4.bam_5mc.new.methyl1_filtered_1.CpG_amp.txt,$p/AD_11_0562/filtered/filtered_CpG/filtered_Methy5mC_AD_11_0562.bam_5mc.new.methyl1_filtered_1.CpG_amp.txt,$p/AD_13_0582/filtered/filtered_CpG/filtered_Methy5mC_AD_13_0582.bam_5mc.new.methyl1_filtered_1.CpG_amp.txt,$p/AD_12_0575/filtered/filtered_CpG/filtered_Methy5mC_AD_12_0575.bam_5mc.new.methyl1_filtered_1.CpG_amp.txt,$p/AD_17_0604/filtered/filtered_CpG/filtered_Methy5mC_AD_17_0604.bam_5mc.new.methyl1_filtered_1.CpG_amp.txt,$p/AD_19_0608/filtered/filtered_CpG/filtered_Methy5mC_AD_19_0608.bam_5mc.new.methyl1_filtered_1.CpG_amp.txt --case $p/AD_ad0405/filtered/filtered_CpG/filtered_Methy5mC_AD_ad0405.bam_5mc.new.methyl1_filtered_1.CpG_amp.txt,$p/AD_ad0487/filtered/filtered_CpG/filtered_Methy5mC_AD_ad0487.bam_5mc.new.methyl1_filtered_1.CpG_amp.txt,$p/AD_7_0515/filtered/filtered_CpG/filtered_Methy5mC_AD_7_0515.bam_5mc.new.methyl1_filtered_1.CpG_amp.txt,$p/AD_8_0519/filtered/filtered_CpG/filtered_Methy5mC_AD_8_0519.bam_5mc.new.methyl1_filtered_1.CpG_amp.txt,$p/AD_9_0223/filtered/filtered_CpG/filtered_Methy5mC_AD_9_0223.bam_5mc.new.methyl1_filtered_1.CpG_amp.txt,$p/AD_10_0525/filtered/filtered_CpG/filtered_Methy5mC_AD_10_0525.bam_5mc.new.methyl1_filtered_1.CpG_amp.txt,$p/AD_15_0572/filtered/filtered_CpG/filtered_Methy5mC_AD_15_0572.bam_5mc.new.methyl1_filtered_1.CpG_amp.txt,$p/AD_16_0594/filtered/filtered_CpG/filtered_Methy5mC_AD_16_0594.bam_5mc.new.methyl1_filtered_1.CpG_amp.txt,$p/AD_18_ad0595/filtered/filtered_CpG/filtered_Methy5mC_AD_18_ad0595.bam_5mc.new.methyl1_filtered_1.CpG_amp.txt,$p/AD_20_0637/filtered/filtered_CpG/filtered_Methy5mC_AD_20_0637.bam_5mc.new.methyl1_filtered_1.CpG_amp.txt,$p/AD_21_0379/filtered/filtered_CpG/filtered_Methy5mC_AD_21_0379.bam_5mc.new.methyl1_filtered_1.CpG_amp.txt,$p/AD_adSample1/filtered/filtered_CpG/filtered_Methy5mC_AD_adSample1.bam_5mc.new.methyl1_filtered_1.CpG_amp.txt -m 6 --output $p/methylkit/AD -w FALSE




#### position based MO ND


p="."

Rscript run_methylkit.r --control $p/MO1/filtered/filtered_CpG/MO1_5mc.new.methyl1_filtered_1.CpG_amp.txt,$p/MO2/filtered/filtered_CpG/MO2_5mc.new.methyl1_filtered_1.CpG_amp.txt,$p/MO3/filtered/filtered_CpG/MO3_5mc.new.methyl1_filtered_1.CpG_amp.txt --case $p/ND1/filtered/filtered_CpG/ND1_5mc.new.methyl1_filtered_1.CpG_amp.txt,$p/ND2/filtered/filtered_CpG/ND2_5mc.new.methyl1_filtered_1.CpG_amp.txt,$p/ND3/filtered/filtered_CpG/ND3_5mc.new.methyl1_filtered_1.CpG_amp.txt -m 2 --output $p/methylkit/MO_ND -w FALSE --destrand TRUE --assembly hg38 --context CpG

## position based WTA KOB

p="."


Rscript run_methylkit.r --control $p/final_R10_bc_5mc_5hmc_sup/analysis/filtered_bam/WTA/filtered/filtered_Methy5mC_WTA.bam_5mc.methyl1_filtered_1_amp.txt,$p/final_R9_bc_5mc_5hmc_sup/analysis/filtered_bam/WTA/filtered/filtered_Methy5mC_WTA.bam_5mc.methyl1_filtered_1_amp.txt --case $p/final_R10_bc_5mc_5hmc_sup/analysis/filtered_bam/KOB/filtered/filtered_Methy5mC_KOB.bam_5mc.methyl1_filtered_1_amp.txt,$p/final_R9_bc_5mc_5hmc_sup/analysis/filtered_bam/KOB/filtered/filtered_Methy5mC_KOB.bam_5mc.methyl1_filtered_1_amp.txt -w FALSE --destrand TRUE --assembly hg38 --context CpG --output $p/methylkit/KOB_WTA


## position based liver cancer

p="."

Rscript run_methylkit.r --control $p/trimmed_SRR2074687/new_destranded_filtered_marked_output_forward_paired_bismark_bt2_pe.CpG_report.txt,$p/trimmed_SRR2074683/new_destranded_filtered_marked_output_forward_paired_bismark_bt2_pe.CpG_report.txt,$p/trimmed_SRR2074679/new_destranded_filtered_marked_output_forward_paired_bismark_bt2_pe.CpG_report.txt,$p/trimmed_SRR2074675/new_destranded_filtered_marked_output_forward_paired_bismark_bt2_pe.CpG_report.txt --case $p/trimmed_SRR2074689/new_destranded_filtered_marked_output_forward_paired_bismark_bt2_pe.CpG_report.txt,$p/trimmed_SRR2074685/new_destranded_filtered_marked_output_forward_paired_bismark_bt2_pe.CpG_report.txt,$p/trimmed_SRR2074681/new_destranded_filtered_marked_output_forward_paired_bismark_bt2_pe.CpG_report.txt,$p/trimmed_SRR2074677/new_destranded_filtered_marked_output_forward_paired_bismark_bt2_pe.CpG_report.txt -b 1 -c 5 -w FALSE --destrand TRUE --assembly hg38 --context CpG --output $p/methylkit/bismark_liver


# position b and NK

case=`ls ./Blood-NK/GSM*/data_CpG_report.txt | paste -sd,`
ctr=`ls ./Blood-B/GSM*/data_CpG_report.txt | paste -sd,`

Rscript run_methylkit.r --case "$case" --control "$ctr" -b 1 -c 5 -w FALSE --destrand TRUE --assembly hg38 --context CpG --output $p/methylkit/NK_B/



######################################################## NEW AS 10/06/2025




### window based AD data

p="."


Rscript run_methylkit.r --control $p/AD_ctrlSample2/filtered/filtered_CpG/filtered_Methy5mC_AD_ctrlSample2.bam_5mc.new.methyl1_filtered_1.CpG_amp.txt,$p/AD_ctrlSample3/filtered/filtered_CpG/filtered_Methy5mC_AD_ctrlSample3.bam_5mc.new.methyl1_filtered_1.CpG_amp.txt,$p/AD_ctrlSample4/filtered/filtered_CpG/filtered_Methy5mC_AD_ctrlSample4.bam_5mc.new.methyl1_filtered_1.CpG_amp.txt,$p/AD_11_0562/filtered/filtered_CpG/filtered_Methy5mC_AD_11_0562.bam_5mc.new.methyl1_filtered_1.CpG_amp.txt,$p/AD_13_0582/filtered/filtered_CpG/filtered_Methy5mC_AD_13_0582.bam_5mc.new.methyl1_filtered_1.CpG_amp.txt,$p/AD_12_0575/filtered/filtered_CpG/filtered_Methy5mC_AD_12_0575.bam_5mc.new.methyl1_filtered_1.CpG_amp.txt,$p/AD_17_0604/filtered/filtered_CpG/filtered_Methy5mC_AD_17_0604.bam_5mc.new.methyl1_filtered_1.CpG_amp.txt,$p/AD_19_0608/filtered/filtered_CpG/filtered_Methy5mC_AD_19_0608.bam_5mc.new.methyl1_filtered_1.CpG_amp.txt --case $p/AD_ad0405/filtered/filtered_CpG/filtered_Methy5mC_AD_ad0405.bam_5mc.new.methyl1_filtered_1.CpG_amp.txt,$p/AD_ad0487/filtered/filtered_CpG/filtered_Methy5mC_AD_ad0487.bam_5mc.new.methyl1_filtered_1.CpG_amp.txt,$p/AD_7_0515/filtered/filtered_CpG/filtered_Methy5mC_AD_7_0515.bam_5mc.new.methyl1_filtered_1.CpG_amp.txt,$p/AD_8_0519/filtered/filtered_CpG/filtered_Methy5mC_AD_8_0519.bam_5mc.new.methyl1_filtered_1.CpG_amp.txt,$p/AD_9_0223/filtered/filtered_CpG/filtered_Methy5mC_AD_9_0223.bam_5mc.new.methyl1_filtered_1.CpG_amp.txt,$p/AD_10_0525/filtered/filtered_CpG/filtered_Methy5mC_AD_10_0525.bam_5mc.new.methyl1_filtered_1.CpG_amp.txt,$p/AD_15_0572/filtered/filtered_CpG/filtered_Methy5mC_AD_15_0572.bam_5mc.new.methyl1_filtered_1.CpG_amp.txt,$p/AD_16_0594/filtered/filtered_CpG/filtered_Methy5mC_AD_16_0594.bam_5mc.new.methyl1_filtered_1.CpG_amp.txt,$p/AD_18_ad0595/filtered/filtered_CpG/filtered_Methy5mC_AD_18_ad0595.bam_5mc.new.methyl1_filtered_1.CpG_amp.txt,$p/AD_20_0637/filtered/filtered_CpG/filtered_Methy5mC_AD_20_0637.bam_5mc.new.methyl1_filtered_1.CpG_amp.txt,$p/AD_21_0379/filtered/filtered_CpG/filtered_Methy5mC_AD_21_0379.bam_5mc.new.methyl1_filtered_1.CpG_amp.txt,$p/AD_adSample1/filtered/filtered_CpG/filtered_Methy5mC_AD_adSample1.bam_5mc.new.methyl1_filtered_1.CpG_amp.txt -m 6 --output $p/methylkit/AD -w TRUE -c 10 -b 2 --destrand TRUE


#### window based MO ND


p="."

Rscript run_methylkit.r --control $p/MO1/filtered/filtered_CpG/MO1_5mc.new.methyl1_filtered_1.CpG_amp.txt,$p/MO2/filtered/filtered_CpG/MO2_5mc.new.methyl1_filtered_1.CpG_amp.txt,$p/MO3/filtered/filtered_CpG/MO3_5mc.new.methyl1_filtered_1.CpG_amp.txt --case $p/ND1/filtered/filtered_CpG/ND1_5mc.new.methyl1_filtered_1.CpG_amp.txt,$p/ND2/filtered/filtered_CpG/ND2_5mc.new.methyl1_filtered_1.CpG_amp.txt,$p/ND3/filtered/filtered_CpG/ND3_5mc.new.methyl1_filtered_1.CpG_amp.txt -m 2 --output $p/methylkit/MO_ND -w TRUE --destrand TRUE --assembly hg38 --context CpG -c 10 -b 2

## window based WTA KOB

p="."


Rscript run_methylkit.r --control $p/final_R10_bc_5mc_5hmc_sup/analysis/filtered_bam/WTA/filtered/filtered_Methy5mC_WTA.bam_5mc.methyl1_filtered_1_amp.txt,$p/final_R9_bc_5mc_5hmc_sup/analysis/filtered_bam/WTA/filtered/filtered_Methy5mC_WTA.bam_5mc.methyl1_filtered_1_amp.txt --case $p/final_R10_bc_5mc_5hmc_sup/analysis/filtered_bam/KOB/filtered/filtered_Methy5mC_KOB.bam_5mc.methyl1_filtered_1_amp.txt,$p/final_R9_bc_5mc_5hmc_sup/analysis/filtered_bam/KOB/filtered/filtered_Methy5mC_KOB.bam_5mc.methyl1_filtered_1_amp.txt -w TRUE --destrand TRUE --assembly hg38 --context CpG --output $p/methylkit/KOB_WTA -c 10 -b 2


## window based liver cancer

p="."

Rscript run_methylkit.r --control $p/trimmed_SRR2074687/new_destranded_filtered_marked_output_forward_paired_bismark_bt2_pe.CpG_report.txt,$p/trimmed_SRR2074683/new_destranded_filtered_marked_output_forward_paired_bismark_bt2_pe.CpG_report.txt,$p/trimmed_SRR2074679/new_destranded_filtered_marked_output_forward_paired_bismark_bt2_pe.CpG_report.txt,$p/trimmed_SRR2074675/new_destranded_filtered_marked_output_forward_paired_bismark_bt2_pe.CpG_report.txt --case $p/trimmed_SRR2074689/new_destranded_filtered_marked_output_forward_paired_bismark_bt2_pe.CpG_report.txt,$p/trimmed_SRR2074685/new_destranded_filtered_marked_output_forward_paired_bismark_bt2_pe.CpG_report.txt,$p/trimmed_SRR2074681/new_destranded_filtered_marked_output_forward_paired_bismark_bt2_pe.CpG_report.txt,$p/trimmed_SRR2074677/new_destranded_filtered_marked_output_forward_paired_bismark_bt2_pe.CpG_report.txt -b 1 -c 5 -w TRUE --destrand TRUE --assembly hg38 --context CpG --output $p/methylkit/bismark_liver


## window b and NK

case=`ls ./Blood-NK/GSM*/data_CpG_report.txt | paste -sd,`
ctr=`ls ./Blood-B/GSM*/data_CpG_report.txt | paste -sd,`

Rscript run_methylkit.r --case "$case" --control "$ctr" -b 1 -c 5 -w TRUE --destrand TRUE --assembly hg38 --context CpG --output $p/methylkit/NK_B/

