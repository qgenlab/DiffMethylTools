

########################## Liver


Rscript run_methylSig.r --control "./liver_GSE70090/trimmed_SRR2074689/new_destranded_filtered_marked_output_forward_paired_bismark_bt2_pe.CpG_report.txt,./liver_GSE70090/trimmed_SRR2074685/new_destranded_filtered_marked_output_forward_paired_bismark_bt2_pe.CpG_report.txt,./liver_GSE70090/trimmed_SRR2074681/new_destranded_filtered_marked_output_forward_paired_bismark_bt2_pe.CpG_report.txt,./liver_GSE70090/trimmed_SRR2074677/new_destranded_filtered_marked_output_forward_paired_bismark_bt2_pe.CpG_report.txt" --case "./liver_GSE70090/trimmed_SRR2074687/new_destranded_filtered_marked_output_forward_paired_bismark_bt2_pe.CpG_report.txt,./liver_GSE70090/trimmed_SRR2074683/new_destranded_filtered_marked_output_forward_paired_bismark_bt2_pe.CpG_report.txt,./liver_GSE70090/trimmed_SRR2074679/new_destranded_filtered_marked_output_forward_paired_bismark_bt2_pe.CpG_report.txt,./liver_GSE70090/trimmed_SRR2074675/new_destranded_filtered_marked_output_forward_paired_bismark_bt2_pe.CpG_report.txt" --output ./bismark_liver/ --cov 5 --spcov1 2 --spcov2 2 -w FALSE


########################## WTA_KOB

Rscript run_methylSig.r --control "./final_R10_bc_5mc_5hmc_sup/analysis/filtered_bam/WTA/filtered/filtered_CpG/new_destranded_filtered_Methy5mC_WTA.bam_5mc.methyl1_filtered_1.CpG_bsseq_report.txt,./final_R9_bc_5mc_5hmc_sup/analysis/filtered_bam/WTA/filtered/filtered_CpG/new_destranded_filtered_Methy5mC_WTA.bam_5mc.methyl1_filtered_1.CpG_bsseq_report.txt" --case "./final_R10_bc_5mc_5hmc_sup/analysis/filtered_bam/KOB/filtered/filtered_CpG/new_destranded_filtered_Methy5mC_KOB.bam_5mc.methyl1_filtered_1.CpG_bsseq_report.txt,./final_R9_bc_5mc_5hmc_sup/analysis/filtered_bam/KOB/filtered/filtered_CpG/new_destranded_filtered_Methy5mC_KOB.bam_5mc.methyl1_filtered_1.CpG_bsseq_report.txt" --output ./KOB_WTA/ --cov 10 --spcov1 1 --spcov2 1 -w FALSE


########################## AD


p="."


Rscript run_methylSig.r --control "$p/AD_ctrlSample2/filtered/filtered_CpG/new_destranded_filtered_Methy5mC_AD_ctrlSample2.bam_5mc.new.methyl1_filtered_1.CpG_bsseq_report.txt,$p/AD_ctrlSample3/filtered/filtered_CpG/new_destranded_filtered_Methy5mC_AD_ctrlSample3.bam_5mc.new.methyl1_filtered_1.CpG_bsseq_report.txt,$p/AD_ctrlSample4/filtered/filtered_CpG/new_destranded_filtered_Methy5mC_AD_ctrlSample4.bam_5mc.new.methyl1_filtered_1.CpG_bsseq_report.txt,$p/AD_11_0562/filtered/filtered_CpG/new_destranded_filtered_Methy5mC_AD_11_0562.bam_5mc.new.methyl1_filtered_1.CpG_bsseq_report.txt,$p/AD_13_0582/filtered/filtered_CpG/new_destranded_filtered_Methy5mC_AD_13_0582.bam_5mc.new.methyl1_filtered_1.CpG_bsseq_report.txt,$p/AD_12_0575/filtered/filtered_CpG/new_destranded_filtered_Methy5mC_AD_12_0575.bam_5mc.new.methyl1_filtered_1.CpG_bsseq_report.txt,$p/AD_17_0604/filtered/filtered_CpG/new_destranded_filtered_Methy5mC_AD_17_0604.bam_5mc.new.methyl1_filtered_1.CpG_bsseq_report.txt,$p/AD_19_0608/filtered/filtered_CpG/new_destranded_filtered_Methy5mC_AD_19_0608.bam_5mc.new.methyl1_filtered_1.CpG_bsseq_report.txt" --case "$p/AD_ad0405/filtered/filtered_CpG/new_destranded_filtered_Methy5mC_AD_ad0405.bam_5mc.new.methyl1_filtered_1.CpG_bsseq_report.txt,$p/AD_ad0487/filtered/filtered_CpG/new_destranded_filtered_Methy5mC_AD_ad0487.bam_5mc.new.methyl1_filtered_1.CpG_bsseq_report.txt,$p/AD_7_0515/filtered/filtered_CpG/new_destranded_filtered_Methy5mC_AD_7_0515.bam_5mc.new.methyl1_filtered_1.CpG_bsseq_report.txt,$p/AD_8_0519/filtered/filtered_CpG/new_destranded_filtered_Methy5mC_AD_8_0519.bam_5mc.new.methyl1_filtered_1.CpG_bsseq_report.txt,$p/AD_9_0223/filtered/filtered_CpG/new_destranded_filtered_Methy5mC_AD_9_0223.bam_5mc.new.methyl1_filtered_1.CpG_bsseq_report.txt,$p/AD_10_0525/filtered/filtered_CpG/new_destranded_filtered_Methy5mC_AD_10_0525.bam_5mc.new.methyl1_filtered_1.CpG_bsseq_report.txt,$p/AD_15_0572/filtered/filtered_CpG/new_destranded_filtered_Methy5mC_AD_15_0572.bam_5mc.new.methyl1_filtered_1.CpG_bsseq_report.txt,$p/AD_16_0594/filtered/filtered_CpG/new_destranded_filtered_Methy5mC_AD_16_0594.bam_5mc.new.methyl1_filtered_1.CpG_bsseq_report.txt,$p/AD_18_ad0595/filtered/filtered_CpG/new_destranded_filtered_Methy5mC_AD_18_ad0595.bam_5mc.new.methyl1_filtered_1.CpG_bsseq_report.txt,$p/AD_20_0637/filtered/filtered_CpG/new_destranded_filtered_Methy5mC_AD_20_0637.bam_5mc.new.methyl1_filtered_1.CpG_bsseq_report.txt,$p/AD_21_0379/filtered/filtered_CpG/new_destranded_filtered_Methy5mC_AD_21_0379.bam_5mc.new.methyl1_filtered_1.CpG_bsseq_report.txt,$p/AD_adSample1/filtered/filtered_CpG/new_destranded_filtered_Methy5mC_AD_adSample1.bam_5mc.new.methyl1_filtered_1.CpG_bsseq_report.txt" --output ./AD/ --cov 10 --spcov1 6 --spcov2 8 -w FALSE


########################## MO_ND

p="."

Rscript run_methylSig.r --control "$p/MO1/filtered/filtered_CpG/new_destranded_MO1_5mc.new.methyl1_filtered_1.CpG_bsseq_report.txt,$p/MO2/filtered/filtered_CpG/new_destranded_MO2_5mc.new.methyl1_filtered_1.CpG_bsseq_report.txt,$p/MO3/filtered/filtered_CpG/new_destranded_MO3_5mc.new.methyl1_filtered_1.CpG_bsseq_report.txt" --case "$p/ND1/filtered/filtered_CpG/new_destranded_ND1_5mc.new.methyl1_filtered_1.CpG_bsseq_report.txt,$p/ND2/filtered/filtered_CpG/new_destranded_ND2_5mc.new.methyl1_filtered_1.CpG_bsseq_report.txt,$p/ND3/filtered/filtered_CpG/new_destranded_ND3_5mc.new.methyl1_filtered_1.CpG_bsseq_report.txt" --output ./MO_ND/ --cov 10 --spcov1 2 --spcov2 2

############################ NK B

case=`ls ./Blood-NK/GSM*/data_CpG_report.txt | paste -sd,`
ctr=`ls ./Blood-B/GSM*/data_CpG_report.txt | paste -sd,`


Rscript run_methylSig.r --control "$ctr" --case "$case" --output ./NK_B/ --cov1 5 --cov2 5 --spcov1 2 --spcov2 2
