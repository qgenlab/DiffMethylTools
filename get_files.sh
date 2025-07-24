# for rmsk.txt file
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/rmsk.txt.gz
gunzip rmsk.txt.gz

# for gencode.v41.chr_patch_hapl_scaff.annotation.gtf file
wget wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_41/gencode.v41.chr_patch_hapl_scaff.annotation.gtf.gz
gunzip gencode.v41.chr_patch_hapl_scaff.annotation.gtf.gz

# for encodeCcreCombined
wget https://hgdownload.soe.ucsc.edu/gbdb/hg38/encode3/ccre/encodeCcreCombined.bb
wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/bigBedToBed
chmod +x bigBedToBed
./bigBedToBed encodeCcreCombined.bb encodeCcreCombined.bed
rm ./bigBedToBed

# For easier use
mv CpG_gencode_annotation_2025July08.bed CpG_gencodev42ccrenb_repeat_epic1v2hm450.bed
