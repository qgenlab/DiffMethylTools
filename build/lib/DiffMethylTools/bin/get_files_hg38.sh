echo "Creating hg38 repository..."
f="hg38"
mkdir -p $f

# for rmsk.txt file
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/rmsk.txt.gz -P $f
gunzip $f/rmsk.txt.gz

# for gencode.v41.chr_patch_hapl_scaff.annotation.gtf file
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_41/gencode.v41.chr_patch_hapl_scaff.annotation.gtf.gz -P $f
gunzip $f/gencode.v41.chr_patch_hapl_scaff.annotation.gtf.gz
mv $f/gencode.v41.chr_patch_hapl_scaff.annotation.gtf $f/gencode.chr_patch_hapl_scaff.annotation.gtf

# for encodeCcreCombined
wget https://hgdownload.soe.ucsc.edu/gbdb/hg38/encode3/ccre/encodeCcreCombined.bb -P $f
chmod +x bin/bigBedToBed
bin/bigBedToBed $f/encodeCcreCombined.bb $f/encodeCcreCombined.bed
rm $f/encodeCcreCombined.bb

# For easier use
wget -P $f https://github.com/qgenlab/DiffMethylTools/releases/download/v0.1/CpG_gencode_annotation_2025July08.bed.gz
gzip -d $f/CpG_gencode_annotation_2025July08.bed.gz
mv $f/CpG_gencode_annotation_2025July08.bed $f/CpG_gencode_annotation.bed
