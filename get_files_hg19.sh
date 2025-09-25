echo "Creating hg19 repository..."
f="hg19"
mkdir -p $f


# for rmsk.txt file (download hg19 version done)
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/rmsk.txt.gz -P $f
gunzip $f/rmsk.txt.gz

# for gencode.v41.chr_patch_hapl_scaff.annotation.gtf file (download hg19 version)
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.chr_patch_hapl_scaff.annotation.gtf.gz -P $f
gunzip $f/gencode.v19.chr_patch_hapl_scaff.annotation.gtf.gz
mv $f/gencode.v19.chr_patch_hapl_scaff.annotation.gtf $f/gencode.chr_patch_hapl_scaff.annotation.gtf


# for encodeCcreCombined
wget https://hgdownload.soe.ucsc.edu/gbdb/hg38/encode3/ccre/encodeCcreCombined.bb -P $f
chmod +x bin/bigBedToBed
bin/bigBedToBed $f/encodeCcreCombined.bb $f/encodeCcreCombined.bed


# listover step for encodeCcreCombined
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz -P $f
chmod +x bin/liftOver
python bin/liftover.py $f/encodeCcreCombined.bed $f/hg38ToHg19.over.chain.gz $f/encodeCcreCombined.hg19
./bin/liftOver $f/encodeCcreCombined.bed $f/hg38ToHg19.over.chain.gz $f/encodeCcreCombined.hg19.bed $f/encodeCcreCombined.unmapped.bed -bedPlus=3 -tab
mv $f/encodeCcreCombined.hg19.bed $f/encodeCcreCombined.bed

# 
wget -P $f https://github.com/qgenlab/DiffMethylTools/releases/download/v0.1/CpG_gencode_annotation_2025July08.bed.gz
gzip -d $f/CpG_gencode_annotation_2025July08.bed.gz
./bin/liftOver $f/CpG_gencode_annotation_2025July08.bed $f/hg38ToHg19.over.chain.gz $f/CpG_gencode_annotation_2025July08.hg19.bed $f/CpG_gencode_annotation_2025July08.unmapped.bed -bedPlus=3 -tab
mv $f/CpG_gencode_annotation_2025July08.hg19.bed $f/CpG_gencode_annotation.bed
