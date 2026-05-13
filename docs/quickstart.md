# Quick Start Tutorial

`DiffMethylTools` supports both **default methylation input formats** and **fully customizable formats**. This guide covers standard analyses using the default formats.

## 1. Core Analysis (Generate DML/DMR)

You can run the entire core analysis pipeline using the `all_analysis` command. 

### Option A: BED Format
If your data is in BED format (containing chromosome, position, coverage, and methylation percentage):

```bash
DiffMethylTools all_analysis \
  --case_data_file case1.bed case2.bed \
  --ctr_data_file ctr1.bed ctr2.bed \
  --input_format BED \
  --ref_folder hg38
```

### Option B: Bismark CpG Report (CR) Format

If your data is from Bismark (containing positive/negative methylation counts):

```bash
DiffMethylTools all_analysis \
  --case_data_file case1.txt case2.txt \
  --ctr_data_file ctr1.txt ctr2.txt \
  --input_format CR \
  --ref_folder hg38
```

## 2. Visualizing Results
Once your data is processed, you can easily generate comprehensive plots.

### Generate All Plots
To generate Volcano plots, Manhattan plots, and gene region mappings all at once:

```bash
DiffMethylTools all_plots \
  --data_file position_based.csv --data_has_header \
  --window_data_file generate_DMR_0.csv --window_data_has_header \
  --gene_file data/map_positions_to_genes_genes.csv --gene_has_header \
  --ccre_file data/map_positions_to_genes_CCRE.csv --ccre_has_header \
  --ref_folder hg38
```

### Annotation pie chart
Region-based DMR annotation pie chart:


```bash
python ../DiffMethylTools.py match_region_annotation \
  --regions_df_file generate_DMR_0.csv \
  --regions_df_has_header \
  --ref_folder (hg19 or hg38)
```

Annotation-based DMR annotation pie chart:


```bash
python ../DiffMethylTools.py match_region_annotation \
  --regions_df_file generate_DMR_0.csv \
  --regions_df_has_header \
  --annotation_or_region annotation \
  --ref_folder (hg19 or hg38)
```

### Plot Specific Methylation Curves

To plot DMR regions on a specific chromosome (e.g., chr1 between positions 3,664,000 and 3,668,000):

```bash
DiffMethylTools plot_methylation_curve \
  --region_data_file data/generate_DMR_0.csv --region_data_has_header \
  --position_data_file data/position_based.csv --position_data_has_header \
  --chr_filter chr1 \
  --start_filter 3664000 \
  --end_filter 3668000 \
  --ref_folder hg38
```

(Note: Omitting the --chr_filter, --start_filter, and --end_filter options will generate plots for all DMRs).
