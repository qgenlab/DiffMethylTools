# DiffMethylTools
# DiffMethylTools
DiffMethylTools is a Python-based toolkit for the comprehensive analysis of DNA methylation differences between two groups of samples. Designed for both short-read (e.g., WGBS, RRBS) and long-read (e.g., Nanopore) methylome data, DiffMethylTools enables accurate and streamlined detection of differentially methylated loci (DMLs) and regions (DMRs). The tool accepts flexible input formats including Bismark reports and generic BED-style methylation calls, making it compatible with most upstream methylation profiling workflows.

The package integrates statistical testing, biological annotation, and high-quality visualization into a single-command pipeline. Users can merge and filter candidate regions, map methylation changes to gene features and cis-regulatory elements (CCREs), and generate summary plots such as volcano plots, Manhattan plots, heatmaps, and gene-region profiles. DiffMethylTools also includes a module for generating annotation pie charts that quantify overlap of DMRs with genomic and functional features.

By combining flexibility with usability, DiffMethylTools provides researchers with a practical and efficient platform for epigenomic analysis. It supports high-throughput and reproducible workflows and is especially well-suited for studies investigating the role of DNA methylation in development, differentiation, and disease progression.

*This tutorial is still a work in progress.*

## Installation
### Install python dependencies
To create DiffMethylTools environment, run:
```
conda env create -f environment.yml
```
Then activate the environment with
```
conda activate DiffMethylTools
```

### Download annotation databases
#### For hg19
```
chmod +x get_files_hg19.sh
./get_files_hg19.sh
```
#### For hg38
```
chmod +x get_files_hg38.sh
./get_files_hg38.sh
```

## Usage


### Generate DML/DMR and map positions to genes
`DiffMethylTools` supports both **default methylation input formats** and **fully customizable formats**. Users can either specify a standard format via `--input_format` or manually define column indices.

### Default Input Formats
#### Input Format 1: BED Format with Methylation Percentage

Supported via: ```--input_format BED```

Expected BED-like format (example):
```
chr1    10468   10469   5mC  743   +   10468  10469  0,0,0   52   95.21
chr1    10470   10471   5mC  850   +   10470  10471  0,0,0   58   94.26
```

Run all analysis with:
```
python ../DiffMethylTools/DiffMethylTools.py all_analysis \
  --case_data_file case1.bed case2.bed \
  --ctr_data_file ctr1.bed ctr2.bed \
  --input_format BED \
  --ref_folder hg38 (or hg19)
```
For BED input, DiffMethylTools automatically interprets:
- Chromosome
- Position
- Coverage
- Methylation percentage

No manual column specification is required.

#### Input Format 2: Bismark CpG Report Format (CR)

Supported via: ```--input_format CR```

Expected Bismark CpG report format:

```
chr1    10468   +    12   4   CG   CGC
chr1    10470   +     8   4   CG   CGC
chr1    10483   +     7   2   CG   CGC
```


Run full analysis:

```
python ../DiffMethylTools/DiffMethylTools.py all_analysis \
  --case_data_file case1_CpG_report.txt case2_CpG_report.txt case3_CpG_report.txt \
  --ctr_data_file ctr1_CpG_report.txt ctr2_CpG_report.txt ctr3_CpG_report.txt \
  --input_format CR \
  --ref_folder hg38 (or hg19)
```
For CR input, DiffMethylTools automatically detects:
- Chromosome
- CpG position
- Methylated counts
- Unmethylated counts


#### Flexible / Custom Input Format
If the input files do not conform to standard BED or Bismark CpG report formats, users can manually specify column indices. For both case and control files, the user must define:
- Field separator
- Chromosome column index (0-based)
- Start position column index
- Either:
  - Methylation percentage + coverage column indices
  or
  - Methylated + unmethylated count column indices



##### Examples:
Custom format with CpG report file as input
```
python ../DiffMethylTools.py all_analysis \
  --case_data_file space_separated_case_file_paths \
  --ctr_data_file space_separated_ctr_file_paths \
  --case_data_chromosome_column_index 0 \
  --ctr_data_chromosome_column_index 0 \
  --case_data_position_start_column_index 1 \
  --ctr_data_position_start_column_index 1 \
  --case_data_positive_methylation_count_column_index 3 \
  --case_data_negative_methylation_count_column_index 4 \
  --ctr_data_positive_methylation_count_column_index 3 \
  --ctr_data_negative_methylation_count_column_index 4 \
  --case_data_separator 't' \
  --ctr_data_separator 't' \
  --ref_folder (hg19 or hg38)
```

Custom format with Bed file as input

```
python ../DiffMethylTools.py all_analysis \
  --case_data_file space_separated_case_file_paths \
  --ctr_data_file space_separated_ctr_file_paths \
  --case_data_separator 't' \
  --ctr_data_separator 't' \
  --case_data_chromosome_column_index 0 \
  --ctr_data_chromosome_column_index 0 \
  --case_data_position_start_column_index 1 \
  --ctr_data_position_start_column_index 1 \
  --case_data_methylation_percentage_column_index 10 \
  --case_data_coverage_column_index 9 \
  --ctr_data_methylation_percentage_column_index 10 \
  --ctr_data_coverage_column_index 9 \
  --ref_folder (hg19 or hg38)
```

Running *all_analysis* will generate two folders:

##### `plot/`

- This folder is initially empty.
- It will be populated by output figures and plots generated by `DiffMethylTools` (e.g., volcano plots, methylation curves, pie charts).

##### `data/`

This folder contains the required input and intermediate files used during the analysis. Below is a description of each file:

- `merge_tables.csv`: Merged and coverage-filtered methylation data for both case and control samples.
- `position_based.csv`: Methylation data summarized at the individual CpG or position level, including q-value and difference information.
- `filters.csv`: Filtered data based on methylation difference and statistical significance (q-value).
- `generate_DMR_0.csv`: Differentially methylated regions (DMRs), aggregated from position-level data.
- `generate_DMR_2.csv`: Differentially methylated loci (DMLs) located within identified DMRs.
- `generate_DMR_1.csv`: DMLs that do not fall within any DMR (isolated differential sites).
- `map_positions_to_genes_genes.csv`: Mapping of methylation regions to annotated gene features.
- `map_positions_to_genes_CCRE.csv`: Mapping of methylation regions to candidate cis-regulatory elements (CCREs).
- `map_win_2_pos.csv`: Maps DMR regions to all underlying CpG positions (includes both DML and non-DML positions, unlike `generate_DMR_1.csv`).
- `state.yaml`: YAML configuration file that tracks tool state, parameters, and progress.

### Plot Generation
Generate **volcano plot**, **Manhattan plot**, **upstream clustering**, and **region-based gene plots**:
```
python ../DiffMethylTools.py all_plots \
  --data_file data/position_based.csv \
  --data_has_header \
  --window_data_file data/generate_DMR_0.csv \
  --window_data_has_header \
  --data_separator ',' \
  --window_data_separator ',' \
  --gene_file data/map_positions_to_genes_genes.csv \
  --gene_has_header \
  --ccre_file data/map_positions_to_genes_CCRE.csv \
  --ccre_has_header \
  --ref_folder (hg19 or hg38)
```
To plot all DMR regions on chr1 between positions 3,664,000 and 3,668,000:
```
python ../DiffMethylTools.py plot_methylation_curve \
  --region_data_file data/generate_DMR_0.csv \
  --region_data_has_header \
  --position_data_file data/position_based.csv \
  --position_data_has_header \
  --chr_filter chr1 \
  --start_filter 3664000 \
  --end_filter 3668000 \
  --ref_folder (hg19 or hg38)
```
Omitting the --chr_filter, --start_filter, and --end_filter options will generate plots for all DMRs.

Region-based DMR annotation pie chart:
```
python ../DiffMethylTools.py match_region_annotation \
  --regions_df_file generate_DMR_0.csv \
  --regions_df_has_header \
  --ref_folder (hg19 or hg38)
```

Annotation-based DMR annotation pie chart:

```
python ../DiffMethylTools.py match_region_annotation \
  --regions_df_file generate_DMR_0.csv \
  --regions_df_has_header \
  --annotation_or_region annotation \
  --ref_folder (hg19 or hg38)
```

## Citing DiffMethylTools
If you used DiffMethylTools please cite:

Derbel, Houssemeddine, Evan Kinnear, Justin J-L. Wong, and Qian Liu. "DiffMethylTools: a toolbox of the detection, annotation and visualization of differential DNA methylation." bioRxiv (2025): 2025-07.

*(Manuscript currently under peer review)*
