# DiffMethylTools 

[![Documentation](https://img.shields.io/badge/docs-live-blue.svg)](https://qgenlab.github.io/DiffMethylTools/)
[![PyPI version](https://badge.fury.io/py/diffmethyltools.svg)](https://badge.fury.io/py/diffmethyltools)

**DiffMethylTools** is a Python-based toolkit for the comprehensive analysis of DNA methylation differences between two groups of samples. Designed for both short-read (e.g., WGBS, RRBS) and long-read (e.g., Nanopore) methylome data, it enables the accurate and streamlined detection of differentially methylated loci (DMLs) and regions (DMRs). 


## Key Features
* **Flexible Input Compatibility:** Accepts Bismark reports and generic BED-style methylation calls, making it compatible with most upstream profiling workflows.
* **Integrated, Single-Command Pipeline:** Seamlessly connects statistical testing, biological annotation, and high-quality visualization, reducing the need for custom scripting.
* **Biological Context & Annotation:** Easily map methylation changes to specific gene features and candidate cis-regulatory elements (cCREs).
* **Rich Visualization Module:** Automatically generate publication-ready volcano plots, Manhattan plots, heatmaps, gene-region profiles, and annotation pie charts quantifying DMR overlaps.


---

## Full Documentation
For advanced custom genome setups, deep-dives into the pipeline parameters, and full usage tutorials, please visit our complete documentation website:

👉 **[DiffMethylTools Documentation](https://qgenlab.github.io/DiffMethylTools/)**

---

## Installation & Setup

### 1. Install the Package

**Option A: Install via pip (Recommended)**

The easiest way to install DiffMethylTools is directly using pip:
```bash
pip install diffmethyltools
```

**Option B: Install from Source (GitHub)**

If you prefer to install the latest development version from source, clone the repository and install:

```Bash
git clone https://github.com/qgenlab/DiffMethylTools.git
cd DiffMethylTools
pip install .
```

### 2. Setup Reference Genome
Before running any downstream analysis or annotations, you must initialize the reference files for your target genome build (requires an active internet connection to pull default files from UCSC).

```Bash
DiffMethylTools --setup hg38  # or hg19
```
(To use custom reference genomes, GTFs, or array manifests, please refer to the [Advanced Setup Documentation](https://qgenlab.github.io/DiffMethylTools/installation/#advanced-setup-using-custom-reference-files)).

## Quick Start Tutorial
DiffMethylTools supports both default methylation input formats and fully customizable formats. This guide covers core analyses, understanding your output, and generating visualizations.

### 1. Core Analysis (Generate DML/DMR)
You can run the entire core analysis and annotation pipeline using the all_analysis command.

**Option A: BED Format**

Use this if your data contains chromosome, position, coverage, and methylation percentage.
Supported via: `--input_format BED`

```Bash
DiffMethylTools all_analysis \
  --case_data_file case1.bed case2.bed \
  --ctr_data_file ctr1.bed ctr2.bed \
  --input_format BED \
  --ref_folder hg38
```

**Option B: Bismark CpG Report (CR) Format**

Use this if your data is from Bismark and contains positive/negative methylation counts.
Supported via: --input_format CR

```Bash
DiffMethylTools all_analysis \
  --case_data_file case1.txt case2.txt \
  --ctr_data_file ctr1.txt ctr2.txt \
  --input_format CR \
  --ref_folder hg38
```
(Note: For custom column-separated data formats, see the [flexible input guide](https://qgenlab.github.io/DiffMethylTools/quickstart/#1-core-analysis-generate-dmldmr)).

### 2. Understanding the Output

Running `all_analysis` will generate two directories in your working folder:

- `plot/`: Initially empty; populated by output figures when you run plotting commands.

- `data/`: Contains results, intermediate files, and annotations. Key files include:

- `position_based.csv`: CpG/position-level stats, including q-values and differences.

- `generate_DMR_0.csv`: Differentially methylated regions (DMRs).

- `generate_DMR_1.csv`: Isolated DMLs not in any DMR.

- `generate_DMR_2.csv`: DMLs located within identified DMRs.

- `map_positions_to_genes_genes.csv`: Regions mapped to gene features.

- `map_positions_to_genes_CCRE.csv`: Regions mapped to cCREs.

### 3. Visualizing Results
Once your data is processed, use the generated data/ files to create comprehensive plots.

Generate All Standard Plots (Volcano, Manhattan, gene region mappings):

```Bash
DiffMethylTools all_plots \
  --data_file data/position_based.csv --data_has_header \
  --window_data_file data/generate_DMR_0.csv --window_data_has_header \
  --gene_file data/map_positions_to_genes_genes.csv --gene_has_header \
  --ccre_file data/map_positions_to_genes_CCRE.csv --ccre_has_header \
  --ref_folder hg38
```

Plot Specific Methylation Curves (e.g., chr1 between 3,664,000 and 3,668,000):

```Bash
DiffMethylTools plot_methylation_curve \
  --region_data_file data/generate_DMR_0.csv --region_data_has_header \
  --position_data_file data/position_based.csv --position_data_has_header \
  --chr_filter chr1 \
  --start_filter 3664000 \
  --end_filter 3668000 \
  --ref_folder hg38
```

(Omitting the filter options will generate plots for all DMRs).

Region-based Annotation Pie Charts:

```Bash
DiffMethylTools match_region_annotation \
  --regions_df_file data/generate_DMR_0.csv \
  --regions_df_has_header \
  --ref_folder hg38
```

Feature-based Annotation Pie Charts:

```Bash
DiffMethylTools match_region_annotation \
  --regions_df_file data/generate_DMR_0.csv \
  --regions_df_has_header \
  --annotation_or_region annotation \
  --ref_folder hg38
```

## Citation
If you use DiffMethylTools in your research, please cite our paper:

Derbel, Houssemeddine, Evan Kinnear, Justin J-L. Wong, and Qian Liu. "DiffMethylTools: a toolbox of the detection, annotation and visualization of differential DNA methylation." bioRxiv (2025): 2025-07.
(Manuscript currently under peer review)

## Issues and Support
If you encounter any bugs or have feature requests, please open an issue on our GitHub Issues page.
