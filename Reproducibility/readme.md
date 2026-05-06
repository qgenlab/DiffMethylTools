## Reproducibility Folder Documentation
This folder contains all the necessary scripts, pipelines, and notebooks required to reproduce the differential methylation analysis, tool benchmarking, and synthetic data simulations presented in our study.

The folder is divided into specific subdirectories for each established Differentially Methylated Region (DMR) calling tool, alongside directories for our custom evaluation pipelines and simulation scripts.

### 1. DSS (/DSS)
This folder contains the scripts to format methylation data and execute the DSS (Dispersion Shrinkage for Sequencing data) tool.

#### Data Preparation
DSS requires a specific input format. Before running the main R script, BED files must be destranded and converted:

```
# 1. Destrand the input BED files
for f in *.bed; do echo $f; python new_destrand_bed.py $f; done

# 2. Convert to DSS-compatible format
for f in new_destranded_filtered*.bed; do echo $f; python bed2dss.py $f; done
```
#### Execution (run_dss.r)
This R script reads the formatted text files, constructs a BSseq object, and performs a Wald test using Beta-Binomial distributions (DMLtest). It then outputs both Differentially Methylated Loci (DMLs) and DMRs into CSV format.

#### Usage Example:

```
Rscript run_dss.r \
  --control "./MO1_dss_stranded.txt,./MO2_dss_stranded.txt,./MO3_dss_stranded.txt" \
  --case "./ND1_dss_stranded.txt,./ND2_dss_stranded.txt,./ND3_dss_stranded.txt" \
  --output ./dss/MO_ND/ \
  --destranded TRUE
```

### 2. methylKit (/MethylKit)
This folder contains run_methylkit.r for executing methylKit.

#### Execution
The script handles reading Bismark Cytosine Reports or generic BED files (methRead), uniting samples into a single object (unite), and calculating differential methylation (calculateDiffMeth). It supports both standard base-resolution analysis and tiled window-based analysis (e.g., 1000bp windows).

#### Usage Example:
```
Rscript run_methylkit.r \
  --control "./MO1_bismark.txt,./MO2_bismark.txt" \
  --case "./ND1_bismark.txt,./ND2_bismark.txt" \
  --output ./methylkit_out/ \
  --window TRUE --window_size 1000 --step 500 \
  --coverage 10 --min_group 2
```

### 3. methylSig (/MethylSig)
This folder contains run_methylSig.r for executing methylSig.
The script utilizes bsseq for initial data loading and applies stringent minimum coverage and group coverage filters (filter_loci_by_coverage, filter_loci_by_group_coverage). It computes differential methylation using `diff_methylsig()`, taking advantage of local correction and optionally tiling the genome into windows.

#### Usage Example:
```
Rscript run_methylSig.r \
  --control "./MO1_bsseq_report.txt,./MO2_bsseq_report.txt,./MO3_bsseq_report.txt" \
  --case "./ND1_bsseq_report.txt,./ND2_bsseq_report.txt,./ND3_bsseq_report.txt" \
  --output ./MO_ND/ \
  --cov 10 --spcov1 2 --spcov2 2
```

### 4. bsseq (/bsseq)
This folder contains `run_bsseq.r` for executing the standard bsseq pipeline.

#### Execution
Unlike the purely statistical models, bsseq applies a local smoothing algorithm (BSmooth) across the genome to estimate methylation boundaries. The script calculates local t-statistics (BSmooth.tstat) and identifies DMRs based on specified statistical cutoff thresholds (dmrFinder).

#### Usage Example:
```
Rscript run_bsseq.r \
  --control "./MO1_bsseq_report.txt,./MO2_bsseq_report.txt,./MO3_bsseq_report.txt" \
  --case "./ND1_bsseq_report.txt,./ND2_bsseq_report.txt,./ND3_bsseq_report.txt" \
  --output ./MO_ND/ \
  --cov1 10 --cov2 10 --spcov1 2 --spcov2 2
```
### 5. Benchmarking Pipeline (/Test DMR)
Because there is rarely a definitive biological "ground truth" for experimental data, these bash scripts establish a consensus-based ground truth to evaluate tool performance (TP, FP, FN metrics).

`test_3_res.sh`: Generates all 3-way tool intersections (a region is a "true" DMR if 3+ tools agree).

`CpG_intersect.sh`: Filters intersections to ensure they contain actual CpG sites.

`rename_files_3_res.sh`: Standardizes file naming conventions.

`generate_true_3_res.sh`: Builds the unbiased consensus ground truth per tool.

`test_results_3_res.sh`: Computes True Positives, False Positives, and False Negatives against the consensus.

`test_results_CpG.sh`: Counts exact CpG sites within the resulting TP/FP/FN regions.

### 6. Simulation Pipeline (/simulation)
This folder contains the data simulator used to create synthetic "Case" and "Control" methylation profiles with known effect sizes for precise benchmarking.

`simulate_methylation.py`: Reads real CpG distributions, clusters them, and applies LOWESS smoothing to create a control baseline. It then mathematically injects Gaussian-smoothed methylation changes to generate synthetic Cases, simulating realistic coverage noise and variance.

Outputs: Generates synthetic CSV profiles (`simulation_case.csv`, `simulation_ctr.csv`) alongside a master sheet (`res_sim.csv`) reporting standard deviations and Cohen's d effect sizes.

### 7. Analysis Notebooks (/notebooks)
This directory holds the final Jupyter notebooks used to aggregate the outputs from the empirical benchmarks and the simulations, calculating final performance metrics (Precision, Recall, F1-Scores) and generating the plots used in the manuscript.

`Simulation results analysis.ipynb`: Parses and plots the performance of all tools on the synthetic data generated by the simulation folder.

`Tools results analysis.ipynb`: Parses and plots the performance of all tools on the real-world biological datasets via the Test DMR consensus benchmarking framework.
