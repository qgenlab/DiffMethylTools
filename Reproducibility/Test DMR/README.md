## DMR Tool Benchmarking & Evaluation Pipeline
This repository contains a suite of Bash scripts designed to evaluate the performance of various DMR (Differentially Methylated Region) calling tools.

Because there is rarely a definitive biological "ground truth" for experimental DMRs, these scripts use a consensus approach: a region is considered a "true" DMR if at least 3 distinct tools agree on it. The pipeline then evaluates each individual tool's performance against this consensus by calculating True Positives (TP), False Positives (FP), and False Negatives (FN), and finally counting the valid CpG sites within those regions.

### Pipeline Execution Order
The `run_all.sh` execution flow runs as follows:

1. `./test_3_res.sh`: Generates all 3-way tool intersections.
2. `./CpG_intersect.sh`: Filters intersections to keep only those containing CpGs.
3. `./rename_files_3_res.sh`: Standardizes file naming (mapping generate_DMR to DiffMethylTools).
4. `./generate_true_3_res.sh`: Builds the consensus "ground truth" for each dataset.
5. `./test_results_3_res.sh`: Computes TP, FP, and FN for each tool.
6. `./test_results_CpG.sh`: Counts CpG sites in the final TP/FP/FN files.

### Script Details
#### 1. test_3_res.sh (Data Preparation & 3-Way Intersection)
This script is the starting point of the pipeline. It defines the file paths for the outputs of 5 tools across 5 datasets (e.g., AD, bismark_liver, MO_ND, etc.).

Function: It creates directories named triple_0 through triple_4 (one for each dataset).

Processing: It uses nested loops to find every unique combination of 3 tools. It runs bedtools intersect to find genomic regions where 3 tools overlap, saving these 3-way consensus regions as .bed files.

#### 2. CpG_intersect.sh (CpG Filtering)
DMRs are only biologically meaningful if they actually contain CpG sites.

Function: Iterates through the generated .bed files in the triple_ directories.

Processing: Uses bedtools intersect against a reference CpG position file (CpG_With_Strand.bed) with the -c flag to count the number of CpGs in each region.

Filtering: Uses awk ($NF > 0) to drop any regions that contain 0 CpGs, outputting the filtered data to *.filter.CpG.not.null.bed.

#### 3. rename_files_3_res.sh (File Standardization)
A brief utility script to clean up file names for downstream processing.

Function: Identifies the tools used to create each 3-way intersection file based on file names.

Processing: Specifically replaces the generic term generate_DMR with the specific tool name DiffMethylTools so the naming convention matches the rest of the scripts.

#### 4. generate_true_3_res.sh (Ground Truth Generation)
This script defines the "truth" for evaluation. To evaluate a specific tool (e.g., methylkit), it needs a reference that doesn't artificially inflate that tool's score.

Function: Generates a specific [tool].test.bed ground-truth file for each of the 5 tools within a given triple_ directory.

Processing: For a target tool, it gathers all 3-way consensus files that do not include the target tool. It merges these arrays using cat and bedtools merge to create a robust, independent ground-truth baseline.

#### 5. test_results_3_res.sh (Performance Metrics Evaluation)
This script calculates the actual benchmarking metrics.

Function: Compares the original output of each tool against its corresponding .test.bed ground truth.

Processing: Strips stray quotes (sed -i 's/"//g'), and uses bedtools intersect to generate three output files per tool:

True Positives (TP): Regions called by the tool that overlap with the ground truth.

False Positives (FP): Regions called by the tool that do not overlap with the ground truth (-v flag).

False Negatives (FN): Regions in the ground truth that the tool missed (-v flag, reversing -a and -b).

Bonus: It uses awk to calculate and append the length of each region ($3 - $2) to the output.

#### 6. test_results_CpG.sh (Final Feature Counting)
Function: Analyzes the TP, FP, and FN files generated in the previous step.

Processing: Intersects these files with an unstranded CpG reference (CpG_Without_Strand.bed) using bedtools -c to count the exact number of CpGs falling into every true positive, false positive, and false negative region. Output is saved as *.count.CpG.bed.
