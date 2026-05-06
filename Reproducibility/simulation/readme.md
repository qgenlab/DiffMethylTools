## Methylation Data Simulator
This script is a data simulation pipeline designed to generate synthetic "Case" and "Control" DNA methylation datasets based on real-world genomic inputs. It is highly useful for benchmarking Differentially Methylated Region (DMR) calling tools by creating a ground truth with known effect sizes.

### What It Does
**CpG Clustering**: Reads a given BED file of CpG sites, groups sites that are within 100bp of each other into regions, and filters out short regions.

**Data Intersection**: Maps real methylation data (BED12 format) to these clustered regions, retaining only regions with high data density (>20 points).

**Control Baseline Generation**: Applies LOWESS smoothing to the real methylation data to establish a clean, continuous "Control" profile.

**Case Simulation**: Generates synthetic "Case" profiles by randomly assigning a biological effect (an increase or decrease in methylation) and applying a Gaussian-smoothed offset to the Control baseline.

**Replicate Sampling**: Uses a truncated normal distribution to generate biological replicates (default: 3 cases, 3 controls) around the simulated profiles.

**Coverage Simulation**: Synthesizes realistic sequencing coverage (read depth) for the generated samples based on real-world distributions and variance.

**Statistical Summary**: Calculates evaluation metrics, including Cohen's d (effect size) and standard deviations across the simulated datasets.

### How to Use
Prerequisites: Ensure you have the required Python libraries installed: pandas, numpy, scipy, and statsmodels.

#### Execution:
Run the script from the command line by providing two arguments: a CpG positions BED file and a Methylation BED file.

```
python simulate_methylation.py <path_to_CpG_file.bed> <path_to_Methylation_file.bed>
```

### Example:

```
python simulate_methylation.py ./CpG_With_Strand.bed ./wgbs_methylation_data.bed
```

### Output Files
The script will generate several files in the directory where it is run:

1. `CpG_regions_long_hg38.bed`: The clustered intermediate regions.

2. `simulation_case_[1-3].csv`: The simulated methylation levels and coverage for the Case replicates.

3. `simulation_ctr_[1-3].csv`: The simulated methylation levels and coverage for the Control replicates.

4. `res_sim.csv`: A summary table containing the baseline values, simulated differences, standard deviations, and Cohen's d for each region.

5. `sim.pos.cohen.csv`: A simplified mapping of chromosome positions to the calculated Cohen's d effect sizes.
