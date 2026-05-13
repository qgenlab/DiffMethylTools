# Installation

Installing `DiffMethylTools` and its dependencies is quick and straightforward.

## 1. Install the Package

### Option 1 — Install from GitHub

Clone the repository:

```bash
git clone https://github.com/qgenlab/DiffMethylTools.git
cd DiffMethylTools
```
Ensure that Python and pip are installed, then install the package from the project root directory:

```bash
pip install .
```
### Option 2 — Install directly from PyPI

Alternativly, you can install the package using pip

```bash
pip install diffmethyltools
```
## 2. Setup Reference Genome
Before running any downstream analysis or annotations, you must initialize the reference files for your target genome build.

### Default Setup
By default, DiffMethylTools comes pre-packaged with the necessary reference files for both hg38 and hg19. To initialize the tool with the default files, simply run:

```bash
DiffMethylTools --setup hg38 # or hg19
```
Once this command completes, the tool is ready to use!

### Advanced Setup: Using Custom Reference Files

If you prefer to use your own custom reference files (e.g., a specific GTF version, custom CpG tracks, or custom array manifests), you can do so by editing the shared configuration file (hg38_min.yml).

In this YAML file, you must define the path, the column separator (sep), and for array manifests, the specific columns to extract (using 0-based indexing). If an optional file path is left empty or is incorrect, the script will safely generate dummy files to prevent pipeline failure.

Here is the template for the yml file:
```yml
# Mandatory Information
ref_name: "hg38_custom"

gtf:
  path: "../gencode.v42.chr_patch_hapl_scaff.annotation.gtf"
  sep: "\t"

CpG:
  path: "../CpG.bed"
  sep: "\t"

# Optional Files (Leave path empty "" if not using)
ccre:
  path: "../encodeCcreCombined.bed"
  sep: "\t"
  
rmsk:
  path: "../rmsk.txt"
  sep: "\t"

epic_v2:
  path: "../EPIC-8v2-0_A1.csv"
  sep: ","
  columns: [1, 15, 16]

epic_v1:
  path: "../infinium-methylationepic-v-1-0-b5-manifest-file-csv"
  sep: ","
  columns: [1, 48, 49, 50]

hm450:
  path: "../HM450.hg38.manifest.tsv"
  sep: "\t"
  columns: [0, 1, 2, 8]
```

#### Required File Formats
To ensure your custom files are parsed correctly, they must match the following structures:

* **Reference Name (`ref_name`) - *Mandatory***
    A unique identifier for your custom reference. The script will generate a directory with this exact name to store all processed files. You must use this name as the input for the `ref_folder` parameter in all subsequent downstream tasks.
* **GTF (`gtf`) - *Mandatory*** <br/>
    Standard Ensembl/GENCODE GTF format (tab-separated). Must contain standard feature columns (`chr`, `source`, `feature`, `start`, `end`, `score`, `strand`, `frame`, `attribute`).

* **CpG Base (`CpG`) - *Mandatory*** <br/>
    A simple 3-column, tab-separated BED-like file. Note: The index must start from 1. <br/>
    Format: `chr` | `pos` | `index` <br/>
    Example: `chr1 \t 10469 \t 1`

* **cCRE (`ccre`) - *Optional***  <br/>
    Standard ENCODE candidate cis-Regulatory Elements BED file (tab-separated). <br/>
    Example: `chr1 \t 181251 \t 181601 \t EH38E1310153 ...`

* **RepeatMasker (`rmsk`) - *Optional***  <br/>
    Tab-separated repeats file. Must contain basic coordinates and repeat categories. <br/>
    Format: `chr` | `start` | `end` | `strand` | `repeat_cat` | `repeat_subcat` ...

* **Illumina Array Manifests (epic_v2, epic_v1, hm450) - *Optional***  <br/>
    When providing Illumina manifest files (CSV or TSV), you must specify the 0-based column indices that   correspond to the Probe ID, Chromosome, and Position.

    - EPIC v2: columns: `[1, 15, 16]` extracts Name, CHR, and MAPINFO.

    - HM450: columns: `[0, 1, 2, 8]` extracts CpG_chrm, CpG_beg, CpG_end, and Probe_ID.

    - EPIC v1: columns: `[1, 48, 49, 50]` extracts the corresponding Probe ID and coordinate columns based on the v1 manifest structure.
