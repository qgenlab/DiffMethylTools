# Installation

Installing `DiffMethylTools` and its dependencies is quick and straightforward.

### 1. Install the Package
Ensure you have Python installed, then run the following command from the root directory of the project:

```bash
pip install .
```

### 2. Setup Reference Genome
Before running any downstream analysis or annotations, you must download the required reference files for your target genome build (e.g., hg19 or hg38).

Note: This command requires an active internet connection to pull the files from UCSC.

```bash
DiffMethylTools --setup hg38 #(or hg19)
```
