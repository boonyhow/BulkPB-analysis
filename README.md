# Bulk & PseudoBulk RNA-seq Analysis Pipeline

This repository contains the full set of code used to analyze transcriptomic differences between matched bulk RNA-seq and pseudo-bulk (PB) RNA-seq data. The pipeline is modular and flexible — while it is designed for comprehensive end-to-end analysis, users are encouraged to extract and adapt individual components for their own scientific workflows.

---

## Getting Started

The pipeline is structured into four main stages:

### **Step 0: Environment Setup**

Create and configure the Conda environment and install all required R packages.

```bash
bash 0_env_setup.sh
```

This script will:

- Create a Conda environment based on `env_config_files/environment.yml`
- Install required R packages listed in `env_config_files/r_packages_with_source.csv` using the `code_modules/2_install_r_packages.R` function

---

### **Step 1: Data Procurement & Preprocessing**

Run this script to download all required datasets and preprocess them.

```bash
bash 1_data_modules.sh
```

This script executes all Python scripts under `code_modules/` in the correct order:

- Downloads datasets and GENCODE annotations
- Preprocesses individual datasets according to current standards

The following folder structure will be automatically generated under the `data/` directory:

```
data/
├── GSEXXXXXX/
│   ├── raw_tars/                 # Raw compressed files downloaded
│   ├── extracted/                # Extracted contents from raw files
│   └── data_for_study/
│       ├── bulk_rna_seq/        # Cleaned, averaged bulk expression matrices
│       ├── single_cell/         # Combined scRNA-seq matrix, features, barcodes
│       └── sample_subtype.csv   # Sample-subtype mapping
```

---

### **Step 2: Data Harmonization**

This script loads, normalizes, and integrates the datasets using R.

```bash
Rscript 2_data_preprocessing.r [config_file]
```

Current default config file is 2_config.yaml file, where the `exp_id` is a parameter set to the experiments downloaded earlier. Ensure that the necessary changes are made in the config file, and environment has been activated prior to running.

Output of this step is saved into: `data/GSEXXXXX/data_for_study/intermediate_data/`

This folder stores filtered, harmonized, and normalized expression matrices and metadata used in downstream analysis.

---

### **Step 3: Main Analysis Script**

This script runs the main analysis and visualization modules. All plots generated will be stored under the `plots` folder

You can either:

**Option A: Run as a full R script**

```bash
Rscript 3_main_file.r
```

**Option B: Run interactively in R / RStudio**

Open `3_main_file.r` and run each section manually to explore results and figures line by line.

This is recommended for users who want to investigate intermediate outputs, adjust parameters, or inspect visualizations in detail.

---

## Adding New Datasets

To integrate additional datasets into this pipeline:

1. Navigate to the `code_modules/` folder.
2. Create a new Python script named in the format:  
   ```
   0_2_X_GSEXXXXX_preprocessing.py
   ```
   where `X` is the next available number and `GSEXXXXX` is the dataset accession.
3. Follow the existing preprocessing structure:
   - Download raw files
   - Extract metadata
   - Convert to bulk/pseudo-bulk formats
   - Output to appropriate folders

Once added, remember to modify `1_data_modules.sh` to include the new script.

---

## Project Structure

```
├── 0_env_setup.sh                # Conda + R package setup
├── 1_data_modules.sh             # Python data download and preprocessing
├── 2_data_preprocessing.r        # R-based data loading and integration
├── 3_main_file.r                 # Full analysis pipeline
├── code_modules/                 # All dataset-specific processing code
├── env_config_files/             # environment.yml and R package CSV
├── 2_config.yaml                 # Custom config file for experiments
├── data/                         # Automatically generated data directory
├── plots/                        # Automatically generated results directory
```

---

## Requirements

- Conda (24.11)
- R (≥ 4.4.3; managed via Conda)
- Internet access to download packages and datasets

---

## Authors

B.H. Low, M.M. Rashid, K. Selvarajoo  
(Refer to manuscript for full author contributions.)

---
