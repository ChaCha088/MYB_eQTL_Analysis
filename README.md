# Analysis Code for: [Cis–trans divergence in the regulatory architectures of MYBL2 and MYB underlies distinct transcriptional networks in immuno-metabolic traits]

This repository contains the source code used for the cis-eQTL, colocalization, and Mendelian Randomization (MR) analyses presented in the manuscript **"[Cis–trans divergence in the regulatory architectures of MYBL2 and MYB underlies distinct transcriptional networks in immuno-metabolic traits]"**.

## 📂 Repository Structure

The analysis pipeline consists of three main steps:

1.  **`scripts/01_gcta_eqtl.sh`**:
    * Performs cis-eQTL analysis using **GCTA-MLMA**.
    * Includes primary and conditional analyses for *MYB*, *MYBL2*, etc.
2.  **`scripts/02_coloc_susie.R`**:
    * Conducts fine-mapping and colocalization using **SuSiE** and **coloc**.
    * Requires summary statistics within a 1 Mb window.
3.  **`scripts/03_mr_analysis.R`**:
    * Performs Two-sample MR using `TwoSampleMR` and `MRPRESSO`.
    * Includes automated pipelines for local clumping, Steiger filtering, and visualization.

## 🛠️ Requirements & Dependencies

To reproduce the analysis, the following software and R packages are required:

* **System Tools:**
    * PLINK (v1.9)
    * GCTA (v1.94.1 or later)
* **R Packages (R >= 4.0.0):**
    * `TwoSampleMR`
    * `coloc`
    * `susieR`
    * `MRPRESSO`
    * `data.table`, `dplyr`, `ggplot2`

## 📂 Data Availability

Due to file size limitations, the raw GWAS summary statistics and Reference Panel files are not included in this repository.

To reproduce the analysis:

Download the GWAS summary statistics from the GWAS Catalog (Accession IDs provided in metadata/metadata_*_mr.csv).

Download the 1000 Genomes Phase 3 Reference Panel.

Place them in the data/ directory following the structure described in the metadata files.
