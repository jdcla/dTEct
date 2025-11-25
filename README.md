# Differential Expression and Translation Efficiency Analysis Pipeline

## Overview

This R pipeline is designed for the rigorous statistical analysis of RNA-sequencing (RNA-seq) and Ribosome Profiling (Ribo-seq) data. It automates the process of quality control, differential gene expression (DE), and differential translation efficiency (dTE) analysis.

The script is built upon `edgeR`'s generalized linear models (GLMs) and solves specific statistical challenges inherent to integrating multi-omics data, such as handling the different variance profiles of Ribo-seq and RNA-seq data.

### Core Features

* **Dual-Mode Analysis:**
    * **Gene Level:** Standard intersection of RNA and Ribo counts.
    * **Transcript/Translon Level:** Implements "Expansion Logic" to map single RNA transcripts to multiple Open Reading Frames (ORFs/Translons), allowing specific analysis of uORFs and isoforms.
* **Dual-Model Statistics:** Automatically fits two parallel GLMs:
    * **Paired Model:** For **dTE** analysis (accounts for RNA-Ribo interaction).
    * **Independent Model:** For **RNA-seq** analysis (estimates variance using only RNA samples for strict transcriptional statistics).
* **Automated Contrasts:** Intelligently identifies relevant comparisons from metadata, including hierarchical (Subgroup vs Supergroup) and combinatorial groups.
* **Smart Deduplication:** When running in Transcript mode, the pipeline automatically deduplicates RNA results by selecting the best-fit statistical result per transcript to avoid redundancy.
* **Rich Visualization:** Generates hierarchical QC plots and publication-ready Volcano plots.

## System Requirements

This script is designed to run in an R environment.

### Installation

The easiest use of this script is by running the singularity container:
```bash
singularity exec /nfs/turbo/umms-prensnerturbo/shared/workflows/singularities/DE_tools.sif dTEct.R <function flags>
```

Alternatively, you can install all required packages in R:

```R
install.packages(c("optparse", "tidyverse", "pheatmap", "RColorBrewer", "ggplot2", "svglite"))

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("edgeR", "DESeq2", "EnhancedVolcano", "DEFormats", "BiocParallel"))
```

## Input Files

1.  **Metadata File (`--metadata`)**
    * **Format:** CSV.
    * **Required Columns:** `smart_id`, `data_type` (`RNA_seq` or `Ribo_seq`), `replicate_num`.
    * **Contrast Columns:** Columns defining experimental groups (e.g., `treatment`, `disease`).

2.  **RNA Counts File (`--rna_counts`)**
    * **Format:** CSV matrix (Rows: Features, Cols: Samples).

3.  **Ribo Counts File (`--ribo_counts`)**
    * **Format:** CSV matrix (Rows: Features, Cols: Samples).

4.  **Transcript-to-Gene Map (`--tx_table_path`)**
    * **Required for Transcript Mode.**
    * **Format:** CSV.
    * **Columns:** Must contain `translon_id`, `transcript_id`, `gene_id`, `gene_name`.

---

## Quick Start Guide

```bash
Rscript dTEct.R \
    --metadata example/metadata.csv \
    --rna_counts example/rna_numreads.csv \
    --ribo_counts example/ribo_numreads.csv \
    --tx_table_path example/tx_table.csv \
    --feature_level "gene" \
    --contrast_cols "source_id,treatment_id" \
    --outdir example/out/ \
    --cores 4
```

---

## Advanced Usage

### Command-Line Options

| Option | Argument | Default | Description |
| :--- | :--- | :--- | :--- |
| `-m`, `--metadata` | path | `NULL` | **Required.** Path to metadata CSV. |
| `-i`, `--rna_counts` | path | `NULL` | Path to RNA-seq count matrix. |
| `-j`, `--ribo_counts` | path | `NULL` | Path to Ribo-seq count matrix. |
| `-t`, `--tx_table_path` | path | `''` | Path to table linking IDs (Required for `--feature_level transcript`). |
| `-v`, `--feature_level` | string | `gene` | Analysis Mode: `gene` (Simple intersection) or `transcript` (Expansion logic). |
| `-a`, `--contrast_cols` | string | `treatment_id` | Comma-separated metadata columns to group samples by. |
| `-f`, `--tx_table_col` | string | `gene_id` | Column name in `tx_table` to use as key (primary use in `gene` mode). |
| `-p`, `--plot_ids` | bool | `FALSE` | If TRUE, labels PCA/MDS with Sample IDs instead of replicate numbers. |
| `-o`, `--outdir` | path | `out` | Output directory. |
| `-l`, `--cores` | int | `1` | Number of cores for parallel processing. |
| `-e`, `--no_batch_factor` | bool | `FALSE` | If TRUE, prevents automatic inclusion of `batch_date` in the model. |

### Analysis Types & Logic

The pipeline automatically performs three distinct types of analysis:

1.  **Transcriptional Output (RNA-seq DE):**
    * **Model:** Uses the **Independent Model** (`fit_rna`). This model estimates dispersion using *only* RNA samples. This ensures that RNA statistics are not penalized by the typically higher variance of Ribo-seq data.
    * **Output:** Results are saved in `RNA/` with the suffix `_RNA_full`.

2.  **Translational Output (Ribo-seq DE):**
    * **Model:** Uses the **Paired Model** (`fit_paired`). Calculates the total change in ribosome occupancy.
    * **Output:** Saved in `Ribo/` with the suffix `_Ribo`.

3.  **Differential Translation Efficiency (dTE):**
    * **Model:** Uses the **Paired Model** (`fit_paired`). Identifies genes where the change in ribosome loading is significantly different from the change in mRNA abundance (Interaction term).
    * **Output:** Saved in `dTE/`.

### Deduplication Logic (Transcript Mode)
When running in `transcript` mode, multiple Translons (ORFs) may map to a single RNA transcript.
* **Problem:** This creates duplicate rows in the RNA analysis, where the counts are identical but P-values vary due to the influence of different Ribo-seq profiles on the row's variance.
* **Solution:** The pipeline automatically post-processes the paired RNA results. It groups by `transcript_id` and retains only the entry with the **lowest P-value**. This represents the "best-case" scenario where Ribo-seq noise had the least negative impact on the statistic.
