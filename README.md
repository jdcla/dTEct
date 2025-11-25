# Differential Expression and Translation Efficiency Analysis Pipeline

## Overview

This R pipeline is designed for comprehensive analysis of RNA-sequencing (RNA-seq) and Ribosome Profiling (Ribo-seq) data. It automates the process of quality control, differential gene expression (DE), and differential translation efficiency (dTE) analysis. The script is built upon robust Bioconductor packages like `edgeR` and leverages `ggplot2` and `pheatmap` for generating publication-quality visualizations.

A key feature of this pipeline is its ability to automatically discover and evaluate complex comparisons (contrasts) based on a structured metadata file, making it highly adaptable to various experimental designs.

### Core Features

* **Integrated Analysis:** Simultaneously analyzes RNA-seq and Ribo-seq data to distinguish between transcriptional and translational regulation.
* **Dual-Mode Analysis:** Supports both **Gene-Level** (standard intersection) and **Transcript-Level** analysis (expanding RNA transcripts to match multiple Open Reading Frames/Translons).
* **Robust Variance Modeling:** Fits two separate statistical models (Paired and Independent) to ensure RNA-seq statistics are not penalized by the higher variance typically found in Ribo-seq data.
* **Flexible Design:** Supports single-factor and multi-factorial experimental designs.
* **Automated Contrast Generation:** Intelligently identifies relevant comparisons from the metadata, including hierarchical and combinatorial groups.
* **Comprehensive QC:** Generates a suite of quality control plots, including MDS, PCA, sample distance heatmaps, and gene clustering heatmaps.
* **Rich Visualization:** Produces volcano plots for every differential analysis to easily identify significant genes.
* **Parallel Processing:** Utilizes the `BiocParallel` package to speed up computationally intensive steps.

## System Requirements

This script is designed to run in an R environment. You will need the following packages from CRAN and Bioconductor.

### Installation

The easiest use of this script is by running the singularity container present on the turbo share, on any UM Linux system:
|||bash
singularity exec /nfs/turbo/umms-prensnerturbo/shared/workflows/singularities/DE_tools.sif dTEct.R <function flags>
|||

Alternatively, you can install all required packages by running the following commands in your R console:

|||R
# Install packages from CRAN
install.packages(c("optparse", "tidyverse", "pheatmap", "RColorBrewer", "ggplot2", "svglite"))

# Install packages from Bioconductor
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("edgeR", "DESeq2", "EnhancedVolcano", "DEFormats", "BiocParallel"))
|||

## Input Files

The pipeline requires up to four types of input files, specified via command-line arguments.

1.  **Metadata File (`--metadata`)**
    A Prensner Lab Sample sheet or CSV metadata sheet with the following requirements:

    * **Format:** CSV.
    * **Header:** Must contain a header row. Lines starting with `#` are ignored.
    * **Required Columns:** The script is flexible, but key columns used in the example are:
        * `smart_id`: A unique identifier for each biological sample.
        * `data_type`: The sequencing type, e.g., `Ribo_seq` or `RNA_seq`.
        * `replicate_num`: The replicate number for the sample.
        * `test_or_control`: A column indicating if a sample is a "test" or "control" group, used for ordering contrasts.
        * **Contrast Columns:** Any column that defines an experimental variable (e.g., `source_id`, `treatment_id`, `disease_id`). These are specified with the `--contrast_cols` argument.

    * **Optional Columns**
        * `batch_date`: Can be used in combination with `--no_batch_factor` to automatically control for batch effects.

2.  **RNA Counts File (`--rna_counts`)**
    A matrix of raw gene counts from an RNA-seq experiment.

    * **Format:** CSV/TSV. (see `--sep` argument)
    * **Rows:** Each row represents a feature (Gene or Transcript), with the ID in the first column.
    * **Columns:** Each column represents a sample. The column headers must match the `smart_id` values from the metadata file.

3.  **Ribo Counts File (`--ribo_counts`)**
    A matrix of raw gene counts from a Ribo-seq experiment.

    * **Format:** Identical to the RNA Counts File.

4.  **Transcript-to-Gene Map (`--tx_table_path`)**
    A lookup table to map feature IDs to readable gene names and link Translons to Transcripts. **Required for Transcript Mode.**

    * **Format:** CSV.
    * **Required Columns:**
        * `translon_id` (Required for `--feature_level transcript`)
        * `transcript_id` (Required for `--feature_level transcript`)
        * `gene_id`
        * `gene_name`

---

## Quick Start Guide

Follow these steps to run the analysis on the provided example data.

### 1. Prepare Your Workspace

Clone or download the repository. The example files should be located in the `example/` directory.

### 2. Run the Script

Open your terminal or command prompt, navigate to the script's directory, and execute the following command:

|||bash
Rscript dTEct.R \
    --metadata example/metadata.csv \
    --rna_counts example/rna_numreads.csv \
    --ribo_counts example/ribo_numreads.csv \
    --tx_table_path example/tx_table.csv \
    --feature_level "transcript" \
    --outdir example/out/ \
    --count_col 5 \
    --contrast_cols "source_id,treatment_id" \
    --cores 4
|||

### 3. Check the Output

Make sure to create an output directory first (e.g., `example_out/`). Within this directory, the script will generate the following structure:

* **`example_out/`**
    * `RNA/`: Results for differential expression on RNA-seq data.
        * `*_RNA_full.csv`: Results from the **Independent Model** (Recommended for transcriptional analysis).
        * `*_RNA.csv`: Results from the **Paired Model** (Shared RNAs, deduplicated).
    * `Ribo/`: Results for differential expression on Ribo-seq data (Paired Model).
    * `TE/`: Results for Translation Efficiency (TE) within a single group.
    * `dTE/`: Results for differential Translation Efficiency (dTE) between groups.
    * `QC/`: Directory containing MDS, PCA, Heatmaps, and Gene Clustering plots.
    * `run_info.txt`: A log file containing details about the run, including the design formula and contrast groups.

Each analysis subfolder (`RNA`, `Ribo`, `dTE`, etc.) will contain:
* `*.csv`: A table of differential expression results for a specific contrast.
* `*_Volcano.png`: A volcano plot visualizing the results from the corresponding CSV file.

---

## Advanced Usage

This section details the command-line options and the logic behind the automated contrast generation, allowing you to tailor the pipeline to your specific experimental design.

### Command-Line Options

| Option              | Argument      | Default            | Description                                                                                                                               |
| ------------------- | ------------- | ------------------ | ----------------------------------------------------------------------------------------------------------------------------------------- |
| `-m`, `--metadata`    | path          | `NULL`             | **Required.** Path to the project sample sheet (metadata CSV).                                                                            |
| `-i`, `--rna_counts`  | path          | `NULL`             | Path to the RNA-seq count matrix.                                                                                                         |
| `-j`, `--ribo_counts` | path          | `NULL`             | Path to the Ribo-seq count matrix.                                                                                                        |
| `-t`, `--tx_table_path`    | path          | `''`               | Path to the table linking transcript and gene IDs/names (Required for `--feature_level transcript`).                                      |
| `-v`, `--feature_level` | string        | `"gene"`           | Analysis mode: `"gene"` (Simple intersection) or `"transcript"` (Expands RNA rows to match Ribo Translons).                               |
| `-f`, `--tx_table_col`  | character     | `"gene_id"`        | The column name of `tx_table` to use as key (mostly used in `gene` mode).                                                                 |
| `-a`, `--contrast_cols`| character    | `"treatment_id"`   | Comma-separated metadata column names from which contrasts are derived.                                                                   |
| `-c`, `--count_col`    | integer       | `5`                | First column in the count files containing sample count data.                                                                             |
| `-d`, `--id_col`       | integer       | `1`                | Column in the count files containing feature identifiers.                                                                                 |
| `-s`, `--sep`          | character     | `,`                | Separator used in input tables.                                                                                                           |
| `-o`, `--outdir`       | path          | `out`              | Output directory for results.                                                                                                             |
| `-l`, `--cores`        | integer       | `1`                | Number of cores for parallel processing.                                                                                                  |
| `-p`, `--plot_ids`     | boolean       | `FALSE`            | If `TRUE`, labels PCA/MDS plots with Sample IDs. If `FALSE`, uses replicate numbers.                                                      |
| `-e`, `--no_batch_factor` | boolean | `FALSE` | If `TRUE`, prevents the script from automatically adding `batch_date` to the design model. |

### Designing Contrasts (`--contrast_cols`)

This is the most powerful feature of the pipeline. The `--contrast_cols` argument tells the script which columns in the metadata define the experimental groups to be compared.

1.  **Multiple Factors:** The script is designed to take into account multiple factors. By default, it will factor in an effect for Sequencing type.
    - **Multiple Columns:** You can provide multiple column names separated by commas (e.g., `"source_id,treatment_id"`). The script will create a new, combined factor.
    - **Combinatorial Groups:** The script can parse group names with combinatorial factors separated by `_and_`.
    - **Batch date:** Batch dates are automatically detected and added to the model unless `--no_batch_factor` is used.

2.  **Hierarchical Groups:** The script understands hierarchical relationships defined by a dot (`.`) in group names (e.g., `Tumor.TypeA`). It generates contrasts between subgroups and supergroups. **The code automatically weighs the contribution of subgroups** based on sample size when calculating supergroup averages.

### Analysis Types & Statistical Logic

The pipeline automatically performs several types of analyses:

* **Transcriptional Output (RNA-seq DE):**
    * **Model:** Uses an **Independent Model**. This estimates dispersion using *only* RNA samples to ensure RNA statistics are not penalized by the higher variance of Ribo-seq data.
    * **Output:** `RNA/*_RNA_full.csv` (Recommended) and `RNA/*_RNA.csv` (Paired/Shared model).

* **Translational Output (Ribo-seq DE):**
    * **Model:** Uses the **Paired Model**. Measures changes in total ribosome occupancy.
    * **Output:** `Ribo/`.

* **Differential Translation Efficiency (dTE):**
    * **Model:** Uses the **Paired Model** (Interaction term). Identifies genes where the *change* in ribosome loading differs significantly from the *change* in transcription.
    * **Output:** `dTE/`.

### Deduplication (Transcript Mode)
When running in `transcript` mode, multiple Translons (ORFs) can map to a single Transcript. This requires expanding the RNA matrix for dTE analysis.
* For RNA-only results (`_RNA.csv`), this creates duplicate rows.
* **Logic:** The pipeline automatically deduplicates these results by grouping by `transcript_id` and retaining only the result with the **lowest P-value** (representing the best-fit statistical scenario for that transcript).