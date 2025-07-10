# Differential Expression and Translation Efficiency Analysis Pipeline

## Overview

This R pipeline is designed for comprehensive analysis of RNA-sequencing (RNA-seq) and Ribosome Profiling (Ribo-seq) data. It automates the process of quality control, differential gene expression (DE), and differential translation efficiency (dTE) analysis. The script is built upon robust Bioconductor packages like `edgeR` and leverages `ggplot2` and `pheatmap` for generating publication-quality visualizations.

A key feature of this pipeline is its ability to automatically discover and evaluate complex comparisons (contrasts) based on a structured metadata file, making it highly adaptable to various experimental designs.

### Core Features

* **Integrated Analysis:** Simultaneously analyzes RNA-seq and Ribo-seq data to distinguish between transcriptional and translational regulation.
* **Flexible Design:** Supports single-factor and multi-factorial experimental designs.
* **Automated Contrast Generation:** Intelligently identifies relevant comparisons from the metadata, including hierarchical and combinatorial groups.
* **Comprehensive QC:** Generates a suite of quality control plots, including MDS, PCA, sample distance heatmaps, and gene clustering heatmaps.
* **Rich Visualization:** Produces volcano plots for every differential analysis to easily identify significant genes.
* **Parallel Processing:** Utilizes the `BiocParallel` package to speed up computationally intensive steps.

## System Requirements

This script is designed to run in an R environment. You will need the following packages from CRAN and Bioconductor.

### Installation

The easiest use of this script is by running the singularity container present on the turbo share, on any UM Linux system:
```bash
singularity exec /nfs/turbo/umms-prensnerturbo/shared/workflows/singularities/DE_tools.sif dTE.R <function flags>
```

Alternatively, you can install all required packages by running the following commands in your R console:

```R
# Install packages from CRAN
install.packages(c("optparse", "tidyverse", "pheatmap", "RColorBrewer", "ggplot2"))

# Install packages from Bioconductor
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("edgeR", "DESeq2", "EnhancedVolcano", "DEFormats", "BiocParallel"))
```

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
        * `batch_date`: Can be used in combination with `--batch_factor` to take into account batch dates when applying sample normalization. **CAN CAUSE THE SCRIPT TO FAIL** due to the resulting matrix not being full rank and therefore unsolvable.

2.  **RNA Counts File (`--rna_counts`)**
    A matrix of raw gene counts from an RNA-seq experiment.

    * **Format:** CSV/TSV. (see `--sep` argument)
    * **Rows:** Each row represents a gene, with the gene ID in the first column.
    * **Columns:** Each column represents a sample. The column headers must match the `smart_id` values from the metadata file.

3.  **Ribo Counts File (`--ribo_counts`)**
    A matrix of raw gene counts from a Ribo-seq experiment.

    * **Format:** Identical to the RNA Counts File. The same gene IDs should be used for direct comparison.


4.  **OPTIONAL: Transcript-to-Gene Map (`--tx_table`)**
    A lookup table to map feature IDs (e.g., gene IDs or transcript IDs) to more readable gene names for plotting.

    * **Format:** CSV.
    * **Columns:** Must contain a column with the feature ID used in the count matrices (e.g., `gene_id`) and a column with the desired name (e.g., `gene_name`). The column containing the ID should be named according to the feature level, e.g., `gene_id` or `transcript_id`.

---

## Quick Start Guide

Follow these steps to run the analysis on the provided example data.

### 1. Prepare Your Workspace

Clone or download the repository. The example files should be located in the `example/` directory.

### 2. Run the Script

Open your terminal or command prompt, navigate to the script's directory, and execute the following command:

```bash
Rscript dTE.R \
    --metadata example/metadata.csv \
    --rna_counts example/rna_numreads.csv \
    --ribo_counts example/ribo_numreads.csv \
    --tx_table example/tx_table.csv \
    --outdir example_out/ \
    --contrast_cols "source_id,treatment_id" \
    --cores 4
```

### 3. Check the Output

make sure to create an output directory first (e.g., `example_out/`). Within this directory, the script will generate the following structure:

* **`example_out/`**
    * `RNA/`: Results for differential expression on RNA-seq data.
    * `Ribo/`: Results for differential expression on Ribo-seq data.
    * `TE/`: Results for Translation Efficiency (TE) within a single group.
    * `dTE/`: Results for differential Translation Efficiency (dTE) between groups.
    * `run_info.txt`: A log file containing details about the run, including the design formula and contrast groups.
    * `formula.txt`: The final formula used for the GLM model.
    * `*_MDS_*.png`, `*_PCA_*.png`, `*_Heatmap_*.png`, `*_GeneClusters_*.png`: Various QC plots for RNA, Ribo, and combined data.

Each analysis subfolder (`RNA`, `Ribo`, `dTE`, etc.) will contain:
* `*.csv`: A table of differential expression results for a specific contrast.
* `*_Volcano.png`: A volcano plot visualizing the results from the corresponding CSV file.

---

## Advanced Usage

This section details the command-line options and the logic behind the automated contrast generation, allowing you to tailor the pipeline to your specific experimental design.

### Command-Line Options

| Option              | Argument      | Default            | Description                                                                                             |
| ------------------- | ------------- | ------------------ | ------------------------------------------------------------------------------------------------------- |
| `-m`, `--metadata`    | path          | `NULL`             | **Required.** Path to the project sample sheet (metadata CSV).                                          |
| `-i`, `--rna_counts`  | path          | `NULL`             | Path to the RNA-seq count matrix.                                                                       |
| `-j`, `--ribo_counts` | path          | `NULL`             | Path to the Ribo-seq count matrix.                                                                      |
| `-t`, `--tx_table`    | path          | `''`               | Path to the table linking transcript and gene IDs/names.                                                |
| `-a`, `--contrast_cols`| character   | `"treatment_id"`   | Comma-separated metadata column names from which contrasts are derived.                                 |
| `-c`, `--count_col`   | integer       | `5`                | First column in the count files containing sample count data. *(Note: script logic seems to override this).* |
| `-d`, `--id_col`      | integer       | `1`                | Column in the count files containing feature identifiers. *(Note: script logic seems to override this).* |
| `-s`, `--sep`         | character     | `,`                | Separator used in input tables.                                                                         |
| `-o`, `--outdir`      | path          | `out`              | Output directory for results.                                                                           |
| `-f`, `--tag`         | character     | `"gene"`     | The feature level of the counts (e.g., "gene", "transcript"). Used to find the correct ID in `tx_table`.  |
| `-l`, `--cores`       | integer       | `1`                | Number of cores for parallel processing.                                                                |
| `-e`, `--batch_factor` | boolean | `FALSE` | If `TRUE`, include `batch_date` as a factor in the design model. |

### Designing Contrasts (`--contrast_cols`)

This is the most powerful feature of the pipeline. The `--contrast_cols` argument tells the script which columns in the metadata define the experimental groups to be compared.

1.  **Multiple Factors:** The script is desinged to take into account multiple factors. By default, it will factor in an effect for Sequencing type.
    - **Multiple Columns:** You can provide multiple column names separated by commas (e.g., `"source_id,treatment_id"`). The script will create a new, combined factor by pasting these columns together (e.g., `CD8T__ACU_AGST`).
    - **Combinatorial Groups:** The script can parse group names with combinatorial factors separated by `_and_`. For instance, it can evaluate a group of samples with combined factors `DrugA_and_Time24h` with `DrugB_and_Time24h` and find the effect of the individual effects as well as the combined effects.
    - **Batch date:** Batch dates can be taken in as a factor as well, see `--batch_factor`. However, based on the number of batch dates across individual samples, this might make it impossible to factor our individual effects due to correlated factors (for which the code will return an error).  

2.  **Hierarchical Groups:** The script understands hierarchical relationships defined by a dot (`.`) in group names. For example, if you have groups `Tumor.TypeA` and `Tumor.TypeB`, the script recognizes they both belong to the `Tumor` supergroup and will generate a contrast between them. It will not compare them to an unrelated group like `Normal.Tissue`. By default, the script will evaluate contrasts between all supergroups if present. **The code will automatically weigh the contribution of the subgroups to the number of samples within this group**. For example, in a setting with 2 Sonic Hedgehog Medulloblastoma (MBL.SHH) and 3 Group 4 Medulloblastoma (MBL) samples, the effect of the supergroup medulloblastoma is ($\mu_{MBL} = 0.4 \mu_{MBL.SHH} + 0.6 \mu_{MBL.MYC}$.

### Analysis Types Explained

The pipeline automatically performs several types of analyses based on the provided data and contrasts:

* **Differential Expression (DE) (RNA or Ribo):** When comparing two groups (e.g., `GroupA` vs `GroupB`), the script calculates the differential expression at the transcriptional level (if RNA-seq data is provided) and the translational level (if Ribo-seq data is provided). The results are saved in the `RNA/` and `Ribo/` directories, respectively.
    * **Contrast:** `GroupA - GroupB`

* **Translation Efficiency (TE):** This measures the ribosome loading onto a transcript, calculated as Ribo-seq counts divided by RNA-seq counts. The script evaluates TE for unique experimental groups that have both RNA and Ribo data. This analysis identifies genes that are efficiently or inefficiently translated *within* a specific condition. The results are saved in the `TE/` directory.
    * **Contrast:** `Ribo.GroupA - RNA.GroupA`

* **Differential Translation Efficiency (dTE):** This is the most insightful analysis when both RNA-seq and Ribo-seq are available. It identifies genes where the *change* in ribosome loading between two conditions is greater than the *change* in transcription. This pinpoints genes that are specifically regulated at the translational level. Results are saved in the `dTE/` directory.
    * **Contrast:** `(Ribo.GroupA - RNA.GroupA) - (Ribo.GroupB - RNA.GroupB)`

### Example Walkthrough

Let's trace the logic using the Quick Start command: `--contrast_cols "source_id,treatment_id"`.

1.  **Metadata Processing:** The script identifies the unique values in `source_id` (`CD8T`, `CD8T_Eif4g2_KO`) and `treatment_id` (`ACU_AGST`, `CHR_AGST`).
2.  **Contrast Generation:** It combines these factors to create interaction terms (e.g., `CD8T__ACU_AGST`). The script then finds all valid pairs to compare. For example, it will generate a comparison between `CD8T__ACU_AGST` and `CD8T_Eif4g2_KO__ACU_AGST`.
3.  **Analysis Execution:** For this specific pair, the pipeline will perform three analyses:
    * **RNA DE:** Is gene expression different between Knockout and Control cells under acute stimulation?
    * **Ribo DE:** Is ribosome occupancy different between Knockout and Control cells under acute stimulation?
    * **dTE:** Is the translation efficiency of genes regulated differently between Knockout and Control cells under acute stimulation?

The script will do this for all other valid pairs it identifies, such as comparing acute vs. chronic stimulation within the control cells (`CD8T__ACU_AGST` vs. `CD8T__CHR_AGST`), providing a fully automated, comprehensive analysis of the entire experiment.
