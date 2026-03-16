# Plan: Neurolincs Test Data Capture
## Objective
Capture checkpoint test data in `.rds` format while processing the Neurolincs dataset to test:
1. The `.parquet` data import parser
2. The `duplicateCorrelation` feature handling technical replicates.

## Data Source
- Dataset: `/Users/ignatiuspang/Workings/2026/neurolincs_book_chapter_analysis/data/proteomics/report.parquet`
- Contains 21 samples. One of the samples (`iMN_5`) has three technical replicates (`3A`, `3B`, `3C`).

## Execution Plan
1.  **Configure Environment Options:**
    Start an R session (or the MultiScholaR app) and set global options to ensure test data capturing is active and targets the correct directory:
    ```r
    options(multischolar.capture_checkpoints = TRUE)
    options(multischolar.checkpoint_dataset = "neurolincs_book_chapter")
    options(multischolar.checkpoint_omics_layer = "proteomics")
    ```

2.  **Dataset Import:**
    - File: `report.parquet`
    - Let the new `.parquet` reader handle the input.
    - Set columns as per the standard output configurations.

3.  **Define Sample Metadata & Technical Replicates:**
    Upload or use the provided design matrix at `docs/neurolincs_design.csv`.
    -   *Sample Column:* `Run`
    -   *Grouping:* `Group` (Control, C9_ALS, SMA_Type1)
    -   *Technical Replicate Grouping:* `BioReplicate`
        -   This ensures the columns ending in 3A, 3B, 3C share the same `BioReplicate` ID (`iMN_5`), while others have unique IDs (e.g. `iMN_1`, `iMN_10`, etc.).

4.  **Run Pipeline Modules to DA:**
    - Proceed through Peptide QC, Imputation, Protein Rollup, Protein QC, and Normalization.
    - Go to **Differential Analysis (DA)**. Set the design to use the standard grouping and, critically, configure the duplicate correlation to use the `BioReplicate` column as the blocking/random effect variable.
    - Execute DA.

5.  **Validation of Capture:**
    - Verify that the `tests/testdata/neurolincs_book_chapter/proteomics/` directory is populated.
    - Verify `.rds` files exist for each pipeline checkpoint.
