# Implementation Plan: Group-Aware Intensity Filter & Test Data Reorg

This implementation plan details the strategy to introduce a group-aware intensity filter for peptides, enable data capturing, and reorganize the captured test data.

## Proposed Changes

### Filtering Logic
#### [MODIFY] func_prot_qc_peptide.R 
Update `peptideIntensityFilteringHelper` to incorporate group-aware logic. We will:
1.  Require `design_matrix`, `sample_id_column`, and `grouping_variable` arguments.
2.  Calculate the threshold (or use provided `min_peptide_intensity_threshold`).
3.  Calculate the percentage of samples below the threshold *per group*.
4.  Remove peptides if the percentage of groups failing the (`groupwise_percentage_cutoff`) is higher than `max_groups_percentage_cutoff`.
**Note:** This logic will be closely mirrored from the optimized `removeRowsWithMissingValuesPercentHelper`.

#### [MODIFY] peptideVsSamplesS4Objects.R / mod_prot_qc_*.R
Update callers to `peptideIntensityFilteringHelper` to pass the necessary arguments for grouping metadata.

---

### Data Checkpoint Capturing & Reorganization
#### [MODIFY] utils_shiny_logging.R
Update `.capture_checkpoint` to construct its save path dynamically.
- Read `options("multischolar.checkpoint_dataset")` (default: "sepsis").
- Read `options("multischolar.checkpoint_omics_layer")` (default: "proteomics").
- Define `base_dir <- file.path("tests", "testdata", dataset, omics_layer)`.

#### [MODIFY] Modules with commented captures
Uncomment `.capture_checkpoint()` in the following files:
-   `R/mod_prot_import.R`
-   `R/mod_prot_qc_peptide_impute.R`
-   `R/mod_prot_qc_protein_rollup.R`
-   `R/mod_prot_design.R`
-   `R/mod_prot_norm.R`
-   `R/mod_prot_da.R`
-   `R/mod_prot_enrich.R`

---

## Storage Restructuring
Move existing test data from old structures into the new structure:
1.  Create dir `tests/testdata/sepsis/proteomics`.
2.  Move all `.rds` objects from `tests/testdata/prot_checkpoints` to the new sepsis folder.

## Verification Plan
### Automated Verification
None currently defined via unit tests, as we rely on the captured data first.

### Manual Verification
1.  Start the app and load the **Sepsis dataset**.
2.  Toggle "Capture Test Checkpoints (Developer)" in the UI (or set the options manually in the R console).
3.  Proceed through the proteomics pipeline. Confirm that `.rds` files are newly captured in `tests/testdata/sepsis/proteomics`.
4.  Change the dataset to **Neurolincs**, input the `.parquet` files, and run the pipeline making sure to utilize `duplicateCorrelation`. Verify the data is saved in `tests/testdata/neurolincs_book_chapter/proteomics`.
