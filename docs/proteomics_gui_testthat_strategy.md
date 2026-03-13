# Proteomics GUI testthat Strategy

> High-coverage test harness for MultiScholaR proteomics pipeline

## 1. Overview

This strategy creates a comprehensive `testthat` test suite for the proteomics module by:

1. **Inserting temporary checkpoint code** into the Shiny server pipeline to capture RDS snapshots of data at each major stage
2. **Running the app once** with real data to generate all test fixtures
3. **Removing checkpoint code** and leaving reference comments
4. **Writing testthat scripts** that use the captured snapshots as input/output pairs

### Data Flow & Checkpoint Map

```
                         PROTEOMICS PIPELINE
                         ═══════════════════

  ┌─────────────┐   CP01    ┌─────────────┐   CP02    ┌──────────────┐
  │  Raw File   │──────────▶│  Imported    │──────────▶│  QC Filtered │
  │  (Upload)   │           │  (long fmt)  │           │  (peptide)   │
  └─────────────┘           └─────────────┘           └──────────────┘
                                                             │
                             CP03                            │
  ┌─────────────┐◀───────────────────────────────────────────┘
  │  Rolled Up  │
  │  (protein)  │
  └─────────────┘
        │
        │  CP04                    CP05                    CP06
        ▼                          ▼                        ▼
  ┌──────────┐              ┌────────────┐           ┌───────────┐
  │  Design  │              │ Normalised │           │ RUV-III   │
  │  Matrix  │              │  (loess)   │           │ Corrected │
  └──────────┘              └────────────┘           └───────────┘
                                                          │
                                                          │  CP07
                                                          ▼
                                                   ┌────────────┐
                                                   │ DA Results │
                                                   │  (limma)   │
                                                   └────────────┘
                                                     │        │
                                           CP08      │        │  CP09
                                           ▼         │        ▼
                                    ┌───────────┐    │  ┌──────────┐
                                    │  Volcano  │    │  │ Heatmap  │
                                    │  (Glimma) │    │  │ (Complex │
                                    └───────────┘    │  │  Heatmap)│
                                                     │  └──────────┘
                                           CP10      │
                                           ▼         │
                                    ┌───────────┐    │
                                    │ Enrichment│    │
                                    └───────────┘    │
```

**CP** = Checkpoint. Each checkpoint captures a `.rds` snapshot of the data state.

---

## 2. Naming & Numbering System

### Test File Naming Convention

```
tests/testthat/test-prot-{NN}-{stage_name}.R
```

| File | Stage | Tests |
|------|-------|-------|
| `test-prot-01-import.R` | Data Import | Format detection, all 5 importers |
| `test-prot-02-qc-filtering.R` | QC & Filtering | Empty row removal, replicate filter, missing values |
| `test-prot-03-rollup.R` | Peptide→Protein Rollup | Precursor→peptide, counting functions |
| `test-prot-04-design.R` | Design Matrix | S4 object creation, validation |
| `test-prot-05-normalisation.R` | Normalisation | Between-sample norm, scaling, RUV helpers |
| `test-prot-06-ruv.R` | RUV-III-C | Neg ctrl selection, best k, cancor |
| `test-prot-07-da-analysis.R` | DA Analysis | limma contrasts, result formatting |
| `test-prot-08-volcano.R` | Volcano Plots | Glimma and static volcano |
| `test-prot-09-heatmap.R` | Heatmap | ComplexHeatmap generation |
| `test-prot-10-annotation.R` | Annotation | ID cleaning, isoform handling, UniProt |

### Testdata Directory Structure

```
tests/testdata/
├── prot_checkpoints/
│   ├── cp01_raw_imported.rds         # Raw imported data (list)
│   ├── cp02_qc_filtered_peptide.rds  # After peptide QC
│   ├── cp03_rolled_up_protein.rds    # After rollup to protein
│   ├── cp04_design_matrix.rds        # Design matrix + S4 object
│   ├── cp05_normalised.rds           # After normalisation
│   ├── cp06_ruv_corrected.rds        # After RUV-III
│   ├── cp07_da_results.rds           # DA results list
│   ├── cp08_volcano_input.rds        # Volcano plot input data
│   ├── cp09_heatmap_input.rds        # Heatmap input data
│   └── cp10_enrichment_input.rds     # Enrichment input data
├── mock_data/
│   ├── mock_diann_report.tsv         # Small synthetic DIA-NN file
│   ├── mock_fragpipe_protein.tsv     # Small synthetic FragPipe file
│   ├── mock_maxquant_pg.txt          # Small synthetic MaxQuant file
│   └── mock_pd_tmt.xlsx              # Small synthetic PD-TMT file
└── glimma_fixtures/                  # (existing)
    └── glimma_snapshot_*.rds
```

---

## 3. Checkpoint Capture Code (Temporary)

### 3.1 Capture Helper Function

A reusable helper inserted at the top of the server module:

```r
# --- TESTTHAT CHECKPOINT CAPTURE (TEMPORARY) ---
.capture_checkpoint <- function(data, checkpoint_id, label) {
  if (getOption("multischolar.capture_test_checkpoints", FALSE)) {
    tryCatch({
      base_dir <- getOption("multischolar.checkpoint_dir",
                            file.path("tests", "testdata", "prot_checkpoints"))
      if (!dir.exists(base_dir)) dir.create(base_dir, recursive = TRUE)
      filename <- file.path(base_dir, paste0(checkpoint_id, "_", label, ".rds"))
      saveRDS(data, file = filename)
      logger::log_info(sprintf("CHECKPOINT %s saved: %s", checkpoint_id, filename))
    }, error = function(e) {
      logger::log_warn(sprintf("Checkpoint %s failed: %s", checkpoint_id, e$message))
    })
  }
}
# --- END CHECKPOINT CAPTURE ---
```

### 3.2 Insertion Points

Checkpoints are inserted at these locations in the module server files:

| CP | File | After Line/Section | Data Captured |
|----|------|---------------------|---------------|
| CP01 | `mod_prot_import.R` | After `importDIANNData()` / `importFragPipeData()` etc. returns | The full import result list |
| CP02 | `mod_prot_qc_peptide.R` | After peptide-level QC filtering completes | Filtered peptide data |
| CP03 | `mod_prot_qc_protein_rollup.R` | After `rollUpPrecursorToPeptide()` and protein table created | Protein quant table |
| CP04 | `mod_prot_design_builder.R` | After design matrix assignment and S4 object creation | `ProteinQuantitativeData` S4 object |
| CP05 | `mod_prot_norm.R` | After `normaliseBetweenSamples()` completes | Normalised S4 object |
| CP06 | `mod_prot_norm.R` | After `ruvIII_C_Varying()` completes | RUV-corrected S4 object |
| CP07 | `mod_prot_da.R` | After `differentialAbundanceAnalysis()` returns | Full DA results list |
| CP08 | `mod_prot_da.R` | Before `generateProtDAVolcanoPlotGlimma()` call | Volcano input args |
| CP09 | `mod_prot_da.R` | Before `generateProtDAHeatmap()` call | Heatmap input args |
| CP10 | `mod_prot_enrich.R` | Before enrichment analysis call | Enrichment input data |

### 3.3 Example Checkpoint Insertion (CP01)

In `mod_prot_import.R`, after import completes:

```r
# After the import function returns successfully:
import_result <- importDIANNData(filepath, use_precursor_norm)

# --- TESTTHAT CHECKPOINT CP01 (see test-prot-01-import.R) ---
.capture_checkpoint(import_result, "cp01", "raw_imported")
# --- END CP01 ---
```

After data capture, this gets replaced with just a comment:

```r
import_result <- importDIANNData(filepath, use_precursor_norm)
# [TESTTHAT CP01: raw_imported — see test-prot-01-import.R]
```

---

## 4. Test Script Templates

### 4.1 `test-prot-01-import.R` — Data Import

```r
library(testthat)
library(MultiScholaR)

# === Stage 01: Data Import ===
# Tests the format detection and data import functions.
# Fixtures: tests/testdata/prot_checkpoints/cp01_raw_imported.rds
#           tests/testdata/mock_data/mock_diann_report.tsv

# --- 01.01: Format Detection ---
test_that("detectProteomicsFormat identifies DIA-NN correctly", {
  diann_headers <- c("Protein.Group", "Protein.Ids", "Run",
                     "Precursor.Id", "Modified.Sequence",
                     "Stripped.Sequence", "Precursor.Charge",
                     "Q.Value", "PG.Q.Value", "Precursor.Normalised")
  result <- detectProteomicsFormat(diann_headers, "report.tsv")
  expect_equal(result$format, "diann")
  expect_gt(result$confidence, 0.5)
})

test_that("detectProteomicsFormat identifies FragPipe correctly", {
  fp_headers <- c("Protein ID", "Protein", "Gene", "Description",
                  "WLP530_1 Intensity", "WLP530_2 Intensity",
                  "WLP530_1 MaxLFQ Intensity")
  result <- detectProteomicsFormat(fp_headers, "protein.tsv")
  expect_equal(result$format, "fragpipe")
  expect_gt(result$confidence, 0.3)
})

test_that("detectProteomicsFormat identifies PD-TMT correctly", {
  tmt_headers <- c("Protein FDR Confidence", "Master", "Accession",
                   "Exp. q-value", "Sum PEP Score",
                   "Abundance: F1: 126, Sample1")
  result <- detectProteomicsFormat(tmt_headers, "proteins.xlsx")
  expect_equal(result$format, "pd_tmt")
  expect_gt(result$confidence, 0.5)
})

test_that("detectProteomicsFormat returns unknown for unrecognised data", {
  bad_headers <- c("col_a", "col_b", "col_c")
  result <- detectProteomicsFormat(bad_headers, "random.csv")
  expect_equal(result$format, "unknown")
})

# --- 01.02: DIA-NN Import ---
test_that("importDIANNData returns expected structure", {
  # Uses mock file or captured checkpoint data
  fixture_path <- testthat::test_path("..", "testdata", "prot_checkpoints",
                                      "cp01_raw_imported.rds")
  skip_if(!file.exists(fixture_path), "CP01 fixture not available")
  cp01 <- readRDS(fixture_path)

  expect_type(cp01, "list")
  expect_true("data" %in% names(cp01))
  expect_true("data_type" %in% names(cp01))
  expect_true("column_mapping" %in% names(cp01))
  expect_s3_class(cp01$data, "data.frame")
  expect_true(nrow(cp01$data) > 0)
})

# --- 01.03: Config defaults ---
test_that("getDefaultProteomicsConfig returns valid config list", {
  config <- getDefaultProteomicsConfig()
  expect_type(config, "list")
  expect_true("generalParameters" %in% names(config))
  expect_true("deAnalysisParameters" %in% names(config))
  expect_true("normalizationParameters" %in% names(config))
  expect_true("ruvParameters" %in% names(config))
  expect_equal(config$generalParameters$q_value_threshold, 0.01)
})
```

### 4.2 `test-prot-02-qc-filtering.R` — QC & Filtering

```r
library(testthat)
library(MultiScholaR)

# === Stage 02: QC & Filtering ===
# Tests: removeEmptyRows, removeRowsWithMissingValuesPercentHelper,
#        removeProteinsWithOnlyOneReplicateHelper

# --- 02.01: removeEmptyRows ---
test_that("removeEmptyRows removes all-zero/NA rows", {
  test_df <- data.frame(
    Protein.Ids = c("P1", "P2", "P3"),
    S1 = c(100, 0, 50),
    S2 = c(200, NA, 0),
    S3 = c(300, 0, 75)
  )
  result <- removeEmptyRows(test_df, col_pattern = "^S", row_id = Protein.Ids)
  expect_equal(nrow(result), 2)
  expect_true("P1" %in% result$Protein.Ids)
  expect_false("P2" %in% result$Protein.Ids)
})

# --- 02.02: Checkpoint-based QC validation ---
test_that("QC filtering reduces data appropriately (CP01 -> CP02)", {
  cp01_path <- testthat::test_path("..", "testdata", "prot_checkpoints",
                                   "cp01_raw_imported.rds")
  cp02_path <- testthat::test_path("..", "testdata", "prot_checkpoints",
                                   "cp02_qc_filtered_peptide.rds")
  skip_if(!file.exists(cp01_path) || !file.exists(cp02_path),
          "CP01/CP02 fixtures not available")
  
  cp01 <- readRDS(cp01_path)
  cp02 <- readRDS(cp02_path)

  # QC should reduce or maintain row count, never increase
  expect_lte(nrow(cp02$data), nrow(cp01$data))
  # Column structure should be preserved
  expect_true(all(names(cp01$column_mapping) %in% names(cp02$column_mapping)))
})
```

### 4.3 Summary of Additional Test Files

Each follows the same pattern. Key tests per file:

| File | Key Tests |
|------|-----------|
| `test-prot-03-rollup.R` | `rollUpPrecursorToPeptideHelper` aggregation correctness, `calcPeptidesPerProtein`, `countProteinsPerRun`, `count_num_peptides` |
| `test-prot-04-design.R` | `ProteinQuantitativeData` S4 construction, validation rules (mismatched samples → error), slot accessors |
| `test-prot-05-normalisation.R` | `scaleCenterAndFillMissing` output properties, `updateRuvParameters` config updates, `getRuvIIIReplicateMatrixHelper` matrix shape |
| `test-prot-06-ruv.R` | `getNegCtrlProtAnovaHelper` returns boolean vector, `findBestK` with mock cancor plot, `findBestKForAssayList` handles edge cases |
| `test-prot-07-da-analysis.R` | CP06→CP07 transition, DA result structure validation, contrast matching |
| `test-prot-08-volcano.R` | `generateProtDAVolcanoPlotGlimma` returns `htmlwidget` (existing test enhanced), `generateProtDAVolcanoStatic` returns `ggplot` |
| `test-prot-09-heatmap.R` | `generateProtDAHeatmap` returns `ComplexHeatmap` or NULL for edge cases |
| `test-prot-10-annotation.R` | `cleanIsoformNumber` regex, `detectProteomicsFormat` edge cases, `goIdToTerm` formatting |

---

## 5. Mock Data Generation

For tests that should **not** depend on captured checkpoints (pure unit tests), we create small synthetic datasets:

### `tests/testdata/mock_data/mock_diann_report.tsv`

A 50-row synthetic DIA-NN report with:
- 5 proteins, 10 peptides, 5 runs
- Realistic column names matching DIA-NN output
- Some deliberate NA values for QC filter testing

### `tests/testdata/mock_data/mock_fragpipe_protein.tsv`

A 20-row synthetic FragPipe `protein.tsv` with:
- 20 proteins, 4 samples
- Both regular and MaxLFQ Intensity columns

These will be generated by an R script: `tests/generate_mock_data.R`

---

## 6. Execution Plan

### Phase 1: Setup & Capture Infrastructure
1. Create `tests/testdata/prot_checkpoints/` directory
2. Add `.capture_checkpoint()` helper to a utility file
3. Add UI toggle checkbox (similar to existing Glimma capture toggle)
4. Insert checkpoint calls at all 10 locations in module server files

### Phase 2: Data Capture
1. Launch the app with `options(multischolar.capture_test_checkpoints = TRUE)`
2. Run through a complete proteomics workflow (import → DA → visualisation)
3. Verify all 10 `.rds` checkpoint files are created

### Phase 3: Write Tests
1. Create the 10 `test-prot-{NN}-*.R` files
2. Create `tests/generate_mock_data.R` and generate synthetic files
3. Wire up `tests/testthat.R` (if not already present)

### Phase 4: Cleanup
1. Remove all `.capture_checkpoint()` calls from module server files
2. Replace with comment markers: `# [TESTTHAT CP{NN}: label — see test-prot-{NN}-*.R]`
3. Remove the capture helper function
4. Keep the UI toggle for future re-captures (optional)

### Phase 5: Verification
1. Run `devtools::test()` to verify all tests pass
2. Run `covr::package_coverage()` to measure coverage improvement
3. Document coverage metrics

---

## 7. Verification Plan

### Automated Tests

```bash
# Run all proteomics tests
cd /Users/ignatiuspang/Workings/2025/MultiScholaR-volcano-fix
Rscript -e 'devtools::test(filter = "prot")'

# Run a specific test file
Rscript -e 'testthat::test_file("tests/testthat/test-prot-01-import.R")'

# Check coverage for proteomics functions
Rscript -e 'covr::file_coverage(
  source_files = c("R/func_prot_import.R", "R/func_prot_qc.R",
                    "R/func_prot_rollup.R", "R/func_prot_norm.R",
                    "R/func_prot_da.R", "R/func_prot_annotation.R"),
  test_files = list.files("tests/testthat", pattern = "test-prot-",
                           full.names = TRUE)
)'
```

### Manual Verification

1. **Checkpoint capture**: Launch the app, enable capture, run a workflow, verify `.rds` files appear in `tests/testdata/prot_checkpoints/`
2. **Post-cleanup**: Confirm no `capture_checkpoint` calls remain in R source files (search for the string)
3. **Test isolation**: Each test file should run independently without requiring other tests

---

## 8. Priority Order

Given the recent volcano plot debugging work, priority order is:

1. **`test-prot-01-import.R`** — Pure unit tests, no checkpoint needed
2. **`test-prot-08-volcano.R`** — Enhances existing `test-glimma-plot.R`
3. **`test-prot-10-annotation.R`** — Pure unit tests for ID cleaning
4. **`test-prot-02-qc-filtering.R`** — Data-dependent, needs checkpoints
5. **`test-prot-04-design.R`** — S4 constructor validation
6. Remaining files in numerical order

<!-- APAF Bioinformatics | proteomics_gui_testthat_strategy.md | Approved | 2026-03-13 -->
