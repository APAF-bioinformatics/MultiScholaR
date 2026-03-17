# Metabolomics GUI testthat Strategy

> High-coverage test harness for MultiScholaR metabolomics pipeline

## 1. Overview

This strategy creates a comprehensive `testthat` test suite for the metabolomics module by:

1. **Inserting temporary checkpoint code** into the Shiny server pipeline to capture RDS snapshots of data at each major stage
2. **Running the app once** with real data to generate all test fixtures
3. **Removing checkpoint code** and leaving reference comments
4. **Writing testthat scripts** that use the captured snapshots as input/output pairs

### Data Flow & Checkpoint Map

```
                          METABOLOMICS PIPELINE
                          ═════════════════════

  ┌──────────────┐   CP01    ┌──────────────┐   CP02    ┌──────────────┐
  │   Raw File   │──────────▶│   Imported   │──────────▶│  QC Filtered │
  │   (Upload)   │           │   (tables)   │           │ (Duplicates, │
  └──────────────┘           └──────────────┘           │ Intensity)   │
                                                        └──────────────┘
                                                              │
                              CP03                            │
   ┌──────────────┐◀──────────────────────────────────────────┘
   │    Design    │
   │    Matrix    │
   │  (S4 Object) │
   └──────────────┘
         │
         │   CP04                    CP05
         ▼                           ▼
   ┌────────────┐             ┌────────────┐
   │ Normalised │             │ DA Results │
   │   (ITSD,   │────────────▶│  (limma)   │
   │   loess)   │             └────────────┘
   └────────────┘                  │        │
                                   │        │
                         CP06      │        │  CP07
                         ▼         │        ▼
                  ┌───────────┐    │  ┌──────────┐
                  │  Volcano  │    │  │ Heatmap  │
                  │  (Glimma) │    │  │ (Complex │
                  └───────────┘    │  │  Heatmap)│
                                   │  └──────────┘
```

**CP** = Checkpoint. Each checkpoint captures a `.rds` snapshot of the data state.

---

## 2. Naming & Numbering System

### Test File Naming Convention

```
tests/testthat/test-metab-{NN}-{stage_name}.R
```

| File | Stage | Tests |
|------|-------|-------|
| `test-metab-01-import.R` | Data Import | File format detection, parser functions from `mod_metab_import.R` |
| `test-metab-02-qc-filtering.R` | QC & Filtering | Duplicate removal, ITSD QC, Intensity filtering |
| `test-metab-03-design.R` | Design Matrix | S4 object (`MetaboliteAssayData`) creation, validation |
| `test-metab-04-normalisation.R` | Normalisation | ITSD normalisation, missing value imputation, scaling |
| `test-metab-05-da-analysis.R` | DA Analysis | limma contrasts, result formatting in `mod_metab_da.R` |
| `test-metab-06-volcano.R` | Volcano Plots | Glimma and static volcano rendering |
| `test-metab-07-heatmap.R` | Heatmap | ComplexHeatmap generation for metabolomics data |

### Testdata Directory Structure

```
tests/testdata/
├── sepsis/
│   └── metabolomics/
│       ├── cp01_raw_imported.rds         # Raw imported data (list)
│       ├── cp02_qc_filtered.rds          # After QC filtering stages
│       ├── cp03_design_matrix.rds        # Design matrix + S4 object
│       ├── cp04_normalised.rds           # After normalisation/imputation
│       ├── cp05_da_results.rds           # DA results list
│       ├── cp06_volcano_input.rds        # Volcano plot input data
│       └── cp07_heatmap_input.rds        # Heatmap input data
├── mock_data/
│   ├── mock_metab_quant.tsv          # Small synthetic quant file
│   └── mock_metab_annotation.tsv     # Small synthetic annotation file
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
                            file.path("tests", "testdata", "sepsis", "metabolomics"))
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

| CP | File | After Line/Section | Data Captured |
|----|------|---------------------|---------------|
| CP01 | `mod_metab_import.R` | After file parser function returns | The full import result list |
| CP02 | `mod_metab_qc_s4.R` (or similar) | After all multi-stage QC filtering completes | Filtered metabolomics data |
| CP03 | `mod_metab_design_builder.R`| After design matrix assignment and S4 object creation | `MetaboliteAssayData` S4 object |
| CP04 | `mod_metab_norm.R` | After normalisation and imputation steps complete | Normalised S4 object |
| CP05 | `mod_metab_da.R` | After differential abundance analysis returns | Full DA results list |
| CP06 | `mod_metab_da.R` | Before `generateMetabDAVolcanoPlotGlimma` call | Volcano input args |
| CP07 | `mod_metab_da.R` | Before `generateMetabDAHeatmap` call | Heatmap input args |

---

## 4. Test Script Templates

### `test-metab-01-import.R` — Data Import

```r
library(testthat)
library(MultiScholaR)

# === Stage 01: Data Import ===
# Tests the format detection and data import functions for metabolomics.

test_that("Metabolomics import returning expected structure", {
  fixture_path <- testthat::test_path("..", "testdata", "sepsis", "metabolomics",
                                      "cp01_raw_imported.rds")
  skip_if(!file.exists(fixture_path), "CP01 fixture not available")
  cp01 <- readRDS(fixture_path)

  expect_type(cp01, "list")
  expect_true("data" %in% names(cp01))
  expect_s3_class(cp01$data, "data.frame")
  expect_true(nrow(cp01$data) > 0)
})
```

---

## 5. Execution Plan

### Phase 1: Setup & Capture Infrastructure
1. Create `tests/testdata/sepsis/metabolomics/` directory.
2. Add `.capture_checkpoint()` helper.
3. Add UI toggle checkbox for capturing checkpoints (if not active globally).
4. Insert checkpoint calls at the 7 locations in the `mod_metab_*` server files.

### Phase 2: Data Capture
1. Launch the app with `options(multischolar.capture_test_checkpoints = TRUE)`.
2. Run through a complete metabolomics workflow (import → QC → Norm → DA → Visualisation).
3. Verify all 7 `.rds` checkpoint files are created.

### Phase 3: Write Tests
1. Create the `test-metab-{NN}-*.R` files in `tests/testthat/`.
2. Create `tests/generate_mock_data.R` and generate small synthetic files so some tests can run without the large RDS fixtures.
3. Wire up `tests/testthat.R`.

### Phase 4: Cleanup
1. Remove all `.capture_checkpoint()` calls from module server files.
2. Replace with comment markers: `# [TESTTHAT CP{NN}: label — see test-metab-{NN}-*.R]`.
3. Check `devtools::test(filter="metab")` passes.

<!-- APAF Bioinformatics | metabolomics_gui_testthat_strategy.md | Approved | 2026-03-17 -->
