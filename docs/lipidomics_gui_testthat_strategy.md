# Lipidomics GUI testthat Strategy

> High-coverage test harness for MultiScholaR lipidomics pipeline

## 1. Overview

This strategy creates a comprehensive `testthat` test suite for the lipidomics module by:

1. **Inserting temporary checkpoint code** into the Shiny server pipeline to capture RDS snapshots of data at each major stage
2. **Running the app once** with real data to generate all test fixtures
3. **Removing checkpoint code** and leaving reference comments
4. **Writing testthat scripts** that use the captured snapshots as input/output pairs

### Data Flow & Checkpoint Map

```
                            LIPIDOMICS PIPELINE
                            ═══════════════════

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
tests/testthat/test-lipid-{NN}-{stage_name}.R
```

| File | Stage | Tests |
|------|-------|-------|
| `test-lipid-01-import.R` | Data Import | File format detection, parser functions from `mod_lipid_import.R`, mapping lipid classes |
| `test-lipid-02-qc-filtering.R` | QC & Filtering | Duplicate removal, ITSD QC, Intensity filtering |
| `test-lipid-03-design.R` | Design Matrix | S4 object (`LipidomicsAssayData`) creation, validation |
| `test-lipid-04-normalisation.R` | Normalisation | ITSD normalisation, missing value imputation, scaling |
| `test-lipid-05-da-analysis.R` | DA Analysis | limma contrasts, result formatting in `mod_lipid_da.R` |
| `test-lipid-06-volcano.R` | Volcano Plots | Glimma and static volcano rendering |
| `test-lipid-07-heatmap.R` | Heatmap | ComplexHeatmap generation for lipidomics data |

### Testdata Directory Structure

```
tests/testdata/
├── sepsis/
│   └── lipidomics/
│       ├── cp01_raw_imported.rds         # Raw imported data (list)
│       ├── cp02_qc_filtered.rds          # After QC filtering stages
│       ├── cp03_design_matrix.rds        # Design matrix + S4 object
│       ├── cp04_normalised.rds           # After normalisation/imputation
│       ├── cp05_da_results.rds           # DA results list
│       ├── cp06_volcano_input.rds        # Volcano plot input data
│       └── cp07_heatmap_input.rds        # Heatmap input data
├── mock_data/
│   ├── mock_lipid_quant.tsv          # Small synthetic quant file
│   └── mock_lipid_annotation.tsv     # Small synthetic annotation file
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
                            file.path("tests", "testdata", "sepsis", "lipidomics"))
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
| CP01 | `mod_lipid_import.R` | After file parser function returns | The full import result list |
| CP02 | `mod_lipid_qc_s4.R` (or similar) | After all multi-stage QC filtering completes | Filtered lipidomics data |
| CP03 | `mod_lipid_design_builder.R`| After design matrix assignment and S4 object creation | `LipidomicsAssayData` S4 object |
| CP04 | `mod_lipid_norm.R` | After normalisation and imputation steps complete | Normalised S4 object |
| CP05 | `mod_lipid_da.R` | After differential abundance analysis returns | Full DA results list |
| CP06 | `mod_lipid_da.R` | Before `generateLipidDAVolcanoPlotGlimma` call | Volcano input args |
| CP07 | `mod_lipid_da.R` | Before `generateLipidDAHeatmap` call | Heatmap input args |

---

## 4. Test Script Templates

### `test-lipid-01-import.R` — Data Import

```r
library(testthat)
library(MultiScholaR)

# === Stage 01: Data Import ===
# Tests the format detection and data import functions for lipidomics.

test_that("Lipidomics import returning expected structure", {
  fixture_path <- testthat::test_path("..", "testdata", "sepsis", "lipidomics",
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
1. Create `tests/testdata/sepsis/lipidomics/` directory.
2. Add `.capture_checkpoint()` helper.
3. Add UI toggle checkbox for capturing checkpoints (if not active globally).
4. Insert checkpoint calls at the 7 locations in the `mod_lipid_*` server files.

### Phase 2: Data Capture
1. Launch the app with `options(multischolar.capture_test_checkpoints = TRUE)`.
2. Run through a complete lipidomics workflow (import → QC → Norm → DA → Visualisation).
3. Verify all 7 `.rds` checkpoint files are created.

### Phase 3: Write Tests
1. Create the `test-lipid-{NN}-*.R` files in `tests/testthat/`.
2. Create `tests/generate_mock_data.R` and generate small synthetic files so some tests can run without the large RDS fixtures.
3. Wire up `tests/testthat.R`.

### Phase 4: Cleanup
1. Remove all `.capture_checkpoint()` calls from module server files.
2. Replace with comment markers: `# [TESTTHAT CP{NN}: label — see test-lipid-{NN}-*.R]`.
3. Check `devtools::test(filter="lipid")` passes.

<!-- APAF Bioinformatics | lipidomics_gui_testthat_strategy.md | Approved | 2026-03-17 -->
