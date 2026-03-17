# Metabolomics and Lipidomics GUI Refactor Strategy

> A unified strategy to refactor and modularize the Shiny UI components for the Metabolomics and Lipidomics pipelines, bringing them in line with the recent Proteomics GUI refactoring and applying testable checkpoints.

## 1. Context and Motivation

The Proteomics Quality Control (QC) GUI has recently been aggressively refactored. Rather than a monolithic `mod_prot_qc.R` file, the UI and server logic is now split across many highly specialized modules (e.g., `mod_prot_qc_peptide_intensity.R`, `mod_prot_qc_protein_rollup.R`, `mod_prot_da.R`). This architectural shift provides several benefits:
- **Maintainability:** Smaller, tightly scoped files are easier to navigate and debug.
- **Testability:** State can be captured (`.rds` checkpoints) and tested (`testthat`) at granular boundaries.
- **Reusability:** Common UI patterns can be abstracted.

Currently, the Metabolomics and Lipidomics modules still rely on larger, more monolithic files natively (e.g. `mod_metab_norm.R` is over 2,100 lines long).

### Goal
Refactor the Metabolomics and Lipidomics GUI components to mirror the modular architecture of the Proteomics pipeline, while simultaneously integrating the `testthat` checkpoint strategies defined in:
- `metabolomics_gui_testthat_strategy.md`
- `lipidomics_gui_testthat_strategy.md`

---

## 2. Current Architecture Audit

Taking Metabolomics as the primary example (Lipidomics is structurally identical):

| Module | Lines of Code | Refactoring Need |
| :--- | :--- | :--- |
| `mod_metab_norm.R` | ~2,107 | **HIGH**. Very large monolithic file handling both missing value imputation, ITSD normalization, scaling, and the entire UI layout for these steps. |
| `mod_metab_da.R` | ~1,250 | **MEDIUM**. Handles the Differential Abundance UI, contrasts, and heatmap/volcano plot rendering. |
| `mod_metab_qc_itsd.R` | ~408 | **LOW**. Already partially modularized. |
| `mod_metab_qc_duplicates.R` | ~400 | **LOW**. Already partially modularized. |
| `mod_metab_qc_intensity.R`| ~340 | **LOW**. Already partially modularized. |
| `mod_metab_qc_s4.R` | ~386 | **LOW**. |

*Note: The Lipidomics `mod_lipid_*.R` files have almost identical line counts and structures.*

---

## 3. Proposed Refactoring Strategy

### 3.1 Unifying the `testthat` Checkpoints with Module Boundaries

The checkpoint strategy must map directly to the new file boundaries. 

#### Data Flow & File Map (Metabolomics / Lipidomics)

```
                            PIPELINE
                            ════════
  ┌────────────────┐ 
  │ Import Module  │  ──▶ [CP01: Raw Imported]
  │ (mod_*_import) │
  └────────────────┘
          │
          ▼
  ┌────────────────┐
  │ QC Modules     │  ──▶ [CP02: QC Filtered]
  │ (mod_*_qc_*)   │  (duplicates, intensity, ITSD QC)
  └────────────────┘
          │
          ▼
  ┌────────────────┐
  │ Design Builder │  ──▶ [CP03: Design Matrix S4]
  └────────────────┘
          │
          ▼
  ┌────────────────┐
  │ Normalisation  │  ──▶ [CP04: Normalised]
  │ (mod_*_norm_*) │  **[TARGET FOR REFACTORING]**
  └────────────────┘
          │
          ▼
  ┌────────────────┐
  │   DA Analysis  │  ──▶ [CP05: DA Results]
  │ (mod_*_da_*)   │  **[TARGET FOR REFACTORING]**
  └────────────────┘
          │
      ┌───┴───┐
      ▼       ▼
 [CP06:    [CP07:
 Volcano]  Heatmap]
```

### 3.2 Splitting the Normalisation Module (`mod_*_norm.R`)

The highest priority for refactoring is `mod_metab_norm.R` and `mod_lipid_norm.R`. These should be split into smaller sub-modules based on functional steps:

1.  **`mod_*_norm_impute.R`**: Handing missing value imputation UI and server logic (e.g., KNN, Min, Zero).
2.  **`mod_*_norm_itsd.R`**: Handling Internal Standard (ITSD) assignment, visualization, and application to the dataset.
3.  **`mod_*_norm_scale.R`**: Handling data transformation (e.g., log2) and scaling (e.g., Pareto, Auto).
4.  **`mod_*_norm.R` (Orchestrator)**: Reduced to a lightweight wrapper that calls the sub-modules sequentially, collects their outputs, and manages the `testthat` [CP04] snapshot.

### 3.3 Splitting the DA Module (`mod_*_da.R`)

While smaller than normalization, DA is complex and benefits from UI isolation:

1.  **`mod_*_da_setup.R`**: Contrast builder, formula selection, and running the `limma` model ([CP05] checkpoint).
2.  **`mod_*_da_volcano.R`**: UI and server logic exclusively for the Glimma/Static volcano plot rendering ([CP06] checkpoint).
3.  **`mod_*_da_heatmap.R`**: UI and server logic exclusively for the ComplexHeatmap rendering ([CP07] checkpoint).
4.  **`mod_*_da.R` (Orchestrator)**: Coordinates the sub-modules and manages state sharing.

---

## 4. Execution Plan

This refactoring should be done sequentially to avoid breaking the application.

### Phase 1: Test Harness Preparation (Current State Capture)
1. Inject the `.capture_checkpoint()` functions into the *existing* monolithic `mod_metab_*.R` and `mod_lipid_*.R` files exactly as outlined in the `*_gui_testthat_strategy.md` docs.
2. Run the application and generate the `tests/testdata/sepsis/*/*.rds` snapshot fixtures.
3. Write the `tests/testthat/test-metab-*.R` and `test-lipid-*.R` scripts.
4. Verify all tests pass. **(This ensures we have a safety net before breaking apart the files).**

### Phase 2: Normalisation Refactoring
1. Create `mod_metab_norm_impute.R`, `mod_metab_norm_itsd.R`, and `mod_metab_norm_scale.R`.
2. Extract the corresponding `ui_` and `server_` functions from `mod_metab_norm.R`.
3. Re-wire `mod_metab_norm.R` to call these new sub-modules using `shiny::callModule()` or `shiny::moduleServer()`.
4. Run `devtools::test(filter = "metab")` to ensure the outputs match the Phase 1 checkpoints perfectly.
5. Repeat for Lipidomics (`mod_lipid_norm_*.R`).

### Phase 3: Differential Abundance (DA) Refactoring
1. Create `mod_metab_da_setup.R`, `mod_metab_da_volcano.R`, and `mod_metab_da_heatmap.R`.
2. Extract the corresponding UI and server components from `mod_metab_da.R`.
3. Re-wire the orchestrator.
4. Run tests to verify `[CP05]`, `[CP06]`, and `[CP07]` outputs remain identical.
5. Repeat for Lipidomics.

### Phase 4: Cleanup & Validation
1. Remove `.capture_checkpoint()` code.
2. Verify UI consistency in the running application (ensure no broken CSS or missing inputs).

<!-- APAF Bioinformatics | gui_refactor_strategy.md | Approved | 2026-03-17 -->
