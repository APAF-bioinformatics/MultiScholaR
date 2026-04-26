---
ticket_id: LIP-001
title: Lipidomics S4 objects and helpers
status: completed
priority: P0
depends_on: []
series: LIP
created_at: 2026-01-15
---

# LIP-001: Lipidomics S4 Objects & Helpers

## Description

Establish the data structures for Lipidomics by cloning and refactoring Metabolomics S4 classes.

## Source Files

- `R/func_metab_s4_objects.R`
- `R/utils_workflow_state.R` (Check if updates needed)

## New Files

- `R/func_lipid_s4_objects.R`

## Tasks

- [ ] **Duplicate File**: Copy `R/func_metab_s4_objects.R` to `R/func_lipid_s4_objects.R`.
- [ ] **Rename Class**:
  - Change `setClass("MetaboliteAssayData", ...)` to `setClass("LipidomicsAssayData", ...)`.
  - Update slots if needed (likely keep same structure).
  - Update `@exportClass` and Roxygen comments.
- [ ] **Rename Constructor**:
  - Change `createMetaboliteAssayData` to `createLipidomicsAssayData`.
  - Update `new("MetaboliteAssayData", ...)` to `new("LipidomicsAssayData", ...)`.
- [ ] **Rename Methods**:
  - Update method signatures: `setMethod(f = "plotPca", signature = "LipidomicsAssayData", ...)`.
  - Update `plotRle`, `plotDensity` signatures.
- [ ] **Refactor Internals**:
  - Search for variables named `metabolite_data` and rename to `lipid_data` (optional but recommended for clarity).
  - Ensure error messages refer to "Lipidomics" not "Metabolomics".

## Verification

- Load package: `devtools::load_all()`
- Test:
  ```r
  obj <- createLipidomicsAssayData(
      metabolite_data = list(Pos = data.frame(...)), # Arguments might still be named metabolite_data unless changed in definition
      design_matrix = data.frame(...)
  )
  stopifnot(is(obj, "LipidomicsAssayData"))
  ```
