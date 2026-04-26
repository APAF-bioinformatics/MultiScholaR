---
ticket_id: LIP-002
title: Core data acquisition modules
status: completed
priority: P0
depends_on:
  - LIP-001
series: LIP
created_at: 2026-01-15
---

# LIP-002: Core Data Acquisition Modules

## Description

Implement the initial steps of the workflow: Import, Design Matrix, and QC.

## Source Files

- `R/mod_metab_import.R`
- `R/func_metab_import.R`
- `R/mod_metab_design.R`
- `R/mod_metab_design_builder.R`
- `R/func_general_design.R` (Shared?)

## New Files

- `R/mod_lipid_import.R`
- `R/func_lipid_import.R`
- `R/mod_lipid_design.R`
- `R/mod_lipid_design_builder.R`

## Tasks

### Import Module

- [ ] **Duplicate**: `mod_metab_import.R` -> `mod_lipid_import.R`
- [ ] **Refactor UI**:
  - Rename function: `mod_lipid_import_ui`.
  - Rename ID: `mod_metab_import` -> `mod_lipid_import`.
  - Update labels: "Metabolomics Import" -> "Lipidomics Import".
- [ ] **Refactor Server**:
  - Rename function: `mod_lipid_import_server`.
  - Use `createLipidomicsAssayData` instead of `createMetaboliteAssayData`.
  - Update `omic_type` hardcoding if present.

### Design Modules

- [ ] **Duplicate**: `mod_metab_design.R` -> `mod_lipid_design.R`
- [ ] **Duplicate**: `mod_metab_design_builder.R` -> `mod_lipid_design_builder.R`
- [ ] **Refactor**:
  - Updates IDs and references to S4 objects.
  - Ensure `LipidomicsAssayData` is accepted.

## Verification

- Run `run_lipidomics_app()` (after Orchestrator is done) and test Import tab.
- Verify Design tab populates correctly.
