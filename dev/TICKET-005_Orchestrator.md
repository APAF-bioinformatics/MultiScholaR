---
ticket_id: LIP-005
title: Orchestrator and summary
status: completed
priority: P0
depends_on:
  - LIP-001
  - LIP-002
  - LIP-003
  - LIP-004
series: LIP
created_at: 2026-01-15
---

# LIP-005: Orchestrator & Summary

## Description

Finalize the Lipidomics workflow by creating the top-level orchestrator and the summary/export module.

## Source Files

- `R/mod_metabolomics.R`
- `R/mod_metab_summary.R` (and `mod_prot_summary.R` for reference)

## New Files

- `R/mod_lipidomics.R`
- `R/mod_lipid_summary.R`

## Tasks

### Orchestrator

- [ ] **Duplicate**: `mod_metabolomics.R` -> `mod_lipidomics.R`.
- [ ] **Refactor**:
  - `mod_lipidomics_ui`: Update ID and title ("Lipidomics Workflow").
  - `mod_lipidomics_server`:
    - Update keys: `omic_type`, `paths_key`.
    - Initialize S4 object (`LipidomicsAssayData`).
    - Call lipidomics sub-modules (`mod_lipid_import`, etc.).
    - Update Stepper icons/labels if needed.
  - `run_lipidomics_app`: Update test launcher to use `mod_lipidomics` and `Lipidomics Test`.

### Summary Module

- [ ] **Duplicate**: `mod_metab_summary.R` -> `mod_lipid_summary.R`.
- [ ] **Refactor**:
  - Ensure it saves/exports `LipidomicsAssayData` correctly.
  - Check for any hardcoded "metabolomics" strings in logging or filenames.

## Verification

- **End-to-End Test**:
  1. Run `run_lipidomics_app()`.
  2. Go through all tabs.
  3. Check logs for any "Metabolomics" errors or incorrect paths.
