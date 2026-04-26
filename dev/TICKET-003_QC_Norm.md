---
ticket_id: LIP-003
title: Quality control and normalization modules
status: completed
priority: P0
depends_on:
  - LIP-001
  - LIP-002
series: LIP
created_at: 2026-01-15
---

# LIP-003: Quality Control & Normalization Modules

## Description

Implement QC and Normalization steps.

## Source Files

- `R/mod_metab_qc.R`
- `R/mod_metab_qc_*.R` (duplicates, intensity, itsd, s4)
- `R/func_metab_qc.R`
- `R/mod_metab_norm.R`
- `R/func_metab_norm.R`

## New Files

- `R/mod_lipid_qc.R`
- `R/mod_lipid_qc_duplicates.R`
- `R/mod_lipid_qc_intensity.R`
- `R/mod_lipid_qc_itsd.R`
- `R/mod_lipid_qc_s4.R`
- `R/func_lipid_qc.R`
- `R/mod_lipid_norm.R`
- `R/func_lipid_norm.R`

## Tasks

### QC Module

- [ ] **Duplicate**: `mod_metab_qc.R` and all `mod_metab_qc_*.R` files.
- [ ] **Refactor**:
  - Update `mod_lipid_qc` to call `mod_lipid_qc_*` submodules.
  - Rename S4 manipulation functions in `mod_lipid_qc_s4.R`.
- [ ] **Functions**: Refactor `func_lipid_qc.R` to handle `LipidomicsAssayData`.

### Normalization Module

- [ ] **Duplicate**: `mod_metab_norm.R` and `func_metab_norm.R`.
- [ ] **Refactor**:
  - Ensure normalization methods work on the lipidomics data structure.
  - Update plots labels.

## Verification

- Test interactive QC plots.
- Run normalization and check success message.
