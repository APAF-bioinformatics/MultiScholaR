---
ticket_id: LIP-004
title: Differential expression
status: completed
priority: P1
depends_on:
  - LIP-003
series: LIP
created_at: 2026-01-15
---

# LIP-004: Differential Expression

## Description

Implement Differential Expression Analysis using Limma.

## Source Files

- `R/mod_metab_de.R`
- `R/func_metab_de.R`

## New Files

- `R/mod_lipid_de.R`
- `R/func_lipid_de.R`

## Tasks

- [ ] **Duplicate**: `mod_metab_de.R` -> `mod_lipid_de.R`.
- [ ] **Refactor**:
  - Update UI: `mod_lipid_de_ui`.
  - Update Server: `mod_lipid_de_server`.
  - Update `load_filtered_session` logic to look for lipidomics session data (filename: `lipid_filtered_session_data_latest.rds`).
  - Rename `runMetabolitesDE` to `runLipidomicsDE`.
  - Update `MetabolomicsDifferentialAbundanceResults` references if that class is also renamed (likely yes).

## Verification

- Run DE analysis on dummy data.
- Verify Volcano plot and Heatmap generation.
