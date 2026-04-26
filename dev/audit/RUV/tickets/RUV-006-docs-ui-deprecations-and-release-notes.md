---
ticket_id: RUV-006
title: Docs, UI deprecations, and release notes
status: pending
priority: P1
depends_on:
  - RUV-001
  - RUV-002
  - RUV-003
  - RUV-004
  - RUV-005
verified_at: 2026-04-26
source_plan: dev/audit/RUV/fix-plan.md
---

# RUV-006 Docs, UI deprecations, and release notes

## Problem

User-facing text still advertises behavior the code does not currently implement. The UIs still present `weighted_difference` as a normal option, and the proteomics workbooks still describe the optimizer as testing different percentages and finding an “elbow” even though the live code is not yet doing that consistently across omics.

## Current repo verification

- `R/mod_prot_norm_ui.R`, `R/mod_metab_norm_ui.R`, and `R/mod_lipid_norm_ui_helpers.R` still list `"Weighted Difference"` without a deprecation label or help text.
- proteomics workbook text in:
  - `inst/workbooks/proteomics/DIA_workflow_starter.rmd`
  - `inst/workbooks/proteomics/DIA_workflow_limpa_starter.rmd`
  - `inst/workbooks/proteomics/DIA_workflow_limpa_experienced.rmd`
  - `inst/workbooks/proteomics/DIA_workflow_limpa.rmd`
  - `inst/workbooks/proteomics/DIA_workflow_experienced.rmd`
  still claims the optimizer “tests different percentages” and uses the “inflection point” or “elbow”.
- `R/func_pept_s4_objects.R` and `R/func_prot_norm_optimization_helpers.R` still document `weighted_difference` as a normal metric.
- package version is still `0.4.1.2` in `DESCRIPTION`.
- `NEWS.md` does not yet describe the upcoming automatic-RUV behavior changes.

## Files in scope

- `R/mod_prot_norm_ui.R`
- `R/mod_metab_norm_ui.R`
- `R/mod_lipid_norm_ui_helpers.R`
- `R/func_pept_s4_objects.R`
- `R/func_prot_norm_optimization_helpers.R`
- `DESCRIPTION`
- `NEWS.md`
- `inst/workbooks/proteomics/*.rmd` files listed above

## Required changes

- mark `weighted_difference` deprecated in UI labels and help text
- update roxygen/user-facing metric descriptions to match the new contract
- rewrite workbook language so it describes the first-plateau rule truthfully
- add release notes that call out:
  - first-plateau K selection
  - true metabolomics/lipidomics percentage search
  - peptide optimization now using the full objective
  - `weighted_difference` deprecation
  - automatic-mode result changes
- bump the package minor version in `DESCRIPTION`

## Acceptance criteria

- all listed workbook locations use wording consistent with the fixed algorithm
- UI labels/help make `weighted_difference` clearly deprecated
- roxygen/docs no longer describe broken or obsolete behavior as current truth
- `NEWS.md` includes explicit user-facing migration notes for automatic-RUV behavior changes
- `DESCRIPTION` minor version is bumped
- docs state that manual mode semantics are unchanged

## Tests to update or add

- add focused assertions where practical for once-per-call deprecation wording
- ensure documentation-related tests or snapshots, if present, are updated to the new text

## Non-goals

- do not introduce new user-facing tuning parameters for the plateau rule
