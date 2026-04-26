---
ticket_id: RUV-005
title: Peptide optimizer reconciliation
status: pending
priority: P0
depends_on:
  - RUV-001
verified_at: 2026-04-26
source_plan: dev/audit/RUV/fix-plan.md
---

# RUV-005 Peptide optimizer reconciliation

## Problem

Peptide automatic optimization is still solving the wrong objective and still has a conditional ANOVA design bug. The optimizer runs through `ruvCancorFast()`, not the full `ruvCancor()` path later used for peptide QC and analysis, and `getNegCtrlProtAnovaPeptides()` still removes `group_id` even when that column is the requested grouping variable.

## Current repo verification

- `R/func_pept_s4_norm_methods.R::findBestNegCtrlPercentagePeptides()` still calls `ruvCancorFast()` inside the percentage loop.
- the peptide optimizer still uses the old trace schema: `percentage`, `num_controls`, no realized-percentage fields, no status/error diagnostics.
- winner selection still relies on `which.max()` over `composite_score`.
- `getNegCtrlProtAnovaPeptides()` still builds the ANOVA design with `dplyr::select(-!!sym(group_id))`, which drops the grouping column when `ruv_grouping_variable == group_id`.
- `ruvCancorFast()` and `ruvCancor()` still resolve `ctrl` through `checkParamsObjectFunctionSimplify(...)`, but the audit re-check confirmed the earlier “caller ctrl is silently overridden” claim is a retracted non-defect and should not be reopened without new evidence.
- `tests/testthat/test-prot-04-design.R` still treats `weighted_difference` as a normal successful path and still mocks `ruvCancorFast()`.
- `tests/testthat/test-prot-04d-peptide-s4-normalization-characterization.R` still asserts the old trace columns and the old helper math.

## Files in scope

- `R/func_pept_s4_norm_methods.R`
- `tests/testthat/test-prot-04-design.R`
- `tests/testthat/test-prot-04d-peptide-s4-normalization-characterization.R`

## Required changes

- guard ANOVA design-column removal so `group_id` is only dropped when it is not also the active grouping variable
- switch peptide automatic optimization from `ruvCancorFast()` to `ruvCancor()`
- mark `ruvCancorFast()` deprecated, but retain it for one release cycle
- move peptide automatic K selection, scoring, and deterministic winner selection onto the hardened contract from `RUV-001`
- emit the `weighted_difference` deprecation warning once per public optimizer call
- expand peptide `optimization_results` to the same richer trace schema used in proteomics

## Acceptance criteria

- `getNegCtrlProtAnovaPeptides()` keeps the grouping column when `ruv_grouping_variable == group_id`
- peptide optimization no longer selects candidates using `ruvCancorFast()`
- `ruvCancorFast()` warns on direct use
- peptide automatic optimization preserves one trace row per tested percentage, including failures
- peptide results expose requested and realized control counts and percentages
- ties are resolved with the same stable ordering as the proteomics optimizer
- top-level `best_k` is always scalar

## Tests to update or add

- add a regression proving the grouping column is preserved when `group_id` and `ruv_grouping_variable` are the same
- replace fast-path mocks with full-objective mocks in peptide optimizer tests
- add a once-per-call deprecation-warning assertion for `weighted_difference`
- update characterization expectations that currently pin the old trace columns and old helper math

## Non-goals

- do not implement a new peptide fast-path approximation in this ticket
- do not treat the retracted `ctrl`-override claim as an acceptance criterion
