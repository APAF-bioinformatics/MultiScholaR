---
ticket_id: RUV-002
title: Proteomics automatic optimizer hardening
status: pending
priority: P0
depends_on:
  - RUV-001
verified_at: 2026-04-26
source_plan: dev/audit/RUV/fix-plan.md
---

# RUV-002 Proteomics automatic optimizer hardening

## Problem

Proteomics is the only omics path that actually iterates percentages today, but it still sits on the old K-selection contract, uses implicit `which.max()` selection, and returns a thin optimization trace that does not record realized control counts or failure reasons.

## Current repo verification

- `R/func_prot_norm_optimization_helpers.R::findBestNegCtrlPercentage()` still uses `findBestK()` rather than a hardened scalar elbow helper.
- winner selection still relies on `which.max()` over `composite_score`.
- `optimization_results` still uses the old shape: `percentage`, `separation_score`, `best_k`, `composite_score`, `num_controls`, `valid_plot`.
- `weighted_difference` still passes through as a normal metric with no once-per-call deprecation warning.
- downstream tests still read old columns such as `optimization_results$num_controls` and `optimization_results$percentage`.

## Files in scope

- `R/func_prot_norm_optimization_helpers.R`
- `R/mod_prot_norm_ruv_helpers.R` if compatibility shims or richer display tables are needed
- `tests/testthat/test-prot-05b-norm-module-contracts.R`
- `tests/testthat/test-prot-05d-norm-module-characterization.R`
- `tests/testthat/test-prot-norm-optimization-direct-shared.R`

## Required changes

The `findBestK()` → `findBestKElbow()` call-site migration in `findBestNegCtrlPercentage()` is handled by RUV-001. This ticket focuses on the optimizer's consumer behavior: deterministic tie-breaking, trace schema, and deprecation warnings.

- replace implicit `which.max()` selection with explicit stable ordering:
  1. highest `composite_score`
  2. lowest `best_k`
  3. highest `separation_score`
  4. lowest `percentage_requested`
- emit the `weighted_difference` deprecation warning once per public optimizer call
- expand `optimization_results` so every tested percentage produces one row with status and error diagnostics
- record both requested and realized control quantities
- preserve enough compatibility for current consumers that still read `best_k`, control vectors, and top-level summary values

## Acceptance criteria

- `findBestNegCtrlPercentage()` always returns scalar `best_k`
- `optimization_results` includes one row per tested percentage, including skipped or failed candidates
- each row includes:
  - `percentage_requested`
  - `candidate_feature_count`
  - `realized_num_controls`
  - `realized_percentage`
  - `sample_size`
  - `best_k`
  - `separation_score`
  - `composite_score`
  - `status`
  - `error_reason`
- ties are broken deterministically by the stable ordering from the fix plan
- `weighted_difference` warns once per call and not once per percentage
- manual proteomics RUV behavior is unchanged
- if a deprecated compatibility alias such as `num_controls` is retained for one cycle, it must be clearly secondary to `realized_num_controls`

## Tests to update or add

- update proteomics optimizer tests so they assert one-row-per-percentage trace behavior
- add explicit tie-break tests
- add a once-per-call deprecation-warning assertion
- update contract tests that currently assume only the old trace columns exist

## Non-goals

- do not change manual-mode control selection or application semantics
- do not add UI parameters for elbow tuning
