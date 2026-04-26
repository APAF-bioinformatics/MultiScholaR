---
ticket_id: RUV-004
title: Lipidomics automatic search rewrite
status: pending
priority: P0
depends_on:
  - RUV-001
verified_at: 2026-04-26
source_plan: dev/audit/RUV/fix-plan.md
---

# RUV-004 Lipidomics automatic search rewrite

## Problem

Lipidomics automatic RUV mirrors the metabolomics defect: the live code allocates a range, ignores it, always uses `percentage_max`, and returns one-row automatic results.

## Current repo verification

- `R/func_lipid_norm_ruv_helpers.R::runLipidPerAssayRuvOptimization()` still creates `percentage_range` but never iterates it.
- automatic mode still calls `getNegCtrlMetabAnova(... percentage_as_neg_ctrl = percentage_max)` and returns `best_percentage = percentage_max`.
- automatic mode still uses `findBestK()` and `calculateSeparationScore()` directly.
- `tests/testthat/test-lipid-12b-norm-ruv-helpers-shared.R` currently asserts automatic `best_percentage == 15`, which encodes the bug.
- lipid helper/observer code still expects the result objects to be renderable in combined tables, so compatibility of existing summary columns matters.

## Files in scope

- `R/func_lipid_norm_ruv_helpers.R`
- `R/mod_lipid_norm_observer_helpers.R` if compatibility shims are needed for richer traces
- `R/mod_lipid_norm_workflow_helpers.R` if result-shape assumptions need adjustment
- `tests/testthat/test-lipid-12b-norm-ruv-helpers-shared.R`
- additional focused lipidomics RUV test file if needed

## Required changes

- replace the one-shot automatic path with a percentage-first whole-object search
- evaluate whole-object helpers once per percentage
- derive per-assay result rows from shared percentage evaluations
- compute sample size from assay columns aligned to design-matrix samples
- record requested versus realized control counts
- isolate failed percentages and failed assays into retained trace rows
- move automatic K selection onto `findBestKElbow()` and stable winner ordering
- warn once at the public entry boundary when `weighted_difference` is selected

## Acceptance criteria

- automatic mode tests the full requested percentage range
- `optimization_results` has one row per tested percentage for each assay
- skipped and failed candidate rows are preserved with `status` and `error_reason`
- whole-object failures at one percentage do not abort later percentages
- final assay result objects expose the same top-level contract defined for metabolomics
- combined table helpers keep their existing primary columns and names
- manual mode remains unchanged

## Tests to update or add

- replace assertions that automatic `best_percentage` equals `percentage_max`
- add call-count coverage for once-per-percentage whole-object helper calls
- add per-assay trace tests on a multi-assay lipid object
- add at least one real-path automatic smoke test

## Non-goals

- do not redesign lipid UI workflow steps beyond compatibility updates required by the new result schema
