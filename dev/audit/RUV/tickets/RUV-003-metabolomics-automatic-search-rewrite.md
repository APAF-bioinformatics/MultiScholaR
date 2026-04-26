---
ticket_id: RUV-003
title: Metabolomics automatic search rewrite
status: pending
priority: P0
depends_on:
  - RUV-001
verified_at: 2026-04-26
source_plan: dev/audit/RUV/fix-plan.md
---

# RUV-003 Metabolomics automatic search rewrite

## Problem

Metabolomics automatic RUV is still not a search. The live code allocates `percentage_range`, ignores it, always evaluates `percentage_max`, and returns `best_percentage = percentage_max` with a one-row trace.

## Current repo verification

- `R/func_metab_norm.R::runPerAssayRuvOptimization()` still creates `percentage_range <- seq(percentage_min, percentage_max, by = 1)` and never iterates it.
- automatic mode still calls `getNegCtrlMetabAnova(... percentage_as_neg_ctrl = percentage_max)`.
- automatic mode still calls `findBestK()` and `calculateSeparationScore()` exactly once per assay result path.
- return payload still hardcodes `best_percentage = percentage_max` and a one-row `optimization_results` table.
- `R/mod_metab_norm_server_helpers.R` still passes a full optimization parameter set, confirming the UI intent is a real search loop.
- `tests/testthat/test-metab-norm-standalone-shared.R` currently asserts `best_percentage == 12`, which encodes the broken behavior rather than the intended contract.

## Files in scope

- `R/func_metab_norm.R`
- `R/mod_metab_norm_server_helpers.R`
- `tests/testthat/test-metab-norm-standalone-shared.R`
- `tests/testthat/test-metab-03d-norm-module-characterization.R`
- additional focused metabolomics RUV test file if the existing characterization file is too stub-heavy

## Required changes

- replace the assay-local one-shot path with the percentage-first whole-object loop described in the fix plan
- call `getNegCtrlMetabAnova()` once per tested percentage, not once per assay-percentage
- call `ruvCancor()` once per tested percentage, not once per assay-percentage
- derive per-assay rows from the whole-object percentage results
- compute sample size from assay columns intersected with design-matrix sample IDs
- compute requested versus realized control counts explicitly
- isolate percentage-level and assay-level failures into trace rows instead of aborting the whole run
- switch K selection to `findBestKElbow()` and stable winner ordering
- warn once at the public entry boundary if `weighted_difference` is selected

## Acceptance criteria

- automatic mode evaluates the full requested percentage range
- each assay receives a full `optimization_results` trace with one row per tested percentage
- one failed percentage does not abort later percentages
- one failed assay at one percentage does not abort other assays
- whole-object helper failures at a given percentage become failed rows for each assay and the run continues
- each assay result exposes:
  - `best_k`
  - `best_percentage`
  - `best_realized_num_controls`
  - `best_realized_percentage`
  - `control_genes_index`
  - `cancor_plot`
  - `separation_score`
  - `composite_score`
  - `optimization_results`
- manual mode remains behaviorally unchanged

## Tests to update or add

- replace characterization assertions that `best_percentage` equals `percentage_max`
- add call-count tests proving whole-object helpers run once per percentage
- add per-assay trace-shape tests on a two-assay object
- add at least one end-to-end smoke test that runs the real metabolomics automatic path without stubbing away the optimizer

## Non-goals

- do not redesign the server-module orchestration beyond what is needed for the new result contract
- do not change manual-mode output semantics
