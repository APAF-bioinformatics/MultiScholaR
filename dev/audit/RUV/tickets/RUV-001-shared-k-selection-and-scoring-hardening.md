---
ticket_id: RUV-001
title: Shared K-selection and scoring hardening
status: pending
priority: P0
depends_on: []
verified_at: 2026-04-26
source_plan: dev/audit/RUV/fix-plan.md
---

# RUV-001 Shared K-selection and scoring hardening

## Problem

The shared automatic-RUV math contract is still wrong in the live repo. `findBestK()` is argmax rather than plateau detection, can return multiple K values on ties, assumes row-order alignment, and feeds a composite-score helper that still divides by zero when `max_acceptable_k == 1`.

## Current repo verification

- `R/func_prot_norm_optimization_helpers.R` still defines `findBestK()` as positional `All - Control` subtraction plus `which(...)`.
- `findBestKForAssayList()` still masks vector-return behavior by silently coercing non-scalar results to `NA_integer_`.
- `calculateSeparationScore()` still subtracts by row position and still treats `weighted_difference` as a normal metric.
- `calculateCompositeScore()` still uses `(max_acceptable_k - 1)` without a dedicated `== 1` guard.
- `R/func_pept_s4_norm_methods.R` still duplicates the same separation/composite logic in `.peptide_calculateSeparationScore()` and `.peptide_calculateCompositeScore()`.
- `tests/testthat/test-prot-norm-optimization-direct-shared.R` still only covers the easy current-path behavior and does not enforce plateau, tie, invalid-input, or shuffled-row semantics.

## Files in scope

- `R/func_prot_norm_optimization_helpers.R`
- `R/func_pept_s4_norm_methods.R`
- `tests/testthat/test-prot-06-ruv.R`
- `tests/testthat/test-prot-norm-optimization-direct-shared.R`
- generated `NAMESPACE`
- generated `man/` entries affected by roxygen

## Required changes

- Add a strict internal curve-extraction helper that aligns `"All"` and `"Control"` by `K`, rejects malformed input, and treats duplicate `K` rows within a `featureset` as invalid.
- Add exported `findBestKElbow()` with the first-plateau contract from the fix plan.
- Convert exported `findBestK()` into a compatibility wrapper that warns and delegates to `findBestKElbow()`.
- Rewrite `calculateSeparationScore()` to consume aligned curve data instead of positional subtraction.
- Harden `calculateCompositeScore()` and the peptide mirror/helper path for non-finite input and `max_acceptable_k == 1`.
- Update internal package call sites to use `findBestKElbow()` directly so package code does not trigger its own deprecation warnings.
- Simplify `findBestKForAssayList()` so invalid input stays `NA_integer_`, but valid curves never degrade into silent vector-length handling.

## Acceptance criteria

- `findBestKElbow()` returns:
  - smallest plateau `K`
  - smallest tied `K`
  - `1L` for valid weak-effect curves
  - `NA_integer_` for malformed or empty aligned curves
- duplicate `K` rows within a single `featureset` are treated as invalid input
- row order no longer changes `best_k` or score results
- `calculateCompositeScore()` never divides by zero and returns `NA_real_` for non-finite inputs
- peptide helper math delegates to or exactly mirrors the shared contract
- internal code paths in proteomics and peptide stop calling the deprecated wrapper directly (metabolomics and lipidomics call sites are migrated as part of their respective automatic-mode rewrites in RUV-003 and RUV-004)
- roxygen regeneration updates `NAMESPACE` and any affected `man/` files without manual `NAMESPACE` edits

## Tests to update or add

- extend `tests/testthat/test-prot-06-ruv.R` with the Phase 1 matrix from the fix plan
- update `tests/testthat/test-prot-norm-optimization-direct-shared.R` so it no longer encodes the old helper semantics as the contract

## Non-goals

- do not change manual-mode semantics
- do not expose `epsilon` or `min_effect` in the UI
- do not re-open the retracted peptide `ctrl`-override claim
