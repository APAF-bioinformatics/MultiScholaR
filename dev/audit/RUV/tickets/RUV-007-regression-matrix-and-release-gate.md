---
ticket_id: RUV-007
title: Regression matrix and release gate
status: pending
priority: P0
depends_on:
  - RUV-001
  - RUV-002
  - RUV-003
  - RUV-004
  - RUV-005
  - RUV-006
verified_at: 2026-04-26
source_plan: dev/audit/RUV/fix-plan.md
---

# RUV-007 Regression matrix and release gate

## Problem

The current test suite still contains characterization that encodes the known broken behavior, especially in metabolomics and lipidomics automatic mode. The RUV fixes are not complete until those tests are rewritten around the corrected contract and the full regression matrix passes.

## Current repo verification

- `tests/testthat/test-metab-norm-standalone-shared.R` currently asserts automatic `best_percentage == 12`.
- `tests/testthat/test-lipid-12b-norm-ruv-helpers-shared.R` currently asserts automatic `best_percentage == 15`.
- `tests/testthat/test-prot-04-design.R` still uses `weighted_difference` as a normal success-path expectation.
- `tests/testthat/test-prot-norm-optimization-direct-shared.R` still lacks plateau, tie, malformed-curve, and shuffled-row coverage.
- `tests/testthat/test-metab-03d-norm-module-characterization.R` is still largely wiring/stub coverage for automatic RUV rather than real-path optimizer coverage.

## Files in scope

- `tests/testthat/test-prot-06-ruv.R`
- `tests/testthat/test-prot-norm-optimization-direct-shared.R`
- `tests/testthat/test-prot-05b-norm-module-contracts.R`
- `tests/testthat/test-prot-05d-norm-module-characterization.R`
- `tests/testthat/test-prot-04-design.R`
- `tests/testthat/test-prot-04d-peptide-s4-normalization-characterization.R`
- `tests/testthat/test-metab-norm-standalone-shared.R`
- `tests/testthat/test-metab-03d-norm-module-characterization.R`
- `tests/testthat/test-lipid-12b-norm-ruv-helpers-shared.R`
- any new focused metabolomics/lipidomics RUV test files created to avoid overloading characterization suites

## Required changes

- implement the Phase 1 through Phase 7 test matrix from `fix-plan.md`
- replace bug-encoding assertions with corrected-contract assertions
- add deprecation-warning coverage for:
  - `findBestK()`
  - `ruvCancorFast()`
  - `weighted_difference` at public entry boundaries
- add real automatic-path smoke tests for metabolomics and lipidomics
- verify combined summary helpers still preserve their current primary columns
- rerun the normalization test matrix across proteomics, peptide, metabolomics, and lipidomics

## Acceptance criteria

- every required case listed in `fix-plan.md` Phase 7 is covered or explicitly superseded by an equal-or-stronger test
- characterization files no longer encode `percentage_max` as the automatic winner contract
- manual-mode tests stay green across all omics
- deprecation warnings fire once per public call, not once per inner iteration
- release-gate checklist items from the fix plan are satisfied before merge

## Release gate

- all new tests pass
- existing normalization suites pass
- invalid cancor input is not silently converted into `K = 1`
- metabolomics and lipidomics automatic mode no longer hardcode `percentage_max`
- peptide optimization no longer uses `ruvCancorFast()`
- roxygen-generated files, `NEWS.md`, `DESCRIPTION`, and workbook updates are present

## Non-goals

- do not preserve tests whose only purpose is to pin known-bug behavior
