# KIT-101 W4: the per-file test baseline

## Method

- **Runner:** `tools/testing/test-baseline.R` (KIT-101), extending the series' one-process-per-file runner
  (`kit/scripts-078/run1.R`). Each test file runs in its own R process (`timeout 1800` per file) with a list reporter and a
  second reporter that records every result raised inside or outside `test_that()`. testthat reports a skip raised outside
  `test_that()` (a file-level guard) to the reporter only: the list reporter alone shows such a file as empty, with no skip.
- **Per file** (`test-baseline.tsv`, one row per file and load mode): static `test_that` blocks (grep), blocks run, hidden
  blocks (static minus run, when a file-level skip or error fired), expectations, passed, failed, error blocks, skipped
  blocks, warnings, the file-level skip with its reason and category, a file-level error, the skip results by category,
  the seconds, and the status (`no_result` when the process left no row: timeout, ceiling or crash).
- **Skip categories** (regex on the reason, first match): `binding_guard` (helper(s) / binding(s) not present or
  unavailable), `file_guard` (file(s) not present, Target-only, legacy duplicates present), `quarantine` (quarantined,
  Q-nn), `environment` (Linux, installed, capable, browser, On CRAN, not available, interactive), `fixture` (LFS pointer,
  snapshot, not found, corpus generator, source tree, fixture), `other`.
- **Header** (`header.tsv` of each run, repeated below): git SHA, R version, load mode, the number and sha256 of the sorted
  `R/` file list, library paths, `NOT_CRAN`, `MULTISCHOLAR_STRICT_SKIPS`.
- **Load modes.** `load_all`: `pkgload::load_all(export_all = TRUE, helpers = TRUE)` on the source tree. `check`: `R CMD
  build` and `R CMD check --no-tests --no-manual --no-build-vignettes` of a tarball of the tree (the regression check
  against the 3 WARNING / 3 NOTE baseline); then the tarball's `tests/` copied into `MultiScholaR.Rcheck/tests`, and each
  file run from `MultiScholaR.Rcheck/tests/testthat` with the check-installed package (`library(MultiScholaR)`,
  `test_file(package = "MultiScholaR", load_package = "installed")`, as `test_check` runs it), so `../..` is the check
  directory. `R CMD check` itself runs only `quality-runtime-contracts` (`tests/testthat.R`); the per-file run in check's
  layout is how every file's behaviour under check is recorded. Both with `NOT_CRAN=true`, as the series checkpoints run.
- **Settings:** 6 files in parallel (the ticket's number; see the note on KIT-099 in `findings.md`), ceiling 3 h per load
  mode, kit_100 library, R 4.6.1, no other heavy job alongside, a sleep lock held.

## Header (both runs)

| Field | Value |
|---|---|
| tree | `7591288f` (`bb87a463` plus KIT-101's test commits; `R/`, `inst/`, `DESCRIPTION` unchanged) |
| R | 4.6.1 (2026-06-24) |
| `R/` files | 385; sha256 of the sorted list `04e3c0f4c16b7025bd0a11f68ee7c09fc67cd4b949eefb37d4aaa8c360b24614` |
| libraries | load_all: kit_100 (`core/browser-smoke.2jjdqX/library`), user, system; check: `MultiScholaR.Rcheck`, then the same |
| `NOT_CRAN` / `MULTISCHOLAR_STRICT_SKIPS` | `true` / unset |
| files | 381 test files, 3,528 static `test_that` blocks |

## Results

| | load_all | check's test mode |
|---|---:|---:|
| run (wall) | 08:20 to 08:48 (1,688 s) | 08:50 to 09:07 (check about 3 min, files 900 s) |
| files with a row / `no_result` | 381 / 0 | 381 / 0 |
| blocks run / hidden | 3,555 / 4 | 2,844 / 715 |
| expectations passed | **33,672** | 21,748 |
| failed / error blocks | 3 / 5 | 152 / 410 |
| skipped blocks | 315 | 352 |
| file-level skips | 1 (Q-21, `test-prot-02o`, quarantine) | 11 (10 `R/`-path guards, Q-21) |
| file-level errors | 0 | 13 |
| skip results by category | file_guard 0, binding_guard 0, environment 52, quarantine 190, fixture 9, other 65 | file_guard 10, binding_guard 0, environment 51, quarantine 186, fixture 9, other 107 |
| `R CMD check --no-tests` | | **3 WARNINGs, 3 NOTEs**, the same six sections as KIT-099's `check.log` |

(Blocks run can exceed the static count: a few files define tests in loops.)

**load_all, the checkpoint.** 33,672 pass; the two known diagnostics error (`test-workflow-perf101-diagnostic.R`,
`test-workflow-ruv-pearson-diagnostic.R`); `test-gui-entry-proteomics.R` fails 3 and errors 3: the three protein-route
lanes have no frozen oracle (`entry-suite.md`; the helper now stops after the missing-oracle failure, so a later run shows
3 failures and no errors). Skips 315 = KIT-099's 314 plus the new trx e2e lane (no browser in the suite environment).
Against KIT-099's checkpoint (`kit-099/suite.txt`), every file's counts are unchanged except the files KIT-100 and
KIT-101 extended (`test-e2e-manifest.R`, `test-e2e-cross-omic.R`, `test-workflow-kit-harness.R`,
`test-workflow-omics-binding-steps.R`, `test-workflow-partitions.R`, `test-workflow-trx-declarations.R`) and the five new
files. Every guarded file runs in full (hidden 0): the sentinel edits changed nothing by default.

**Check's test mode: the `R/`-path guards skip by construction, with two exceptions.** From `MultiScholaR.Rcheck/tests/
testthat`, `../..` has no `R/`, so 10 of the 12 path-guarded files skip at file level (`file_guard`, 295 hidden blocks:
`lipid-00b` 101, `lipid-12` 69, `metab-00b` 4 + 45, `metab-01p` 26, `metab-01s` 35, `metab-02aa-shared` 7, `prot-05c`
4, `prot-07` 2, `prot-09` 2). `test-general-filemgmt-duplicate-resolution-contracts.R` (13 blocks) and
`test-prot-04-design.R` (182 blocks) do not skip: their first line, `repoRoot <- pkgload::pkg_path()`, errors outside a
source package before the guard runs, so under check they are file-level errors, not skips (`findings.md`). The 8
namespace-guarded files run in full under check (the internals are in the installed namespace). The other check-mode
failures are files that read the source tree (`R/` sources, `tools/`, `.github`, `docs/`) or the source package
(`pkgload`): 113 files in all (metab 31, workflow 29, lipid 13, prot 9, ci-impact 7, perf* 14, others). `R CMD check`
itself runs only `quality-runtime-contracts` (pass 3). The entry-point suite's frozen lanes pass in check mode as in
load_all: the oracles hold under both load modes.

## The rule for KIT-102 to 109

Compare each ticket's per-file baseline, per load mode, with this TSV. The acceptance fails on a new file-level skip, an
increase in hidden blocks, a drop in tests run (blocks run or expectations passed) for a file, or a new skip category in
a file, unless a quarantine entry (Q-23 onwards in `quarantined-tests.md`) names the internal the guard or seam test
covered and what replaced it. A family ticket runs its targeted files with `MULTISCHOLAR_STRICT_SKIPS=1`.
