# KIT-101 outcome: the tests work (W2 to W5)

The test half of KIT-101, on `kit-101-tests` from `bb87a463`. `design.md` (W1) and W6 are the other half. Nothing under
`R/`, `inst/` or `DESCRIPTION` changed (`git diff --stat bb87a463 -- R inst DESCRIPTION` is empty). Library: kit_100
(`core/browser-smoke.2jjdqX`), R 4.6.1.

## Test first (the pins)

- `test-e2e-manifest.R` pinned to ten lanes, the `trx_screen` filter `^e2e-transcriptomics-screen$`, its fixture
  directory and its per-file sha256, before the lane existed: **pass 105, fail 16, error 1** (four tests). After the
  fixture and the manifest entry: **pass 131, fail 0**.
- `test-module-ci-manifest.R` pinned to MCI-026.1 to MCI-032.1 on the transcriptomics lane before the registration:
  **pass 32, fail 1** (the transcriptomics pin). After registering MCI-030.1, MCI-031.1 and MCI-032.1: **pass 33, fail 0**.
- `test-e2e-fixture-integrity.R` failed once the trx lane named `transcriptomics_report.rmd` (no stub); after the stub
  (`tests/testdata/e2e/report_templates/transcriptomics_report.rmd`, the other lanes' stub form): pass.

## Simulation Setup: the strict-skip sentinel

A scratch copy of `R/` and `tests/testthat/` (scratchpad `kit101/strict-sim`), the package loaded from the worktree,
each file run in its own process with the baseline runner:

| Case | `MULTISCHOLAR_STRICT_SKIPS` | Result |
|---|---|---|
| worktree, `test-metab-02aa` (path guard intact) | unset | pass 43 |
| worktree, `test-metab-02aa` | 1 | pass 43 |
| scratch, `R/mod_metab_da_display_helpers.R` moved | unset | pass 0, skip 0: the file vanishes (testthat reports a top-level skip only to the reporter; the baseline runner records it as a file-level skip, 7 hidden blocks) |
| scratch, same | 1 | **error 1**: `strict skip (MULTISCHOLAR_STRICT_SKIPS=1): Target-only metab DA display file(s) not present` |
| scratch, partial load, `R/mod_metab_qc_s4.R` moved | unset | pass 163, error 6: the file is silently not sourced |
| scratch, same | 1 | **error 1**: `strict skip ...: partial load dropped mod_metab_qc_s4.R` |
| worktree, `test-lipid-01a` (namespace guard intact) | 1 | pass 5 |
| scratch, `test-lipid-01a`'s guard naming an internal not in the namespace | unset | file-level skip ("binding(s) unavailable") |
| scratch, same | 1 | **error 1**: `strict skip ...` |

Default behaviour is unchanged: with the variable unset each guard calls `testthat::skip()` with the same message, and
`strictSkipPresent()` returns `paths[file.exists(paths)]`.

## W2: fixtures and the UniProt double

- `trx_screen` is written by `generate_fixtures.R` through the KIT-040 corpus tool (`kitOmicsRecipe("trx_screen", 3,
  300)`); its six files are byte-identical to the module-CI fixture `tests/testdata/mci/fx/transcriptomics/trx_screen/`
  (same sha256), and the manifest entry carries each file's sha256. Running the generator also showed that it does not
  reproduce the committed metab and lipid fixtures (`findings.md`); those files were restored, not regenerated.
- The offline UniProt double is a cache file, not a mock: `getUniprotAnnotations()` reads
  `<results_dir>/cache/uniprot_annotations/uniprot_annotations.RDS` before any request, so the lane places the committed
  cache there. The design step logged "UniProt annotations retrieved successfully. Found 6 annotations" in every
  proteomics lane, with the HTTP(S) proxy pointed at a closed port.

## W3: the entry-point suite

`entry-suite.md`: 10 lanes, 7 frozen (determinism over runs 3, 4 and the freeze run after the one correction pass),
prot_tmt, prot_lfq and prot_lfq_fragpipe not frozen (the store file count varies); mutation fails at step `da`; an
erroring step fails, not skips; the suite passes under the strict sentinel.

## W4: the baseline

`test-baseline.md` and `test-baseline.tsv` (762 rows: 381 files x 2 load modes). load_all: 33,672 pass, the two known
diagnostics' errors, and `test-gui-entry-proteomics.R`'s three unfrozen lanes (3 fail, 3 error; after the helper's
missing-oracle fix, 3 fail); skips 315. Check's test mode: 10 of the 12 `R/`-path-guarded files skip by construction,
2 error first (`pkgload::pkg_path()`); `R CMD check --no-tests`: 3 WARNINGs, 3 NOTEs (baseline).

## W5: the trx e2e lane, once

`test-e2e-transcriptomics-screen.R` (manifest lane `trx_screen`, filter `^e2e-transcriptomics-screen$`), run once
the KIT-097 way (`CHROMOTE_CHROME=/usr/bin/brave`, the app process on the kit_100 library, the test process from source,
`timeout 1800`): **pass 14, fail 0, error 0, skip 0, 75 s** (START 07:52:40, tree `7f3ec8f9`, `R/` unchanged). The
harness and its failure artifacts are the KIT-097-proven ones; the lane's own assertions are the analysis record's
steps, the DA's source state and model, the RUVIII `M`, the manifest's three expected exports, the report and the
session export.
