# Merge train 28 of 2026-11-26: KIT-278 and KIT-269

**This train is NOT stacked.** Branched from `origin/v0.9-staging` at `623f7550fc6efcd33ddbb28ebeddfe1c9da29b2c` (fetched
before building; the newest staging head, the ledger commit VAL-344 corrected). It does not contain trains 25 to 27 (KIT-255,
KIT-271, KIT-270/KIT-246c/KIT-221: other agents', out), so KIT-255's test file is NOT in the scope below.

## Train

Base: `623f7550fc6efcd33ddbb28ebeddfe1c9da29b2c` (staging).

| Order | Ticket | Branch head | Merge commit | Ledger row | Audit |
|---|---|---|---|---|---|
| 1 | KIT-278 | `origin/kit-278` 993f41b9 | 98d19950 | `VAL-360` | AUDIT-KIT-278-impl-20261126, Closure check 1 (ready) |
| 2 | KIT-269 | `origin/kit-269` ac5606d8 | 9bfe6c4f | `VAL-361` | AUDIT-KIT-269-impl-20261126, Closure check 1, C5 (ready) |

Head (bound by the run): `9bfe6c4fe` (full hash in the request). This file is in the declaration commit after it, under `docs/` only.

### Conflicts met and their resolutions

- **KIT-278**: clean (the ledger row auto-merged at the table's end).
- **KIT-269**: conflicts in `NEWS.md` and `docs/methods-decision-ledger.md` only.
  - `NEWS.md`: staging's Fixes bullets first, then KIT-269's bullet; no marker left (+9 lines against the base).
  - Ledger: every staging row kept (to VAL-344); directly after it KIT-278's row (moved from the table's end, where its branch put
    it) as `VAL-360`, then KIT-269's as `VAL-361`. The diff of the ledger against the base is exactly those two added lines.
- No conflict under `R/`, `tests/`, `inst/`, `man/`, `DESCRIPTION` or `NAMESPACE`.

### The merged code against the base (git diff 623f7550 HEAD, outside docs/qa)

Exactly the two members' changes, and nothing else:
- `R/func_general_filemgmt_config_io_helpers.R` (KIT-269: `formatConfigList`)
- `R/func_general_filemgmt_study_parameter_writer.R` (KIT-269: the handler in the writer)
- `tests/testthat/test-kit269-config-vector-values.R` (KIT-269, new)
- `tests/testthat/test-omics-module-globals.R` (KIT-278)
- `NEWS.md`, `docs/methods-decision-ledger.md` (the two merges' text), and files under `docs/qa/workflow-execution/kit/kit-269/` and `kit-278/`.

No file under `inst/`, `man/`, `DESCRIPTION` or `NAMESPACE` differs. Code line counts against the base (`git diff --numstat`): `R/` +9/-5 over two files, tests +129/-1 over two files.

**The merged vector branch of `formatConfigList`** (`R/func_general_filemgmt_config_io_helpers.R:482-487`):
```
if (length(value) > 1) {
    # KIT-269: every vector value in full; a named vector as name=value, the record path's form (REC-072, REC-001).
    text <- as.character(value)
    if (!is.null(names(value))) text <- paste0(names(value), "=", text)
    paste(text, collapse = ", ")
} else {
```
**The merged allowance** (`tests/testthat/test-omics-module-globals.R:60`):
```
hits <- omics_global_hits(sources, "<<-", allow = "shown_idle <<- idle")
```

### Path lengths

The train's own files are short (`docs/qa/workflow-execution/kit/merge-train-20261126u*`); the members' record paths are theirs, as audited.

## Declaration (written before the run is sent)

**Purpose.** The merge check of KIT-278 and KIT-269 on the merged tree, base against head; the full suite is not repeated (next checkpoint).

**Scope recount: 49 test files**, as AUDIT-KIT-269 Closure check 1, C5 asks.
- KIT-269's grep (`git grep -lE "createWorkflowArgsFromConfig|createWorkflowArgsFn|saveWorkflowSummaryParameters|study_parameters|Parameters from Final S4 Object|Additional Configuration Parameters|formatConfigList|createStudyParametersFile" -- tests/testthat ':!tests/testthat/helper-*' ':!data/benchmark'`)
  on the merged head: 47 files, that is 46 on the base (including `test-kit213-study-parameters-labels.R`) plus the new KIT-269 file.
- Added: `test-gui-entry-proteomics.R` (not in that grep) and `test-omics-module-globals.R`. 47 + 2 = **49**.
- No train carrying KIT-255 is on this base, so `test-kit255-omics-correlation-assays.R` is not in scope and the base values of
  `test-kit254-omics-correlation-skip.R` (379/0) and `test-kit236-report-sentences.R` (289/0) are train 24's/q's head values as committed.

**Run settings.** 3 shards per side; `file_timeout_minutes` 45, `file_timeouts` `test-gui-entry-proteomics.R` 50; `blas` reference; `check` true
on both sides; no `env`; no `extra_packages`; no scripts. Base `623f7550`, head `9bfe6c4f`; one run.

**Ceiling.** 600 job-minutes = 40.0 core-hours (core-hours = job-minutes x 4 / 60). Expected about 315 job-minutes (21.0 core-hours): train 24's
307.8 for 47 files plus about 0.25 for the two short files.

### Predictions (AUDIT-KIT-269 C5; base counts checked against the newest committed tables, trains p and q)

Counts are passed/failed/errors/skipped. The train's own base run is the reference where a base count differs from a figure here; that is not a
result, a head that differs from its base is. For `test-gui-entry-proteomics.R` the committed value (92/101, train 23's head, merge-train-20261126p) predates KIT-265's landing: the
base run of this train is its reference, and the prediction is `same` with identical failures.

| File | Base | Head | Verdict |
|---|---|---|---|
| test-kit269-config-vector-values.R | absent | 34/0/0/0 | only-on-head |
| test-omics-module-globals.R | 13/1/0/0 | 14/0/0/0 | fixed |
| test-kit213-study-parameters-labels.R | 56/0/0/0 | 56/0/0/0 | same (the G3 log line prints a different digest on head: a log line, not an assertion) |
| test-integration-record.R | 200/0/0/0 | same | same |
| test-workflow-figure-paths.R | 181/0/0/0 | same | same (a base-only failure quoting a folder name containing 259 is the known flake: recorded, not re-run) |
| test-gui-entry-proteomics.R | 92/101/0/0 (see above) | as base | same, failures identical |
| test-omic-summary-module.R | 1675/28/0/0 | as base | same, failures identical |
| test-omics-summary-core.R | 676/25/0/2 | as base | same, failures identical |
| test-module-ci-trx-summary.R | 87/1/0/0 | as base | same, failure identical |
| test-kit254-omics-correlation-skip.R 379/0; test-kit236-report-sentences.R 289/0 | as listed | as base | same |
| test-kit263-prot-study-parameters-sections.R 11/0; test-kit235-lipid-metab-report-parser.R 54/0; test-general-filemgmt-config-direct-shared.R 33/0; test-general-filemgmt-config-roundtrip-contracts.R 24/0; test-general-filemgmt-report-helpers-shared.R 85/0/0/1; test-general-filemgmt-s4-contracts.R 30/0 | as listed | as base | same |
| test-workflow-native-session.R 734/0/0/2; test-workflow-analysis-record.R 351/0; test-kit212-da-current.R 309/0; test-kit207-summary-head-guard.R 292/0; test-kit243-summary-current-outputs.R 246/0; test-workflow-integration-input.R 241/0; test-workflow-omics-application-reads.R 230/0; test-integration-input-lipid-metab.R 208/0; test-kit248-correlation-skip-record.R 206/0; test-kit227-report-ruv-sentences.R 162/0; test-module-ci-metab-da.R 138/0; test-kit206-prot-norm-export-guard.R 116/0; test-kit172-downstream.R 84/0; test-integration-input-transcriptomics.R 79/0; test-kit224-ruv-feature-levels.R 69/0 | as listed | as base | same |
| test-module-ci-prot-summary.R 20/0/0/5; test-module-ci-lipid-summary.R 14/0/0/6; test-module-ci-metab-summary.R 12/0/0/6; test-workflow-memory-mechanism-member.R 0/0/0/13; test-integration-record-mofa.R 0/0/0/1; test-integration-workbook-mofa.R 0/0/0/3; the six test-e2e-*.R (skips 2, 3, 1, 1, 2, 3); the five 0/0/0/0 files (lipid-14a, lipid-15, metab-05, prot-12, prot-12a) | as listed | as base | same |
| every other file in scope | as base | as base | same |
| R CMD check | 3 WARNING, 3 NOTE | | identical, nothing only on head |

### Scope (49 files)

- test-e2e-lipidomics-canonical.R
- test-e2e-metabolomics-lc-gc.R
- test-e2e-proteomics-dia-limpa.R
- test-e2e-proteomics-dia.R
- test-e2e-proteomics-enrichment.R
- test-e2e-proteomics-tmt-lfq.R
- test-general-filemgmt-config-direct-shared.R
- test-general-filemgmt-config-roundtrip-contracts.R
- test-general-filemgmt-report-helpers-shared.R
- test-general-filemgmt-s4-contracts.R
- test-gui-entry-proteomics.R
- test-integration-input-lipid-metab.R
- test-integration-input-transcriptomics.R
- test-integration-record-mofa.R
- test-integration-record.R
- test-integration-workbook-mofa.R
- test-kit172-downstream.R
- test-kit206-prot-norm-export-guard.R
- test-kit207-summary-head-guard.R
- test-kit212-da-current.R
- test-kit213-study-parameters-labels.R
- test-kit224-ruv-feature-levels.R
- test-kit227-report-ruv-sentences.R
- test-kit235-lipid-metab-report-parser.R
- test-kit236-report-sentences.R
- test-kit243-summary-current-outputs.R
- test-kit248-correlation-skip-record.R
- test-kit254-omics-correlation-skip.R
- test-kit263-prot-study-parameters-sections.R
- test-kit269-config-vector-values.R
- test-lipid-14a-summary-module-characterization.R
- test-lipid-15-summary-module-seams.R
- test-metab-05-summary-module-characterization.R
- test-module-ci-lipid-summary.R
- test-module-ci-metab-da.R
- test-module-ci-metab-summary.R
- test-module-ci-prot-summary.R
- test-module-ci-trx-summary.R
- test-omics-module-globals.R
- test-omics-summary-core.R
- test-omic-summary-module.R
- test-prot-12a-summary-module-characterization.R
- test-prot-12-summary-module-contracts.R
- test-workflow-analysis-record.R
- test-workflow-figure-paths.R
- test-workflow-integration-input.R
- test-workflow-memory-mechanism-member.R
- test-workflow-native-session.R
- test-workflow-omics-application-reads.R

### Decision and attempts

- Every prediction holds: TRAIN GREEN for both members. Any other outcome is a result: recorded with its files and counts, not repaired, not re-run.
- One run. A re-send only for an instrument defect with no result judged; an HTTP 500/502 at dispatch with no run created is started once by hand.

## What the coordinator resolves at landing

- Ledger: `VAL-360` (KIT-278) and `VAL-361` (KIT-269) take their numbers (`docs/methods-decision-ledger.md`, the two lines after VAL-344); status cells "merged (<commit>), train 28".
- Placeholders in the members' own records are not edited here (listed in the hand-back).

## Results (run 45730995940, head 9bfe6c4f, base 623f7550; request commit e399a132)

Run success; every file ran to completion (no timeout, no crash). `compare.tsv` and `check-compare.txt` are committed beside this file
(`merge-train-20261126u-compare.tsv`, `merge-train-20261126u-check-compare.txt`). Base versus head: same 47, new-failure 0, message-changed 0,
count-changed 0, fixed 1, only-on-head 1. **Clean under the rule.** Head: passed 7411, failed 160; base: passed 7376, failed 161.

### Verdict table against the predictions (passed/failed/errors/skipped)

| File | Base | Head | Verdict | Predicted |
|---|---|---|---|---|
| test-kit269-config-vector-values.R | absent | 34/0/0/0 | only-on-head | held |
| test-omics-module-globals.R | 13/1/0/0 | 14/0/0/0 | fixed | held (base failure `"mod_omic_norm.R"`) |
| test-kit213-study-parameters-labels.R | 56/0/0/0 | 56/0/0/0 | same | held |
| test-integration-record.R | 200/0/0/0 | 200/0/0/0 | same | held |
| test-workflow-figure-paths.R | 181/0/0/0 | 181/0/0/0 | same | held; no 259-folder flake appeared |
| test-gui-entry-proteomics.R | 92/101/0/0 | 92/101/0/0 | same, failures identical | held (also at the carried count) |
| test-omic-summary-module.R | 1675/28/0/0 | same | same | held |
| test-omics-summary-core.R | 676/25/0/2 | same | same | held |
| test-module-ci-trx-summary.R | **83/2/0/0** | 83/2/0/0 | same | verdict held; **base count missed** (predicted 87/1) |
| test-integration-input-transcriptomics.R | **67/3/0/0** | 67/3/0/0 | same | verdict held; **base count missed** (predicted 79/0, status ok) |
| test-kit224-ruv-feature-levels.R | **65/1/0/0** | 65/1/0/0 | same | verdict held; **base count missed** (predicted 69/0, status ok) |
| every other file in scope (43 files) | as predicted | as base | same | held |

Missed predictions (recorded, not repaired, not re-run): the three base counts above. They are identical on both sides, so they are not caused by
the two members. Their failures are one cause: `trxDeclReference("trx_edge_12")` in `tests/testthat/helper-workflow-trx-declarations.R:34`, expected
manifest digest `e17d1ed5...4fc4344`, actual `faa19db7...a914fa797f` (the `corpus_manifest_sha256` check), and in
`test-module-ci-trx-summary.R` also the `ruviii_newY` product digest (`1c02bfdd...` against `a1f43e2d...`). Train q's committed tables hold the same three
files passing (79/0, 69/0, 87/1), so the stored trx reference digest changed between those runs and this one on the same base code or the runner's
snapshot library; this run did not investigate. Belongs to the transcriptomics declarations (helper-workflow-trx-declarations), not to KIT-269 or KIT-278.

### R CMD check

Base and head both 0 ERROR, 3 WARNING, 3 NOTE; findings identical (6); only on head: none; only on base: none. Prediction held.

### Budget

Job-minutes (jobs API, summed): 301.9 (tests head 46.08 + 60.63 + 27.82; tests base 35.25 + 41.85 + 66.30; check base 13.28, head 9.53; plan 0.37; results 0.78).
Core-hours = 301.9 x 4 / 60 = **20.13**. Ceiling 600 job-minutes (40.0 core-hours); expected about 315 (21.0).

### Decision

Every prediction about the two members holds: TRAIN GREEN for KIT-278 and KIT-269. The three base-count misses are a recorded result outside the members.

## Merge record (coordinator)

Landed on v0.9-staging on 2026-11-26 by a no-fast-forward merge of this branch at f09083cf onto ec27f69e (the train was built on 623f7550; since then staging gained merge trains 26, 27 and 25 (KIT-271; KIT-270, KIT-246c, KIT-221; KIT-255), KIT-211's documents and ledger rows).

- Verdict: green. Run 45730995940 (49 files and R CMD check; 301.9 of 600 job-minutes, 20.13 core-hours): the new KIT-269 file 34/0 only on head; test-omics-module-globals.R 13/1 to 14/0 fixed; 47 files same; check identical.
- Three base counts of the declaration were missed, identically on both sides (test-module-ci-trx-summary.R 83/2 for 87/1; test-integration-input-transcriptomics.R 67/3 for 79/0; test-kit224-ruv-feature-levels.R 65/1 for 69/0), all on the transcriptomics corpus-manifest digest. Cause, by comparing requests: this train's request set no environment, while the trains whose tables show those files passing set MULTISCHOLAR_TRX_MANIFEST_NOTE=1 (the same effect was seen in KIT-242's pre-run 45723354055). A property of the request, not of the code; verdict same on all three. Later train requests that reach transcriptomics files set that variable.
- Conflicts at landing: NEWS (both kept: staging's bullets first, then KIT-269's) and the ledger. R/func_general_filemgmt_study_parameter_writer.R merged by itself: staging's side (KIT-255) adds lines near 371 and 571, this train's change is at 991; read at landing.
- Ledger at landing: staging's rows kept through VAL-359; KIT-278 = VAL-360; KIT-269 = VAL-361. Placeholders in this record replaced.
- Not covered by the train run: the combination with merge trains 25 to 27 (KIT-255's test of the study parameters among them). The next checkpoint covers it.
- With KIT-278 on staging, test-omics-module-globals.R reads 14/0 there; checkpoint 2's one unpredicted difference is closed.
