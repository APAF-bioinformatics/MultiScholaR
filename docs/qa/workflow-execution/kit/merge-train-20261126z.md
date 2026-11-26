# Merge train 33 of 2026-11-26: KIT-274

**Not stacked.** Branched from `origin/v0.9-staging` at `221af09e87ef23849a717da8cea63b57907f687b` (fetched before building; the
newest staging head, ledger VAL-362). It contains trains 25 to 28 and KIT-278 (so `test-omics-module-globals.R` reads 14/0 on the base).

## Train

Base: `221af09e87ef23849a717da8cea63b57907f687b` (staging).

| Order | Ticket | Branch head | Merge commit | Ledger row | Audit |
|---|---|---|---|---|---|
| 1 | KIT-274 | `origin/kit-274` 4c5ab27f | 6861b213 | `VAL-370` (the coordinator numbers it at landing) | AUDIT-KIT-274-impl-20261126 (no Closure check section was present when the train was built) |

Head (bound by the run): `b6a3f6b798e8843caafc2f9cc04e0495ef848a3f` = the merge 6861b213 plus one commit that moves the ledger row. This file
and the later commits are under `docs/` only.

### The tip contains the fix (the two tree checks, as asked, quoted)

`git diff 4c5ab27f edb15841 --stat -- R inst` shows **25 files, 80 insertions, 636 deletions**, not only the fix's four files. The 80 insertions are
the fix; the 636 deletions and 21 further files are the difference between edb15841 (the revert, whose tree was the older staging 31a4d684 at
`ba2ad567`) and the tip, which has merged newer staging (6bbdf6da; KIT-255, KIT-269, KIT-279 and others). So the check as worded does not
come out to "exactly the fix's files"; it carries the staging movement. The fix itself is in the tip: `git diff ec27f69e 4c5ab27f --stat -- R inst`
(ec27f69e = the staging commit the branch last merged) gives exactly `R/utils_workflow_analysis_record.R` +28, `R/utils_workflow_summary.R` +28/-,
`R/utils_workflow_trx_declarations.R` +31/-, `inst/reports/transcriptomics/transcriptomics_report.rmd` 3 lines: 4 files, 80 insertions, 10 deletions;
`trxAssertRuviiiExportApplied` is defined in the tip's R files.

`git diff origin/v0.9-staging 4c5ab27f --stat -- R inst` (staging 221af09e) shows **6 files**: the four above plus
`R/func_general_filemgmt_config_io_helpers.R` (9 lines) and `R/func_general_filemgmt_study_parameter_writer.R` (5 lines). The two extra files are
KIT-269's change, which landed on staging after the tip's last staging merge (ec27f69e); in this two-dot diff they appear reversed. The merged
tree (below) does not lose them.

Decision: the tip does contain the fix; the wording's "nothing else" is met on the merge, not on the two-dot diffs. I proceeded and say so here.

### Conflicts

None. `git merge --no-ff origin/kit-274` completed with no conflict (NEWS and the ledger auto-merged). One edit after the merge: the ledger row had landed
at the end of the file (`:2571`, after the KIT-157 row); it is moved directly after VAL-362 (`:2535`) and its id written `VAL-370`
(commit b6a3f6b7). The ledger diff against the base is exactly that one added row.

### The merged tree against the base (`git diff 221af09e HEAD`, outside docs/qa)

Exactly KIT-274's change: `NEWS.md` (+8), `R/utils_workflow_analysis_record.R`, `R/utils_workflow_summary.R`, `R/utils_workflow_trx_declarations.R`,
`inst/reports/transcriptomics/transcriptomics_report.rmd`, `tests/testthat/test-kit274-trx-no-ruv-export.R` (new, 512 lines), `docs/methods-decision-ledger.md` (+1 row), and
files under `docs/qa/workflow-execution/kit/kit-274/`. Nothing under `man/`, `NAMESPACE` or `DESCRIPTION`. `R/`, `inst/`, `tests/`: 5 files, +592/-10.
The two files KIT-269 changed are not in the diff (the merge kept staging's version).

## Declaration (written before the run is sent)

**Purpose.** The merge check of KIT-274 on the merged tree, base against head. The full suite is not repeated.

**Scope: 37 files.** The ticket's grep (`git grep -lE "trx_ruviii_export|runTrxReferenceRuviiiExport|acceptTrxRuviiiExport|trxRuviiiExportParameters|run_ruviii_export|trxIntegrationRuvState|ruviii-text|trxSummaryIntegrationRequires|workflowRecordExportIsCurrent|replayWorkflowAnalysisRecord|replayWorkflowRecord|transcriptomics_report" -- tests/testthat ':!tests/testthat/helper-*' ':!data/benchmark'`) recounted on the merged tree: **36 files** (the audit's number), plus `test-omics-module-globals.R` = 37. The list equals R2's 37 exactly.

**Settings** (copied from R2's request, `kit274-r2-87106bc3.json`, which DID set the variable): 3 shards per side; `file_timeout_minutes` 45;
`file_timeouts` `test-prot-norm-core.R` 120 and `test-workflow-saved-da.R` 90 (as R2; neither file is in scope); `extra_packages` limpa, mixOmics, sva;
`env` `MULTISCHOLAR_TRX_MANIFEST_NOTE=1`; `check` true on both sides; `blas` reference; `os` ubuntu-22.04; no scripts. Base `221af09e`, head `b6a3f6b7`; one run.
Because R2 already set the variable, the base counts below (from R2's tables) are with it set, and nothing changes against the audit's counts.

**Ceiling.** 660 job-minutes = 44.0 core-hours (core-hours = job-minutes x 4 / 60); expected about 350 job-minutes (R2: 352.32 for the same files).

### Predictions (passed/failed/errors/skipped; the run's own base side is the reference wherever a base count differs from a figure here)

Counts WITH the variable set, from R2's committed base table (`kit-274/r2-tests-base.tsv`, variable set) and the newest train tables that set it
(train r: kit-harness, operation-registry, kit129, binding-steps; train q: ci-trx-summary).

| File | Base | Head | Verdict |
|---|---|---|---|
| test-kit274-trx-no-ruv-export.R | absent | 84/0/0/0 (T11 restored; R2b 45739818528: base 55/25, head 84/0) | only-on-head |
| test-kit264-trx-export-current.R | 85/0 | 85/0 | same |
| test-omics-module-globals.R | 14/0 (KIT-278 is on staging) | 14/0 | same |
| test-omic-summary-module.R | 1675/28 | same, failures identical | same |
| test-omics-summary-core.R | 676/25/0/2 | same, failures identical | same |
| test-workflow-trx-declarations.R | 313/31 (R2 base, variable set) | same, failures identical | same |
| test-module-ci-trx-summary.R | 87/1 (train q, variable set) | same | same |
| test-workflow-kit-harness.R | 734/1 (train r) | same | same |
| test-workflow-omics-binding-steps.R | 57/2 (train r; its message may change between runs, as train r) | same | same (or message-changed on that known line only) |
| test-workflow-operation-registry.R | 1107/2 (train r) | same | same |
| test-kit129-run-path.R | 711/0/0/1 or 687/6/0/1 (the runner decides), equal on both sides | same | same |
| test-integration-input-transcriptomics.R 79/0; test-kit224-ruv-feature-levels.R 69/0 | as listed (variable set) | same | same |
| every other file (21 files) | R2's base counts (`r2-tests-base.tsv`) | same | same |
| R CMD check | 0 ERROR, 3 WARNING, 3 NOTE | identical | nothing only on head |

New failures: none. A head that differs from its base is a result.

### Scope (37 files)

- test-e2e-manifest.R
- test-e2e-transcriptomics-screen.R
- test-e2e-trx-journey.R
- test-integration-input-transcriptomics.R
- test-integration-model.R
- test-integration-mofa-model.R
- test-kit129-run-path.R
- test-kit200-prot-da-no-value-proteins.R
- test-kit201-treat-moderation.R
- test-kit205-enrichment-background.R
- test-kit207-summary-head-guard.R
- test-kit212-da-current.R
- test-kit213-study-parameters-labels.R
- test-kit215-not-tested-lists.R
- test-kit224-ruv-feature-levels.R
- test-kit227-report-ruv-sentences.R
- test-kit243-summary-current-outputs.R
- test-kit248-correlation-skip-record.R
- test-kit254-omics-correlation-skip.R
- test-kit263-prot-study-parameters-sections.R
- test-kit264-trx-export-current.R
- test-kit274-trx-no-ruv-export.R
- test-module-ci-trx-summary.R
- test-omic-gui-config.R
- test-omic-summary-module.R
- test-omics-module-globals.R
- test-omics-summary-core.R
- test-prot-annotation-uniprot-empty.R
- test-prot-missing-value-counts.R
- test-workflow-analysis-record.R
- test-workflow-kit-harness.R
- test-workflow-lipid-declarations.R
- test-workflow-metab-declarations.R
- test-workflow-omics-application-reads.R
- test-workflow-omics-binding-steps.R
- test-workflow-operation-registry.R
- test-workflow-trx-declarations.R

### Decision and attempts

- Every prediction holds: TRAIN GREEN for KIT-274. Anything else is a result, recorded and not re-run. One run; a re-send only for an instrument defect with
  no result judged; an HTTP 500/502 at dispatch with no run created is started once by hand.

## What the coordinator resolves at landing

- Ledger: `VAL-370` (`docs/methods-decision-ledger.md:2535`) takes its number. Status cells: VAL-318 (`:2490`) "implemented by KIT-274 (VAL-n), merged by train 33";
  VAL-332 (`:2504`) "KIT-274 part implemented (VAL-n)"; VAL-336 (`:2508`) append "the k = 0 export (X1) and its replay closed by KIT-274 (VAL-n)" (AUDIT-KIT-274 I04).
- Placeholders in the member's own records (not edited here): `docs/qa/workflow-execution/kit/kit-274/outcome.md:63` and `:88` read `VAL-<merge>`.

## Results (run 45741997862, head b6a3f6b7, base 221af09e; request file kit274-train33-b6a3f6b7.json)

Run success; every file ran to completion. Tables beside this file: `merge-train-20261126z-compare.tsv`, `merge-train-20261126z-check-compare.txt`.
Base versus head: same 36, only-on-head 1, new-failure 0, message-changed 0, count-changed 0, fixed 0. **Every prediction held: TRAIN GREEN for KIT-274.**

| File | Base | Head | Verdict | Predicted |
|---|---|---|---|---|
| test-kit274-trx-no-ruv-export.R | absent | 84/0/0/0 | only-on-head | held (T11 restored) |
| test-kit264-trx-export-current.R | 85/0 | 85/0 | same | held |
| test-omics-module-globals.R | 14/0 | 14/0 | same | held |
| test-omic-summary-module.R | 1675/28 | same | same | held |
| test-omics-summary-core.R | 676/25/0/2 | same | same | held |
| test-workflow-trx-declarations.R | 313/31 | same | same | held |
| test-module-ci-trx-summary.R | 87/1 | same | same | held |
| test-workflow-kit-harness.R | 734/1 | same | same | held |
| test-workflow-omics-binding-steps.R | 57/2 | same | same (no message change) | held |
| test-workflow-operation-registry.R | 1107/2 | same | same | held |
| test-kit129-run-path.R | 711/0/0/1 | same | same | held (the 711 variant) |
| test-integration-input-transcriptomics.R 79/0; test-kit224-ruv-feature-levels.R 69/0 | as listed | same | same | held |
| the other 24 files (the declaration's "21" miscounted them) | R2's base counts | same | same | held |

R CMD check: base and head both 0 ERROR, 3 WARNING, 3 NOTE; findings identical (6); only on head none. Held.

Budget (jobs API, summed, skipped jobs 0): **339.98 job-minutes = 22.67 core-hours** (tests head 60.42 + 71.03 + 41.10; tests base 45.25 + 64.63 + 37.30; check head 9.13, base 10.23; plan 0.38; results 0.50). Ceiling 660 (44.0); expected about 350.
R2b (run 45739818528, which the branch did not read): 23.52 job-minutes = 1.57 core-hours.

Decision: KIT-274 is green on the merged tree.

## Merge record (coordinator)

Landed on v0.9-staging on 2026-11-26 by a no-fast-forward merge of this branch at 3cf5ce33 onto a5de1c2e (the train was built on 221af09e; since then staging gained merge trains 30, 29, 31 and 32 (the KIT-258 texts; KIT-247; KIT-272; KIT-268) and ledger rows).

- Verdict: green. Run 45741997862 (37 files and R CMD check; 339.98 of 660 job-minutes, 22.67 core-hours; the request set MULTISCHOLAR_TRX_MANIFEST_NOTE=1): the new file 84/0 only on head; 36 files same with identical failures; no message change; check identical.
- The member branch holds two revert-and-restore pairs (comparison bases built by reverting the fix). Checked at landing and by the auditor's closure check: against the staging it last merged, the branch's R/ and inst/ differ by exactly the fix (three R files and the transcriptomics report template); after this merge the tree differs from staging by exactly KIT-274's change.
- Conflicts at landing: NEWS (both kept: staging's first Fixes bullets, then KIT-274's) and the ledger. R/utils_workflow_analysis_record.R and R/utils_workflow_summary.R merged by themselves: staging's side changed line 420 of the first (KIT-268) and lines 232 to 249 of the second; this train's changes are at 385 to 396, 593 and 1136 to 1165, and at 623 to 727; read at landing.
- Ledger at landing: staging's rows kept through VAL-369; KIT-274 = VAL-370. Status cells of VAL-318, VAL-332 and VAL-336 updated. Placeholders in this record and in kit-274/outcome.md replaced.
- Recorded with it (AUDIT-KIT-274-impl-20261126 I01): a test repair commit on the member branch deleted the test T11 without saying so; the coordinator did not accept the deletion; T11 was restored byte for byte and run (45739818528).
- Not covered by the train run: the combination with merge trains 29 to 32 (KIT-268's change to the same record file among them). KIT-216's train, a full suite, covers it.
- Stated limits, not changed here: a project whose report was generated before this fix keeps its old copy until Save is pressed again; a leftover RUVIII figure stays in the project and Copy delivers it until the normalisation is run again; in a run without RUV the report does not say that the integration file holds the TMM logCPM (a candidate for the next report-wording ticket).
