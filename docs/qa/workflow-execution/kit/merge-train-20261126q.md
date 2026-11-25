# Merge train 24 of 2026-11-26: KIT-213

**This train is NOT stacked.** It was branched from `origin/v0.9-staging` at `176f9d3ffa100f0480fd5100a32c002c5a5cdbe5`
(fetched before building; staging holds trains 15 to 18 and VAL-329; no code merge since the audit). It does not
contain merge trains 19 to 23 (`merge-train-20261126l` to `p`, other agents'), which are still out.

**Order at landing.** KIT-213 merges first. KIT-216's and KIT-269's merges must run
`tests/testthat/test-kit213-study-parameters-labels.R` (the writer's record block and its labels are KIT-213's; KIT-269's
handler in the same writer is a separate hunk).

## Train

Base: `176f9d3ffa100f0480fd5100a32c002c5a5cdbe5` (staging).

| Order | Ticket | Branch head | Merge commit | Ledger row | Audit |
|---|---|---|---|---|---|
| 1 | KIT-213 | `origin/kit-213` 4ec3f4b2 (contains staging 4d588ba3, not the newer merges) | 4915d322 | `VAL-341` (replaces `VAL-<merge>`) | AUDIT-KIT-213-impl-20261126 (verdict ready) |

Head (bound by the run): `4915d322790831d78c0d6e43d4f06b87ef3cac4d`. This file is in the declaration commit after it,
under `docs/` only.

### Conflicts met and their resolutions

- **`docs/methods-decision-ledger.md`** only: staging's rows VAL-314 to VAL-329 against KIT-213's row. Kept staging's rows in
  order, then KIT-213's row after VAL-329 as `VAL-341` (the placeholder is left for the coordinator to number at
  landing); its status cell now reads "merged (<commit>), train 24". Nothing else in the ledger differs from staging.
- No conflict under `R/`, `tests/`, `inst/`, `man/`, `DESCRIPTION` or `NAMESPACE`. `git diff --stat` against staging for those
  paths and `NEWS.md`: `NEWS.md` +12, the writer +37/-6 lines in two hunks, the new test file +448, `test-omic-summary-module.R` 3 lines.
- `docs/qa/workflow-execution/kit/kit-213/outcome.md` still names `VAL-<merge>` (lines 13, 16, 83, 104): not edited here.

### Auto-merged files, read in the merged result

- **`NEWS.md`**: KIT-213's one bullet added (12 lines), nothing removed.
- **`R/func_general_filemgmt_study_parameter_writer.R`**: KIT-263's hunks (staging) at the top and KIT-213's two hunks
  (1006-1022 and the helper at 1065-1096) are separate; both present. Function head (line 19):
  `createWorkflowArgsFromConfig <- function(workflow_name, description = "", organism_name = NULL, taxon_id = NULL, source_dir_path = NULL, final_s4_object = NULL, contrasts_tbl = NULL, workflow_data = NULL, omic_type = NULL) {`
  (as laid out over lines 19 to 25), with KIT-263's filter `if (identical(omic_type, "proteomics")) {` at line 103. The record
  block's three print calls (lines 1006 to 1022):
  ```
  record_lines <- c(record_lines, sprintf("  %s.%s [%s] = %s", previous$operation, parameter$name,
      parameter$provenance, studyParameterRecordValueText(parameter, step$order)))
  for (parameter in step$parameters) record_lines <- c(record_lines, sprintf("  %s [%s] = %s", parameter$name,
      parameter$provenance, studyParameterRecordValueText(parameter, step$order)))
  ...
  for (parameter in analysis$parameters) record_lines <- c(record_lines, sprintf("  %s [%s] = %s", parameter$name,
      parameter$provenance, studyParameterRecordValueText(parameter, source_order %||% NA_integer_)))
  ```
- **`tests/testthat/test-omic-summary-module.R`**: against staging the diff is exactly three lines: 2533
  (`trx_main/save_before_export`, `47 lines 96b3fdbde9b1bb30` to `387b8047e011288e`), 2619 (`trx_main/save`, `50 lines
  6309448be78f095c` to `bbadf5a52d70875f`) and 3060 (`trx_managed/save`, the same pair). KIT-263's recorder (`...` passed
  through, lines 452 to 459) is staging's and unchanged.

### Path lengths

Checked before sending: 43 paths are new against staging; the longest are KIT-213's own records, 100 bytes with the
`MultiScholaR/` prefix (`docs/qa/workflow-execution/kit/kit-213/control/text-head-1/write-0N-transcriptomics.txt`), exactly the
paths R3 (45700894419) checked with identical findings. This train's own files are shorter.

## Declaration (written before the run is sent)

**Purpose.** KIT-213's merge check on the merged tree; the full suite is NOT repeated (the ticket's R3 45700894419 ran
it; the next checkpoint follows).

**Scope recount.** `git grep -lE "createWorkflowArgsFromConfig|study_parameters|workflowRecordValueText|
saveWorkflowSummaryParameters|Analysis Record" -- tests/testthat` gives 50 files on the head (47 test files and 3 helpers,
`helper-module-ci-{lipid,metab,prot}-summary.R`, not run as files) and 49 on the base (46 test files); the difference
is the new file `test-kit213-study-parameters-labels.R`. **47 test files**, as the audit expected, including that file and
`test-kit263-prot-study-parameters-sections.R`. `test-prot-norm-core.R` is not in scope, so `MULTISCHOLAR_KIT105_W3_DUMP` is
not set. This list is not train 18's (that used KIT-263's wider grep): it adds `test-omic-enrich-module.R` and
omits `test-gui-entry-proteomics.R` and `test-kit217-intensity-filter-missing.R`.

**Run settings.** Copied from KIT-213's R2 (`kit213-r2-a68fd69c.json`): 3 shards per side, `file_timeout_minutes` 45,
no extra packages, `blas` reference, `check` true, `fetch_refs` empty. Changed: no scripts (R2's control C is not
repeated); the env switch `MULTISCHOLAR_TRX_MANIFEST_NOTE=1` is added as in the later trains' requests (R2 set no env).
Base = `176f9d3f`, head = `4915d322`; one run.

**Ceiling.** 600 job-minutes = 40.0 core-hours. Expected about 300 job-minutes (about 20 core-hours): train 18's 47-file
run with the same shard count took 294.1 job-minutes.

### Predictions

| File | Base | Head | Verdict |
|---|---|---|---|
| test-e2e-lipidomics-canonical.R | 0/0/2/ok as train 18's head | same as base | same |
| test-e2e-metabolomics-lc-gc.R | 0/0/3/ok as train 18's head | same as base | same |
| test-e2e-proteomics-dia-limpa.R | 0/0/1/ok as train 18's head | same as base | same |
| test-e2e-proteomics-dia.R | 0/0/1/ok as train 18's head | same as base | same |
| test-e2e-proteomics-enrichment.R | 0/0/2/ok as train 18's head | same as base | same |
| test-e2e-proteomics-tmt-lfq.R | 0/0/3/ok as train 18's head | same as base | same |
| test-general-filemgmt-config-direct-shared.R | 0/0/0/ok as train 18's head | same as base | same |
| test-general-filemgmt-config-roundtrip-contracts.R | 0/0/0/ok as train 18's head | same as base | same |
| test-general-filemgmt-report-helpers-shared.R | 0/0/1/ok as train 18's head | same as base | same |
| test-general-filemgmt-s4-contracts.R | 0/0/0/ok as train 18's head | same as base | same |
| test-integration-input-lipid-metab.R | 0/0/0/ok as train 18's head | same as base | same |
| test-integration-input-transcriptomics.R | 0/0/0/ok as train 18's head | same as base | same |
| test-integration-record-mofa.R | 0/0/1/ok as train 18's head | same as base | same |
| test-integration-record.R | 0/0/0/ok as train 18's head | same as base | same |
| test-integration-workbook-mofa.R | 0/0/3/ok as train 18's head | same as base | same |
| test-kit172-downstream.R | 0/0/0/ok as train 18's head | same as base | same |
| test-kit206-prot-norm-export-guard.R | 0/0/0/ok as train 18's head | same as base | same |
| test-kit207-summary-head-guard.R | 0/0/0/ok as train 18's head | same as base | same |
| test-kit212-da-current.R | 0/0/0/ok as train 18's head | same as base | same |
| test-kit213-study-parameters-labels.R | absent | 56/0/0/0 | only-on-head |
| test-kit224-ruv-feature-levels.R | 0/0/0/ok as train 18's head | same as base | same |
| test-kit227-report-ruv-sentences.R | 0/0/0/ok as train 18's head | same as base | same |
| test-kit235-lipid-metab-report-parser.R | 0/0/0/ok as train 18's head | same as base | same |
| test-kit236-report-sentences.R | 0/0/0/ok as train 18's head | same as base | same |
| test-kit243-summary-current-outputs.R | 0/0/0/ok as train 18's head | same as base | same |
| test-kit248-correlation-skip-record.R | 0/0/0/ok as train 18's head | same as base | same |
| test-kit254-omics-correlation-skip.R | read from run (376/3 at train 18; VAL-325 may have moved it) | same as base | same |
| test-kit263-prot-study-parameters-sections.R | 11/0 | same as base | same |
| test-lipid-14a-summary-module-characterization.R | 0/0/0/ok as train 18's head | same as base | same |
| test-lipid-15-summary-module-seams.R | 0/0/0/ok as train 18's head | same as base | same |
| test-metab-05-summary-module-characterization.R | 0/0/0/ok as train 18's head | same as base | same |
| test-module-ci-lipid-summary.R | 0/0/6/ok as train 18's head | same as base | same |
| test-module-ci-metab-da.R | 0/0/0/ok as train 18's head | same as base | same |
| test-module-ci-metab-summary.R | 0/0/6/ok as train 18's head | same as base | same |
| test-module-ci-prot-summary.R | 0/0/5/ok as train 18's head | same as base | same |
| test-module-ci-trx-summary.R | 1/0/0/failed as train 18's head | same as base | same |
| test-omic-enrich-module.R | read from run | same as base | same |
| test-omic-summary-module.R | 1675/28 (train 18 head; ticket R3) | same as base | same, identical failure list |
| test-omics-summary-core.R | 676/25 (train 18 head; ticket R3) | same as base | same, identical failure list |
| test-prot-12-summary-module-contracts.R | 0/0/0/ok as train 18's head | same as base | same |
| test-prot-12a-summary-module-characterization.R | 0/0/0/ok as train 18's head | same as base | same |
| test-workflow-analysis-record.R | 0/0/0/ok as train 18's head | same as base | same |
| test-workflow-figure-paths.R | 0/0/0/ok as train 18's head | same as base | same |
| test-workflow-integration-input.R | 0/0/0/ok as train 18's head | same as base | same |
| test-workflow-memory-mechanism-member.R | 0/0/13/ok as train 18's head | same as base | same |
| test-workflow-native-session.R | 0/0/2/ok as train 18's head | same as base | same |
| test-workflow-omics-application-reads.R | 0/0/0/ok as train 18's head | same as base | same |
| R CMD check | | | identical findings, nothing only on head |

Counts are passed/failed/errors/skipped as the train 18 head table (`merge-train-20261126k-compare.tsv`) shows them; staging
now holds trains 15 to 18, so those are the expected base counts; a base count that differs is read from the run, and the
prediction is the verdict (`same`), not the count, as in train 18. `test-omic-enrich-module.R` was not in train 18's scope.

### Decision and attempts

- Every prediction holds: TRAIN GREEN for KIT-213.
- Any other outcome stops the train for this member: recorded with its files and counts, not repaired, not re-sent.
- One run; a re-send only for an instrument defect with no result judged. An HTTP 500/502 at dispatch with no run
  created: started once by hand with the README's command, recorded.

## What the coordinator resolves at landing

- The ledger: `VAL-341` takes its number (next free after the rows landed before it); status cell "merged (<commit>), train 24"; VAL-169's
  status cell: "implemented by KIT-213 (VAL-341); merged <commit>".
- `NEWS.md`: KIT-213's bullet in landing order.
- `kit-213/outcome.md` lines 13, 16, 83, 104 name `VAL-<merge>`.
- Records owed at landing (audit I13): plan section 11 row, section 2 rows, ticket and index, DISPATCH, E-7, KIT-216 partner line, H3.
- KIT-216 and KIT-269 merges run `tests/testthat/test-kit213-study-parameters-labels.R`.

## Result

Run 45716990567 (request `merge-train-20261126q-4915d322`, request commit 9d5ca71d), created by the dispatcher without error
(no hand dispatch); green (infrastructure). Head 4915d322 vs base 176f9d3f. Tables: `merge-train-20261126q-compare.tsv`,
`merge-train-20261126q-check-compare.txt` (same directory).

**Job-minutes 307.8 = 20.52 core-hours** (job-minutes x 4 / 60) against the ceiling of 600 job-minutes (40.0 core-hours)
and the expected about 300: tests head 43.7/52.7/41.9, tests base 27.4/55.3/57.4, check head 15.7, check base 12.9,
plan 0.3, results 0.5 (job start to end from the run's job times).

| File | Predicted | Observed base | Observed head | Verdict |
|---|---|---|---|---|
| test-kit213-study-parameters-labels.R | only-on-head, 56/0 | absent | 56/0/0/0 | only-on-head: HELD |
| test-omic-summary-module.R | same, 1675/28, identical failures | 1675/28 | 1675/28 | same: HELD |
| test-omics-summary-core.R | same, 676/25, identical failures | 676/25 (2 skipped) | 676/25 (2 skipped) | same: HELD |
| test-kit263-prot-study-parameters-sections.R | same, 11/0 | 11/0 | 11/0 | same: HELD |
| test-kit254-omics-correlation-skip.R | same (count read from run) | 379/0 | 379/0 | same: HELD (VAL-325's edit is on staging) |
| test-module-ci-trx-summary.R | same | 87/1 | 87/1 | same: HELD (the known failure, identical message) |
| test-omic-enrich-module.R | same (count read from run) | 1011/59 | 1011/59 | same: HELD (identical failures apart from digest values) |
| the other 39 files | same | same | same | same: HELD |
| R CMD check | identical, nothing only on head | 0 E / 3 W / 3 N | 0 E / 3 W / 3 N | identical (6 findings on both; only on head: none; only on base: none): HELD |

Totals: 47 files, 46 `same`, 1 `only-on-head` (the new file), no `new-failure`, no `message-changed`, no `count-changed`.
Tests head: passed 8302, failed 113; base: passed 8246, failed 113 (the 56 are the new file). Four files fail on both sides
with identical failure lists: test-module-ci-trx-summary.R, test-omic-enrich-module.R, test-omic-summary-module.R,
test-omics-summary-core.R.

**Decision.** Every prediction held. TRAIN GREEN for KIT-213: the new file passes 56/0 on head, the re-frozen leaves of
test-omic-summary-module.R pass against the fix with the failure list unchanged, and no file or check finding moved.

## Merge record (coordinator)

Landed on v0.9-staging on 2026-11-26 by a no-fast-forward merge of this branch at 417305aa onto 75acd421 (merge trains 19 to 22 landed after this train was built on 176f9d3f; the train is not stacked on them).

- Verdict: green. Run 45716990567 (47 files and R CMD check; 307.8 of 600 job-minutes, 20.52 core-hours): the new file 56/0 only on head; 46 files same; check identical (0 errors, 3 warnings, 3 notes on both sides).
- Conflicts at landing: the ledger only (NEWS merged by itself). No code conflict: staging did not change the study-parameter writer or the summary module's test file after 176f9d3f.
- Ledger at landing: staging's rows kept through VAL-340; KIT-213 = VAL-341; the status cell of VAL-169 now names VAL-341 and this train. The placeholder in kit-213/outcome.md (four places) and in this record is replaced by VAL-341.
- Not covered by the train run: the combination with merge trains 19 to 22 (KIT-264, KIT-258, KIT-265, KIT-267), none of which touches the files this train changes. The next checkpoint run covers it.
- Owed by the next merges: KIT-216 and KIT-269 run tests/testthat/test-kit213-study-parameters-labels.R in their trains.
