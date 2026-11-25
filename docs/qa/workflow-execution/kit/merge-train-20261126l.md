# Merge train 19 of 2026-11-26: KIT-264

**This train is STACKED on merge train 18 (KIT-263), itself stacked on train 17 (KIT-246), train 16 (KIT-261) and
train 15 (KIT-190, KIT-204, KIT-217).** It was branched from `origin/merge-train-20261126k` at its tip
`4bcff128bbe016d67fc5a6fa805b035bf87799dc` (fetched before building; it carries the KIT-263 merge 0090ae3d and train
18's declaration). Staging is 4d588ba3. **Train 19 lands only if trains 15, 16, 17 and 18 all land with their code
unchanged.** If any is rebuilt or its code changes, this train is rebuilt and its check repeated. Train 15's branch is
receiving one tests-only commit that re-derives constants in `test-kit254-omics-correlation-skip.R`; the base does not
contain it (expected; see the prediction for that file).

Train 18's own run 45699156973 was `in_progress` when this train was built; its result is not part of this record.

## Train

Base: `4bcff128bbe016d67fc5a6fa805b035bf87799dc` (train 18's tip).

| Order | Ticket | Branch head | Merge commit | Ledger row | Audit |
|---|---|---|---|---|---|
| 1 | KIT-264 | `kit-264` 3e299387 | 9f65a58c | `VAL-<train19-1>` (replaces `VAL-<merge>`) | AUDIT-KIT-264-impl-20261126 (verdict ready) |

Head (bound by the run): `9f65a58c9a46d63a1a43f966ded0acc032363fc6`. This file and the predictions are in the
declaration commit after it, under `docs/` only; the records commit after the run likewise.

### Conflicts met and their resolutions

- **`NEWS.md`**: the base's bullets (KIT-261, KIT-246, KIT-263) against KIT-264's, at the top of `### Fixes`. Both
  kept, the base's first, then KIT-264's; no text changed.
- **`docs/methods-decision-ledger.md`**: the base's appended rows against KIT-264's. Kept the base's rows in order,
  then KIT-264's as `VAL-<train19-1>`. The partly numbered rows inherited from the base are NOT repaired here; the
  coordinator resolves the ledger at landing.
- `docs/qa/workflow-execution/kit/kit-264/outcome.md` line 55 names the placeholder; it now also names
  `VAL-<train19-1>`.
- No conflict under `R/`, `tests/`, `inst/`, `man/`, `DESCRIPTION` or `NAMESPACE` (auto-merged without conflict:
  `R/mod_omic_norm.R`, `R/utils_workflow_analysis_record.R`). Code files changed against the base:
  `R/mod_omic_norm.R`, `R/utils_workflow_analysis_record.R`, `inst/reports/transcriptomics/transcriptomics_report.rmd`,
  `tests/testthat/test-kit264-trx-export-current.R` (new). KIT-264 adds no replay marker.

### Path lengths

Checked before sending: every path new or changed since the base is at most 72 bytes with the `MultiScholaR/` prefix
(longest: `docs/qa/workflow-execution/kit/kit-264/r2-check-compare.txt`, 72; the new test file, 63).

## Declaration (written before the run is sent)

**Purpose.** KIT-264's merge check on the merged tree, stacked on train 18.

**Scope recount.** The ticket's grep (`git grep -lE "buildWorkflowAnalysisRecord|workflowRecordExportAnalysis|
trx_ruviii_export|transcriptomics_report|omicNormOwnFiles|runTrxNormalisation|replayWorkflowAnalysisRecord|
exportWorkflowAnalysisRecordBundle" <commit> -- 'tests/testthat/test-*.R'`) gives **70 files on the base 4bcff128** and
**71 on the head 9f65a58c** (the new `test-kit264-trx-export-current.R`; no other difference). That is the register's
count (I12: 70 + the new file = 71). No difference. `test-kit217-intensity-filter-missing.R`,
`test-kit248-correlation-skip-record.R`, `test-kit263-prot-study-parameters-sections.R`,
`test-ruv-control-anova-grouping.R`, `test-prot-norm-core.R`, `test-kit254-omics-correlation-skip.R`,
`test-integration-record.R`, `test-workflow-operation-registry.R` and `test-workflow-omics-binding-steps.R` are among
the 71. The optional 72nd file of the register (`test-gui-entry-transcriptomics.R`) is not added.

**Run settings.** `MULTISCHOLAR_TRX_MANIFEST_NOTE=1` and `MULTISCHOLAR_KIT105_W3_DUMP=1` on both sides; extra packages
limpa, mixOmics, sva; reference BLAS; 6 shards per side; default file timeout (45 min); check on both sides; no
scripts.

**Ceiling.** 980 job-minutes = 65.3 core-hours (job-minutes x 4 / 60): KIT-264's R2 row. Expected about 680
job-minutes (about 45 core-hours).

### Predictions

| File | Base | Head | Verdict |
|---|---|---|---|
| test-kit264-trx-export-current.R | absent | 85/0/0/0 | only-on-head |
| test-kit254-omics-correlation-skip.R | 376/3 (the known train-15 state before its tests-only edit) | 376/3 | same |
| test-integration-record.R | 200/0/0/0 (train 17) | 200/0/0/0 | same |
| test-workflow-operation-registry.R | 1107/2 (train 15) | 1107/2 | same |
| test-workflow-omics-binding-steps.R | 57/2 (train 15) | 57/2 | same |
| test-kit217-intensity-filter-missing.R | 344/0 | 344/0 | same |
| test-ruv-control-anova-grouping.R | 139/0 | 139/0 | same |
| test-kit248-correlation-skip-record.R | 206/0 | 206/0 | same |
| test-kit263-prot-study-parameters-sections.R | 11/0 | 11/0 | same |
| test-prot-norm-core.R 1548/306; test-workflow-saved-da.R 755/196; test-omic-enrich-module.R 1011/59; test-omic-summary-module.R 1675/28; test-omics-summary-core.R 676/25/0/2; test-omics-enrich-core.R 221/20/0/1; test-workflow-trx-declarations.R 313/31; test-module-ci-trx-da.R 52/1; -norm 63/1; -summary 87/1; -design 66/0/1; -qc 38/0/1; test-workflow-kit-harness.R 734/1; test-workflow-omics-owner.R 78/1 | KIT-264's R2 counts (known red) | as base | same |
| every other file of the 71 | count read from the run (not in an earlier declaration) | as base | same |
| R CMD check | | | identical findings, nothing only on head |

Where "read from the run" is shown, the prediction is `same` and the count is read from the run (stated as such).

### Decision and attempts

- Every prediction holds: TRAIN GREEN, landing only after trains 15 to 18 land with their code unchanged.
- Any other verdict is a result: recorded with its files and counts, not repaired or re-sent.
- One run; one re-send only for an instrument defect with no result judged. An HTTP 500/502 at dispatch with no run
  created: started once by hand with the README's command, recorded.

## What the coordinator resolves at landing

- The ledger: `VAL-<train19-1>` takes its number; its `<id>` and `<t>` fill; the earlier trains' placeholders and
  staging's VAL-298 to VAL-313 rows are reconciled as the trains left them.
- `NEWS.md`: KIT-261's, KIT-246's, KIT-263's and KIT-264's bullets in that order if landed in sequence.
- `kit-264/outcome.md` line 55 placeholders.
- Train 15's tests-only edit of `test-kit254-omics-correlation-skip.R` (376/3 becomes 379/0 on both sides).

## Result

Run 45707095842 (request commit 0b5ad4b6; conclusion success; jobs tests success, check success; started
2026-11-25 21:44 UTC). Base 4bcff128, head 9f65a58c; 71 files, 6 shards per side. No re-send; dispatch was automatic.

**Verdict: every prediction held. TRAIN GREEN** (landing only after trains 15 to 18 land with their code unchanged;
train 18 has since landed on staging).

| File | Base | Head | Verdict |
|---|---|---|---|
| test-kit264-trx-export-current.R | absent | 85/0/0/0 | only-on-head (as predicted) |
| test-kit254-omics-correlation-skip.R | 376/3 | 376/3 | same (as predicted; staging has 379/0 through train 15's tests-only edit, VAL-325) |
| test-integration-record.R | 200/0/0/0 | 200/0/0/0 | same |
| test-workflow-operation-registry.R | 1107/2 | 1107/2 | same |
| test-workflow-omics-binding-steps.R | 57/2 | 57/2 | same |
| test-kit217-intensity-filter-missing.R | 344/0 | 344/0 | same |
| test-ruv-control-anova-grouping.R | 139/0 | 139/0 | same |
| test-kit248-correlation-skip-record.R | 206/0 | 206/0 | same |
| test-kit263-prot-study-parameters-sections.R | 11/0 | 11/0 | same |
| known red files (prot-norm-core 1548/306, saved-da 755/196, omic-enrich-module 1011/59, omic-summary-module 1675/28, omics-summary-core 676/25/0/2, omics-enrich-core 221/20/0/1, trx-declarations 313/31, module-ci-trx-da 52/1, -norm 63/1, -summary 87/1, -design 66/0/1, -qc 38/0/1, kit-harness 734/1, omics-owner 78/1) | as predicted | identical | same |
| every other file of the 71 | read from the run | identical | same |

compare.tsv: 71 files; the only non-`same` verdict is the new file (`only-on-head`).

**R CMD check** (both sides exit 0): base ERRORs 0, WARNINGs 3, NOTEs 3; head ERRORs 0, WARNINGs 3, NOTEs 3; 6
findings identical on both; only on head: none; only on base: none.

**Cost.** 40822 job-seconds = **680.4 job-minutes = 45.36 core-hours** (job-minutes x 4 / 60); ceiling 980 (65.3
core-hours), expected about 680. Check jobs 685 s (base) and 934 s (head).
