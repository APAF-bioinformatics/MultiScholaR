# Merge train 18 of 2026-11-26: KIT-263

**This train is STACKED on merge train 17 (KIT-246), itself stacked on train 16 (KIT-261) and train 15 (KIT-190,
KIT-204, KIT-217).** It was branched from `origin/merge-train-20261126j` at its tip
`2b5fa8435d326b3c4d76cdaa639e0911361cdcf9` (fetched before building). Staging is 4d588ba3. **Train 18 lands only if
trains 15, 16 and 17 all land with their code unchanged.** If any is rebuilt or its code changes, this train is rebuilt
and its check repeated.

## Train

Base: `2b5fa8435d326b3c4d76cdaa639e0911361cdcf9` (train 17's tip).

| Order | Ticket | Branch head | Merge commit | Ledger row | Audit |
|---|---|---|---|---|---|
| 1 | KIT-263 | `origin/kit-263` 1cf78e20 | 0090ae3d | `VAL-<train18-1>` (replaces `VAL-<merge>`) | AUDIT-KIT-263-impl-20261126 (verdict ready) |

Head (bound by the run): `0090ae3d453e207cd8061790dee960ae8cb1bbf7`. This file and the predictions are in the
declaration commit after it, under `docs/` only; the records commit after the run likewise.

### Conflicts met and their resolutions

- **`NEWS.md`**: KIT-261's and KIT-246's bullets (base) against KIT-263's, at the top of `### Fixes`. Both kept, the
  base's first, then KIT-263's; no text changed.
- **`docs/methods-decision-ledger.md`**: the base's appended rows against KIT-263's. Kept the base's rows in order,
  then KIT-263's as `VAL-<train18-1>`. Trains 16 and 17 hold staging's later VAL rows only partly: NOT repaired here;
  the coordinator resolves the ledger at landing.
- `docs/qa/workflow-execution/kit/kit-263/outcome.md` line 26 names the placeholder; it now also names
  `VAL-<train18-1>`.
- No conflict under `R/`, `tests/`, `inst/`, `man/`, `DESCRIPTION` or `NAMESPACE`. Code files changed against the
  base: `R/func_general_filemgmt_study_parameter_writer.R`, `R/utils_workflow_summary.R`,
  `man/createWorkflowArgsFromConfig.Rd`, `tests/testthat/test-omic-summary-module.R`,
  `tests/testthat/test-kit263-prot-study-parameters-sections.R`.

### Path lengths

Checked before sending: every path new since the base is at most 72 bytes with the `MultiScholaR/` prefix
(`docs/qa/workflow-execution/kit/kit-263/r2-check-compare.txt`, `tests/testthat/test-kit263-prot-study-parameters-sections.R`).

## Declaration (written before the run is sent)

**Purpose.** KIT-263's merge check on the merged tree.

**Scope recount.** The ticket's grep (`git grep -lE "createWorkflowArgsFromConfig|createWorkflowArgsFn|
saveWorkflowSummaryParameters|study_parameters|Parameters from Final S4 Object|Additional Configuration Parameters"
-- tests/testthat`) gives 47 hits on the base 2b5fa843 (44 test files and 3 helpers) and 48 on the head (the new KIT-263
file). Scope = the 44 base test files, plus `test-kit263-prot-study-parameters-sections.R`,
`test-gui-entry-proteomics.R` and `test-kit217-intensity-filter-missing.R` = **47 test files**, as the audit's I09
expected (47). No difference from 47. `test-prot-norm-core.R` is not in scope, so `MULTISCHOLAR_KIT105_W3_DUMP` is not set.

**Run settings.** `MULTISCHOLAR_TRX_MANIFEST_NOTE=1` on both sides; no extra packages (as KIT-263's R2); reference
BLAS; 3 shards per side (the ticket's R2); `test-gui-entry-proteomics.R` file timeout 120 min; check on both sides; no
scripts.

**Ceiling.** 600 job-minutes = 40.0 core-hours: KIT-263's R2 row (ticket section 10: 6 x 70 + 2 x 90), which is its
merge-check ceiling. Expected about 350 job-minutes (about 23 core-hours).

### Predictions

| File | Base | Head | Verdict |
|---|---|---|---|
| test-e2e-lipidomics-canonical.R | read from run | same as base | same |
| test-e2e-metabolomics-lc-gc.R | read from run | same as base | same |
| test-e2e-proteomics-dia-limpa.R | read from run | same as base | same |
| test-e2e-proteomics-dia.R | read from run | same as base | same |
| test-e2e-proteomics-enrichment.R | read from run | same as base | same |
| test-e2e-proteomics-tmt-lfq.R | read from run | same as base | same |
| test-general-filemgmt-config-direct-shared.R | read from run | same as base | same |
| test-general-filemgmt-config-roundtrip-contracts.R | read from run | same as base | same |
| test-general-filemgmt-report-helpers-shared.R | read from run | same as base | same |
| test-general-filemgmt-s4-contracts.R | read from run | same as base | same |
| test-gui-entry-proteomics.R | 92/101 (train 15's prediction) | same as base | same |
| test-integration-input-lipid-metab.R | read from run | same as base | same |
| test-integration-input-transcriptomics.R | read from run | same as base | same |
| test-integration-record-mofa.R | read from run | same as base | same |
| test-integration-record.R | 200/0/0/0 (train 17's head) | 200/0/0/0 | same |
| test-integration-workbook-mofa.R | read from run | same as base | same |
| test-kit172-downstream.R | read from run | same as base | same |
| test-kit206-prot-norm-export-guard.R | read from run | same as base | same |
| test-kit207-summary-head-guard.R | read from run | same as base | same |
| test-kit212-da-current.R | read from run | same as base | same |
| test-kit217-intensity-filter-missing.R | 344/0/0/0 | 344/0/0/0 | same |
| test-kit224-ruv-feature-levels.R | read from run | same as base | same |
| test-kit227-report-ruv-sentences.R | read from run | same as base | same |
| test-kit235-lipid-metab-report-parser.R | read from run | same as base | same |
| test-kit236-report-sentences.R | 289/0 | 289/0 | same |
| test-kit243-summary-current-outputs.R | read from run | same as base | same |
| test-kit248-correlation-skip-record.R | 206/0 | 206/0 | same |
| test-kit254-omics-correlation-skip.R | 379/0 | 379/0 | same |
| test-kit263-prot-study-parameters-sections.R | absent | 11/0/0/0 | only-on-head |
| test-lipid-14a-summary-module-characterization.R | read from run | same as base | same |
| test-lipid-15-summary-module-seams.R | read from run | same as base | same |
| test-metab-05-summary-module-characterization.R | read from run | same as base | same |
| test-module-ci-lipid-summary.R | read from run | same as base | same |
| test-module-ci-metab-da.R | read from run | same as base | same |
| test-module-ci-metab-summary.R | read from run | same as base | same |
| test-module-ci-prot-summary.R | read from run | same as base | same |
| test-module-ci-trx-summary.R | read from run | same as base | same |
| test-omics-summary-core.R | read from run | same as base | same |
| test-omic-summary-module.R | 1675/28 | 1675/28 | same |
| test-prot-12a-summary-module-characterization.R | read from run | same as base | same |
| test-prot-12-summary-module-contracts.R | read from run | same as base | same |
| test-workflow-analysis-record.R | read from run | same as base | same |
| test-workflow-figure-paths.R | read from run | same as base | same |
| test-workflow-integration-input.R | read from run | same as base | same |
| test-workflow-memory-mechanism-member.R | read from run | same as base | same |
| test-workflow-native-session.R | read from run | same as base | same |
| test-workflow-omics-application-reads.R | read from run | same as base | same |
| R CMD check | | | identical findings, nothing only on head |

Where "read from run" is shown the base count was not found in the earlier trains' declarations; the prediction is
`same` and the count is read from the run (stated as such).

### Decision and attempts

- Every prediction holds: TRAIN GREEN, landing only after trains 15, 16 and 17 land with their code unchanged.
- Any other verdict is a result: recorded with its files and counts, not repaired or re-sent.
- One run; one re-send only for an instrument defect with no result judged. An HTTP 500/502 at dispatch with no run
  created: started once by hand with the README's command, recorded.

## What the coordinator resolves at landing

- The ledger: `VAL-<train18-1>` takes its number; its `<id>` and `<t>` fill; the earlier trains' placeholders and
  staging's VAL-298 to VAL-301 rows are reconciled as trains 16 and 17 left them.
- `NEWS.md`: KIT-261's, KIT-246's and KIT-263's bullets in that order if landed in sequence.
- `kit-263/outcome.md` line 26 placeholder.

## Result

(filled after the run)
