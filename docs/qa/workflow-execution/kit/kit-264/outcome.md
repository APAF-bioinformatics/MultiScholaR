# KIT-264 outcome: the transcriptomics record and report no longer state a leftover RUVIII export

Branch `kit-264`, from origin/v0.9-staging 988512e6 (trains 12 and 13 landed). Commits: b374e2c5 (tests only), d68ec438
(tests only: instrument fix), 317c3630 (the fix: `R/utils_workflow_analysis_record.R`, `R/mod_omic_norm.R`,
`inst/reports/transcriptomics/transcriptomics_report.rmd`, `NEWS.md`). Ticket r1.2 (closure check C1, C2, N-a, N-c, N-d
written in before R1). Predictions: the header of `tests/testthat/test-kit264-trx-export-current.R` (before R1) and
`predictions.md` (before R2).

## What was changed

1. `R/utils_workflow_analysis_record.R`: new `workflowRecordExportIsCurrent(snapshot, export)`; the export loop of
   `buildWorkflowAnalysisRecord` lists every export except `INTEGRATION_MODEL_OPERATION` only when
   `export$source$entity_id` is the `entity_id` of `snapshot$states[[state]]` for a state in `snapshot$state_history`
   (by entity, never by name; a NULL source entity is not current). The state head is not rewritten.
2. `inst/reports/transcriptomics/transcriptomics_report.rmd`, chunk `ruviii-figure`: prints the existing
   "Figure not available: ruviii_before_and_after.png" text when the record lists no export, else `show_figure` as before.
3. `R/mod_omic_norm.R`, `omicNormOwnFiles` (transcriptomics): `ruviii_before_and_after.png` and `.pdf` added, so a run's
   start removes them from `publication_graphs_dir` and Copy no longer delivers the earlier export's figure.
4. `NEWS.md` (disclosure that a project saved before the fix keeps the earlier statement until Save is pressed again).
No existing test edited; no `man/` or NAMESPACE change; no stored reference moved.

## Observed (labels: OBSERVED = in the run named)

- R1, run 45684412917 on b374e2c5: E1, E2, E3, E5, E6, E7 failed at the predicted base values (entry listed, text "k = 2",
  figure value an `include_graphics` object, `omicNormOwnFiles` returns three names, both figure files remain); L1, L2, G1,
  G2, G5, G6, G7, G8 passed; X1 observed. E8 and its control could not run: **instrument defect 1** (the helpers
  `module_ci_trx_summary_steps` and `module_ci_trx_docx_text` are defined in `test-module-ci-trx-summary.R`, not in a
  helper; "could not find function"); no expected value, threshold or stop rule changed by the fix (VAL-54 conditions
  i to iv; rerun 1 of 3).
- R1-2, run 45687693688 on d68ec438 (the file with its own copies of the two helpers): 70 passed, 15 failed, 0 errors. Every
  E row fails at its base value, E8 included: the never-exported sentence is absent, "RUVIII (ruv::RUVIII) was applied" is
  present, the figure is in the results summary; E8's preconditions held (no step failed, export in the state, last
  `trx_ruvs` step `ruv_applied` FALSE, exactly one docx); E8 control passes (RUVIII sentence and figure present); L and G
  rows pass. Stop rule 2 did not fire (no G row failed, every E row failed on base).
- R2, run 45690457065, head 317c3630 against base d68ec438, 67 files, 5 shards per side, check both: 66 files `same`
  (equal counts, identical failure messages), the new file `fixed` (base 70/15, head 85 passed, 0 failed, 0 errors, 0 skipped).
  The known-failing files kept their counts: test-prot-norm-core 1548/306, test-workflow-saved-da 755/196,
  test-omic-enrich-module 1011/59, test-omic-summary-module 1675/28, test-omics-summary-core 676/25, test-omics-enrich-core
  221/20, test-workflow-trx-declarations 313/31, test-workflow-operation-registry 1095/2, test-module-ci-trx-da 52/1,
  -norm 63/1, -summary 87/1, -design 66/0+1 error, test-workflow-kit-harness 734/1, test-workflow-omics-binding-steps 58/1;
  also `same` and red on both sides, not named in the ticket: test-module-ci-trx-qc 38/0+1 error,
  test-workflow-omics-owner 78/1. `test-kit257-norm-own-files.R` 133/0, `test-omics-norm-state-guards.R` 651/0,
  `test-kit254-omics-correlation-skip.R` 379/0, `test-kit207-summary-head-guard.R` 292/0 (KIT-207 E6 and E9 = G4), the
  integration-record file (G3) and the three KIT-186/192 renders in `test-module-ci-trx-summary.R`: same.
  `R CMD check`: 0 ERRORs, 3 WARNINGs, 3 NOTEs on both sides, 6 findings identical, none only on head.

## Acceptance criteria

| # | Criterion | Evidence |
|---|---|---|
| A1 | R1: every E row fails at its base value, G and L pass, X1 observed | R1-2 (45687693688); R1 (45684412917) for the rows it could run; instrument defect 1 recorded |
| A2 | The fix: every E, G and L row passes | R2 head: 85 passed, 0 failed, 0 errors, 0 skipped |
| A3 | R2: every file `same` except the new file; check identical | R2: 66 same, new file fixed; check identical |
| A4 | NEWS, ledger row, wording rows, this file with "For the change report" and "Not checked", X1 reported | NEWS.md; ledger `VAL-<merge>`; wording.md; below |

## X1 (finding for the coordinator, to be put to the owner)

X1 is a characterisation; it did not and cannot pass or fail the ticket. The same line printed on both sides (the code
under test does not touch this case), verbatim:

`KIT-264 X1 ruv_applied=FALSE export=ok entry=listed k=0 text=RUVIII (ruv::RUVIII) was applied to the TMM logCPM, with the RUVs inferred negative controls and k = 0, for the integration export only, with the replicate structure M taken from the group. The differential expression results do not use it.`

Meaning: after a run with no RUV (design `none`), the RUVIII export is **accepted** (`ruv::RUVIII` returns with k = 0 and no
controls), the record **lists** it (its source is the current `trx_ruvs` state, so the fix leaves it), and the report says
RUVIII was applied with k = 0 beside the caption "No correction for unwanted variation was applied". This is out of scope
(ticket 2.3, R-5); whether such an export should be refused or described is the owner's decision.

## Not checked

- Whether a project saved before the fix is repaired by Save (by reading: the record is rebuilt at each Save; no run).
- The replay of a record written before the fix that lists a stale export (pre-existing; candidate).
- The Summary tab's in-session RUVIII status and figure after a later normalisation (candidate, ticket section 3).
- The export's two other figures (`trx_ruviii_pca.png`, `trx_ruviii_rle.png`): delivered by nothing, left.
- No browser; no real data; no full suite (the merge check rides a merge train).
- Core-hours: 15.78 + 14.53 + 629.15 = 659.46 job-minutes = 43.96 core-hours (expected 51.27, ceiling 145.67), summed from
  the runs' job start and completion times.

## For the change report

New to 0.9 (transcriptomics, its analysis record and report are not on `main`; the tag v0.9.0 has the defect).

- **What was wrong, who was affected.** A transcriptomics user who exported the RUVIII matrix (Summary tab) and then ran
  the normalisation again (another k by default; without RUV when the RUV design section is shown, which is off by default
  and on under the installer's test features) and then saved and made the report: the analysis record, the record workbook,
  `study_parameters.txt`, a replay and the report still listed the earlier export, with its earlier k and its figure,
  beside a stage caption saying no correction had been applied (or a different k). No number was affected: with RUV
  applied the integration input refused the stale export (KIT-207 SD-4); without RUV it ignored it.
- **Code.** Base 988512e6: the loop at `R/utils_workflow_analysis_record.R:532` listed every export of the head;
  `inst/reports/transcriptomics/transcriptomics_report.rmd:172-174` showed the figure file whenever it existed;
  `omicNormOwnFiles` (`R/mod_omic_norm.R:1058-1062`) did not name the figure. Fix 317c3630: `workflowRecordExportIsCurrent`
  and the loop condition; the figure chunk follows the record's export; two names added to `omicNormOwnFiles`.
- **What changed in results.** None. A report of an analysis saved earlier keeps the earlier statement until Save is
  pressed again.
- **Evidence.** Runs 45684412917, 45687693688, 45690457065 (tables in this folder). Synthetic data only (MCI-030
  transcriptomics corpus). Limits: no real data; not changed: an export made after a run without RUV (k = 0), see X1.
