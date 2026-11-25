# KIT-254 outcome

Ticket `tickets/kit/KIT-254.md` r1.2 (Re-check 1 of AUDIT-KIT-254-plan-20261125 applied). Branch `kit-254` from
**v0.9-staging e0747fc1** (merge train 4 of 2026-11-25; staging moved from 16112917 while step 0 was written).
Coordinator's ruling: OD-1 and OD-2 are the owner's and open; this implementation makes no change to replay or to
`R/utils_workflow_analysis_record.R` (E8 is a G row, the (a0) form); R2 is the 45-file scope, counted again at the
branch point (below). **Whatever the owner chooses under OD-1 ((a), (b) or (c)) is a later, separate commit with its
own failing-first run (R1b) and its own comparison over the replay consumers; nothing here presumes the answer.**

Labels: OBSERVED = a runner's output (run id given). READ = by reading the code at the commit named. INFERRED = follows
from what was read or observed, not executed itself.

## Step 0 (by reading, before G0)

**Branch point.** `git rev-parse origin/v0.9-staging` = e0747fc1 (78 commits after 16112917). Among the files the
ticket cites, these did not change between 16112917 and e0747fc1 (`git diff --stat`): `R/mod_omic_norm.R`,
`R/utils_workflow_omics_actions.R`, `R/utils_workflow_omics_declarations.R`, `R/utils_workflow_execution.R`,
`R/func_lipid_s4_correlation_methods.R`, `R/func_metab_s4_qc_methods.R`, `R/utils_workflow_state.R`,
`R/utils_omic_gui_config.R`, `R/utils_workflow_summary.R`, `R/utils_workflow_saved_exports.R`,
`R/func_general_filemgmt_study_parameter_writer.R`, `tests/testthat/test-omics-norm-state-guards.R`,
`tests/testthat/test-kit207-summary-head-guard.R`, `tests/testthat/test-kit212-da-current.R` and the helpers the ticket
names. Every line the ticket cites in them therefore reads as r1.1 states. Changed (READ, `git diff 16112917 e0747fc1`):
- `R/mod_omic_da.R`: three lines added at `:263-269` (KIT-220); the DA's source read is now `:419` (was `:416`); the
  head-forgetting block `:158-162` unchanged.
- `R/mod_omic_norm_panels.R`: one argument at `:861` (KIT-222); the summary `:1297-1322` and the final PCA `:1336`
  unchanged.
- `R/utils_workflow_dia_state.R`: KIT-222 refusals added in the proteomics design builder; `selectWorkflowStateReference`
  is `:1625-1637` (truncation `:1632-1634`), `truncateDIAStepProvenance` `:846`.
- `R/utils_workflow_dia_da.R`: KIT-220's fold-change cutoff for the interactive volcano; `workflowDaResultStatus`
  `:125-143` unchanged.
- `R/utils_workflow_analysis_record.R`: KIT-222's replay statement (`design_statement` in replay's value,
  `workflowReplayKeyCollisionStatement`). Not touched by this ticket (coordinator's ruling); the mechanism OD-1 (a) would
  use now exists on staging.

**R2's scope, counted again (READ).** The ticket's grep over `tests/testthat` without helpers finds **43** files at
e0747fc1 (41 at 16112917, plus `test-kit220-da-threshold-consumers.R` and `test-kit222-replicate-keys.R`, merge train 4).
With `test-gui-entry-lipidomics.R`, `test-gui-entry-metabolomics.R`, `test-metab-03c-norm-module-seams.R` and the new
file: **47 files** (the ticket's "45" counted at 16112917). The two new files ran 51.7 s and 53.5 s
(`merge-train-20261125d/tests.tsv`): about 1.8 minutes per side, inside R2's expected figure's margin.

**Tests that press Apply or Skip at a correlation-step head (READ).** `test-omics-norm-state-guards.R`: F1 (Apply,
Apply again, Skip at the correlation state, `:420-426`), T1 (`:543-560`) and T8 (`:1028-1040`) (the declared set-up
edit), T6 (e) (`:908-913`, head names only). `test-kit212-da-current.R`: E4's Q3 (`:724`, Apply again after a DA) and
G5 (`:1077`, the store step through the wrapper on the filtered head; F4). `test-kit207-summary-head-guard.R`,
`test-module-ci-lipid-norm.R`, `test-module-ci-metab-norm.R`, `test-omic-norm-module.R`: Apply or Skip at the run's
result only (no cut).

**E12's queued press in the worker-journey harness (READ).** `omicsManagedChainMocks` binds `pollWorkflowAction` to
"complete" (`helper-workflow-omics-worker.R:155`), and its `runWorkflowAction` runs the worker and accepts at submission
(`:151-154`); the step delivers on its poll observer (`R/utils_workflow_execution.R:949-971`). A second binding of
`pollWorkflowAction` returning "running" while a flag is set keeps the first Apply's pairs action active, so a second
press is queued (`R/utils_workflow_omics_actions.R:84-87`); releasing the flag and advancing the mock clock
(`session$elapse`) delivers it. G0 runs this (case E12).

**Stored records (READ, `git ls-files`, `data/benchmark/` excluded).** The committed JSON files naming
`<prefix>_norm_complete` are the KIT-001 step maps and `tests/testdata/mci/manifest.json` (descriptive lists, not
analysis records). No committed analysis record has a `<prefix>_norm_complete` step after a
`<prefix>_correlation_filtered` step. The committed lipid/metab records naming `<prefix>_correlation_filtered` (KIT-031,
KIT-032 journeys) come from one press each; they were not checked value by value for a correlation step whose input is
its own result.

## G0 (run 45599722106; script `g0-base-values.R` on branch `kit-254-g0` at 671aba2c = e0747fc1 + the script)

Predictions were written in the script and in `g0-predictions.tsv` before any value was read. Tables:
`g0-base-values.tsv` (665 rows), `g0-best-partner.tsv` (supplementary, chose nothing).

**Thresholds (OBSERVED; the rule fixed in the script).** F172 in process, the design's samples removed per listed
threshold: lipid and metab remove none up to 0.9, A_5 at 0.95, ten samples at 0.99; lipid with RUV (manual, k 2, 5 %)
none up to 0.95, B_4 at 0.99. Hence **t_hi = 0.95, t_lo = 0.9** (lipid, metab) and **t_hi = 0.99, t_lo = 0.95** (lipid
with RUV). Stop rule 2 not triggered for any configuration.

**Phase 0 (R1-01; OBSERVED): passed** for lipid and metab. A managed Apply at 0.95 completed at
`<prefix>_correlation_filtered` with the in-process sample list (A_5 removed); a managed Apply then Skip completed at
`<prefix>_norm_complete`; every managed action `delivered`, no error in the tab's log.

**What base does (OBSERVED, in process; lipid, metab and lipid with RUV alike unless stated; every prediction held).**

| Case | Observed on base |
|---|---|
| I: run, Apply(t_hi), Skip | head `<prefix>_norm_complete`, `<prefix>_correlation_filtered` before it in the active history; head design 11 samples (A_5; with RUV B_4 removed); head components not the run result's; summary "original 12, after 12, removed 0"; final PCA 12 samples; DA source `<prefix>_norm_complete`, 11 samples, status current; Save done, integration object 11 samples, `study_parameters.txt` lists the selection with `min_pearson_correlation_threshold = 0.95` (0.99); record: `<prefix>_pearson_selection` at t_hi then the store step; replay completes, no statement, replayed head 11 samples = live |
| II: run, Apply(t_hi), Apply(t_lo) | the second correlation's source = the t_hi state (not the run's result); design 11 samples (A_5, B_4 out); summary "removed 1"; record: one selection at t_lo whose input is its own result; replay completes with 12 samples: **replay differs from live** |
| II0: run, Apply(t_lo) alone | 12 samples |
| C1: run, Skip | 12 samples; components = run result's; no selection step; replay 12 = live |
| C2 vs C3 | Apply(t_lo) then Apply(t_hi) = Apply(t_hi) alone (A_5; B_4): one assay, as read |
| E13: run, Apply(t_hi), Skip, Apply(t_lo) | head `<prefix>_correlation_filtered` with `<prefix>_norm_complete` after it in the active history; 11 samples; replay ends at `<prefix>_norm_complete` with 12: differs |
| E14: run, Apply(t_hi), DA, Skip, Save | DA status after the Skip "current"; Save done; integration object 11 samples |
| E10: Apply(t_hi), then an Apply with an absent grouping column | refused (FALSE, notice); head the t_hi state; fields of t_hi |
| E11: Apply(t_hi), Skip with the commit refused | head the t_hi state; fields of t_hi (lipid put back, metab untouched) |
| E9: two assays, Apply(0.5) twice | first removes A_1; second removes **A_1 and A_2** (the cascade), no error |
| E12 (managed): Apply(t_hi) held on the worker, Apply(t_lo) queued, released | the queue held (the first press active, the second queued); after both, 11 samples (A_5 out): the queued press read the t_hi head |
| I, II managed (head and design only) | as in process: I `<prefix>_norm_complete` with A_5 out; II A_5 out |

**Instrument limits (OBSERVED; VAL-54, not rerun).** (1) Managed I and II: the script read the tab's reactive values
after the test server had ended (`shiny.destroyed.error`), six errors; head, history and design were recorded first.
Not rerun: R1 runs every managed E row on base and records its value. (2) Lipid with RUV, managed: the worker-journey
harness dispatched no RUV action for these inputs (`worker inputs`: no `lipid_ruv`), so the managed run ended at
`lipid_normalized` and Apply(0.99) removed ten samples. E7's managed half is therefore not run (Not checked).

Stop rule 1 not triggered: after Apply then Skip the final state, DA source and integration input hold the filtered
samples on every route observed, and a re-Apply reads the filtered head.

## R1 (run 45604366935, tests-only commit 9dc99037)

`r1/tests.tsv`, `r1/failures.txt`. `test-kit254-omics-correlation-skip.R`: 191 passed, 101 failed, 0 errors. Every
failure is an E-row assertion at G0's base value (OBSERVED, classified by row): E1 to E4 (lipid, metab; in process and
managed), E7 (lipid with RUV, in process), E5, E6 (both routes), E9 (A_1, A_2 removed), E10, E11, E12 (source the t_hi
state, A_5 out), E13, E14 (status current, Save done; the DA tab's results kept). No G row failed (G1, G2, G3, E8,
G7, G8) and no precondition failed. `test-omics-norm-state-guards.R` with the T1/T8 set-up edit: 651 passed, 0 failed
(647 expectations at merge train 3: T1's split precondition adds them, A4).

## The fix (bff1eb38)

`R/mod_omic_norm.R`: `omicNormCorrelationStart` (beside `omicNormApplyRunStart`) returns the active history to the state
holding the run's result when the head is not the run's result (after the KIT-194 guard accepted it); it logs "Correlation
filter starts again from '<state>', the table the normalisation produced". `applyOmicNormCorrelation` calls it in process at
the click and, managed, from `start` when the actions start; `skipOmicNormCorrelation` at the click. A per-press flag
records whether it returned; a failure after it calls `clearOmicNormCorrelation` (the three fields at their initial
values). `R/utils_workflow_omics_actions.R`: `runWorkflowOmicsManagedCorrelation(..., start = NULL)` passes the actions as a
function that calls `start` unless `runner$retrying()`. No change to replay, records, kernels, texts, `man/`, `NAMESPACE`.

## R2

Run 45609913680: head bff1eb38 (the fix) against base 9dc99037 (the tests-only commit), the 47 files longest first
(`r2/compare.tsv`, `r2/check-compare.txt`, `r2/tests*.tsv`). OBSERVED: **46 files `same`**, including
`test-omics-norm-state-guards.R` (the T1/T8 edit runs on both sides), `test-kit212-da-current.R`,
`test-kit207-summary-head-guard.R`, `test-omic-norm-module.R` and `test-prot-norm-core.R` (identical failure messages on
both sides, as on every runner); **the new file `fixed`**: base 191 passed 101 failed, head 292 passed 0 failed, every
E row fixed, every G row passing on both. R CMD check: ERRORs 0, WARNINGs 3, NOTEs 3 on each side, findings identical, 0
only on head. Stop rule 3 not triggered. 385 job-minutes (25.67 core-hours).

## Acceptance criteria

| # | Criterion | Evidence | State |
|---|---|---|---|
| A1 | G0 before any code | run 45599722106; tables above | met, with the two instrument limits stated (managed I/II tab fields observed in R1 instead; managed RUV not run) |
| A2 | each E row fails at G0's value on the tests-only commit; G rows pass; edited T1/T8 pass on base; all pass on the fix | R1 45604366935; R2 head side | met (R1; head side of R2: 292 passed, 0 failed) |
| A3 | R2 base against head: every file `same` except the new file; check nothing only on head | R2 45609913680 | met |
| A4 | the merge train's full suite; `test-omics-norm-state-guards.R` count-changed by the T1/T8 edit | the train | not this branch's run (the coordinator's train) |
| A5 | NEWS, ledger, outcome with "For the change report"; `man/`, `NAMESPACE` unchanged | 964d6f65; this file | met (ledger `VAL-<merge>`, numbered at merge) |

**OD-1 (owner, decided after this implementation started: (a), VAL-247).** This change makes no replay change (the
coordinator's ruling). The replay statement under (a) is a later, separate commit with its own failing-first run and its
own comparison over the replay consumers (the 54-file scope), merged after KIT-222 (on staging since merge train 4).

## Not checked

- E7 managed (lipid with RUV through the worker): the harness ran no RUV action for these inputs (G0); why was not traced.
- The managed I/II tab fields at base in G0 (observed in R1 instead).
- G0 was not rerun after its instrument errors.
- A Skip pressed while a managed Apply runs (in process at the click; traced by reading in Re-check 1 only).
- The full suite (the merge train's); Windows; real data; how often either sequence happens in practice.
- Lipid DA refusal after a Skip on multi-assay data where a sample failed in only some assays (R1-07 (b)): by reading
  only (E9 uses constructed pairs and does not run a DA).

## For the change report

### Appendix: a defect of code new to 0.9 (regression from `main`)

- **What was wrong.** Lipidomics and metabolomics, correlation filter. (i) Pressing "Skip Filtering & Proceed to DE"
  after "Apply Correlation Filter" committed the head's record, which after an Apply is the filtered table. The tab then
  said "Removed: 0 samples" and drew every sample in the final PCA, but the differential analysis, the exported
  normalised session, the integration input, `study_parameters.txt` and `analysis_record.json` used or named the filtered
  table. Observed on F172: A_5 left out of all of them while the screen showed 12 of 12 (G0, case I). (ii) Applying the
  filter again filtered the already filtered table: a lower threshold could not bring a sample back (G0 II: 0.9 after 0.95
  kept A_5 out) while the record named the lower threshold and a replay brought it back; with two or more assays a second
  press, even at the same threshold, could remove a further sample (G0 E9: A_2, through A_1's removal in the other assay).
  A DA run between Apply and Skip stayed "current" and was saved filtered (G0 E14).
- **Who was affected.** Users of lipidomics and metabolomics development builds from 31 October 2026 (6506af27 lipid,
  4819f030 metab; the only path from a93d4c14, 1 November 2026), in tag v0.9.0; no release shipped. Every route.
- **Code.** Base e0747fc1: `R/mod_omic_norm.R:616-653` (the Skip commits the head's record, `:636-638`), `:564-609` (the
  Apply reads the head: `R/utils_workflow_omics_declarations.R:1563`; managed `R/utils_workflow_execution.R:942`). Fix
  bff1eb38: `R/mod_omic_norm.R` `omicNormCorrelationStart`, `clearOmicNormCorrelation`, and their calls in
  `applyOmicNormCorrelation` and `skipOmicNormCorrelation`; `R/utils_workflow_omics_actions.R`
  `runWorkflowOmicsManagedCorrelation(start)`.
- **What changes in results.** After Apply then Skip the analysis uses every sample; a re-applied filter at threshold t
  gives the run's result filtered once at t; record and replay agree with live. A DA run between Apply and Skip becomes
  stale and the Summary asks for it again. The record no longer shows a filter that was tried and undone (as for a Revert).
- **How to recognise affected outputs.** Records: a `<prefix>_norm_complete` step directly after a
  `<prefix>_correlation_filtered` step (case i); a `<prefix>_pearson_selection` step whose input entity is its own result
  (case ii). Outputs: the DA tables, the exported session, the integration input and final S4 file, and
  `study_parameters.txt` (a `min_pearson_correlation_threshold` line although Skip was pressed) of such sessions.
- **How to recompute.** Run the normalisation again, then Skip (or Apply once), then the DA and the Summary.
- **Evidence.** G0 45599722106 (base values), R1 45604366935 (each E row failing at them), R2 45609913680. F172 fixtures (12
  samples, one assay; E9 two assays with constructed pair correlations); limits under Not checked.

### Main part: legacy defect confirmed on `main` (436cdbe8; `git show`, not checked out)

On `main` (and v0.5.0-rc1) the Apply read the run's result (`main:R/mod_lipid_norm.R:1697-1701`) and the Skip saved the run's
result as `lipid_norm_complete` (`:1752-1781`; metabolomics `main:R/mod_metab_norm.R:1752-1768`), and the DA read the state
manager's current state (`main:R/mod_lipid_da.R:811-814`), so the analysis after Apply then Skip used every sample. But the
Skip did not reset `correlation_filtered_obj`: the tab's summary (`main:R/mod_lipid_norm.R:1792` on) kept reporting the
Apply's removals and the final PCA (`:1842` on) drew the filtered samples. A user who reported "k samples removed" from that
screen reported a filter the analysis did not apply. The development line corrected that display on 2026-05-05 (b5f26546,
ac052181) and broke the analysis on 2026-10-31; after this fix screen and analysis agree.
