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
