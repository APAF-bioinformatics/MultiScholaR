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
(`r2/compare.tsv`, `r2/check-compare.txt`, `r2/tests*.tsv`). OBSERVED: **46 files `same`** (seven of them fail on both sides, with the same 446 failure messages byte for byte: test-prot-norm-core.R 306, test-omic-da-module.R 53, test-workflow-trx-declarations.R 33, test-omic-summary-module.R 28, test-workflow-omics-declarations.R 23, test-workflow-operation-registry.R 2, test-workflow-kit-harness.R 1; all fail identically in staging's merge trains), including
`test-omics-norm-state-guards.R` (the T1/T8 edit runs on both sides), `test-kit212-da-current.R`,
`test-kit207-summary-head-guard.R`, `test-omic-norm-module.R` and `test-prot-norm-core.R` (identical failure messages on
both sides, as on every runner); **the new file `fixed`**: base 191 passed 101 failed, head 292 passed 0 failed, every
E row fixed, every G row passing on both. R CMD check: ERRORs 0, WARNINGs 3, NOTEs 3 on each side, findings identical, 0
only on head. Stop rule 3 not triggered. 385 job-minutes (25.67 core-hours).

## Correction pass (AUDIT-KIT-254-impl-20261126) and the follow-up under OD-1 (a) (VAL-247; ticket r1.3 section 18)

- **Merge.** `origin/v0.9-staging` c4505c7f merged into `kit-254` (1a5291d3): conflicts only in `NEWS.md` (both entries
  kept) and the ledger (staging's rows VAL-230 to VAL-261 kept, KIT-254's after them); no code conflict.
- **First statement superseded.** 6c43b26b (per-assay columns, Skip after Apply only) was replaced after the audit (I08,
  point 6; coordinator's ruling): its R1b (45619709256) failed as predicted, its R2b (45623518474) was cancelled with no
  result read; 993c454b reverts it.
- **I01 (E7 after KIT-224).** R1c (45624689427, merged tree without statement code): the new precondition failed for
  lipid with RUV (Apply(0.99) now kept B_2 and B_3 only; lipid and metab held). G0b (45628705535, G0's rule unchanged, on
  the merged tree): no listed threshold removes some samples and keeps two per group (none up to 0.95, ten at 0.99):
  **stop rule 2** for E7. Declared test-only edit (44e3e910): lipid with RUV at 0.5 / 0.3, nothing removed; E7 keeps
  its state and record assertions; its sample difference is not shown on the fixture after KIT-224 (it was, at e0747fc1:
  B_4). Under KIT-248's G0b precedent (m + m2) / 2 from `g0b/g0b-best-partner.tsv` would be about 0.960 (A_1 0.9573, B_4
  0.9628); not chosen here (the coordinator's call).
- **R1d (45631887861; staging c4505c7f plus the test file only; predictions committed first), OBSERVED as predicted:**
  257 passed, 114 failed: every E row fails at base, E7 on its state and record assertions; E8, E15 to E18 at "no
  statement"; no G row or precondition failed. E17's first run (R1c) was an instrument defect (the replay's inputs
  matched by assay order; fixed by SHA-256 in 44e3e910).
- **Statements (904a5514).** `replayWorkflowAnalysisRecord` returns `correlation_statement` and raises it as a warning of
  class `workflow_replay_correlation_skip`; the replay is unchanged. Per lipid/metab prefix: a selection whose input is not
  the result of the step before it (the filter applied more than once, including Apply, Skip, Apply, whose selection's
  input is the skip state, G0 E13) gets the second statement and no other; otherwise the skip step directly after the
  selection gets the first, naming the samples left out of the design with the assays they failed in, and none when no
  sample left the design. The register's point 6 said Apply, Skip, Apply records carry "input equals result"; G0's E13
  record shows `input_is_result=FALSE`, hence the definition above. Texts are drafts (`wording.md`).
- **The re-apply signature (stands for the closure check to judge).** "A selection whose input is not the result
  of the step before it", not the register's "input equals result": G0's E13 record (Apply, Skip, Apply on base) has
  `input=b30f32faacad;result=1a202f6dc924;input_is_result=FALSE` (lipid; metab the same), its input being the skip
  state; under the register's wording that record would get the first statement and a replay that differs from the
  analysis that was run (G0 E13: replay 12 samples, live 11).
- **E7 (coordinator's ruling; a declared-comparison change, for the owner):** t_hi 0.960024 = (0.957297 + 0.962751) / 2,
  the two smallest best within-group partner correlations G0b observed on the merged tree; t_lo 0.95 by G0's rule;
  PREDICTED to remove A_1 only; E7's sample assertions restored (db727407).
- **E17 both ways (coordinator's ruling):** the opposite error (an assay wrongly omitted) is caught by a second case in
  which A_1 fails in both assays and must be named with both (db727407).
- **R2c:** predictions in `r2c-predictions.md` (committed before the request). **R2c (45635516160, base 44e3e910, head 904a5514; 506 job-minutes), OBSERVED:**
  57 files `same` (the seven with pre-existing failures at staging's counts; `test-omics-norm-state-guards.R` 651/0 and
  `test-workflow-analysis-record.R` 351/0 on both sides; the replay consumers unchanged); R CMD check identical, none only
  on head. The new file: base 351 passed 20 failed, head 369 passed **2 failed: E17 lipid and metab**, against the
  prediction of 0. E8, E15, E16, E18 and every guard pass on head (statements in the head log). E17's statement names
  the right sample, left out of the design ("leaves out A_1"), but its assay detail is false: "(below the threshold in
  LCMS_Pos, LCMS_Neg)", while A_1 failed in LCMS_Pos only. **Cause (read):** `workflowReplayDesignSamples` intersects each
  assay's columns with the state's own design, and A_1 has left the after-state's design, so A_1 seems to leave every
  assay. **Repair (not made):** compare the after-state's assay columns as stored (`colnames(assay)`), intersected with
  the before-state's design only. **Reading of the budget (coordinator's ruling, recorded here):** R2c was the statement code's first comparison
  run; VAL-247 made the statement its own later commit ("the fix ships first ... and the statement follows"), so the
  correction pass already used belonged to the main implementation's audit, and the statement keeps its own one
  correction: the one-place fix and R2d below. No further correction after it.

- **R2d (45658805131, base db727407, head 267c9507; 3 files, no check, 57 job-minutes; dispatched by hand once after
  the dispatcher run 45649854076 failed on a GitHub HTTP 500 and started none), OBSERVED, every prediction held
  (`r2d-predictions.md`):** the new file base 377 passed / 2 failed (exactly E17 "pos", lipid and metab: "below the
  threshold in LCMS_Pos, LCMS_Neg") to head 379 / 0 (`fixed`; "pos" names LCMS_Pos only, "both" names both); E7's
  restored sample assertions pass on both sides, so Apply(0.960024) on lipid with RUV removed A_1 only, as predicted;
  `test-omics-norm-state-guards.R` 651/0 and `test-workflow-analysis-record.R` 351/0 `same`. The statement's one
  correction is closed; the merge train repeats the full scope and R CMD check.
- **Budget:** 1,108 job-minutes in all = 73.87 core-hours against the ticket's 115.00 ceiling.

## Closure check (AUDIT-KIT-254-impl-20261126, "Closure check (2026-11-26)"): ready for a merge train

No code or test repair. Record repairs in this commit: C1 (the re-apply recognition signature, here and in ledger row
`VAL-298`), C2 (NEWS's replay sentence), C3 (the second statement's draft, `wording.md`). VAL-247's correction is
written by the coordinator at landing from the register's text.

Notes:
- **C4 (KIT-255).** Statement 1's assay detail is read from the columns the replayed filter removed. If KIT-255 (option
  (a)) removes a failing sample from every assay without a replay rule for older records, the detail would name every
  assay; E17 "pos" will fail in KIT-255's R2 and catch it. KIT-255 declares the move and derives the detail from the
  pairs (per assay, no pair at or above t), or its replay marker keeps the earlier rule for older records.
- **C5.** `workflowReplayDesignSamples` hydrates the two full states only to read column names; acceptable on a path
  reached only by the validation harness and tests (no measured cost).
- **C6 (a line for KIT-224's change report).** On the lipid-with-RUV fixture (F172, manual RUV, k 2, 5 %), threshold
  0.99 removed only B_4 before KIT-224 (G0 45599722106: every other sample's best within-group partner at least 0.99) and
  removes 10 of 12 samples after it (G0b 45628705535: the smallest best partner 0.957), because the old per-feature level
  shift had raised replicate correlations; a threshold chosen under it filters differently.
- **C7 (an inferred link).** E7's sample assertions at 0.960024 were never run on base code (R1d ran E7 at 0.5, nothing
  removed, failing on its state and record assertions; R2d has the fix on both sides). Failing first at the final
  constants follows from R1d's observed state and record failures, the precondition (the first Apply is the same code on
  both commits) and base's Skip committing the filtered head (R1, R1d). No run owed.

**The merge train's scope and predictions** (run-tests; base = the train's base, head = the train's head; `check: true`;
4 shards; longest first; `file_timeouts` 120 for `test-prot-norm-core.R`; the train's env):
- On staging abfc99ac, **61 files**: the ticket's grep (46: R2c's 44 plus `test-kit243-summary-current-outputs.R` and
  `test-kit252-managed-failure-text.R`), the three frame files, the new file, `test-kit214-sample-columns.R`, and the ten
  replay files (`test-integration-model.R`, `test-integration-mofa-model.R`, `test-kit200-prot-da-no-value-proteins.R`,
  `test-kit201-treat-moderation.R`, `test-kit205-enrichment-background.R`, `test-kit215-not-tested-lists.R`,
  `test-prot-annotation-uniprot-empty.R`, `test-prot-missing-value-counts.R`, `test-workflow-analysis-record.R`,
  `test-workflow-omics-binding-steps.R`). With KIT-257 merged first, **62** (`test-kit257-norm-own-files.R`).
- Predictions: the new file absent on base, head 379 passed, 0 failed; `test-omics-norm-state-guards.R` count-changed
  647/0 to 651/0 (the T1/T8 set-up edit); every other file `same` at the base's counts (the files with pre-existing
  failures, `test-workflow-analysis-record.R` 351/0 with 4 warnings, `test-kit243-summary-current-outputs.R` 246/0,
  `test-kit257-norm-own-files.R` when KIT-257 is on the base, the replay consumers with no
  `workflow_replay_correlation_skip` warning added); R CMD check nothing only on head.
- **M2 against KIT-257** (`R/mod_omic_norm.R`): KIT-257's `omicNormOwnFiles` and KIT-254's `omicNormCorrelationStart`
  and `clearOmicNormCorrelation` are inserted at the same place, after `omicNormApplyRunStart`. Resolution: keep
  KIT-257's `omicNormOwnFiles` whole, then a closing `}` and a blank line, then KIT-254's two functions (the shared `}`
  that follows closes `clearOmicNormCorrelation`). No semantic overlap: KIT-257 changes `omicNormApplyRunStart` and the
  proteomics correlation, not `applyOmicNormCorrelation` or `skipOmicNormCorrelation`.

## Acceptance criteria

| # | Criterion | Evidence | State |
|---|---|---|---|
| A1 | G0 before any code (E7's thresholds re-derived on the merged tree: G0b, stop rule 2) | run 45599722106; tables above | met, with the two instrument limits stated (managed I and II: head, history and design recorded in G0; the DA source, Save, record and replay observed on base by R1's managed rows; the tab's correlation summary and PCA on the managed route observed on no run: the Skip sets them in process on every session, `R/mod_omic_norm.R:629-673`, as G0 observed in process; managed RUV not run) |
| A2 | each E row fails at G0's value on the tests-only commit; G rows pass; edited T1/T8 pass on base; all pass on the fix | R1 45604366935; R2 head side | met (R1; head side of R2: 292 passed, 0 failed) |
| A3 | R2 base against head: every file `same` except the new file; check nothing only on head | R2 45609913680 | met |
| A4 | the merge train's full suite; `test-omics-norm-state-guards.R` count-changed by the T1/T8 edit | the train | not this branch's run (the coordinator's train) |
| A5 | NEWS, ledger, outcome with "For the change report"; `man/`, `NAMESPACE` unchanged | 964d6f65; this file | met (ledger `VAL-298`, numbered at merge) |

**OD-1 (owner: (a), VAL-247).** The main fix makes no replay change (the coordinator's ruling); the replay statement
under (a) is the separate follow-up above, with its own failing-first run (R1b) and comparison (R2b).

Conventions: 20 new lines of the fix exceed R.md's 100 characters, following the file's idiom (230 such lines); a
recorded divergence under R.md's existing-code scope rule (I13).

## Not checked

- E7 on the merged tree runs at a constructed threshold (0.960024, the coordinator's ruling after G0b's stop rule 2), not at a listed one.
- The follow-up statement's assay detail on multi-assay records: wrong at 904a5514 (R2c, E17); repaired in 267c9507 (R2d).

- E7 managed (lipid with RUV through the worker): the harness ran no RUV action for these inputs (G0); why was not traced.
- The managed route's correlation summary and PCA reference (base and fix): read only.
- G0 was not rerun after its instrument errors.
- A publication failure of the managed selection after its science was committed: the tab's correlation fields are
  cleared while the head holds the new filtered table, until the retry (any managed action) delivers the Apply's result
  (`R/utils_workflow_execution.R:2072-2083`, `:845`, `:853`; `R/mod_omic_norm.R:599`). Transient; read, not run (I02).
- G7 presses at the run's result; an Apply pressed at a correlation-step head while a display runs is shown only through
  E12's shared mechanism (start-time resolution) (I03).
- After the fix a re-applied step's input is the run's result (its previous state is the state the cut returns to);
  read, and implied by E6 and E13 (replay equals live), not asserted on the record (I12).
- Why the managed harness ran no RUV step (phase 0 copied lipid's pass for lipid with RUV; I11).
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
  `<prefix>_correlation_filtered` step (case i); a `<prefix>_pearson_selection` step whose input is not the result of the step before it (its own result after a second
  press; the skip state after Apply, Skip, Apply) (case ii). Outputs: the DA tables, the exported session, the integration input and final S4 file, and
  `study_parameters.txt` (a `min_pearson_correlation_threshold` line although Skip was pressed) of such sessions. The
  first signature also matches a record a script makes, after the fix too, by calling the wrappers in that order; it
  describes the record, not the version (I12).
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

## Addendum 2026-11-26 (merge train 15): the `lipid_ruv` thresholds re-derived on the train head

On merge-train-20261126h (4308c5ae, which holds KIT-190's decided RUV control rule, VAL-110 and VAL-126) the E7 `lipid_ruv`
expected sample lists no longer held (focused audit `AUDIT-TRAIN15-kit254-e7-20261126`; run 45681065939: three failures,
the head removing A_5, A_6, B_1 at 0.960024 where 11 of 12 samples were expected). By the coordinator's ruling the
constants were re-derived from G0b's script, run unchanged on 4308c5ae (run 45704341824, 12.7 job-minutes, exit 0,
tables in `g0c/`): G0's rule gives `t_hi` = 0.8 (Apply removes B_1) and `t_lo` = 0.7 (removes none; stop rule 2 did not
apply). `test-kit254-omics-correlation-skip.R`: `kit254T$lipid_ruv` = c(hi = 0.8, lo = 0.7), `kit254HiRuns$lipid_ruv` =
all but B_1, `kit254LoRuns$lipid_ruv` = all. Tests-only edit; declared prediction: the file 379 passed, 0 failed.
Not checked: the head's other lipid RUV samples by hand; the metab RUV case does not exist in this file.
