# KIT-174 part (a): outcome

Ticket revision 1.4. Branch `kit-174` from `v0.9-staging` 0cf0e440. Part (b) is held by the coordinator (ticket
section 19) and is not part of this merge.

## Commits

| Commit | Content |
|---|---|
| b2e7c4c2 | tests only: `tests/testthat/test-kit174-ruv-grouping.R` (T1, T2, T9, T10); `feasibility.md` (F2 for part (a)) |
| 002957c2 | R1a records |
| 2b89963c | the fix: help text (`R/utils_omic_gui_config.R`); the warning and its helpers (`R/mod_omic_norm.R`); the two proteomics call sites (`R/mod_prot_norm_workflow_helpers.R`); NEWS; ledger row `VAL-<merge>`, UVR-021 closed, its open-decisions entry removed; parked row 20 done |
| 9f9a16b1 | A3 script `norm-panel.R` |
| f121b057 | `wording.md` |
| 53332f55 | R2a and panel records; the RUV-III-C name check |
| 24960027 | merge of `v0.9-staging` cf8af245 (staging changed only documents since 0cf0e440; the ledger's VAL-191 to VAL-196 first, this ticket's row after) |

## Runs (observed)

See `runs.tsv` for ids, commits, minutes and core-hours.

- **R1a** (45494891087, tests-only b2e7c4c2): 37 passed, 13 failed, 0 errors. T1 failed at its three expectations
  (text not identical; "select 'batch'" present; no Model Formula pointer). T10 failed first at "the function exists"
  in each of its four blocks, then at 0 warnings where 1 was expected (classic crossed, partly crossed, gate-on
  `groups_all_features` with declared technical replicates; DIA with limpa crossed; shared frame crossed and its log
  line). T2 and T9 passed. T2's values: before correction median |group| 0.9779, |batch| 0.6890; `group` as replicate
  sets 0.9771 and 0.0015; `batch` as replicate sets 0.0013 and 0.6893 (30 controls each). The thresholds were met on
  the first run.
- **R2a** (45496986018, fix 9f9a16b1 against 0cf0e440, 21 files, 3 shards per side, check): 20 files `same`; the new
  file `only-on-head`, 56 passed, 0 failed. `test-prot-norm-core.R` (306) and `test-gui-entry-proteomics.R` (4) fail
  identically on base and head (stored laptop references, as in KIT-156's record). R CMD check: 0 ERROR, 3 WARNING,
  3 NOTE on both; findings only on head: none.
- **R2a-panel** (45496979801): the options column of the proteomics normalisation tab, 158 lines on each side; the
  only changed line is the help-block span of the RUV grouping help (`norm-panel.diff`). The script checked that the
  renderer file is identical on base and head and sourced the base commit's configuration over the head namespace.
- **R3a** (45507740071, merge 24960027 against staging cf8af245, all 479 files, 20 shards per side): 478 files
  `same`, the new file `only-on-head` (56 passed), 0 `new-failure`, 0 `message-changed`; 35 files fail identically on
  base and head (`r3a-full-suite-compare.tsv`; `r3a-full-suite-compare-failures.txt`: no non-identical failure).
  `test-gui-entry-proteomics.R`, `-lipidomics.R`, `-metabolomics.R` are `same`. The longest job took 76.8 minutes
  (head shard 11), above the 70-minute per-job ceiling of ticket section 11; the step's total, 80.56 core-hours, is
  inside its ceiling of 186.7. Recorded as a departure; the result is valid and was not re-run.

## Budget

Part (a): R1a 0.68 + R2a 19.10 + R2a-panel 0.82 + R3a 80.56 = 101.16 core-hours, against 76.0 expected and 232.7 ceiling.
Above expected: R2a took 19.10 against 12.0 expected (8 jobs, inside 34.0); R3a 80.56 against 61.0 (the readiness note
of 2026-11-25 had re-priced a full suite at 67 to 80). No correction pass was used.

## Acceptance criteria (part (a))

| # | Criterion | Evidence | Status |
|---|---|---|---|
| A1 | F2 answered before (a)'s code | `feasibility.md` (committed with the tests, before the fix) | met |
| A2 | Strings as proposed; `wording.md` lists them | T1 and T10 assert the exact strings; `wording.md` | met |
| A3 | Panel diff shows only the help text | `norm-panel.diff` (R2a-panel) | met |
| A4 | T1, T10 fail on the tests-only commit at a stated expectation and pass on the fix; T10 first asserts the function exists; T2, T9 pass on both | R1a; R2a (new file 56 passed) | met |
| A6 | Targeted files and dependents, check (no finding only on head), full suite base vs head; GUI-entry files unchanged | R2a; R3a (every difference is the new test file; the three GUI-entry files `same`) | met |
| A7 | NEWS under the unreleased heading, no version change; no Rd change (no exported function) | diff | met |

## Observed against inferred

Observed: the run results above. Inferred, not run in the app: that the warning reads the right design in a live
session. The fix reads the head state's design and group column (`omicNormRuvHeadDesign`, the same metadata read as the
existing `omicNormRuvHeadTechnicalReplicates`) and falls back to the session's design table only without a session
state. The session's saved design (`workflowData$design_matrix`) would have been wrong whenever quality control removed
samples: a removed sample would still be counted, so a value could be reported as spanning groups (or the count "n of
its m values" could be off) because of a sample that is no longer in the table RUV-III reads. The tests exercise the
fallback only (their managers are not session adapters); the head-state read is exercised by the dependents that run
real session managers with `group` (no warning expected, none raised: `test-prot-design-log2.R`,
`test-workflow-saved-qc.R`, both `same`).

## Not checked

- The warning in a live application session (the application is not started on this machine; no GUI-entry lane chooses
  a column other than `group`).
- The head-state read with a column that spans groups in a real session store (only the fallback path is tested with a
  spanning column).
- The managed (worker) path of the classic route with a spanning column: the warning is raised in the controller before
  dispatch, on both paths, by reading the code; the test uses the in-process path.
- Part (b): not started (held).

## For the change report

- **What was wrong.** The proteomics normalisation tab's help for the RUV grouping variable called it the factor "to
  preserve" and advised "For TMT workflows with batches, select 'batch' to remove batch effects". In every RUV design
  the package offers, choosing `batch` makes the samples of one batch replicates of each other, so RUV-III builds its
  correction from the differences between the groups and keeps the batch differences; the controls, the automatic
  search and the sample-correlation filter also follow the batches. Users of proteomics (any route; the advice named
  TMT) who followed the advice were affected; nothing warned (lipidomics and metabolomics have the same selector with a
  neutral help text and also had no warning).
- **What changed.** The help states what the column does and points a known batch to the model formula. At the run
  action of the normalisation tab (classic proteomics, DIA with limpa, lipidomics, metabolomics) a warning
  notification and log line appear when a value of the chosen column holds samples of more than one group and the
  column gives RUV-III's replicate sets in that run. The run proceeds; RUV's calculation, defaults and choices are
  unchanged, so no result changes and earlier records replay as before.
- **Code.** Base 0cf0e440: help `R/utils_omic_gui_config.R:2039-2041`; run actions `R/mod_prot_norm_workflow_helpers.R`
  `runProtNormNormalizationWorkflow` `:503-507`, `runProtNormRuvCovariatesObserver` `:93-95`; `R/mod_omic_norm.R`
  `runOmicNormPipeline` `:355-361`. Fix 2b89963c: help `R/utils_omic_gui_config.R:2039-2058`; new functions
  `omicNormRuvGroupingSpansGroups`, `omicNormRuvGroupingWarningText`, `omicNormRuvGroupingIsReplicateSets`,
  `omicNormRuvHeadDesign`, `omicNormWarnRuvGrouping` in `R/mod_omic_norm.R` (before `runOmicNormPipeline`); calls at
  the three sites.
- **Evidence.** T2's simulation (seed 174; 2 groups x 2 batches x 3 runs; 300 proteins; batch loadings N(0, 1); group
  effect 1 log2 with random sign; noise sd 0.2; k = 1; controls from the ANOVA at 10 %): with `group` as replicate sets
  the median batch coefficient fell from 0.689 to 0.0015 and the group effect stayed 0.977; with `batch` the batch
  coefficient stayed 0.689 and the group effect fell to 0.0013 (run 45494891087). Limits: one simulated design, k = 1,
  complete data; it demonstrates the help's sentence and does not validate RUV. T9 shows the classic route's
  differential analysis runs with `~ 0 + group + batch`. T10 shows the warning's conditions at each site.
