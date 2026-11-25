# KIT-261 outcome: the normalisation Run button shows that a step is running

Ticket KIT-261 r1.2 (r1.1 plus the implementer's line-number note). Branch `kit-261` from v0.9-staging **0a774e0b**
(merge train 11, KIT-257). Commits: tests-only **4de95bbb**; fix **8ee87839** (`R/mod_omic_norm.R` only); records after.
No R on the laptop; both runs on GitHub-hosted runners (`run-tests`). Nothing under `data/benchmark/` read.

## What changed

- `R/mod_omic_norm.R:114` (8ee87839): in `runOmicNormServer`, after the run observer, one call
  `registerOmicNormRunButtonState(session, <omic type>, <run label>)` when the tab has the `actions` panel; no new local
  in the frame.
- `:263` `omicNormRunButtonTexts()`: the busy label, draft "Please wait: a step is still running..." (VAL-41).
- `:275-289` `registerOmicNormRunButtonState`: resolves `workflowOmicsManagedActions(<omic>)` at server start and
  registers nothing without a runner (in process, every GUI-entry and e2e lane). Otherwise one observer of
  `runner$step$active()`: on a change of state only (first state idle) it calls `shinyjs::toggleState("run_normalization",
  condition = idle)` and `shiny::updateActionButton(session, "run_normalization", label = <own label or busy label>,
  icon = shiny::icon("play"))`. A runner without a step (the test double of `test-omic-norm-module.R`) counts as idle.
- KIT-194's refusal (`omicNormRunStart`, `:984-987` at the base, text `:789`) is unchanged and remains the guard.
- Not changed: `R/mod_omic_norm_panels.R`, `R/utils_workflow_omics_actions.R`, `R/utils_workflow_execution.R`,
  `NAMESPACE`, `man/`, DESCRIPTION, any stored reference. The log line "Starting normalization pipeline..." before a
  refused press is out of scope (E-2, the coordinator's ruling).

## Runs (OBSERVED)

| Run | Binds | Result | Job-min |
|---|---|---|---|
| R1 (= G0) 45671049216 | 4de95bbb (tests-only), new file only | 85 passed, 20 failed, 0 errors; every prediction in the file header held | 11.5 |
| R2 45673434696 | head 8ee87839 vs base 4de95bbb, 43 files x 4 shards per side, check both | 42 files `same`; new file `fixed` (85/20 -> 105/0); check ERRORs 0, WARNINGs 3, NOTEs 3 on both, 6 findings identical, none only on one side | 402.3 |

Core-hours: 413.8 job-minutes x 4 / 60 = **27.59** (expected 34.08, ceiling 124.33). No instrument rerun, no
correction pass. Env of R2: `MULTISCHOLAR_TRX_MANIFEST_NOTE=1`, `MULTISCHOLAR_KIT105_W3_DUMP=1` (test-prot-norm-core.R
in scope).

R1 in detail (`r1-tests.tsv`, `r1-failures.txt`): L1, L2, L3 passed (the runner holds and releases an action; a failed
poll is told "Error applying lipidomics step: ..." and the request's `failed` called; the managed harness dispatches
`lipid_qc_figures` at tab entry, holds it, delivers it, and one press completes a run). G1 to G5 passed, including G2:
a press while the entry figures were held gave exactly one error notification with KIT-194's text, no new reservation,
head, history and `normalization_complete` unchanged. Every precondition expectation passed. The 20 failures were
exactly the predicted E-row button expectations (E1/E2 lipid 4, metab 4; E3 3; E5 2; E6 trx 4; E7 3): no
`toggleState` or `updateActionButton` call for `run_normalization` on the base.

R2 scope: the ticket's grep at the branch point gives 34 files (the 32 of the ticket plus `test-kit257-norm-own-files.R`
and `test-kit259-pearson-histogram-figure-failure.R`, as section 10.3 foresaw), plus the 8 helper-driven files and the
new file: 43 (the ticket's "41" at abe2199c). Files that fail on both sides (laptop-recorded references:
test-prot-norm-core.R 306, test-omic-da-module.R 80, test-omic-summary-module.R 28, test-gui-entry-transcriptomics.R 13,
test-module-ci-trx-{summary,da,norm}.R 1 each, test-workflow-kit-harness.R 1, test-workflow-ruv-pearson-diagnostic.R
2 errors) fail identically: `same`, as predicted (`r2-predictions.tsv`, committed c81e0bf1 before the run).

## Acceptance criteria

| # | Criterion | Met | Evidence |
|---|---|---|---|
| A1 | R1: L1-L3 pass; G2 shows the refusal on base with its precondition met; every G row passes; every E row fails at its base value | met | R1 45671049216 |
| A2 | The fix: every E, G and L row passes | met | R2 head: new file 105 passed, 0 failed |
| A3 | R2 base against head: every file `same` except the new file; check identical | met | R2 45673434696 `r2-compare.tsv`, `r2-check-compare.txt` |
| A4 | NEWS, the short ledger row, the wording rows; this file with the divergence, "For the change report", "Not checked"; KIT-237 told of H-1/H-2 | met | `NEWS.md` (### Fixes), `docs/methods-decision-ledger.md` VAL-<train16-1>, `wording.md`; hand-off note in `tickets/kit/KIT-237.md` (2026-11-26, read, not edited) |

## Recorded divergence

`conventions/R-shiny.md` asks for `shiny::ExtendedTask` with `bslib::input_task_button` for long work. The module keeps
its managed actions (its own recorded divergence, `R/mod_omic_norm.R:28-29`, under the existing-code scope rule of
`conventions/R.md`); this change only reflects their state on the existing `actionButton`. `bslib::input_task_button`
would change the button's markup and the frozen UI digests.

## Stated limit (for users)

A user who presses Run within about a second of opening the tab, before the disable reaches the page, still sees
KIT-194's refusal "A step is still running. Press "Run Normalization Pipeline" again when it has finished."; by then the
button is already disabled and says a step is running. Nothing else changes for that user. Not fixed here (the
coordinator's ruling E-7: no disable before the entry dispatch).

## Observed vs inferred

- OBSERVED (testServer, R1/R2): the calls the server makes (`toggleState` condition, `updateActionButton` label) in
  every E case; the refusal kept (G2); the in-process tab makes no call (G1).
- INFERRED, not observed: that shinyjs and the action-button binding render the greyed state and label in a browser,
  and the race in which a click sent before the disable arrives is refused (not testable under testServer).

## Not checked

- No browser, no Windows: the greyed appearance, the label rendering, the race of ticket 2.6, the Windows timing and the
  Windows bundle's shiny version (KIT-237's WJ' with H-1 and H-2 observes these).
- The DOM id the two calls address. The new tests replace `shinyjs::toggleState` and `shiny::updateActionButton` with
  recorders, so the module prefix is inferred, not observed: shinyjs adds it inside a module (`asis = FALSE`), as
  `R/mod_prot_import_path_helpers.R:20-28` already relies on, and `updateActionButton` adds it through the module
  session. The first browser observation is KIT-237 WJ' (H-2 records the seconds the button was disabled at S6.1). A
  press refused there would show a defect.
- Whether Shiny writes a custom message at once or at the end of the flush.
- A UI rebuilt by `output$dynamic_tabs` while a step runs shows an enabled button until the next change of state (today's
  behaviour; the guard refuses); not tested.
- The per-step deadline under the Windows plain launcher (without it the button stays disabled as long as today's press
  would be refused).
- Partners: KIT-254 landed on v0.9-staging in merge train 13 (988512e6), after this branch point. Its hunks in
  `R/mod_omic_norm.R` (`applyOmicNormCorrelation`, `skipOmicNormCorrelation`, two functions after
  `omicNormOwnFiles`) do not overlap this ticket's. This ticket's merge train runs
  `test-kit254-omics-correlation-skip.R`, which the grep picks up at the train's base. KIT-247, KIT-242 and KIT-255
  are not merged: there is no line overlap, and the second to merge runs the other's new test file.
- Real data.

## Implementation audit (AUDIT-KIT-261-impl-20261126): ready for a merge train

**Train scope (I01).** KIT-261 (`kit-261` d6ab5b1b). The scope is:
- the ticket's grep recounted at the train's base: 35 files at staging 988512e6 and at merge-train-20261126h
  069c2832, that is the 34 of R2 45673434696 plus `test-kit254-omics-correlation-skip.R`;
- the 8 helper-driven files: `test-gui-entry-{lipidomics,metabolomics,transcriptomics}.R`,
  `test-workflow-omics-worker.R`, `test-workflow-omics-residency.R` and `test-module-ci-trx-{norm,da,summary}.R`;
- `test-kit261-norm-run-button.R`.

That is **44 files**, with env `MULTISCHOLAR_TRX_MANIFEST_NOTE=1` and `MULTISCHOLAR_KIT105_W3_DUMP=1` and `check: true`.
Predictions: `test-kit261-norm-run-button.R` is only on head, with 105 passed and 0 failed; every other file is `same`,
with identical failure messages; `R CMD check` is identical.

**Notes, no action (I06 to I08).**
- I06: the default label and icon are restated in the helper (`R/mod_omic_norm.R:285-286`) beside the panel's
  (`R/mod_omic_norm_panels.R:70-72`); no frame omic sets `run_label` today, so the restored label is right on every
  route. If the panel's default ever changes, change both.
- I07: `test-lipid-12-norm-module-contracts.R`, `test-metab-03b-norm-module-helper-characterization.R`,
  `test-metab-03c-norm-module-seams.R` (0 expectations) and `test-prot-05c-norm-module-seams.R` (4 skipped) ran no
  expectation on either side in R2; their `same` carries no weight (pre-existing).
- I08: G1 (in process) covers lipid and metab, not trx; for trx in process the helper returns at `:277` as for the
  others (the runner is NULL); that trx sends no call in process is INFERRED.

**Later mergers (I09).** Any ticket merging after KIT-261 that edits `R/mod_omic_norm.R`,
`R/utils_workflow_omics_actions.R`, `R/utils_workflow_execution.R`, `R/mod_omic_norm_panels.R` or the norm groups of
`R/utils_omic_gui_config.R` includes `test-kit261-norm-run-button.R` (the grep picks it up). Named: if KIT-255, KIT-247
or KIT-242 merge after KIT-261, their merge checks run `test-kit261-norm-run-button.R`.

**For KIT-237 (I10).** H-2 waits on the Run button's test id being visible and enabled (the route presses by test id;
the press script treats `disabled`, class `disabled` or `aria-disabled` as disabled), not on the label text, which is a
draft. H-2's "seconds disabled" is an observation, not a criterion, and the only planned browser observation of the
disable.

## For the change report

Appendix entry only: code new to 0.9, nothing legacy.

On the 0.9 managed application the lipidomics and metabolomics normalisation tab drew its pre-normalisation QC figures
on the worker when opened (KIT-078) and, since KIT-194 (09f05774, 2026-11-24; not in tag v0.9.0, where the press waited
behind the figures), refused a run while any step of the omic was running. `main` has no managed application and drew
the figures in the session. The first press after opening the tab was therefore refused although the button looked
ready (KIT-237 WL, three of three runs on Windows: 45605975875, 45619531044, 45637098269); transcriptomics' Run was
refused likewise while any of its steps ran. No result was affected: the refusal returned before any state change. The
button now shows that a step is running (disabled, "Please wait: a step is still running..."); a press within about a
second of opening the tab can still meet the refusal. References: at the base 4de95bbb (= 0a774e0b for `R/`),
`R/mod_omic_norm.R:82-85` (entry observer), `:206-208` and `:225-234` (entry dispatch), `:456-460` (the run rule's call),
`:984-987` (the refusal); at the fix 8ee87839, `:82-85` (unchanged), `:208-210` and `:227-236` (entry dispatch),
`:488-492` (the run rule's call), `:1016-1019` (the refusal, unchanged in content), and the new `:114`, `:263`,
`:275-289`. Evidence: R1 45671049216, R2 45673434696 (testServer on the real runner, step
and tab with the worker's entry point in process; poll held by the test).
