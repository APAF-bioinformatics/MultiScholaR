# KIT-180 outcome: label-free and TMT go through the protein-level QC

Ticket `tickets/kit/KIT-180.md` revision 1.4 (git-ignored). Branch `kit-180-protein-qc-routing` from `v0.9-staging`
a60c4e9a (version 0.9.1, unchanged). Every run is listed in `runs.tsv`; gates F1 to F4 in `feasibility.md`.

## What changed

- `R/utils_omic_gui_config.R`: `chain_rules = list(qc_bypass = c("TMT", "LFQ"))` removed from the proteomics shell;
  its roxygen (`@noRd`) says TMT and LFQ open the QC tab after the design (VAL-26).
- `R/mod_omic_workflow.R`: comments (`@noRd`) say no shell sets `chain_rules$qc_bypass` since KIT-180; the generic
  branch is kept.
- Tests: T1, T2, T3, T5, T6 (section 10); the old block "bypasses QC for TMT workflows" removed in the fix commit.
- Browser lane: `test-e2e-proteomics-tmt-lfq.R` waits for QC `pending` and Normalisation `disabled` after the design,
  presses the four protein QC steps (`e2e_run_prot_protein_qc`) and then waits for `complete`/`pending`; the four
  steps are one list (`e2e_prot_protein_qc_steps`) read by both the DIA and the protein-route helpers.
- `docs/qa/e2e-proteomics-tmt-lfq.md` sentence corrected.
- `tests/testdata/gui-entry/prot_tmt.json`, `prot_lfq.json`, `prot_lfq_fragpipe.json` re-frozen from R4b.
- NEWS under `## MultiScholaR (unreleased; numbered at merge)`; methods-ledger row `VAL-<merge>` after VAL-49.
- `man/`, `NAMESPACE`, `DESCRIPTION`, lock files, `helper-gui-entry.R`, `helper-workflow-checkpoint-fixtures.R` and
  `tools/profiling/` unchanged.

## Commits

| Commit | What |
|---|---|
| 09f242bb | tests only (T1, T2, T3, T5, T6; S1; F2 script) |
| c6bbd40d | fix (the rule removed, comments, old test removed, browser lane and helper, F4 rows, S2) |
| eea8775b | F2 script correction (R0b) |
| 367851cc | S2 base side inside `test_that()` (after R4 stopped; R4b ratified, VAL-63) |
| c4b11b97 | re-frozen references, NEWS, ledger, run evidence |
| later | F2 instrument fixes (VAL-54), records |

## Acceptance criteria

| # | Observed | Evidence |
|---|---|---|
| A1 | Met. On the tests-only commit T1 fails at its first expectations for TMT and LFQ (actual `complete`, then `pending` for Normalisation) and the warning expectation with 0 calls; DIA (control) passes; T2 fails with actual `c("TMT", "LFQ")`. On the fix commit both pass (R2: `test-workflow-server-shared.R` 226/0, `test-omic-gui-config.R` 1911/0) | `r1-tests.tsv`, `r1-failures.txt`, `r2/compare.tsv` |
| A2 | Met. T3 passes on the tests-only commit (R1, `R/` = base) and on head (R2). T5 passes on both (R1, R2, R5). F2 did not trigger T5's drop (see F2) | `r1-tests.tsv`, `r2/compare.tsv`, `compare.tsv` |
| A3 | Met. Arms B and B-ruv: 168 of 168 state digests identical between R0 (tests-only commit) and R3 (fix). `qc-effect.tsv` has the three lanes and four arms with every field of 5.3 (provenance values came back empty; the step texts show the defaults: "Min replicates per group: 2", "Min groups required: 2", "Intensity cutoff percentile: 1.0%"). Figures in NEWS and the ledger row | `r0/state-digests.tsv`, `r3/state-digests.tsv`, `qc-effect.tsv`, `r3/step-texts.tsv` |
| A4 | Met, with the owner's ratification of R4b (VAL-63; R4 stopped, see below). F4 rows committed (c6bbd40d) before R4. R4b: step 1, base against the frozen files: on each lane only `integration_artifacts` at `summary_save`, `summary_copy`, `summary_report`, `summary_export` (72 rows, 6 keys per step), 0 outside, equal to the reference rows; step 2, head against base: exactly 15 `qc_routing (KIT-180)` rows per lane (`tab_status` quality_control and normalization, and the stepper output, at the five steps), 0 outside; head record equal to its frozen file; `scripts.tsv` S2 `ok`, nothing withheld; only the three `oracles/` files committed | `control/verdict.tsv`, `gui-entry-control.tsv`, `control/gui-entry-steps.tsv`, `control/scripts.tsv` |
| A5 | Met. R5 head: `test-gui-entry-proteomics.R` 183 passed, 10 failed: the three lanes pass every step; `prot_dia` (4) and `prot_dia_limpa` (6) fail exactly as on base | `compare.tsv`, `r5/compare-failures.txt` |
| A6 | Met. R2: every difference is a new test (T1, T2, T3's LFQ case, T5, T6), the removed T1 block, or the GUI-entry file (+15 failing expectations, the declared routing captures, before the re-freeze). T6 errors on the tests-only commit ("could not find function") and passes on the fix. The browser lane skipped on both (F3) | `r2/compare.tsv`, `r1-failures.txt` |
| A7 | Met. R CMD check base and head: 0 ERROR, 3 WARNING, 3 NOTE, findings only on head 0 | `check-compare.txt` |
| A8 | Met. R5 (after rerunning one lost base shard): 441 same, 5 count-changed (the new tests), 2 fixed (`test-gui-entry-proteomics.R` as A5; `test-workflow-figure-paths.R`, see below), 0 new-failure, 0 message-changed | `compare.tsv` |
| A9 | Met: NEWS entry under the unnumbered heading; `DESCRIPTION`, `packaging/electron/package.json`, lock files, `man/`, `NAMESPACE` unchanged | `git diff a60c4e9a` |
| A10 | Met. T2 passes on head for all four omics; "bypassing QC tab" occurs 12 times in R2's base GUI-entry log and 0 times on head | R2 `logs/`, `logs-base/` |

## Size of the change on the fixtures (observed, R3; one synthetic draw per lane, no inference)

Default path (RUV `automatic`, k = 2 on every lane), before (B-ruv) and after (Q-ruv):

| Lane | Tested | Significant | Planted significant (of 40) | Nulls significant (of 360) | Post-RUV clean-up removed (cut-offs) |
|---|---|---|---|---|---|
| prot_tmt | 398 -> 392 | 20 -> 18 | 17 -> 16 | 3 -> 2 | 2 (60/60) -> 5 (60/0) |
| prot_lfq | 398 -> 394 | 16 -> 18 | 15 -> 15 | 1 -> 3 | 2 (60/60) -> 3 (60/0) |
| prot_lfq_fragpipe | 397 -> 393 | 25 -> 25 | 22 -> 22 | 3 -> 3 | 3 (60/60) -> 3 (60/0) |

RUV skipped, before (B) and after (Q): tested 400 -> 397, 400 -> 397, 400 -> 396; planted significant 18 -> 17,
16 -> 15, 24 -> 21; nulls significant 1 -> 1 on each lane. The QC removed 3, 3 and 4 proteins, all at the intensity
filter; the accession clean-up changed no identifier and merged no row. Largest absolute log fold-change difference
for proteins tested in both arms: 0.027, 0.055, 0.156 (RUV skipped) and 0.033, 0.817, 0.205 (RUV default).
Inferred, not measured: the effect on real data (more missing values, abundance-dependent) is larger.

## Observed and inferred

Observed: everything in the tables above, from the runs named. Observed (audit comparison of R4 and R4b records,
`AUDIT-KIT-180-impl-20261123.md` I01): R4's control stopped (base against frozen: `da` and `summary_refusal` beyond
the reference list; base against head: 6 rows outside the routing captures; per lane 12 rows, at 6 steps, order of
names only, 0 value differences), because its base side ran outside `test_that()` and sorted names in the runner's
locale; R4's head records and frozen files equal R4b's in content. Nothing was committed from R4. R4b, the base side
run inside `test_that()` (C collation) and nothing else changed, was sent without the owner's word (before VAL-54);
the owner ratified it on 2026-11-24 (VAL-63, the case of VAL-55). Inferred: that `test-workflow-figure-paths.R`'s
base failure (a tempdir name containing "259") is chance and unrelated.

## Unrelated discoveries (not changed)

- `test-workflow-figure-paths.R` fails when the random temporary path contains "259" (R5 base); flaky test.
- GUI-made label-free and TMT records replay in order on base, but no replayed state is `identical()` to the first
  run's (F2, R0d); cause not determined. Relevant to KIT-177's F3 and T5.

## Not checked

- The rewritten browser lane: skipped on run-tests (no browser), base and head; to be run by the distributable bundle
  lane at the merger's version bump.
- Why GUI-made records' replayed states are not `identical()` (F2).
- The effect on real data and on designs smaller than five per group.
- Windows (plan section 5's trigger list does not apply).
