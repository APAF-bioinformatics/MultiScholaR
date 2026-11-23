# KIT-190 outcome (in progress; stopped before T7, the replay row and E)

Ticket `tickets/kit/KIT-190.md` r1.2 (C01 to C04 of the closure check applied at step 0). Branch
`kit-190-ruv-control-selection` from `v0.9-staging` d4e2eb2c (the branch point). Run table: `runs.tsv`.

## Commits

| Commit | Content |
|---|---|
| 06c08bce | tests first: T1 to T6 and T7 (ii) (`tests/testthat/test-ruv-control-anova-grouping.R`); nothing under `R/` |
| 5c7d984d | `effect.R` (E1, E2, E3; lifecycle and fault modes); nothing under `R/` |
| 75a55cff | tests: instrument repairs after R1a (self-check called through the namespace binding; T6 protein table without row names) |
| 4ec3e7c8 | tests: the pool refusal's message matched across white space (cli wraps it); not reached on base |
| 6b6d733c | fix: 7.1, 7.2, 7.3, 7.5 (control methods keep the grouping column; helper guard and pool class; methods re-signal the pool refusal; sweeps pass the refusal through; panel legacy note and source; roxygen and both Rd pages hand-edited) |

Not yet written: 7.4 (the `fixed` marker on the legacy declaration, the row in KIT-187's `workflowReplayRuleMarkers`,
the managed new-design removal at `R/utils_workflow_execution.R:3327-3332`), T7 (i), (iii) to (vi), NEWS, the ledger
row, the D5 re-dump. They wait for KIT-187 to merge into `v0.9-staging` (task instruction).

## Observed

- **T1, the executed confirmation of the defect (R1a, run 45278118848; repeated in R1a-2, run 45279820283), at the
  tests-only commit.** Instrument checks held on lipid and metab: the pass-through double recorded exactly one call;
  the method returned the assay as a logical vector of 40 named values; the trace recorded 40 p-values. Base values:
  selected = F01, F02, F03, F04 (the first round(40 x 10 / 100) rows); every traced p = 1; the value
  `design_matrix[colnames(data_matrix), "group"]` read in the helper was NULL. This is 2.1's reading, executed. The
  ticket is not void.
- Base, same runs: T2 (grouping `factor1`) passes; T3 (i) to (iv) give no error; T4 and T5 select F01 to F04; T7 (ii)
  passes. T6 passes on base (R1a-2): protein and peptide selections equal the oracle, none of F01 to F10, traced p equal
  to the oracle's (A2, base side).
- Instrument defects in R1a (VAL-54; recorded, repaired, rerun once as R1a-2; no expected value changed): the
  self-check called the helper from the test environment, which found the attached package copy, not the namespace
  binding the trace modifies (0 traced p); T6 built its protein table with row names, which `column_to_rownames`
  refuses. T1's own trace check held in R1a, so T1's p assertions stood.
- **Fix commit on head (run 45282639086):** the new file 66 passed, 0 failed; `test-lipid-12c-ruv-s4-characterization.R`
  59, `test-metab-02y-...R` 26, `test-workflow-lipid-declarations.R` 259, `test-workflow-metab-declarations.R` 259,
  `test-workflow-omics-ruv-complete-controls.R` 28, all passing with the same counts as on base. The journey values of
  6.2 did not move (predicted).
- **E lifecycle (run 45279512656, tests-only tree as both sides), counts:**
  - F1 = 0 on both sides in every cell; F3 = 0; oracle check (N-r with grouping `factor1`) 0 of 4: the oracle agrees
    with the helper's ANOVA where it runs.
  - F2 = every head journey in E-r, E-s, E-p, E-g and E1 (as stated), **but 0 of 2 per class in N-r and N-s**, where the
    ticket states every head journey.
  - Fault run: stopped with `k190_falsifier` after writing `falsifiers.tsv` (as required). The dropped control also
    reached the lipid oracle-check journey (count 1), because the fault keys on class and replicate only.
  - E3 at 520 features x 2 assays: 7.0 s RUV step, 19.9 s journey; `times.tsv` written. On the tests-only tree no
    ANOVA runs, so this measures the mechanics only.
  - Mean journey 16.0 s; projected 86 min per side at R = 50 with two workers per side (within 150).

## Inferred

- Why F2 counts 0 in the null cells: in N-r and N-s every journey's head-side `share_tie_decided` is 1 and the oracle's
  selection equals the first rows (`e-lc/journeys.tsv`; the `factor1` journeys, where the real ANOVA runs, also select
  exactly the first rows). Under the null the BH q-values form one top plateau that holds every chosen control, and the
  stable ranking takes it in table order: the designed selection and the row-order selection coincide. This is 2.6's
  tie dependence in its strongest form, measured, not an instrument fault. It means (i) the lifecycle's stated F2
  count is wrong for the null cells, and (ii) in E2 the N cells will show no change from the repair by construction.
- The lifecycle's pass rule ("exactly the counts above") is therefore not met. Changing the stated count is an
  expected-value change, which VAL-54 does not allow as an instrument repair; it needs the coordinator's or owner's
  ruling (see the report).

## Acceptance criteria (state at this stop)

| # | State | Evidence |
|---|---|---|
| A1 | T1 confirmed on base with instrument checks; T2 to T6, T7 (ii) as stated on base; all pass on head. T7 (i), (iii) to (vi) not yet written | R1a, R1a-2, fix run |
| A2 | base side met (R1a-2); head side met for the targeted run | R1a-2; fix run |
| A3, A4, A6 | not run (R2, R3 after KIT-187 merges) | |
| A5 | lifecycle not passed as stated (F2 in N cells); E, E3 not run | R1-lc |
| A7 | NEWS, ledger: not written; `DESCRIPTION`, `NAMESPACE`, versions unchanged | diff |
| A8 | not written | |

## Not checked

- T7 and 7.4 (need KIT-187's table); R2, R3, E, E3 at 5,000 and 20,000 features.
- The pool refusal's text as the tab and the worker show it (only its class and words, by test).
- The E3 cost of the ANOVA (the lifecycle ran on a tree without it).
