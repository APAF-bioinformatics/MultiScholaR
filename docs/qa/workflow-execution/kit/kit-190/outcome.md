# KIT-190 outcome (in progress; E and E3 done; T7, the replay row, R2 tests and R3 wait for KIT-187)

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

## After VAL-100 (owner's ruling on the lifecycle)

- Ticket text amended (revision 1.2, amendment row citing VAL-100): F2's positive control applies to the planted cells
  and E1; N-r and N-s are reported as "designed selection equals row order". `effect.R` (798066fc, d6d67efe) states
  the counts in a `stated` column, keys the fault on the E2 journey, and records what D2's candidate rules would select
  (reported only).
- Lifecycle re-run (run 45291710599): the fault run stopped as required (oracle check 0); the main run failed in the
  controller's new stated-count code (instrument defect, no table). Repaired (d6d67efe) and re-run under VAL-54 (run
  45294460862): **every count as stated** (`falsifier_counts_as_stated: true`).
- Commits: 798066fc (effect.R after VAL-100), d6d67efe (stated-count repair). Run branch `run-k190-lc2` holds the
  tests-only tree 4ec3e7c8 with each effect.R, for the lifecycle runs.

## E and E3 (observed; base d4e2eb2c against head d6d67efe)

Runs 45296793520 (E1, N-r, N-s, E-r), 45296797285 (E-s, E-p, E-g), 45296800839 (E3). No falsifier triggered. Tables:
`e/a/`, `e/b/` (`cells.tsv`, `journeys.tsv`, `falsifiers.tsv`, `fixtures.tsv`, `provenance.json`), `e3/times.tsv`.
`e/d2-summary.tsv` is a reported-only extract of `cells.tsv` for D2 (made from the committed rows, no new run).

**Lipid and metab give identical figures in every cell**: the same generated values go through the two bindings, and
the two code paths compute the same numbers. They are one measurement, not two.

E2, means over R = 50 (95 % Monte Carlo interval in the tables), base -> head:

| Cell | controls planted | head: controls decided by ties | head: one-group controls (features) | k | FDP | power |
|---|---|---|---|---|---|---|
| N-r | 0 -> 0 | 0.933 | 0 | 2.00 -> 2.26 | 0.160 -> 0.200 | - |
| N-s | 0 -> 0 | 0.955 | 0 | 2.22 -> 2.16 | 0.240 -> 0.220 | - |
| E-r | 0.103 -> 0.032 | 0.867 | 0 | 1.72 -> 2.48 | 0.156 -> 0.180 | 0.502 -> 0.561 |
| E-s | 0.100 -> 0.019 | 0.869 | 0 | 1.46 -> 2.46 | 0.129 -> 0.188 | 0.565 -> 0.620 |
| E-p | 0.868 -> 0.422 | 0.847 | 0 | 1.28 -> 2.34 | 0.280 -> 0.276 | 0.481 -> 0.552 |
| E-g | 0.086 -> 0.024 | 1.000 | 0.052 (0.050 of features; all 26 survive the filter) | 2.12 -> 2.48 | 0.194 -> 0.200 | 0.530 -> 0.523 |

Under the null, 93 % (N-r) and 96 % (N-s) of head's controls are still the first rows. The FDP and power differences
have intervals that include 0 in every cell (largest: E-s FDP +0.059 [-0.007, 0.126]); k rises by 0.3 to 1.1.

Candidate rules for D2, head, at the percentage the shipped sweep chose (selection only; k, FDP and power under them
not measured):

| Cell | (b) by p, id tie-break: overlap with shipped / planted / one-group / complete | (c) no one-group candidates: overlap / planted | (d) both: overlap / planted / one-group / complete |
|---|---|---|---|
| N-r | 0.092 / 0 / 0 / 0.371 | 1.000 / 0 | as (b) |
| E-r | 0.146 / 0.006 / 0 / 0.358 | 1.000 / 0.032 | as (b) |
| E-s | 0.140 / 0.007 / 0 / 0.349 | 1.000 / 0.019 | as (b) |
| E-p | 0.176 / 0.007 / 0 / 0.372 | 1.000 / 0.422 | as (b) |
| E-g | 0.048 / 0.002 / **0.724** / 0.088 | 0.904 / 0.025 | 0.047 / 0.005 / 0 / 0.353 |

Shipped rule, for comparison: planted 0.032 (E-r), 0.019 (E-s), 0.422 (E-p), 0.024 (E-g); complete 0.30 to 0.31.

E1 (KIT-172 fixtures): T1b, 42 journeys: k 1.95 -> 2.86, calls 60.1 -> 61.4, FDP 0.226 -> 0.219, power 0.621 ->
0.630, planted controls 0.118 -> 0.031, control overlap (Jaccard) 0.505. T11 concentration fixture (no planted
effect): base 1 %, 5 controls, k 2, 16 calls; head 14 %, 72 controls, k 7, **209 calls** (all false by construction).
T11 F172 (no planted effect): base 37 calls; head 0. The largest |change in log2 FC| per fixture is large (median 3.5
over T1b, up to 1,066; 31 and 195 in T11), which points at unstable fits at some k; not examined further here.

E3: head Automatic RUV step 201.7 s at 5,000 features x 2 assays and 791.1 s at 20,000 x 2 (base 22.2 s there);
Manual 22.7 s and 86.8 s. Projected Automatic step at 50,000 features: 989 s, under the 1,800 s worker deadline.

## Inferred from E

- The repair changes the selection (planted share of controls falls 3 to 5 fold in E-r, E-s, E-g) but, because the
  q-value ties are wide, 85 to 100 % of head's controls are still chosen by row order inside the top tie block, and
  under the null the selection hardly moves.
- With planted features placed first (E-p), the repaired rule still takes 42 % of its controls from them: the top tie
  block includes truly changed features, taken in table order. Rules (b) and (d) cut that to 0.7 %.
- One-group features were not over-represented after the repair in E-g (5.2 % of controls against 5.0 % of features):
  the table's order inside the all-1 tie block decided. Rule (b) alone would make them 72 % of the controls; (c) and
  (d) exclude them.
- On the T11 concentration fixture the repaired route chose 14 % and k = 7 and called 209 features with no planted
  effect, against 16 on base. One fixture, one draw; a result to put before the owner with KIT-188.

## Acceptance criteria (state at this stop)

| # | State | Evidence |
|---|---|---|
| A1 | T1 confirmed on base with instrument checks; T2 to T6, T7 (ii) as stated on base; all pass on head. T7 (i), (iii) to (vi) not yet written | R1a, R1a-2, fix run |
| A2 | base side met (R1a-2); head side met for the targeted run | R1a-2; fix run |
| A3, A4, A6 | not run (R2, R3 after KIT-187 merges) | |
| A5 | met: lifecycle as stated after VAL-100 (R1-lc3); E: F1 to F3 not triggered; `cells.tsv`, `fixtures.tsv` hold the quantities of 5.3; E3's `times.tsv` with the projection (989 s < 1,800 s) | R1-lc3, E-a, E-b, E3 |
| A7 | NEWS, ledger: not written; `DESCRIPTION`, `NAMESPACE`, versions unchanged | diff |
| A8 | not written | |

## Not checked

- T7 and 7.4 (need KIT-187's table); R2's tests and check, R3.
- k, FDP and power under candidate rules (b) to (d) (only their selections were recorded).
- Why the T11 concentration fixture's head journey reached k = 7 and 209 calls.
- The pool refusal's text as the tab and the worker show it (only its class and words, by test).

## Observed, not examined (named to the owner and to KIT-188)

- T11 concentration fixture (`kit172ConcentrationValues()`, no planted effect), legacy Automatic, lipid and metab
  identical: head (first repair, d6d67efe) chose 14 % (72 controls) and k = 7 and made 209 calls, all false by
  construction, against 1 %, 5 controls, k = 2 and 16 calls on base (`e/a/fixtures.tsv`, rows `t11_concentration`).
- KIT-172 T1b (42 journeys): the largest |change in log2 FC| per journey, head against base, has median 3.5 and maximum
  1,066 (`e/a/fixtures.tsv`, column `max_abs_dlogfc`). Neither has been examined; both point at unstable fits at some k.

## After KIT-187's merge (6416dbd6)

- Merged (true merge); no conflict. 7.4 written (1d9c6847): the row `ruv_control_grouping` in `workflowReplayRuleMarkers`,
  `fixed = list(ruv_control_grouping = "kept")` on the legacy RUV declaration, and no marker on managed new-design steps
  (`completeWorkflowAction`). Ticket citations corrected (omics branch `:878`, wrapper `:785`). NEWS and the ledger row
  `VAL-<merge>` written with the first repair's figures; they change if the rule changes (VAL-110).
- T7 failing first (tests-only 2caff088, run 45307510617): (i) at the marker row; (iii) and (vi) at the marker absent;
  (iv) at (a)'s mismatches and (b)'s marker; (v) passed. (vi)'s new-design half failed in the instrument (the fixture's
  default object has duplicate ids, which the declared design refuses); repaired in 518f29b1. KIT-187's T8 file ok.

## Rule (d) (VAL-110; lipid and metab, Automatic and Manual, VAL-126)

- Ticket r1.4 (2.7a, D8 answered). Tests first (2c843db4, run 45325824390): T5 (lipid, metab) and T8 (i) to (iii)
  failed at the first repair's selections; T8 (iv) (proteomics, peptide and the helper's default unchanged) passed.
  Code a1c12d28 (helper argument `candidate_rule`, default `"q_table_order"`; the lipid and metab methods pass
  `"p_id_two_groups"`, or the default under the earlier rule); run 45327321011: all six targeted files green (new file
  110 passed; journey files 259 and 259; characterisation files 59 and 26; complete-controls 28).
- **Declared before the measurement runs** (k190-ed-a/b-a1c12d28): E (E1 and E2, R = 50, the seeds of runs 45296793520
  and 45296797285), base d4e2eb2c against head a1c12d28, `K190_RULE=d` (F2 checks the head against rule (d)). Compared:
  base (row order), the first repair (d6d67efe, earlier runs, same seeds and base) and rule (d), on k, FDP, power and
  the control shares, per cell, and on E1's fixtures (T1b's 42 journeys, T11 concentration, T11 F172). The figures are
  reported; the ticket states no threshold, so no verdict is drawn from them. The head's `share_tie_decided` is computed
  on q ties and does not describe rule (d).

## Rule (d) measured (observed; runs 45329292905, 45329295530; `ed/`, `ed/three-way.tsv`)

The base side reproduced the first E run's base figures exactly in every cell, so base, the first repair and rule (d)
are comparable journey for journey. Means over R = 50 (95 % Monte Carlo intervals in `ed/three-way.tsv`):

| Cell | k: base / first / (d) | FDP: base / first / (d) | power: base / first / (d) | calls per data set (d) | planted controls (d) |
|---|---|---|---|---|---|
| N-r | 2.00 / 2.26 / 5.46 | 0.160 / 0.200 / **0.680** | - | 96.5 (34 of 50 with any call) | 0 |
| N-s | 2.22 / 2.16 / 5.22 | 0.240 / 0.220 / **0.640** | - | 87.7 (32 of 50) | 0 |
| E-r | 1.72 / 2.48 / 5.32 | 0.156 / 0.180 / 0.444 | 0.502 / 0.561 / 0.800 | 127.8 | 0.007 |
| E-s | 1.46 / 2.46 / 5.78 | 0.129 / 0.188 / 0.526 | 0.565 / 0.620 / 0.786 | 146.9 | 0.008 |
| E-p | 1.28 / 2.34 / 5.60 | 0.280 / 0.276 / 0.506 | 0.481 / 0.552 / 0.720 | 134.8 | 0.008 |
| E-g | 2.12 / 2.48 / 6.02 | 0.194 / 0.200 / 0.580 | 0.530 / 0.523 / 0.802 | 164.1 | 0.006 (one-group 0) |

Head-minus-base (d): k +3.0 to +4.3, FDP +0.23 to +0.52, power +0.22 to +0.30; every interval excludes 0.

E1 under (d): T1b (42 journeys) k 1.95 -> 5.21, FDP 0.226 -> 0.239, power 0.621 -> 0.516, false calls in all 1,166 ->
1,748; the largest |change in log2 FC| per journey median 5.95 (23 of 42 above 5), maximum 1,066. T11 concentration:
base 16 false calls, k 2; first repair 209, k 7; (d) **218, k 7** (12 %, 62 controls). T11 F172 (no planted effect):
base 37 false calls; first repair 0; (d) **330, k 8** (18 %, 93 controls).

~~Inferred, not examined: rule (d) chooses the features least associated with the groups, which also carry the least
unwanted variation relative to their noise.~~ **Withdrawn** (audit `tickets/kit/audits/AUDIT-KIT-190-rule-d-20261125.md`,
K190-M07): in this generator (d)'s controls carry more of the planted factor, not less (inferred by the audit's model;
(d)'s fold-change errors are smaller than base's, shown). The audit's account (K190-M05, M06; inferred from code, the k
shift shown): (d) chooses controls on the same group factor the k diagnostic uses, so the diagnostic's control curve is
forced down and the elbow returns k of 5 to 10 on 12 samples; the test on corrected values then keeps 10 residual df
less nothing for k, and false calls rise with k (shown: at equal k, (d) is no worse than base). The rule meets its selection goals (planted controls under
1 %, no one-group controls, no row-order dependence) and fails the error-rate one in every cell measured.

## R2 (run 45329435006; head 4baddb58 = rule (d) + 7.4 + staging d2179e3b; base d2179e3b)

Observed: 70 files `same`; the new file `only-on-head`; check 0 findings only on head. Not as predicted:
- `test-kit172-downstream.R` new failure (T7, lipid and metab): "unobserved >= 1" is 0. ~~Attributed to (d)'s one-group
  exclusion.~~ **Corrected** (audit K190-M12, inferred from code and E1's k): on this fixture (d) chose k = 7; the five
  features with no group-A value leave 5 residual dimensions < 7, RUV-III-C returns them all NA and the method removes
  them, so no unobserved feature reaches the test. The exclusion only decides whether they are controls. The first repair
  also chose k = 7 here; T7 on the first repair is checked by run k190-split-check-c68c2975.
- `test-workflow-omics-ruv-cancor-once.R` new failure (lipid, metab, sweep with missing values): not every assay
  succeeds (inferred: a larger k exceeds the complete controls at every percentage).
- `test-workflow-omics-binding-steps.R` message-changed: the record snapshot gains one fixed parameter (the marker),
  expected from 7.4; its stored snapshot was not in 6.2's list.
- `test-prot-norm-core.R` message-changed: the lipid and metab RUV pinned digests (failing on base too).
- `test-omic-norm-module.R` new failure: the provenance digests (D5, predicted). The D5 control script failed reading
  the dumps (instrument defect; not re-run).

## Observation: T7 and cancor-once on the first repair (run 45471764300)

Observed: on 1ea2e784 merged with staging e305fab3, `test-kit172-downstream.R` T7 fails for lipid and metab exactly as
under (d) ("unobserved >= 1" gives 0); the test prints no k (E1's committed rows record k = 7 for this fixture under the
first repair, `e/a/fixtures.tsv`). `test-workflow-omics-ruv-cancor-once.R` passes on the first repair (16 of 16, as on
base). Inferred: T7 is a consequence of the unbounded k under either rule (audit K190-M12); cancor-once is specific to
(d)'s higher k on that fixture.

## Next: the bounded-k measurement (VAL-182)

Declared in `bound-declaration.md`, to be audited before the lifecycle run is sent. Nothing of it has run.

## Bounded-k measurement: lifecycle (run 45495259308, commit abbe4525, declaration revision 8; not counted)

Observed: every item as declared (`runs.tsv`, row `bound-lc`; outputs in `bound-lc/`). Blinding item (ix): both
stage-1 runs (plain and planted, the plant acting on the rows before any summary) uploaded only `reestimate.tsv` and
`provenance.json` under `out`; the two `reestimate.tsv` are byte-identical, the consoles identical apart from the
`seconds:` line, the provenance identical apart from the two seconds fields. One limit of this lifecycle: at R = 2 the
null cells' discordant share was 0, so the planted swap in N-r could not change it whatever the code did; the swap's
invariance is shown by `reestimate.R selftest` (item (viii)). (xi-b): a stage-2 job on N-w and E-w alone passed the
stage-1 check. 3.15 core-hours.

Declaration revision 9 (e344c691, SHA-256 4ca0913e...288a) cites the owner's acceptance of R4-01 (VAL-192). Stage 1
waits for the coordinator's word.
