# KIT-190 outcome (implementation at ea5c71ad (code) and 6b465f1e (the VAL-280 test edit); final R2 45619021147; implementation audit AUDIT-KIT-190-impl-20261126)

Ticket `tickets/kit/KIT-190.md` r1.6 (C01 to C04 of the closure check applied at step 0). Branch
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

## Acceptance criteria (final states; AUDIT-KIT-190-impl-20261126 I11)

| Criterion | Final state |
|---|---|
| A1 | met: R1a, d-first, bound-first |
| A2 | met |
| A3, A4 | met: final R2, with the four files of section 5 |
| A5 | met |
| A6 | the staging checkpoint |
| A7 | met after I02 to I04 and I12 |
| A8 | met at the branch point d4e2eb2c (the 2.4 and 2.5 re-read) |

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
  `VAL-319` written with the first repair's figures; they change if the rule changes (VAL-110).
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

## Bounded-k stage 1 (counted, blinded; run 45504138202, commit e344c691)

Observed (only the nuisance quantities stage 1 returns, `bound-s1/reestimate.tsv`; n = 50 pairs per cell): discordant
share N-r 0, N-s 0, N-w 0.06; SD of the paired FDP difference E-r 0.0559, E-s 0.1126, E-p 0.2233, E-g 0.1202, E-w 0.0618.
Falsifiers 0. No mean difference, rate or statistic was computed or read. 8.06 core-hours.

`reestimate.R plan`, run once (run 45521457364; spent 11.21): R per cell N-r 100, N-s 100, N-w 250, E-r 200, E-s 200,
E-p 300, E-g 200, E-w 200 (equal to the auditor's recomputation); target per cell 0.97249, achieved joint probability
under equality 0.987; projected 73.6 core-hours in total (ceiling 100). Files: `bound-s1/plan.tsv`, `plan-summary.tsv`.

## Bounded-k stage 2 (counted; 8 runs at 0413a169)

Observed (mechanics only; no endpoint read and `decide.R` not run): all 8 jobs completed; falsifiers 0 in every job;
every job passed the stage-1 check; journey errors as in each `bound-s2/*/provenance.json`. 733.3 job-minutes, 48.89
core-hours; programme total about 61.4 core-hours (ceiling 100). Returned files committed under `bound-s2/` for the
auditor's check before `decide.R decide` runs.

## The decision run, declared before its result (audit "Stage 2 (counted) and the decision run")

One run of `decide.R decide` (660ebf9a) at 2b04b121 on the eight `bound-s2/*/journeys.tsv`; `out_pattern` admits only
`verdicts.tsv`, `reported.tsv`, `decision.tsv`; a re-run only under VAL-54 (mechanics failure with no result).

Decision table (declaration revision 10, sections 5 and 6; intersection–union, each cell one-sided at α = 0.05):

| `decide.R` output | What follows |
|---|---|
| `no_worse` TRUE in all 8 deciding cells | rule (d) with the k bound ships (VAL-182 (3)); the bound's value goes to the owner with the k and FDP-by-k rows (VAL-182 (1)) |
| `no_worse` FALSE in any deciding cell, inconclusive cells included | the first repair with the bound ships; its row-order dependence stays and is stated in NEWS |
| FALSE only in N-w and/or E-w | as the row above; also reported to KIT-188's default ruling with the B+B reference |
| `owner_asked` TRUE | the owner is asked before shipping |
| `tango_bb_agree` FALSE in a null cell | reported to the owner; the Tango verdict stands |
| `k190_decision_incomplete` | no decision; an instrument defect under VAL-54 |
| E1-T1b, T11, secondary intervals | reported; decide nothing |

The descriptive B+B comparison (D1), fixed now and computed after `decide.R` by `bb-compare.py` (committed with this
entry): N-w, over the 250 seeds, the share of seeds with any call under X+B minus the share under B+B (X = F, D); E-w,
over the 200 seeds, the mean over seeds of FDP(X+B) − FDP(B+B); "far worse" = both differences above 0.10 in that cell.
Descriptive only; reported to the owner and KIT-188; decides nothing.

## The decision (run 45574056931, `decide.R` 660ebf9a at 2b04b121; files in `bound-decide/`, unchanged)

`decision.tsv`: all_cells_no_worse TRUE, ships **D+B**, owner_asked FALSE. Decision table row 1 applies: **rule (d) with
the k bound ships** (VAL-182 (3)); the bound's value goes to the owner with the k and FDP-by-k rows (VAL-182 (1)).

| Cell | n | Estimate (D+B − F+B) | Statistic | p | Berger–Boos p | Skewness g | Johnson | No worse |
|---|---|---|---|---|---|---|---|---|
| N-r | 100 | −0.010 (any call) | Tango Z −3.33 | | 0.0051 | | | TRUE |
| N-s | 100 | −0.020 | Z −3.65 | | 0.0050 | | | TRUE |
| N-w | 250 | +0.004 | Z −4.21 | | 0.0050 | | | TRUE |
| E-r | 200 | −0.0159 (FDP) | t_J −13.88 | sign-flip 1.0e-5 | | −4.45 | yes | TRUE |
| E-s | 200 | −0.0255 | t_J −17.40 | 1.0e-5 | | −3.77 | yes | TRUE |
| E-p | 300 | −0.0929 | t_J −16.56 | 1.0e-5 | | −2.12 | yes | TRUE |
| E-g | 200 | −0.0150 | t_J −13.02 | 1.0e-5 | | −2.70 | yes | TRUE |
| E-w | 200 | −0.0004 | t −6.90 | 1.0e-5 | | −0.81 | no | TRUE |

Tango and Berger–Boos agree in all three null cells. Reported (decides nothing): E1-T1b by seed, n = 20, mean −0.040,
sign-flip p 2.0e-5; T11 concentration 1 false call under both bounded arms (k = 2), T11 F172 0 under both (k = 2); D+B
minus F+B (paired t 95 %): k +0.08 to +0.48 in every cell; power +0.03 to +0.20 in the effect cells; fold-change error
−0.03 to −0.17. The declared per-cell confidence bounds are not among `decide.R`'s outputs (only the statistics are).

Descriptive B+B comparison (D1; `bb-compare.py`; `bound-decide/bb-compare.tsv`): N-w any-call share minus B+B: F+B
+0.112, D+B +0.116, **both above 0.10: "far worse"** (reported to the owner and KIT-188 as the measured form of
K190-M11; changes nothing here). E-w mean FDP minus B+B: F+B +0.050, D+B +0.049 (not far worse).

## After the decision: bounds and the metabolomics check (both decide nothing)

**Per-cell one-sided 95 % bounds** (`bound-decide/bounds.tsv`, run 45577176410; computed after the decision;
descriptive; non-simultaneous; decides nothing): any-call difference upper bound N-r +0.026, N-s +0.007, N-w +0.035
(margin 0.10; inverted Tango score test); mean FDP difference upper bound E-r −0.0017, E-s −0.0123, E-p −0.0721,
E-g −0.0023, E-w +0.0116 (margin 0.05; paired t). These equal the audit's figures. The sign-flip p of every effect cell is
reported as p < 3 × 10⁻⁵ (none of 100,000 sign vectors as extreme); in E-r, E-s, E-p and E-g Johnson's t_J decided and
the sign-flip p is reported beside it; the Berger–Boos p is ≤ 0.0051.

**K190-D01, the declared metabolomics identity check, run late** (run 45577183150; `effect.R` 20c7f9ce, phase
`metab-check`; the counted stages used 7d76467a unchanged): 100 metabolomics rows, falsifiers 0; **100 of 100 identical
to the counted lipidomics rows in all 23 compared fields**. The result may therefore be stated for lipidomics and
metabolomics, within the audit's other limits (12 samples, cap 2 after the sweep, simulated data).

## The bound on k (VAL-182, VAL-235; ticket section 18)

Merged staging e0747fc1 at e6cc8496: the one code conflict, `workflowReplayRuleMarkers`, resolved by keeping all rows
(missing_value_comparison, ruv_control_grouping, enrichment_background, intensity_filter_samples), each with its comment,
no other line changed (coordinator's word). Tests first (ef9cda11, run 45602075177): T9 (i) to (iii) failed at the
unbounded k, the missing notice, warning, record fields and marker; T1 to T8 passed. Code 00434b08: `omicsRuvKBound`
(the sweep untouched; Automatic k = min(k_sweep, `calculateAdaptiveMaxK`) through `boundOmicsRuvSweepK`; Manual kept with a
classed warning `multischolar_ruv_k_notice`), `omicsRuvAssayStatus` (per-assay `corrected` / `uncorrected` / `dropped`),
record outcomes `k_sweep_per_assay`, `k_cap_per_assay`, `ruv_k_notice_per_assay`, `ruv_status_per_assay`, the marker
`ruv_k_bound` (also dropped from managed new-design steps). Notice texts are drafts in `wording.md` (VAL-41).

Scope of every figure here: simulated lipidomics data without internal-standard normalisation (requested, not applied;
VAL-253, K242-P01).

Inferred (for the implementation auditor): no cell of the rule-(d) evidence changes. The measured path (the untouched
sweep, its percentage and controls, k = min(k_sweep, 2), RUV-III-C, limma) is the shipped route at n = 12, where
`calculateAdaptiveMaxK` gives 2. Not covered: caps 3 to 5 at larger n (KIT-188).

## Merge of staging c4505c7f (1cae0c3e, d3318820): the two resolved code hunks (coordinator's word)

`R/utils_workflow_omics_declarations.R`, `omicsRuvOperation`: staging's comment line kept; `fixed` becomes
`list(ruv_feature_levels = "restored", ruv_control_grouping = "kept", ruv_k_bound = "sample_size")` (staging's marker
first), with one added comment line; KIT-190's return fields (`k_sweep`, `k_cap`, `k_notice`, `ruv_status`) and its two
functions (`omicsRuvKBound`, `omicsRuvAssayStatus`) kept; no other line of staging's version of the function changed.

`R/utils_workflow_analysis_record.R`, `workflowReplayRuleMarkers`: staging's rows in staging's order
(missing_value_comparison, enrichment_background, intensity_filter_samples, duplicate_resolution_samples,
ruv_feature_levels), then `ruv_control_grouping`, then `ruv_k_bound`, each with its own comment; nothing else changed.
NEWS and the ledger merged as text; the ledger row `VAL-319` follows VAL-261 and cites VAL-235.

Stated limit. The counted evidence ran on the RUV-III-C correction before KIT-224 (features shifted by a constant).
The sweep's k, percentage and controls do not read the correction (by reading), and in all 3,634 bounded journeys
k = min(k_sweep, 2) at the sweep's percentage, so the bounded arms are the shipped computation at n = 12. The FDP and
any-call figures do depend on the correction through limma's mean-abundance trend (eBayes trend and robust); KIT-224's S1
(unbounded k) measured a mean paired change of at most 0.0007 (per object up to 0.044), common to both arms and small
against the smallest margin left (0.065, N-w). Implementation audit AUDIT-KIT-190-impl-20261126, I04.

**R3: departure by the coordinator.** The next full-suite checkpoint on staging stands in for R3 after this merge (as for
KIT-201, KIT-204, KIT-243).

## R2, final comparison on the merged tree (run 45619021147, head d3318820, base staging c4505c7f)

Observed, every file as predicted before sending: 98 `same`, including KIT-202's, KIT-205's, KIT-214's and KIT-224's
files and every reader of the replay markers; the new file `only-on-head`, 126 passed; `test-kit172-downstream.R` ok
(with k capped at 2); two `new-failure`: `test-omic-norm-module.R` (only the stored provenance digests of the lipid and
metabolite manual and automatic blocks, from the markers; the D5 control shows 12 dump lines differ, all provenance
digests) and `test-workflow-omics-ruv-cancor-once.R` (the sweep, which the bound leaves untouched, under rule (d): its
success assertion was a property of the row-order defect; audit K190-M12); two `message-changed`:
`test-prot-norm-core.R` (lipid and metabolite RUV pinned digests, failing on base too) and
`test-workflow-omics-binding-steps.R` (its stored snapshot gains the markers and outcomes). R CMD check: no finding only
on head. 41.06 core-hours.

Open at that point for the coordinator and the owner: the D5 re-dump of the 12 provenance fields (exception by file,
D5 (a)); the two stored snapshots (`test-workflow-omics-binding-steps.R`, `test-prot-norm-core.R`) and
`test-workflow-omics-ruv-cancor-once.R`'s success assertion. The re-dump is now applied (I07 below); the other three are
not changed here.

## Implementation audit repairs (AUDIT-KIT-190-impl-20261126; 2026-11-26)

**The correction pass (coordinator's ruling).** Before this pass the ticket's correction pass had not been used:
(1) no correction was ever applied to a result of this ticket: the two results of R2 45329435006 (T7 in
`test-kit172-downstream.R`; the cancor-once test) were put to the owner, not repaired; (2) the bound is scope the owner
added by decision (VAL-182, VAL-185, VAL-235), with its own tests-first run and comparison run; (3) the runs that
`runs.tsv` labelled "the correction allowance (252.7)" are R2 repeats after the owner's scope decision (VAL-235). They
are relabelled with a dated note that the old label was wrong. I07's re-dump is an edit the owner authorised by file
(VAL-88 (5)), and the staging merge is mechanics; neither is a correction. **This pass (I01 and the record repairs) is
the ticket's one correction pass; after it none is left.**

**Merge of staging 828dd105 (f5ea9cc2; marker resolution as authorised, N18).** The one code conflict,
`R/utils_workflow_analysis_record.R`. In the header, the KIT-190 clause now says "KIT-190 adds two rows,
`ruv_control_grouping` and `ruv_k_bound`", followed by staging's clauses and KIT-201's `treat_moderation` sentence. In
`workflowReplayRuleMarkers`, staging's `treat_moderation` row ends with a comma (the list's closing parenthesis moved),
then `ruv_control_grouping` and `ruv_k_bound`, each with its comment, close the list. Staging's two functions after
the list (`workflowTreatModerationCurrent`, `workflowTreatModerationRecord`) are kept unchanged. The ledger was merged as
text; `VAL-319` follows VAL-272. No other line changed.

**I07, the D5 re-dump (VAL-88 (5); the D5 control of run 45619021147, 12 of 2,078 dump lines differ).** Commit
09c08b5e changes only these 12 provenance values in `tests/testthat/test-omic-norm-module.R`, each replaced by the head
column of `r2-final/norm-module-control.tsv`:

| Line(s) | Old value | New value |
|---|---|---|
| 1043 | 1ab8a05076a605820dc0d14e237215f9abfaa251d04c7076ca03825877e57dc2 | 70d15badd8b75ecacecaf1a1828f9f98be516a19f6630d390d8e34a85990d090 |
| 1095 | 0b68624dba919ae6e6852528493ba338d672d75d235b8437c2895e8c6f724bf5 | a9b7326e2a270b95a32381ae27ffff05d2ce7e032116b273a7cd983bb72c550b |
| 1132 | 6f03ec2fff1c90a3a1728422205f6ebe8d929a0e97a02b42dc0e053b026e5942 | 56e0f5614543ba4ed0027756f65aeb31be2cf65ef3364290cd24b6c495998dec |
| 1182 | 24204497ce33beb525005c4869fd60e6f560a7c1f6eb2293cc928049e1410ba2 | 9f3dc69aebe34ac437e7da9abd85c2bc560af6548b0fb66c698b423b04c0de69 |
| 1234 | 667cc9562bd90c598ef04d93d5fe30d6ad1818b89c02f50d86dcfec344a11371 | 357815f322e924dde66af6e38355a135eda96220161075443076b7fc176a29dc |
| 1271 | 279b01b63aa09b94e9ccaa8e085154e7b281eb0fd48361548782ceead23456e3 | 54190570168222bd4009ff992fe6585e852580faf7c06f8da9a741234b34a890 |
| 1320, 1358, 1375 | 024a3e8b9582da5526ed27f5f48f518ccda26236c0eeb37c4ddd94ad294e955a | ff4392a8e096ed84963eb0a6c5d2ab8b0d9faf12d4916215a5e38a346f435b30 |
| 1447, 1485, 1502 | 6858408bb39aa5c6a558ccd712b51e5ad5ce9aca6ddf7abd28832b96d50a8848 | e30a4c778dfd8d9cdc90cff9a0776ec84a3a92b58fdcc8b32827c572f98b00ad |

**I01, the k removed and the notices in the tab.** Tests first in 5185b830: T9 (iv) in
`tests/testthat/test-ruv-control-anova-grouping.R`. It sets up a capped Automatic value (found 8, cap 2, 12 samples).
The tab's summary must read "Factors removed: 2 (the search found 8; at most 2 with 12 samples)" and not "Best k: 8";
an uncapped result keeps "Best k: …". Each notice is logged and sent to the notification once. Both routes call the
annotation and the notices. The export line must read "k=2 (found 8, at most 2 with 12 samples)". Code ea5c71ad:
- `omicsRuvKBound` also returns the assay's sample count (`k_samples`), and the step's value carries it.
- `omicsRuvAnnotateResults` adds `k_removed`, `k_found`, `k_cap` and `k_samples` to an assay's sweep result only
  where the cap reduced its k. An assay the cap left alone is unchanged, so the stored `normdata` digests of
  `test-omic-norm-module.R` (where the cap does not bind) do not move.
- `omicNormRuvSummary` (the tab) and `omicsRuvExportSummaryLine` (the session export,
  `R/utils_workflow_saved_exports.R`) print the k removed, and the k found and the cap where they differ.
- `omicsRuvReportNotices` writes each notice to the tab's log and its notification (type "warning"): in process
  through the tab's `notify`, managed through `shiny::showNotification`, as `omicsLog2ReportResult` does.
- The notice texts are O3's: Automatic "{assay}: the automatic search found {k} unwanted factors. With {n} samples the
  application removes at most {cap} (the limit proteomics uses), so {cap} were removed."; Manual "{assay}: you chose to
  remove {k} unwanted factors. With {n} samples the application removes at most {cap} in Automatic mode (the limit
  proteomics uses); your choice was kept." The earlier "leaving enough samples to judge…" is gone from the code and is
  withdrawn in `wording.md`. The final wording remains the owner's (VAL-41).

Every other place k is printed or exported, by reading (not executed with a capped run):

| Place | Source | Shows |
|---|---|---|
| The tab's log line "RUV optimization complete. Best k per assay" (`R/utils_workflow_omics_actions.R`, both routes) | the step's `value$best_k` | the k removed (unchanged) |
| Study parameters, RUV-III per-assay section (`R/func_general_filemgmt_study_parameter_writer.R:404`, `:601`) | the S4 object's `ruv_number_k`, set by the RUV method to the k it applied (`R/func_lipid_s4_ruv_methods.R:396`; metabolomics the same) | the k removed |
| Study parameters, "Automatic RUV Optimization Results" (`buildWorkflowRuvOptimizationLines`, `R/func_general_filemgmt_study_parameter_support.R:256`) | the top level of the result it is given; for lipidomics and metabolomics that is the per-assay list, whose top level has no `best_k` | "N/A" as before; neither the found nor the removed k (unchanged; not repaired here) |
| Lipidomics and metabolomics reports (`inst/reports/lipidomics/lipidomics_report.rmd:379-460`, `:589`; metabolomics `:367` on) | the study parameters' per-assay k and the S4 `ruv_number_k` | the k removed |
| Analysis record and its viewer (`omicsRouteRecordOutcomes`; `R/utils_workflow_analysis_record.R:159`, `:184`, `:240`, `:913`) | the step's saved value: `k_per_assay`, `ruv_k` = `best_k`; also `k_sweep_per_assay`, `k_cap_per_assay`, `ruv_k_notice_per_assay` | the k removed; the found k, cap and notice are in the record's outcomes |
| Integration (`R/utils_integration_gui.R:1185-1199`) | the record's `ruv_k` or per-view `k` | the k removed |

**I05, stated limit.** The cap does not look at the number of groups. After RUV-III-C the test counts n − g residual
degrees of freedom while k of them have been removed. With k at the cap: 2 groups of 2 keep 0 of 2; 3 groups of 2 keep
1 of 3; 2 groups of 3 keep 2 of 4 (standard errors understated about 1.4-fold); 2 groups of 4 keep 4 of 6; the measured
design, 2 groups of 6, keeps 8 of 10. Nothing in the legacy route stops k reaching n − g. This is not measured and is
stated in NEWS and the ledger row. The owner decided O2 (a), VAL-289: the limit stated, the cap's rule kept; the
degrees-of-freedom question goes to KIT-188 and KIT-262's table. O3 (the notice wording) stays a draft.

**I06 / O1.** Decided by the owner (VAL-280, option (a)); the test edit is 6b465f1e; run 45675334261: base 14/2, head 16/0.

**I08.** Draft ledger row `VAL-320` (the KIT-094 snapshot gains two markers and four outcome fields; only the
message of `test-workflow-omics-binding-steps.R` changes; replaced at the end-of-project re-record) for the coordinator.

**I09.** The two RUV result digests of `tests/testthat/test-prot-norm-core.R` (lipid_ruv pinned 0921c3d1…, base
cc4276ca…, head 1cfd7586…; metab_ruv pinned dea4b047…, base 9ada5791…, head f7593cff…; run 45619021147) are already on
the end-of-project re-record list under VAL-213. They move by design (the control rule; the bound is inert at k = 1).
They are not edited.

**I10.** `tests/testthat/test-workflow-omics-declarations.R`, `test-gui-entry-lipidomics.R` and
`test-gui-entry-metabolomics.R` are added to R2's scope, with predictions declared in `r2-audit/predictions.md` before
the run.

## For the change report

- **What was wrong.** Lipidomics and metabolomics RUV-III-C (Automatic, the default, and Manual) did not choose its
  negative-control features by comparing the groups, as documented. The control method removed the design's group
  column before calling the shared ANOVA helper, which then read that column as NULL without an error: no feature was
  tested, every p-value became 1, and the controls were the first rows of the assay table, in table order. The
  automatic number of factors had no upper bound (proteomics has had one since KIT-123).
- **Who was affected.** Every lipidomics or metabolomics analysis that ran RUV-III-C in Automatic or Manual mode with
  the RUV grouping equal to the design's group column (the default), on `main` and on staging until this change.
  - On `main` 436cdbe8, the methods drop the group column (`R/func_lipid_s4_objects.R:1999`, `:2005`;
    `R/func_metab_s4_objects.R:1999`, `:2005`).
  - The helper reads the grouping (`R/func_prot_norm.R:310`), sets missing p-values to 1 (`:322`) and takes the first
    names of a stable sort (`:344`).
  - Automatic, the tab's default (`R/mod_lipid_norm.R:149`, `R/mod_metab_norm.R:149`), selects once at the largest
    percentage, 20 % by default (`R/func_lipid_norm.R:415`, `:456`; metabolomics the same), so its controls were the
    first 20 % of rows. It takes k as the K at the largest canonical-correlation gap, with no upper bound (`:491`;
    `R/func_prot_norm.R:363-374`).
  - Proteomics and peptide kept their grouping and were not affected.
- **Code.** Base c4505c7f.
  - `R/func_lipid_s4_ruv_control_methods.R`, `R/func_metab_s4_design_control_methods.R`: grouping kept; refusals.
  - `R/func_prot_norm_optimization_helpers.R`: `candidate_rule = "p_id_two_groups"` (q > 0.05, observed in at least two
    groups, largest p first, ties by feature name); guard; pool class.
  - `R/utils_workflow_omics_declarations.R`: `omicsRuvKBound` (Automatic k = min(k found,
    `calculateAdaptiveMaxK`(samples))); records.
  - Replay markers `ruv_control_grouping`, `ruv_k_bound`.
- **Evidence.**
  - T1 (run 45278118848): on the base the controls were F01 to F04, every p = 1, and the grouping read NULL.
  - Unbounded, rule (d) raised k to 5.2 to 6.0, and null data sets with a call to 32 to 34 of 50 (VAL-182).
  - With k capped at 2, rule (d) was no worse than the capped first repair in all 8 pre-declared cells (VAL-235;
    decide run 45574056931). Upper bounds: N-r +0.026, N-s +0.007, N-w +0.035 against 0.10; E cells −0.072 to +0.012
    against 0.05.
  - Scope: simulated lipidomics data, 2 groups of 6, without internal-standard normalisation (VAL-253); metabolomics
    identical in 100 of 100 rows.
- **Limits.** Caps of 3 to 5, and designs below 12 samples, are not measured. The measured correction was before
  KIT-224's level restoration (I04 of the audit). Where the unwanted variation lines up with the groups, both
  group-comparing rules gave false calls in about 20 % of null data sets, against 8 % for group-blind controls
  (KIT-188).

## The correction pass's runs (observed)

**I01 failing-first** (run 45666131951, 5185b830): 126 passed, 7 failed, all in T9 (iv) (no annotation, no notices, no
route call, no export line); T1 to T9 (iii) pass. The first dispatch failed with GitHub HTTP 502/500 during an outage
(no run started); the same request was then dispatched with the dispatcher's own command.

**R2 after the repairs** (run 45651851022, head 687ba6f3, base staging 828dd105; 106 files, check, D5): every file as
predicted in `r2-audit/predictions.md`. 101 `same`, among them `test-omic-norm-module.R` (639 passed, after I07), the
I01 closure files, both GUI-entry files and `test-kit172-downstream.R`. The new file is `only-on-head`, 139 passed.
`test-workflow-omics-ruv-cancor-once.R` is a `new-failure` at `:73` (this head did not hold the VAL-280 edit).
Three files are `message-changed` and fail on both sides: `test-prot-norm-core.R` (VAL-213),
`test-workflow-omics-binding-steps.R` (`VAL-320`) and `test-workflow-omics-declarations.R` (base 45/23, head
43/25). Check: 6 findings identical, 0 only on head. D5: 12 of 2,078 dump lines differ, each head value equal to the
stored I07 value. 44.83 core-hours. Evidence: `r2-audit/`.

**O1, VAL-280** (run 45675334261, head 529b2cc0, base 687ba6f3; the cancor-once file alone): base 14/2, head 16/0,
`fixed`, as declared. The observation decides nothing; it is for KIT-262 M1 (`o1/observation.txt`). Assay `two` was
rejected at every percentage as `skipped_k_exceeds_complete_controls`, with k 2, 2, 2, 3, 5 at 2 to 6 % against 0, 1,
1, 1, 2 controls complete in every sample. Lipidomics and metabolomics gave identical lines. Assay `one` was accepted at
2 to 4 % (k 2) and rejected at 5 and 6 % (k 8 > 6, k 9 > 7). At 2 to 4 % assay `two` is rejected at k = 2, so a cap
applied inside the search would not have rescued it (O1 (b)). 1.57 core-hours.

Correction-pass budget: 220.39 before the pass, then + 0.92 + 44.83 + 1.57 = 267.71 of 510.4 core-hours. No correction
pass remains.

## Closure check (2026-11-26): READY FOR A MERGE TRAIN

The audit file's closure check found no remaining code finding. The records below follow its texts.

- **CC1.** NEWS now states VAL-280's consequence. An assay can end with no usable percentage, and is then left out of
  the RUV-corrected data and recorded as 'dropped'. KIT-242 Part B will reword this when it keeps such an assay
  uncorrected with a notice.
- **CC3.** The title, the ticket revision (r1.6) and the I06 / O1 paragraph are brought up to date.
- **CC4 (note).** T9 (iv) checks the two routes only by searching their function bodies for the helper names. No test
  drives a capped run through either route. The auditor accepted this because the delivery code is two lines per route.
- **CC5 (note).** The NEWS limit sentence names "VAL-289" in user text. That is for the owner's wording pass and is not
  changed here.
- **Hand-off to KIT-242 Part B (not KIT-190's).** An assay whose search succeeded but whose correction then failed
  (control-vector alignment, or no finite value) still has its k shown as removed in the summary and the log.
  `ruv_status_per_assay` records "dropped" or "uncorrected", but nothing displays it. This predates KIT-190, KIT-190
  does not make it worse, and Part B owns it (VAL-242, VAL-290).

**`test-workflow-omics-declarations.R`, observed part by part by the auditor (run 45651851022; base 45/23, head 43/25).**
- **Two new failures:** `values` (the RUV step's value) for lipidomics and metabolomics. Each gains `k_sweep`, `k_cap`,
  `k_samples`, `k_notice` and `ruv_status`; its `best_k`, `control_count` and `percentage` are unchanged (legacy Manual,
  k = 1, 20 %).
- **8 of the 23 failures already present on base change their message in substance,** for both omics: the RUV-corrected
  and correlation-filtered states, the sample correlations, and the DA results. The reason is that rule (d) picks other
  controls at 20 %.
- **Unchanged:** every state before the RUV step gives the same message on both sides. The registry's message is
  unchanged; its `fixed` list was already a list on base and gains the two markers.
- **Ledger:** draft row `VAL-321`, in the VAL-181 form, carries the register's text.

**The train's scope and predictions (from the register; base = staging).**

Scope:
- the 106 files of R2 45651851022;
- staging's norm- and report-related test files added or changed since 828dd105: `test-kit257-norm-own-files.R`,
  `test-workflow-saved-exports.R`, `test-kit236-report-sentences.R`, `test-kit235-lipid-metab-report-parser.R` and
  `test-kit227-report-ruv-sentences.R`;
- if stacked on trains 13 and 14, their files as well: `test-kit254-omics-correlation-skip.R`,
  `test-kit248-correlation-skip-record.R`, `test-omics-norm-state-guards.R` and `test-workflow-operation-registry.R`.

| File | Predicted |
|---|---|
| `test-ruv-control-anova-grouping.R` | only-on-head, 139/0 |
| `test-workflow-omics-ruv-cancor-once.R` | `same` 16/0 |
| `test-omic-norm-module.R` | `same` 639 |
| `test-prot-norm-core.R` | message-changed, 1548/306 on both sides |
| `test-workflow-omics-binding-steps.R` | message-changed, 58/1 on both sides |
| `test-workflow-omics-declarations.R` | 45/23 → 43/25, the parts as above |
| both GUI-entry files | `same` |
| every other file | `same` |
| R CMD check | identical |
