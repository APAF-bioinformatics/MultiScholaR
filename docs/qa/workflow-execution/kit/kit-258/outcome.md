# KIT-258 outcome (in progress: tests-only commit, R1, fix as far as the owner's rule, S1)

Ticket `tickets/kit/KIT-258.md` r1.2 (Re-check 1 R1-01 to R1-04 applied). Branch `kit-258` from `v0.9-staging` at
**63143aa3** (KIT-214 merged, VAL-256). The owner's questions OD-1 (does an intensity of exactly 0 count as measured)
and OD-2 (the strict ordering) are OPEN: nothing here assumes an answer. E5, E6, the final fix commit, NEWS, the
ledger row, the user-facing texts and R2 wait for them.

Labels: OBSERVED = read from a run's returned tables; PREDICTED = written before the run was read.

## Commits

| Commit | What |
|---|---|
| ddbedc3f | tests only: new `test-kit258-duplicate-completeness.R`; declared edits of `test-kit214-sample-columns.R` (`k214Earlier`, `k214OracleKept`), `test-metab-01m-...` and `test-metab-02u-...` source pins. Nothing under `R/` |

## R1 (tests-only commit ddbedc3f, base none; 1 job, ceiling 45 job-minutes)

Purpose: failing-first. Binds ddbedc3f (R/ = staging 63143aa3). Decides: whether the fix may be written against these
tests (every failing row at its stated base value, every guard passing).

PREDICTED (written before the run was read):
- `test-kit258-duplicate-completeness.R`:
  - E1, E2, E3 (2 omics x 2 variants = 12 tests): each fails on the kept row only, at KIT-214's value: E1 row 1 'A',
    E2 row 2 'A', E3 row 2 'Q'; the one-row-per-id expectation passes.
  - E4 kernels (1 test): fails on id D for both kernels ('A' for 'B'); id E passes. E4 S4 (2 tests): fail on D ('A').
  - T1 (b) (2 tests): `kept` 'A', `replayed_kept` 'A', the rank marker NULL in the record and the replayed record;
    mismatches empty and the samples marker pass.
  - T2 (1 test): the rank marker NULL for lipid and metab; worker status and samples marker pass.
  - G1, G3 (8 tests), G2 (2 tests), T1 (a) (2 tests): pass.
- `test-kit214-sample-columns.R`: passes in full (its G1 objects hold only equally measured duplicates).
- `test-metab-01m-...`: the source-pin test fails on its three new pins (`n_measured`, `keepDuplicateResolutionRow()`,
  the `select`); every other expectation passes.
- `test-metab-02u-...`: the source-pin test fails on its one new pin (`keepDuplicateResolutionRow()`); every other
  expectation passes.

OBSERVED (run 45617808183, 14.5 job-minutes; `r1/tests.tsv`, `r1/failures.txt`): **every prediction held.**
- `test-kit258-duplicate-completeness.R` 76 passed, 26 failed, 0 errors: E1 to E3 (12 tests) fail only on the kept
  row, at E1 "1" "A", E2 "2" "A", E3 "2" "Q"; E4 kernels "A" for both kernels, E4 S4 "A" in both omics; T1 (b) `kept`
  and `replayed_kept` "A" and the rank marker NULL in the record and the replayed record, both omics; T2 the rank
  marker NULL for lipid and metab. G1, G2, G3, T1 (a) and every other expectation pass.
- `test-kit214-sample-columns.R` 176 passed, 0 failed (with the declared edits).
- `test-metab-01m-...` 15 passed, 3 failed: exactly the three new pins. `test-metab-02u-...` 16 passed, 1 failed:
  the one new pin.

Verdict: R1 met (A1). The fix may be written against these tests.

## Fix as far as the owner's decided rule (pending OD-1 and OD-2)

The ticket runs S1 on the fix commit (section 10: "base = staging after KIT-214, head = the fix"), so the fix is written
only as far as VAL-254 decided: at the four sites the most measured (not `NA`, today's definition, unchanged) values
among the averaged columns first, then the highest mean of the measured values, then the first in the table
(`slice_max(with_ties = FALSE)`, unchanged); one internal helper `keepDuplicateResolutionRow` (`R/func_omics_qc_shared.R`)
reads the replay marker `duplicate_resolution_rank` (`measured_then_mean` / `mean`), appended to
`workflowReplayRuleMarkers`; the step's `fixed` gains it. Not written yet (wait for OD-1 and OD-2): E5, E6, NEWS, the
ledger row, the help, log and exported-function texts, `man/`.

## S1 (control on the fix commit, base 63143aa3; 1 scripts job, ceiling 60 job-minutes)

Script `dup-rank-control.R` (KIT-214's control re-pointed as Re-check 1 R1-02 specifies: inputs and variants kept,
classes by completeness on the base side, criteria A4 (a) to (e)). Decides A4 and stop rules 1 and 2.

PREDICTED (written before the run was read):
- (a), (b): every row equal; non-vacuous for lipid and metabolite; every equally measured group keeps its row.
- Classes `e` (unequally measured duplicates): e1_lipid, e1_metab, lipid_msdial_dup_import; e1's kept row does not
  change (row B is both more measured and higher).
- (c), (d1), (d2): every row equal.
- (e): exactly one id changes its kept row, in both variants: lipid_msdial_dup_import / LCMS_Pos / P001, base row 2
  ('TG 52:3', 3 of 4 measured, mean 177.333333) to head row 1 ('PC 34:1', 4 of 4, mean 150.5). Any other changed id is
  stop rule 1.

## Not checked

- No R on the laptop; every result above comes from a runner.
