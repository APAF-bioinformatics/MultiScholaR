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

## Not checked

- No R on the laptop; every result above comes from a runner.
