# DI10b: the confirmatory cell for the I10 re-pin, declaration written before the request

Owner ruling: methods ledger VAL-137 (1) ("one confirmatory cell with one-directional changes is
run and reported"). Audit register `AUDIT-KIT-159-impl-recheck-20261124.md`, S2. Budget line DI10b
(1.5 core-hours from the ticket's unassigned total; coordinator 2026-11-24).

## What it measures

The I10 diagnostic (request D, run 45316883316) changed two things at once against T2's generator:
it fixed the feature draw per cell (conditioning the target on one feature configuration) and it
split the 100 changed features +1 and -1 (so median normalisation is not shifted). The re-pinned
C4 and C4e `accuracy` and `threshold_count` rows rest on both. The spike-in benchmarks (KIT-160,
161, 162) change in one direction. DI10b separates the two: **cell C4e-147 (delete-one-block, g = 5,
KIT-162's shape) under the same fixed feature draw as D (same baseline seed, same simulation and
target seeds), with every changed feature moved up (+1 per step of the group effects), as T2's
changes.** Generator `fixed_onedir` (msrval `valCovFixedFeatures`); everything else is D's.

It measures, for the three statistics of C4e (`mae`, `bias`, `rej`) at the 8 levels, the coverage
of the delete-one-block jackknife interval of the target conditional on that feature draw, at 2,000
simulations with targets from 20,000 (plan `t2/di10b-plan.tsv`, request DB, four ranges of 500).
Hypothesis: with the features fixed, the one-directional normalisation shift is a fixed bias that
the conditional target absorbs, so the coverage holds without the symmetric split.

## Mechanism checks (the cell reproduces the one-directional shift)

- N1: the target of `bias` is at or below -0.10 (T2 at g = 5: -0.210; D: 0.0008).
- N2: the target of `rej` (the plain arm's null rejection share at p < 0.05) is at least 0.20
  (T2: 0.366; D: 0.0506).

If N1 or N2 fails, the cell did not reproduce T2's shift: it is not a test of the question; it is
recorded, not re-run, and reported.

## Reading (band as T2 and D ran it: a statistic fails at a level when coverage < level - 0.02;
`degenerate` by r1.13's rule, reported as such)

- **Conditioning suffices** (the symmetric split was not needed for C4e g = 5's coverage): N1 and
  N2 hold and no statistic fails at any of the 8 levels.
- **The symmetric split was needed for the coverage figures**: N1 and N2 hold and at least one
  statistic fails at a level where D's `fixed_symmetric` C4e-147 row passed (D: `mae` 0.989 to 1,
  `bias` 0.9605 to 0.9965, `rej` 0.9835 to 0.999). What that would suggest (for the owner, not
  decided here): the re-pinned C4/C4e `accuracy` and `threshold_count` rows license balanced
  changes only, and KIT-160, 161 and 162 need ticket coverage regimes with their own composition.
- **Anything else** (for example a `degenerate` statistic and no failure): reported as observed.

## What it decides

Nothing by itself (VAL-137 (1): "run and reported"). No pinned row changes; the result goes to the
owner with the numbers. One run: an unfavourable valid result is recorded, never re-run; a failure
in the mechanics with no result is VAL-54's instrument rule (at most three further runs).

## Binding, ceiling, charge

One `run-tests` scripts request at the fix commit that adds `fixed_onedir` (the commit is named in
the request and in `runs.tsv`), `val-tests.R <out> --files metrics --coverage-plan
t2/di10b-plan.tsv --request DB --nsim 2000 --chunks 4`, `timeout_minutes` 60. Expected: C4e-147
took about 804 s of process time in D, so about 4 min of chunk time on 4 cores plus set-up, about
0.5 to 0.8 core-hours. Ceiling 1.5 core-hours, charged to line DI10b.
