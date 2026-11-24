# DI10c: the one-directional confirmatory cell for C4 g = 4, declaration written before the request

Owner ruling: methods ledger VAL-155 ("One cell, C4-143, is run under one-directional changes with
2,000 simulations, declared before it is sent ... If it passes, the two tickets use the pinned
rows; if it fails, they need intervals calibrated to their own composition before they run").
Audit register `AUDIT-KIT-159-impl-recheck-20261124.md`, re-audit 2, R2.4. Budget line DI10c
(1.5 core-hours, ticket r1.21, section 11 and section 36).

## What it measures

DI10b (run 45336507597) showed that the C4e g = 5 block-form interval keeps its coverage when the
feature draw is fixed and every change goes one way, but loses its conservative margin. KIT-160 and
KIT-161 use a different regime: C4, jackknife form (ii) `replicate`, g = 4 (cell C4-143: two groups
A and B of 4 replicates each, `valCovLayout("replicate", 4)`; 500 features). Its pinned rows come
from the I10 diagnostic (request D, run 45316883316) under `fixed_symmetric`, where its lowest
coverages at 95 % were `fcs_diff` 0.9685, `fcs` 0.971 and `mae_diff` 0.981. It was never run with
one-directional changes under a fixed feature draw.

DI10c runs **cell C4-143 under generator `fixed_onedir`**: the same fixed feature draw as D and
DI10b's generator (`valCovFixedFeatures`: baselines drawn once with seed 2,100,000,000 + the cell's
index, so the same baselines as D's C4-143; the first fifth of the features changed), every changed
feature moved +1 log2 in group B (as T2 and the spike-in benchmarks), the same simulation and
target seeds as D, residual SD 0.3, per-sample offset SD 0.1, median normalisation, limma. The code
path is DI10b's: `msrval/R/coverage.R` and `val-tests.R` are unchanged between DI10b's commit
647706dd and the commit this request binds (`git diff 647706dd <commit> -- coverage.R
val-tests.R` is empty). Only the cell differs.

It measures, for the six statistics of C4 replicate (`mae`, `mae_diff` (class accuracy); `sens`,
`sens_diff`, `fcs`, `fcs_diff` (class threshold_count)) at the 8 levels, the coverage of the
delete-one-replicate-pair jackknife interval of the target conditional on that feature draw, at
2,000 simulations with targets from 20,000 (plan `t2/di10c-plan.tsv`, request DC, four ranges of 500
simulations and 5,000 target simulations).

## Mechanism checks (the cell carries the one-directional shift), fixed before the run

C4's replicate statistics include neither DI10b's `bias` nor its `rej`, and the code path is held
as DI10b's, so the checks read the two C4 targets that the normalisation shift moves. Why they move:
with 100 of 500 features up by 1 in group B, the median of each B sample rises by about 0.19 log2
(the median of 0.8 N(20, 1) + 0.2 N(21, 1)), so median normalisation shifts every unchanged
feature's B-minus-A estimate by about -0.19 and every changed feature's by the same amount. Under
`fixed_symmetric` the B medians do not move.

- **N1c (the shift of the estimates):** the target of `mae` is at least **0.175**. D
  (`fixed_symmetric`, C4-143): 0.1443; T2 (redraw, one-directional, C4-143): 0.2049. Expected under
  the shift: about 0.20 (median of |N(-0.19, 0.21)|; 0.21 = sqrt(2 x 0.09 / 4), against 0.143 =
  0.674 x 0.21 without it). The threshold is the midpoint of D and T2, rounded.
- **N2c (false calls among unchanged features, the C4 counterpart of `rej`):** the target of `fcs`
  (share of unchanged features called at BH < 0.05 by the moderated arm) is at least **0.03**. D:
  0.0103; T2: 0.0475. The threshold is the midpoint of D and T2, rounded.

The targets are read from `di10c-coverage-detail.tsv` (column `target`), as DI10b's were.

## Pass rule and reading (band as T2, D and DI10b ran it)

A statistic fails at a level when its coverage is below level - 0.02. The rule is on the coverage
values of all six statistics at all 8 levels (48 values), whatever status label a row carries; a
`degenerate` or `conservative` label is reported beside it as observed.

- **Mechanism checks not met** (N1c or N2c fails): the cell did not carry the one-directional
  shift; **the run does not test the question**. Recorded, not re-run, reported.
- **Pass** (N1c and N2c hold, no coverage below level - 0.02 at any of the 8 levels for any of the
  six statistics): **KIT-160 and KIT-161 may use the pinned C4 g = 4 rows** (`accuracy`,
  `threshold_count`) for their one-directional spike-in designs.
- **Fail** (N1c and N2c hold, at least one coverage below level - 0.02): **KIT-160 and KIT-161 need
  intervals calibrated to their own composition before they run** (ticket coverage regimes).

Reported with the result, for the owner: each statistic's coverage at every level beside D's
`fixed_symmetric` C4-143 value and the pinned table's C4 g = 4 row; the misses per level; the
`se_inflation` against D's. No pinned row changes by this run (VAL-155: the result decides what
KIT-160 and KIT-161 may use).

## Binding, ceiling, charge

One `run-tests` scripts request at the commit that adds this declaration and the plan (named in the
request and in `runs.tsv`): `val-tests.R <out> --files metrics --coverage-plan t2/di10c-plan.tsv
--request DC --nsim 2000 --chunks 4`, script `timeout_minutes` 12. Expected: C4-143 took 163 s of
process time in D (59 s simulations, 104 s targets), so under a minute per chunk on 4 cores plus
library set-up (6 to 9.5 min): about 0.15 to 0.25 h of wall, 0.6 to 1.0 core-hours. Ceiling 1.5
core-hours (9.5 min set-up + 12 min script limit = 21.5 job-minutes, 1.43 core-hours), charged to
line DI10c.

One run. The result is a result: recorded whatever it shows, never re-run. A failure in the
mechanics with no result produced is VAL-54's instrument rule (fix the mechanics only, at most three
further runs, each recorded; no threshold, rule or cell changes).
