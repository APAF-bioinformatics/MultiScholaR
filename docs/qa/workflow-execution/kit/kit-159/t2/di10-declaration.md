# T-COV I10 diagnostic and the O4 half-width ratio: declaration, written before the request

Owner rulings: methods ledger VAL-109 (I10 on T2's line, about 2 core-hours; O4 the half-width
ratio only for rows consumers use). Audit register `AUDIT-KIT-159-impl-prerelease-20261124.md`,
I10 and section 6 row 2.

## Request D: the I10 diagnostic (four 500-feature cells)

**Hypothesis under test (I10).** The `fail` statuses of the C4 and C4e `accuracy` and
`threshold_count` rows come from T2's generator, not from the jackknife: T2 draws the feature
baselines anew in every simulation and shifts 20 % of the features by +1 in one direction, so
median normalisation moves every estimate by a shift (about 0.19 log2 in T2's C4e targets) that
is a feature-level quantity common to every unit, which a delete-one-unit jackknife cannot see.

**Generator `fixed_symmetric` (msrval `valCovFixedFeatures`).** Per cell one seeded draw
(seed 2,100,000,000 + the cell's index) of the 500 feature baselines (normal, mean 20, SD 1); the
changed features (the first 100) split alternately +1 and -1 per step of the group effects;
residual SD 0.3, per-sample offset SD 0.1, block effect SD 0.1 (C4e) and median normalisation as
in T2. Simulation and target seeds as in T2 (by cell and index only). Everything else is T2's.

**Cells and sizes.** C4-143 (replicate, g = 4), C4-145 (replicate, g = 6), C4e-147 (block, g = 5),
C4e-149 (block, g = 14); 2,000 simulations and 20,000 target simulations each; plan
`t2/di10-plan.tsv`, request D, one range per cell, 4 concurrent chunks.

**Mechanism checks (must hold for the result to be read as a test of the hypothesis).**
- M1: the C4e targets of `bias` lie in [-0.02, 0.02] (T2: -0.210 at g = 5, -0.186 at g = 14);
- M2: the C4e targets of `rej` (the plain arm's null rejection share at p < 0.05) lie in
  [0.03, 0.08] (T2: 0.366 and 0.633).

If M1 or M2 fails, the generator did not remove the shift: the run is not a test of I10; it is
recorded, not re-run, and reported to the owner.

**Decision rule (the statistics' rows of the four cells, the band as T2 ran it: fail below
level - 0.02).**
- **Confirmed**: M1 and M2 hold and no row of the four cells is `fail` at any of the 8 levels.
- **Refuted**: M1 and M2 hold and none of the 11 rows that failed at 95 % in T2 (C4-143 `mae`
  0.918, `fcs` 0.908, `fcs_diff` 0.9285; C4-145 `mae` 0.898, `fcs` 0.8995; C4e-147 `mae` 0.745,
  `bias` 0.7565, `rej` 0.7065; C4e-149 `mae` 0.7065, `bias` 0.707, `rej` 0.726) rises by more
  than 3 of T2's Monte Carlo SE at 95 % (about 0.015 to 0.03).
- **Partial**: anything else (some rows rise beyond Monte Carlo error, at least one row still
  fails). The four cells are not replaced; reported to the owner.

**What the result decides.** Confirmed: this run's rows replace the four cells' rows in the
pinned table (VAL-109 (2)); ticket 10.2 is revised to the fixed symmetric generator; `USAGE.md`'s
offered-regime table is regenerated; KIT-160, 161, 162, 166 and 170 are told which of their rows
change. Refuted or partial: nothing in the pinned table changes; reported. It is evidence: an
unfavourable valid result is recorded, never re-run (one correction pass only for an instrument
defect, VAL-54).

## Request H: the half-width ratio of the conservative rank rows consumers use (O4)

C4e-148 (block, g = 5, rank: KIT-162 at 95 %, 98.33 % and 99.17 %) and C4-146 (replicate, g = 6,
rank) at 200 simulations and 2,000 target simulations, T2R's generator (`redraw`; the rank cells'
re-pinned effect). Only the column `half_width_ratio` (mean half-width over the simulated SD of the
full statistic times the level's normal quantile; 1.00 for a correctly scaled interval) is read;
the coverage of this small run is not a regime and does not change any row's status. C4e g = 14
rank (no observed miss, no declared consumer) is not run: it licenses only its own statistic (O4).

## Binding, ceiling, charge

One `run-tests` scripts request at the implementation commit named in it, two scripts (D, then H),
`timeout_minutes` 120 each. Expected: D about 20 min of chunk time (T2C measured 2,089 s of process
time for these cells) plus set-up; H about 15 min plus the rank-effect calibration. Ceiling
4 core-hours, charged to T2's line (39 ceiling, 26.93 used).
