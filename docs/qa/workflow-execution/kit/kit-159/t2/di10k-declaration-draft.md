# DI10K: the confirmatory cell in KIT-162's composition, DRAFT declaration (not frozen, not sent)

Status: DRAFT for the harness auditor's check, with the generator change it needs. Nothing here is
frozen, and no coverage request has been sent. The final declaration (`t2/di10k-declaration.md`) is
committed alone, at the commit the request binds, after the audit.

Owner ruling: methods ledger VAL-163 (3) ("One more confirmatory cell is run in KIT-162's
composition (about 4,800 features, about 40 changed, all in one direction), in the form of DI10b
and DI10c, declared before it is sent"). Coordinator's ruling of 2026-11-25: the smallest generator
change (a new cell with its own index), the generator's group effects kept, KIT-162's level
structure not built. Budget line DI10K, 6 core-hours (KIT-159 r1.24, section 39). The prefix is
`di10k` because `t2/di10d-*` holds the I10 diagnostic's request D.

## What it measures

DI10b (run 45336507597) ran cell C4e-147 (C4e, delete-one-block jackknife, g = 5) under
`fixed_onedir`. It passed (lowest coverage at 95 %: `rej` 0.9485), but with 20 % of 500 features
changed: a large median-normalisation shift (target `bias` -0.234, `rej` 0.428). KIT-162 has about
40 changed of about 4,800 (0.8 %). DI10K runs **cell C4e-170, family `c4e:4840acc`, under
generator `fixed_onedir`**. This is DI10b's cell at KIT-162's composition:

| Item | DI10K (C4e-170) | DI10b (C4e-147) | KIT-162 (`tickets/kit/KIT-162.md`) | Same as KIT-162? |
|---|---|---|---|---|
| Units | 5 blocks, delete-one-block jackknife (form (iii)), g = 5 | same | 5 TMT10plex mixtures, leave one mixture out | yes |
| Samples per unit | 10: 5 groups x 2 | same | 10 channels: 4 UPS levels x 2 + `Norm` x 2 | count yes; structure no (below) |
| Features | 4,840 | 500 | about 4,800 proteins | yes |
| Changed | 40 (0.83 %) | 100 (20 %) | 48 UPS1 spiked, about 40 after filters (about 0.8 %) | yes |
| Direction | all up relative to group A | all up | all in one direction within each contrast | yes |
| Group effects (log2) | 0, 0.5, 1, 1.5, 2 (groups A to E) | same | levels 1, 0.667, 0.5, 0.125 = 0, -0.585, -1, -3 log2 relative to level 1; `Norm` a pool (286.5 fmol), no contrast | **no** |
| Contrasts | 4: B, C, D, E each against A (0.5, 1, 1.5, 2 log2) | same | 6 pairwise among the four levels (0.41, 0.59, 1, 2, 2.41, 3 log2); primary `L1 - L05` (1 log2) | **no** |
| Feature-wise variance | one residual SD 0.3 for every feature; baseline N(20, 1) drawn once (seed 2,100,000,170); block-by-feature SD 0.1; per-sample offset SD 0.1 | same (seed 2,100,000,147) | intensity-dependent, heteroscedastic; protein-by-plex effects; TMT ratio compression; missing values (most at the 0.125 level); three technical runs summed | **no** |
| Normalisation and model | median normalisation; limma, group + block terms, moderated (`plain`) and unmoderated arms | same | the package's route (its own normalisation, RUV designs, `~ 0 + group (+ batch)`), declared in KIT-162 | **no** |
| Statistics | `mae` = median \|est - truth\|, `bias` = median (est - truth), over all 4,840 x 4 entries; `rej` = share of unchanged-feature contrasts with plain p < 0.05 (19,200 per simulation) | same (2,000 entries; 1,600 null) | detection on every spiked protein, accuracy on the unique UPS subset, calibration on the unchanged proteins | `rej` corresponds to calibration on the unchanged; `mae` and `bias` are dominated by unchanged features and do **not** correspond to accuracy on the 40 |

What the result therefore says: whether the delete-one-block jackknife interval of C4e g = 5 keeps
its coverage for these three statistics at KIT-162's **proportion** of one-directional changes,
under the harness's homoscedastic generator. It says nothing about KIT-162's level structure,
`Norm` group, pairwise contrasts, heteroscedastic variance, ratio compression or missing values.
**No statistic of this cell counts true positives.** The statistics KIT-162 computes on the
40 spiked proteins (detection, accuracy on the unique subset) are not simulated here. Their
intervals rest on 40 features, and this cell does not speak to them (see "What it cannot show").

Code path: `msrval/R/coverage.R` at the commit this request will bind. The generator, fit and
statistics are DI10b's (`valCovMatrixData`, `valCovFixedFeatures`, `valCovFit`, `valCovStats`'s C4e
branch). The cell differs from DI10b's code path only by family `c4e:4840acc`. Its 4,840 features
give 40 changed (`n_feat == 4840L`), and it keeps the C4e statistics and effect scale 1, not the
rank statistic and the rank-calibrated scale of the 4,840-feature rank cells. The tests that show
this are in `test-di10k.R`, and the diff in `outcome.md` (DI10K).

Size: 2,000 simulations; targets from 20,000 target simulations. Plan `t2/di10k-plan.tsv`, request
DK, four ranges of 500 simulations and 5,000 target simulations, `generator` `fixed_onedir`.

## Mechanism checks, fixed before the run

At 40 of 4,840 the one-directional shift is small, and the checks say what it is. Derivation
(normal theory; per-sample values across features have SD s = sqrt(1 + 0.09 + 0.01) = 1.049).
With a share p of features moved by e in a sample, the sample median moves by delta solving
(1 - p) Phi(delta / s) + p Phi((delta - e) / s) = 1/2. Median normalisation subtracts delta from
every feature of that sample, so each unchanged feature's estimate for group j against A is
centred at -delta_j. The contrast SE is 0.3 x sqrt(2 / 10) = 0.134, because offsets cancel
under normalisation and the block terms absorb the block effects.

| Quantity | p = 0.2 (DI10b's composition) | p = 40 / 4,840 (this cell) |
|---|---|---|
| delta for e = 0.5, 1, 1.5, 2 | 0.098, 0.186, 0.254, 0.298 | 0.0040, 0.0072, 0.0093, 0.0103 |
| expected `bias` target (median of the mixture of N(-delta_j, 0.134)) | -0.210 | **-0.0077** |
| expected `rej` target (normal approximation, at 1.96) | 0.368 | **0.0504** |
| expected `mae` target (median \|N(-delta_j, 0.134)\|, unchanged features) | 0.212 | 0.0907 |
| measured | T2 redraw g = 5: `bias` -0.210, `rej` 0.366; DI10b fixed draw: `bias` -0.234, `rej` 0.428, `mae` 0.235 | (this run) |

The derivation reproduces T2's redraw values at p = 0.2 (bias -0.210 against -0.210; rej 0.368
against 0.366). DI10b's fixed feature draw moved them by a factor of about 1.1 to 1.2.

- **N1k (the shift is present, with its sign and size): the target of `bias` is at most
  -0.0035.** The threshold is the midpoint, rounded, of D's `fixed_symmetric` C4e-147 target
  (0.0008, where there is no shift) and the expected -0.0077. The target's Monte Carlo SE is
  about 0.00002 (a per-simulation SD of about 0.003 over 20,000), so the check decides cleanly
  between "shift present" and "no shift". It **shows** that the cell carries KIT-162's
  one-directional normalisation shift. It **cannot show** coverage under a large shift: here the
  shift is about 6 % of one contrast SE (0.0077 / 0.134), against about 175 % in DI10b. DI10b
  remains the evidence for large shifts.
- **N2k (the composition is KIT-162's, not DI10b's): the target of `rej` lies in [0.045,
  0.060].** The expected value is 0.0504, plus the moderated arm's calibration offset seen in D
  (0.0506 at 500 features), so about 0.051. The band is wider than Monte Carlo noise (target SE
  about 0.00001) to allow for the fixed draw and the moderated degrees of freedom. It is far from
  any 20 % composition (0.366 to 0.428). It **shows** that the cell is at low contamination and
  that null calibration is not visibly inflated by the shift. It **cannot show** an inflation
  below about 0.01: the expected inflation from the shift is 0.0004.

`mae`'s target (expected 0.091; D at 500 features 0.091) does not separate shift from no shift
here. It is reported, not used as a check. The targets are read from
`di10k-coverage-detail.tsv` (column `target`), as DI10b's and DI10c's were.

## Pass rule and readings (band as T2, D, DI10b and DI10c ran it)

A statistic fails at a level when its coverage is below level - 0.02. The rule applies to all
three statistics at all 8 levels (24 values), whatever status label a row carries. A `degenerate`
or `conservative` label is reported beside the value as observed.

- **Mechanism checks not met** (N1k or N2k fails): the cell did not carry KIT-162's composition
  and shift, so **the run does not test the question**. The result is recorded, reported and not
  re-run.
- **Pass** (N1k and N2k hold, and no coverage is below level - 0.02 at any of the 8 levels for
  any of the three statistics): **DI10b's licence for KIT-162's re-pinned C4e g = 5 rows
  (`accuracy`, `threshold_count`) stands at KIT-162's composition**, for statistics of these
  kinds, under the differences listed above.
- **Fail** (N1k and N2k hold, and at least one coverage is below level - 0.02): **KIT-162 needs
  intervals calibrated to its own composition before it runs** (ticket coverage regimes).

## Reported beside the rule (for the owner)

- Each statistic's coverage at each of the 8 levels, beside DI10b's (C4e-147, `fixed_onedir`) and
  D's (`fixed_symmetric`), and the misses per level (out of 2,000).
- For the lowest-covering statistic at each level: the exact (Clopper-Pearson) 95 % interval of
  its coverage from its misses out of 2,000.
- `se_inflation` and `half_width_ratio` per statistic, beside DI10b's (`mae` 0.979, `bias` 0.996,
  `rej` 0.900) and D's.
- The targets of the three statistics against the expected values above.
- Event counts. `rej` counts null rejections among 19,200 unchanged-feature contrasts per
  simulation (about 980 expected). `mae` and `bias` are medians over 19,360 entries, of which 160
  are changed-feature contrasts. None of the three statistics counts true positives, so no event
  count on the 40 changed features exists here, and none is claimed.

## What it cannot show

The coverage of intervals for statistics over the 40 spiked proteins (sensitivity, accuracy on the
unique UPS subset), which rest on about 40 features per contrast. Also KIT-162's real-data features
listed in the table. Covering those would need further statistics in the cell (for example
sensitivity and median absolute error on the changed features). That is outside the coordinator's
ruling of 2026-11-25 and is for the auditor and the owner.

## Binding, ceiling, charge (to be fixed after the cost probe)

One `run-tests` scripts request at the commit that adds the final declaration and the plan:
`val-tests.R <out> --files metrics --coverage-plan t2/di10k-plan.tsv --request DK --nsim 2000
--chunks 4`.

Expected cost, from three measured anchors:
- **Measured on this cell** (run 45403024967, `test-di10k.R`, one process): one simulation (six fits
  of 4,840 x 50) plus one target simulation plus one further fit took 0.995 s. That is about
  0.75 s per simulation and 0.13 s per target: 2,000 x 0.75 + 20,000 x 0.13 = 4,100 s, or
  1.1 process-hours. This is one timing and includes set-up inside the test.
- DI10b at 500 features (same code path): 296.6 s for 2,000 simulations and 516.7 s for 20,000
  targets. Scaled linearly by 4,840 / 500 = 9.68, that is 7,873 s, or 2.2 process-hours.
- T2's C4e-148 (4,840 features, the rank statistic): 11,805 s, or 3.3 process-hours (an upper
  anchor).

DI10b ran 813 process-seconds in 205 s of wall on 4 chunks (factor 3.97). The likely 1.1 to
2.2 process-hours therefore give 0.28 to 0.55 h per chunk, plus about 0.2 h of set-up: 0.48 to
0.75 h of wall, or 1.9 to 3.0 core-hours (4.1 at the rank anchor).

Proposed: script `timeout_minutes` 50. With set-up the job is at most about 62 minutes, a ceiling
of 4.1 core-hours. Line DI10K has 6 - 1.61 = 4.39 core-hours left after the two test runs. The cost
probe (`t2/di10k-probe-plan.tsv`: 20 simulations and 200 targets in four ranges, request DKP, about
0.7 core-hours, almost all of it set-up) would leave 3.6 and lower the script limit to about
42 minutes. That is still above the 0.55 h likely upper bound per chunk, but below the rank
anchor. Recommendation: no separate probe; the cell's own timing above stands in for it. The
coordinator decides.

One run. The result is a result: it is recorded whatever it shows and never re-run. A failure in
the mechanics with no result produced falls under VAL-54's instrument rule: fix the mechanics
only, at most three further runs, each recorded, and no change to a threshold, rule or the cell.
