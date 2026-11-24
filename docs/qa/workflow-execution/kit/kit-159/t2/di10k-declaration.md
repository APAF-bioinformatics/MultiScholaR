# DI10K: the confirmatory cell in KIT-162's composition, declaration written before the request

Owner ruling: methods ledger VAL-163 (3) ("One more confirmatory cell is run in KIT-162's
composition (about 4,800 features, about 40 changed, all in one direction), in the form of DI10b
and DI10c, declared before it is sent"). Coordinator's rulings of 2026-11-25: the smallest generator
change (cell C4e-170, its own index), the generator's group effects kept, KIT-162's level structure
not built, no P1-form rank cell (the owner's decision). Harness auditor: re-audit 4
(`AUDIT-KIT-159-impl-recheck-20261124.md`), repairs R1 to R10 applied here to
`t2/di10k-declaration-draft.md`, which stays as the audited draft. Budget line DI10K, 6 core-hours
(KIT-159 r1.24 and r1.25). The prefix is `di10k` because `t2/di10d-*` holds the I10 diagnostic's
request D.

## What it measures

DI10b (run 45336507597) ran cell C4e-147 (C4e, delete-one-block jackknife, g = 5) under
`fixed_onedir`. It passed (lowest coverage at 95 %: `rej` 0.9485), but with 20 % of 500 features
changed, a large median-normalisation shift (target `bias` -0.234, `rej` 0.428). KIT-162 has about
40 changed of about 4,800 (0.8 %). DI10K runs **cell C4e-170, family `c4e:4840acc`, under
generator `fixed_onedir`**: DI10b's cell at KIT-162's composition.

**Expected reading (R1).** At KIT-162's composition the one-directional shift is about 6 % of a
contrast SE (expected `bias` -0.0077 against SE 0.134; expected `rej` inflation 0.0004). The
expected reading is that coverage equals the symmetric case. What this cell newly tests is the C4e
g = 5 block jackknife for `bias` and `rej` at 4,840 features (the pinned rows come from 500). A
large shift is DI10b's evidence.

| Item | DI10K (C4e-170) | DI10b (C4e-147) | KIT-162 (`tickets/kit/KIT-162.md`) | Same as KIT-162? |
|---|---|---|---|---|
| Units | 5 blocks, delete-one-block jackknife (form (iii)), g = 5 | same | 5 TMT10plex mixtures, leave one mixture out | yes |
| Samples per unit | 10: 5 groups x 2 | same | 10 channels: 4 UPS levels x 2 + `Norm` x 2 | count yes; structure no |
| Features | 4,840 | 500 | about 4,800 proteins | yes |
| Changed | 40 (0.83 %) | 100 (20 %) | 48 UPS1 spiked, about 40 after filters (about 0.8 %) | yes |
| Direction | all up relative to group A | all up | all in one direction within each contrast | yes |
| Group effects (log2) | 0, 0.5, 1, 1.5, 2 (groups A to E) | same | levels 1, 0.667, 0.5, 0.125 = 0, -0.585, -1, -3 log2 relative to level 1; `Norm` a pool (286.5 fmol), no contrast | **no** |
| Contrasts | 4: B, C, D, E each against A (0.5, 1, 1.5, 2 log2) | same | 6 pairwise among the four levels (0.41, 0.59, 1, 2, 2.41, 3 log2); primary `L1 - L05` (1 log2) | **no** |
| Feature-wise variance | one residual SD 0.3 for every feature; baseline N(20, 1) drawn once (seed 2,100,000,170); block-by-feature SD 0.1; per-sample offset SD 0.1 | same (seed 2,100,000,147) | intensity-dependent, heteroscedastic; protein-by-plex effects; TMT ratio compression; missing values (most at the 0.125 level); three technical runs summed | **no** |
| Normalisation and model | median normalisation; limma, group + block terms, moderated (`plain`) and unmoderated arms | same | the package's route (its own normalisation, RUV designs, `~ 0 + group (+ batch)`), declared in KIT-162 | **no** |
| Statistics (R3) | `bias` = median (est - truth) over all 4,840 x 4 entries; `rej` = share of unchanged-feature contrasts with plain p < 0.05 (19,200 per simulation); `mae` = median \|est - truth\| over all entries | same (2,000 entries; 1,600 null) | `bias` corresponds to KIT-162's frame guard (HeLa median bias, deciding, `require: match`); `rej` corresponds to G-cal (deciding; KIT-162 pools six pairwise contrasts, the cell four against A); `mae` corresponds to nothing KIT-162 decides (its UPS M2 is a key secondary, and the cell's `mae` is dominated by unchanged features) | by kind, for `bias` and `rej` |

**Code path.** `msrval/R/coverage.R` at 9adab5fe, unchanged at this commit. The generator, fit and
statistics are DI10b's (`valCovMatrixData`, `valCovFixedFeatures`, `valCovFit`, `valCovStats`'s C4e
branch). The cell differs only by family `c4e:4840acc`. Its 4,840 features give 40 changed
(`n_feat == 4840L`), and it keeps the C4e statistics and effect scale 1, not the rank statistic and
rank-calibrated scale of the 4,840-feature rank cells. Shown by `test-di10k.R` (runs 45400680535 and
45403024967; `outcome.md`, DI10K).

**Size.** 2,000 simulations; targets from 20,000 target simulations. Plan `t2/di10k-plan.tsv`,
request DK, four ranges of 500 simulations and 5,000 target simulations, `generator`
`fixed_onedir`.

## Configuration checks (the cell carries KIT-162's composition), fixed before the run (R2)

Derivation (normal theory). Per-sample values across features have SD s = sqrt(1 + 0.09 + 0.01)
= 1.049. With a share p of features moved by e in a sample, the sample median moves by delta
solving (1 - p) Phi(delta / s) + p Phi((delta - e) / s) = 1/2. Median normalisation subtracts delta
from every feature of that sample, so each unchanged feature's estimate for group j against A is
centred at -delta_j. The contrast SE is 0.3 x sqrt(2 / 10) = 0.134.

| Quantity | p = 0.2 (DI10b's composition) | p = 40 / 4,840 (this cell) |
|---|---|---|
| delta for e = 0.5, 1, 1.5, 2 | 0.098, 0.186, 0.254, 0.298 | 0.0040, 0.0072, 0.0093, 0.0103 |
| expected `bias` target | -0.210 | **-0.0077** |
| expected `rej` target (normal approximation) | 0.368 | **0.0504** |
| expected `mae` target | 0.212 | 0.0907 |
| measured | T2 redraw g = 5: `bias` -0.210, `rej` 0.366; DI10b fixed draw: `bias` -0.234, `rej` 0.428, `mae` 0.235 | (this run) |

- **N1k: the target of `bias` is at most -0.0035.** The threshold is the midpoint, rounded, of
  D's `fixed_symmetric` C4e-147 target (0.0008, no shift) and the expected -0.0077. Robustness:
  over random baseline draws the predicted target is -0.0077 ± 0.0011 (1,500 draws, re-audit 4;
  none above -0.0035), so N1k fails only on a misconfigured cell (no shift, symmetric changes, or
  the 20 % composition).
- **N2k: the target of `rej` lies in [0.045, 0.060].** Expected 0.0504 plus the moderated arm's
  offset seen in D (0.0506 at 500 features), about 0.051. The band separates this composition from
  any 20 % composition (0.366 to 0.428). It cannot see an inflation below about 0.01; the expected
  inflation from the shift is 0.0004.

`mae`'s target (expected 0.091) is reported and not checked. The targets are read from
`di10k-coverage-detail.tsv` (column `target`).

## Pass rule and readings (band as T2, D, DI10b and DI10c ran it)

A statistic fails at a level when its coverage is below level - 0.02. The rule is on the coverage
values of all three statistics at all 8 levels (24 values), whatever status label a row carries.
A `degenerate` or `conservative` label is reported beside it as observed.

- **Configuration checks not met** (N1k or N2k fails): the cell did not carry KIT-162's
  composition; **the run does not test the question**. Recorded, not re-run, reported.
- **Pass** (N1k and N2k hold, no coverage below level - 0.02 at any of the 8 levels for any of the
  three statistics) (R4): **the C4e g = 5 `accuracy` (bias) and `threshold_count` (null share)
  rows license KIT-162's frame guard and G-cal intervals at its composition (feature count,
  proportion, direction), under the listed differences; `mae`'s result is reported and licenses
  nothing KIT-162 decides.**
- **Fail** (N1k and N2k hold, at least one coverage below level - 0.02) (R5): **KIT-162 needs a
  ticket coverage regime at its own composition for the guard whose statistic failed, before it
  runs.** (A failure of `mae` alone is reported; it decides nothing of KIT-162.)

## Reported beside the rule (for the owner)

- Each statistic's coverage at each of the 8 levels, beside DI10b's (C4e-147, `fixed_onedir`) and
  D's (`fixed_symmetric`); the misses per level (of 2,000).
- For the lowest-covering statistic at each level: the exact (Clopper-Pearson) 95 % interval of its
  coverage.
- (R8) For every statistic whose misses exceed nominal: P(X >= misses) under the nominal miss rate
  (binomial, n = 2,000) and the exact 95 % interval. A statistic below nominal at p < 0.01 is
  reported as anti-conservative with its measured coverage, as VAL-163 (1) did for `sens_diff`,
  whatever the band says.
- `se_inflation` and `half_width_ratio` per statistic, beside DI10b's (`mae` 0.979, `bias` 0.996,
  `rej` 0.900).
- The targets of the three statistics against the expected values above.
- Event counts: `rej` counts null rejections among 19,200 unchanged-feature contrasts per
  simulation (about 980 expected); `mae` and `bias` are medians over 19,360 entries, of which 160
  are changed-feature contrasts. None of the three counts true positives.
- (R10) **No pinned row changes by this run; `di10k-coverage-table.tsv` shares the pinned table's
  row key (C4e, jackknife, block, n = 5) and is never merged into it or read in its place.**

## What it cannot show (R6)

- **KIT-162's primary endpoint P1** (pAUC to FPR 0.01 over every UPS key, class `rank`) **is not
  simulated here.** Its licence is the C4e g = 5 rank rows from cell C4e-148, which simulate 40 of
  4,840 features changed in one direction (`redraw`, scale 0.35, four contrasts pooled, 160
  positive entries). P1's single contrast with 40 positives is not simulated. Whether to simulate
  it is the owner's decision (re-audit 4, R4.2).
- The secondaries that count true positives (S1 to S3, sensitivity, accuracy on the UPS subset)
  rest on about 40 features and are not simulated here; they decide nothing in KIT-162.
- KIT-162's real-data features listed in the table above (levels, `Norm`, pairwise contrasts,
  heteroscedastic variance, ratio compression, missing values, normalisation route).

## Binding, ceiling, charge (R7)

One `run-tests` scripts request at the commit that adds this declaration and `t2/di10k-plan.tsv`
(and nothing else), named in the request and in `runs.tsv`: `val-tests.R <out> --files metrics
--coverage-plan t2/di10k-plan.tsv --request DK --nsim 2000 --target-sims 20000 --chunks 4`.

No probe. Script `timeout_minutes` 50. Ceiling 4.13 core-hours (62 job-minutes x 4 cores), within
the 4.39 left on line DI10K. Expected 1.1 to 2.2 process-hours (the cell's own timing in run
45403024967, and DI10b's times scaled by 4,840 / 500), so 0.28 to 0.55 h per chunk, plus set-up:
1.9 to 3.0 core-hours.

One run. The result is a result: recorded whatever it shows, never re-run. A timeout with no result
falls under VAL-54's instrument rule, with the limit raised as the mechanics fix. Any other failure
in the mechanics with no result is VAL-54 too: at most three further runs, each recorded, and no
threshold, rule or cell changes.
