# P1K: the confirmatory rank cell in the form of KIT-162's primary endpoint, declaration written before the request

Status: FINAL, committed alone with `t2/p1k-plan.tsv` at the commit the request binds. The
audited draft `t2/p1k-declaration-draft.md` stays as it is, with repairs R-P1 to R-P4 of the
harness auditor's P1K check (`AUDIT-KIT-159-impl-recheck-20261124.md`, "P1K: change and draft
declaration") applied in both. The request is not sent until the auditor confirms this file.

Owner's decision: methods ledger VAL-176. KIT-162's primary endpoint gets its own confirmatory
coverage cell: one contrast, 40 positives among about 4,800, the effect set so that the partial AUC
is not at its ceiling, a fixed one-directional feature draw, declared before it is sent. The
coordinator's ruling of 2026-11-25 sets the form: C4e five-block layout, delete-one-block, g = 5;
budget line P1K, 6 core-hours (KIT-159 r1.26). Why it cannot be folded into C4e-170, and how
C4e-148 differs, is in re-audit 4 (2), `AUDIT-KIT-159-impl-recheck-20261124.md`.

## What it measures

The pinned C4e g = 5 rank rows license P1. They come from cell C4e-148 under `redraw`, at effect
scale 0.35, with the four contrasts against group A pooled (160 positive entries over 40 features).
At 97.5 % they covered 0.9935, and at 98.33 % 0.996 (T2R, run 45227824000).

P1K runs **cell C4e-171, family `c4e:4840p1`, under generator
`fixed_onedir`**. The cell's statistic is `pauc_diff` on **one contrast** (group C against A): 40
positives and 4,800 nulls. It is the normalised pAUC to FPR 0.01 of the moderated (`plain`) arm
minus that of the ordinary-t arm, with positives ranked by one-sided p in the expected direction
and nulls by two-sided p (`valPartialAuc`, the harness's M7). The interval is the delete-one-block
jackknife at g = 5.

| Item | P1K (C4e-171) | C4e-148 (pinned rank rows) | KIT-162's P1 (`tickets/kit/KIT-162.md` 6.1, 6.7) | Same as P1? |
|---|---|---|---|---|
| Units | 5 blocks, delete-one-block (form (iii)), g = 5 | same | 5 mixtures, leave one mixture out | yes |
| Samples per unit | 10 (5 groups x 2) | same | 10 channels (4 levels x 2 + `Norm` x 2) | count yes; structure no |
| Features | 4,840 | 4,840 | about 4,800 | yes |
| Positives | 40 (1 contrast) | 160 entries (40 features x 4 contrasts) | about 40 UPS keys, one contrast | **yes** (C4e-148: no) |
| Nulls | 4,800; 48 in the FPR-0.01 region | 19,200; 192 in the region | about 4,760 HeLa keys; about 48 in the region | yes |
| Direction | all up (fixed draw, `fixed_onedir`) | all up (`redraw`) | one direction in the contrast | yes |
| Feature draw | fixed (baseline seed 2,100,000,171) | redrawn per simulation | one fixed proteome | yes (C4e-148: no) |
| Contrast and effect | group C against A, true effect 0.40 log2 (scale 0.40 x step 1) for all 40 | 4 contrasts at 0.175, 0.35, 0.525, 0.70 | `L1 - L05`, 1 log2 nominal, compressed by TMT ratio compression; heterogeneous per protein | **no**: one effect size for every positive; size set by the pAUC rule, not by KIT-162's fold change |
| Plain-arm pAUC (calibration mean) | 0.593 (rule below) | 0.602 | not known before data | off the ceiling, by rule |
| Paired difference | moderated minus ordinary t on the same data | same | arm minus comparator (for example against `baseline-plex`) on the same runs | by kind (a paired pAUC difference on one data set); the two rankings differ in what makes them differ |
| Variance, normalisation, model | one residual SD 0.3; block-by-feature SD 0.1; median normalisation; limma with group and block terms | same | heteroscedastic, ratio compression, missing values (union view: an untested UPS key is a miss), the package's route | **no** |
| Decision levels | 8 levels incl. 97.5 % and 98.33 % | same | 97.5 % (primary), 98.33 % (RUV family) | yes |

**What the result will say.** It says whether the delete-one-block jackknife interval covers a
paired pAUC difference computed on one contrast with 40 positives and a fixed feature draw, at
KIT-162's feature count, unit count and decision levels, under the harness's homoscedastic
generator. It says nothing about KIT-162's real-data features listed in the table, nor about
heterogeneous effects across the 40 positives.

**Code path.** `msrval/R/coverage.R` at 3939ffe2 (the cell), with the tests corrected at 70648f6d
(test file only) and the one-sided misses recorded at a303bf3f (R-P1). The fit and generator are C4e-148's and DI10K's (`valCovMatrixData`,
`valCovFixedFeatures`, `valCovFit`). The statistic is `valCovStats`'s rank branch restricted to
contrast column 2 (`contrast = 2L`). The effect is `valRankEffect("block", 5L, contrast = 2L)`.

## The effect rule, stated before any coverage run

`valRankEffect`'s existing rule (KIT-159 r1.13, VAL-69): the generator's group effects are scaled
by the grid value (0.05 to 1 by 0.05) whose mean normalised pAUC at FPR 0.01 of the plain arm lies
in [0.4, 0.8] and is nearest 0.6. The mean is taken over 50 seeded calibration simulations
(`redraw`, seeds 1,900,000,000 + 50,000 + 5,000 + s), the grid is walked upwards, and the walk
stops once the mean passes 0.8.

P1K applies the **same rule, unchanged, to the pAUC of the one contrast** (`contrast = 2L`), on the
same calibration data sets. Nothing in it reads a coverage result.

- **Recorded before any coverage run** (run 45468618081, `test-p1k.R`, console
  `[p1k-calibration]`): one-contrast effect **0.40**, mean plain-arm pAUC **0.5932**. The same run
  reproduced the pooled rule's value: 0.35, mean 0.6017, the effect of the pinned C4e g = 5 rank
  rows.
- **Predicted before that run** (normal theory, contrast SE 0.3 x sqrt(2 / 10) = 0.134; the
  plain-arm pAUC of one contrast at scale s is the mean over FPR f in (0, 0.01) of
  Phi(s / 0.134 - z_{1-f})): 0.343 at 0.30, 0.481 at 0.35, **0.621 at 0.40**, 0.746 at 0.45,
  0.845 at 0.50. So the rule selects 0.40, and the walk stops after 0.50. The same formula for the
  pooled contrasts gives 0.613 at 0.35, the pinned rows' value.
- **Why calibrate under `redraw`** when the run uses the fixed draw: the feature baselines cancel in
  every contrast. They enter only through median normalisation, whose shift at 40 of 4,840 is
  about 8 % of a contrast SE (DI10K). So the draw does not move the pAUC materially. Using the
  existing calibration seeds keeps one rule for every rank cell.
- **Why one contrast at scale 1 cannot be used:** at effect 1 log2 the one-contrast pAUC is 1.000
  (t about 7.5). The statistic would sit at its ceiling and the difference would be near 0 in
  every simulation (re-audit 4 (2)).

## Configuration checks (the cell is the declared one-contrast cell), fixed before the run

- **K1 (the effect):** every coverage row of C4e-171 carries `effect` 0.40, the value recorded
  above. It shows that the run used the rule's value, not a ceiling or a pooled value.
- **K2 (the statistic's size): the target of `pauc_diff` lies in [0.025, 0.075].**
  - Derivation: a model of the two rankings (normal theory; the ordinary arm as t on 41 residual
    df, the moderated arm as known variance; 40 positives at 0.40 / 0.134, 4,800 nulls; positives
    ranked by one-sided p, nulls by two-sided p). A 300-replicate Monte Carlo of the empirical pAUC
    gives an expected 0.044 (MC SE 0.002) and a plain-arm pAUC of 0.627. The population value is
    0.049 by numerical integration (harness auditor).
  - The same model applied to C4e-148's pooled form gives 0.0213, against the measured target
    0.0243 (T2R). Scaling by that ratio, 1.14, gives about 0.050 for P1K.
  - The band is 0.050 ± 0.025, which allows for the model's error on the pooled cell (14 %) with
    margin.
  - What it shows: the paired difference is off the ceiling and of the modelled size; degeneracy is
    read from the zero-width share.
  - What it cannot show: one contrast against several pooled at a similar scale (pooled at 0.35:
    0.024, at the band's edge). That is shown by the tests instead: 40 positives, 4,800 nulls,
    contrast column 2.

## Pass rule and readings (band as T2, D, DI10b, DI10c and DI10K ran it)

A coverage value fails at a level when it is below level - 0.02. The rule is applied to the one
statistic at all 8 levels. A `conservative` label is reported beside it as observed; a
`degenerate` label is a reading of its own (below).

- **Configuration checks not met** (K1 or K2 fails): **the run does not test the question.** It is
  recorded, reported and not re-run.
- **Pass** (K1 and K2 hold, no coverage below level - 0.02 at any of the 8 levels, and the row is
  not `degenerate`) (R-P4): **the C4e g = 5 delete-one-block jackknife covers P1's form (a paired
  pAUC difference to FPR 0.01 on one contrast, 40 positives of 4,840, one direction, fixed draw)
  at the decision levels, under the listed differences.** KIT-162 quotes this cell's coverage and
  one-sided misses at 97.5 % and 98.33 %, not the pooled rows' coverage or label. It does not
  license heterogeneous or compressed fold changes, untested UPS keys counted as misses, the
  package's normalisation route, KIT-162's arm-versus-baseline rankings, or a P1 near its ceiling
  in real data.
- **Degenerate** (R-P2; more than 1 % zero-width intervals, or the target within one Monte Carlo SE
  of ±1): not offered by the harness's rule (r1.13); KIT-162 needs a ticket coverage regime for
  P1.
- **Fail** (K1 and K2 hold, at least one coverage below level - 0.02): **KIT-162 needs a ticket
  coverage regime for P1 at its own composition before it runs.**

## Reported beside the rule (for the owner)

- Coverage at each of the 8 levels, beside the pinned C4e g = 5 rank rows (C4e-148: 0.99, 0.9935,
  0.996, 0.997, 0.998, 0.9985, 0.9995, 0.9995). Misses per level, of 2,000.
- At every level: the exact (Clopper-Pearson) 95 % interval of the coverage. Where misses exceed
  nominal: P(X >= misses) under the nominal miss rate (binomial, n = 2,000). A level below nominal
  at p < 0.01 is reported as anti-conservative with its measured coverage, whatever the band says
  (as VAL-163 (1) did for `sens_diff`).
- (R-P1) At every level, the misses at the lower and at the upper bound separately, each with its
  exact 95 % interval; KIT-162 computes its stated false rates from these (KIT-159 C66; KIT-162
  6.7). They are recorded in the detail as `miss_low` and `miss_high` (counts of 2,000); the
  regime table does not carry them.
- `se_inflation` and `half_width_ratio`; the share of zero-width intervals; the target against
  K2's expectation.
- **The number of positives and the statistic's discreteness.** 40 positives per simulation and
  contrast, and 4,800 nulls, of which 48 fall in the FPR-0.01 region. The true positive rate moves
  in steps of 1/40 = 0.025, and the false positive rate in steps of 1/4,800. In normalised pAUC
  units, one null ranked above every positive costs 1 / (0.01 x 4,800) = 0.021, and one positive
  ranked below the 48th null costs up to 0.025; smaller rearrangements cost proportionally less.
  KIT-162's margin of 0.125 is about five such steps.

## What it cannot show

- Heterogeneous fold changes across the 40 positives (TMT ratio compression, per-protein
  variance), missing values (P1's union view counts an untested UPS key as a miss), and KIT-162's
  normalisation route and arms. The cell's paired difference is moderated against ordinary t,
  which is not KIT-162's arm against a comparator; the match is by kind.
- P1's pairwise-intersection view and the secondaries (S1 to S3, S11, S13), which decide nothing
  in KIT-162.

## Binding, ceiling, charge

One `run-tests` scripts request at the commit that adds the final declaration and `t2/p1k-plan.tsv`
(request PK; C4e-171 in four ranges of 500 simulations and 5,000 target simulations; `generator`
`fixed_onedir`). The command is `val-tests.R <out> --files metrics --coverage-plan t2/p1k-plan.tsv
--request PK --nsim 2000 --target-sims 20000 --chunks 4`.

**Cost**, from measured times:

| Anchor | Source | Process time (2,000 simulations, 20,000 targets) |
|---|---|---|
| The same fits without pAUC | DI10K (run 45407487281) | 2,974 s |
| The same fits with pooled pAUC (19,360 entries) | C4e-148 (T2R) | 10,525 s |
| This cell's pAUC on one contrast (4,840 entries) | pooled pAUC cost (7,551 s) divided by 4 | about 1,890 s more than DI10K |
| Calibration in each chunk process (once) | `test-p1k.R` | about 50 s |

- **Likely:** 4,860 process-seconds, so 1,215 s per chunk plus 50 s, about 21 minutes. The
  in-test timing of one simulation (0.78 s for seven fits and the pAUCs) points lower, near
  DI10K's 13 minutes per chunk.
- **Upper:** the pAUC cost not reduced by the smaller table: 2,680 s per chunk, about 45 minutes.

Script `timeout_minutes` 50; job about 62 minutes with set-up; ceiling 4.13 core-hours.

Coordinator's ruling (ticket r1.27): line P1K is 9 core-hours, and the coverage run takes a
50-minute script limit (ceiling 4.13 core-hours), about 2.4 times the likely per-chunk time and
above the upper anchor.

One run. The result is a result: recorded whatever it shows, never re-run. A timeout with no result
falls under VAL-54's instrument rule (the limit raised as the mechanics fix), and so does any other
failure in the mechanics with no result produced: at most three further runs, each recorded, and no
threshold, rule or cell changes.
