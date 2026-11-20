# KIT-151 d, Q2: per-feature z-score and scale_views. Declaration

Fixed before the counted run. Authority: methods ledger INT-13 (user, 2026-11-21: "go with all recommendations").
Branch `run-mofa-evidence-k151d`; lane `.github/workflows/k151d-evidence.yml`, script `.github/mofa-smoke/k151d.R`
(task `z`), scorer `.github/mofa-smoke/k151d-score.R` (task `z`). Script hashes: end of this file. Lifecycle:
run 44785248037 (commit fc5e8042): all three tasks, 8 jobs; every deliberately failing cell was written as failed and its job marked failed, every other cell trained to the end (0 at the round limit), the score job ran on all of them. One scorer wording defect was fixed after it (35c41ca3, pin history line); the evidence script and the lane are unchanged. Q2 lifecycle (seeds 91, 92 at N = 20, 91 at N = 60; every cell): with scale_views on the scale patterns gave the same iterations and abs r (as predicted); about 3 s per training at N = 20 and 16 s at N = 60. Evidence only: no package code, default or wording changes; the user decides.

## Purpose

In k151b (run 44677751982, `selection-b-results.md`) two effects were measured and not explained:

1. A per-feature z-score after selection gave a large gain under abundance-dependent missing values (raw cap 1,000,
   scenario `full`, worst of three planted factors: +0.347 [+0.250, +0.444]; carried by F2, the weak factor of the two
   small data sets D3 and D4), and none to speak of without them (`trend`: +0.040, inconclusive; `prev`: +0.024).
2. `scale_views` made F2 worse when the data sets are on different scales (none vs scale_views, no selection: trend
   -0.149 [-0.232, -0.077], full -0.086 [-0.166, -0.017]).

This run separates candidate causes, with the selection fixed at the decided default (INT-12).

## Candidate causes and what separates them

- z-score gain: (a) the missing values themselves (the gain appears with them and not without): arms with and without
  missing values; (b) equal weight for the small data sets' features (D3, D4 kept whole, holding F2): z-score on the
  small data sets only; (c) equal weight for the large data sets' selected features (D1, D2): z-score on the large
  ones only; (d) an interaction with `scale_views`: z-score with `scale_views` off.
- scale_views loss: (a) `scale_views` divides each data set by one number, so it undoes whatever scale a data set
  happens to have; in k151b D4 (F2) was multiplied by 4, so with `scale_views` off D4 weighed more and F2 was easier.
  If that is the cause, the sign of the effect follows the scales: off is better when the small data sets are scaled
  up (`ldown`), on is better when they are scaled down (`lup`), and there is little or no effect without scale
  differences (`none`). (b) a loss with no scale differences at all would point at the size balance instead.
  With `scale_views` on, the result should not depend on the scales at all (each data set divided by its own spread;
  the trend selection is unchanged by a scale factor): reported as a check.

## Data (k151b generator, unchanged)

k151b.R's generator (`generate()`; copied as `generate_b` with N and the scales as parameters; N = 20 with scales
`orig` is k151b's arithmetic exactly): D1 20,000 features (F1, F3), D2 2,000 (F1), D3 200 and D4 100 (F1, F2); 10%
loaded per factor and data set; weights +/-1, F2 +/-1.5; offset ~ N(20, 2); noise sd 1.5 x 2^(-(offset - 20) / 4)
(mean-variance trend, always on); missing values (when on): P(missing) = 1 - (1 - plogis((17.5 - x) / 0.5)) x 0.98,
own stream (data seed + 400,000), on the unscaled values (about 19% per data set); scales applied last.

Scale patterns (D1, D2, D3, D4): `none` (1, 1, 1, 1); `orig` (1, 2, 0.25, 4), k151b's; `lup` (4, 2, 0.5, 0.25):
large data sets scaled up, small ones down; `ldown` (0.25, 0.5, 2, 4): large data sets scaled down, small ones up.

Data seeds 151001 to 151030 (30 paired seeds; at N = 20 the same data as k151b) at N = 20 and N = 60 (N = 60 draws
differ; two groups of N / 2).

## Preparation and training (fixed)

Selection as INT-12 decided: features observed in at least half the samples; then (z-score data sets only) the
zero-spread filter; then at most 1,000 features per data set by variance minus the fitted mean-variance trend
(`lowess(mean, variance, f = 0.3)` per data set, the difference form, k151b's arm `trend1000`; the ledger row says
"variance above the fitted mean-variance trend", which this form is); D3 and D4 are kept whole; then the z-score
(centre and divide by the sd of each feature's observed values) on the data sets named by the arm. MOFA2 1.22.0 with
5 factors, double precision, `center_groups` TRUE, convergence mode slow, round limit 10,000, training seed 20261119,
`weight_views` FALSE, `scale_views` as the arm says.

## Cells (code `m<missing>_<scales>_z<which>_sv<scale_views>`)

N = 20, 24 cells per seed (720 trainings):

- z-score grid, scales `orig`, missing on (m1) and off (m0): z none / all / small (D3, D4) / large (D1, D2) with
  `scale_views` on; z none / all with `scale_views` off.
- scale_views grid, no z-score, missing on and off: scales none / orig / lup / ldown, `scale_views` on and off.

N = 60, 13 cells per seed (390 trainings): z none / all / small with `scale_views` on, missing on and off (scales
`orig`); and with missing on, no z-score: scales none / orig / lup / ldown with `scale_views` on and off.

## Measures

Per seed and cell: |r| of each planted factor (F1 group, F2 small data sets, F3 D1) with its greedily matched learned
factor; the worst of the three (primary); variance explained by the matched factors (of the selected features);
features kept per data set; iterations and whether training stopped by the rule or at the round limit.

## Analysis (paired, continuous; no single-seed thresholds)

For each declared comparison: per-seed difference (b minus a), mean with a 95% seed-bootstrap interval (20,000
resamples, `set.seed(151)`), median, seeds better / worse, worst seed. Reading: interval inside [-0.02, +0.02]
"no practical difference"; lower bound > 0 "better"; upper bound < 0 "worse"; otherwise "inconclusive". For the worst
of three, F2, F3 and F1. Unadjusted, descriptive (about 26 comparisons x 4 measures x 2 N).

Declared comparisons (a -> b), at N = 20 (and at N = 60 where both cells exist):

- Z1 missing on: z none -> z all (`scale_views` on); Z2 the same with missing off; Z3 interaction: Z1 minus Z2 per seed.
- Z4 / Z5 missing on: z none -> z small; z none -> z large; Z6 z small -> z all. Z7 / Z8: Z4 / Z5 with missing off.
- Z9 missing on, `scale_views` off: z none -> z all; Z11 the same with missing off; Z10 z all: `scale_views` off -> on.
- S, for each scale pattern (none, orig, lup, ldown) and missing on / off: `scale_views` off -> on (no z-score).
- SX interactions: the `scale_views` effect at orig, lup, ldown minus its effect at none (per seed), missing on / off.

Checks (not gates): with `scale_views` on, the |r| at each scale pattern against `none` (expected near 0); the tie
to k151b (`m1_orig_znone_sv1` at N = 20 is k151b's `full:trend1000_obs50`, `m0_orig_znone_sv1` is
`trend:trend1000`, on the same data, selection and options: expected equal to rounding on the same CPU model).

## What the possible outcomes would mean

- Z1 better and Z2 no practical difference (Z3 better): the gain depends on the missing values. Z4 better and Z5 not:
  it comes from equal weight for the small data sets' features; Z5 better and Z4 not: from the large ones; both:
  from the z-score as such.
- Z1 no practical difference or worse: the k151b gain does not carry over to the decided selection (variance above
  the trend, observed in half the samples); it was specific to the raw-variance cap.
- Z9 better while Z10 is no practical difference: the gain is not an interaction with `scale_views`.
- S: off better at ldown and orig, on better at lup, little at none (SX away from 0 in the predicted directions): the
  k151b loss of `scale_views` is the undoing of arbitrary data-set scales, not a defect of the option; the Methods can
  then say what `scale_views` does and that the result does not depend on how the data sets happen to be scaled.
  A loss at none would mean `scale_views` itself costs the weak factor here.
- No reading picks a default; the user decides.

## Budget

60 jobs (N x 15 blocks of 2 seeds), each one install plus 48 (N = 20) or 26 (N = 60) trainings; job timeout 150
minutes. Ceiling 25 core-hours (job wall time x 4 cores); expected from the lifecycle, about 250 s per N = 20 job and 510 s per N = 60 job, about 13 core-hours. One counted run; one correction
pass at most, for an instrumentation failure only, rerunning only the affected jobs. A failed cell drops that seed
from the comparisons it is in; the number is reported.

## Limits stated in advance

One generator (Gaussian, linear; one trend shape; one missing-value mechanism; four scale patterns); one training
seed; N = 20 and 60 only; the k151b mechanism inference for the z-score was not predicted before k151b and is tested
here, not assumed.

## Script hashes (sha256, fixed for the counted run)

```
8fe56e9ad458fce63443eafe532e9cd47439ca4de73f020f9ae1885024eb6f7e  .github/mofa-smoke/k151d.R
6261544a4ab84361fcfcb6fc7dbc3b0ee19040a61db5520abc3c47cbbc73ba4e  .github/mofa-smoke/k151d-score.R
04c250d48e264f120042f4e3030388fc308480e62c0d9f2ad13a79c577b0eed4  .github/mofa-smoke/install-sva-library.R
405d84a1b78f7cdf2324ef7d85dc309539e88a1264a1ad2a29ca448c0bea2d17  .github/mofa-smoke/install-ruv-library.R
f4caa249375b225ce1d1cf012b6f3c52d4518ee554e465dfc68159274a1051bd  .github/workflows/k151d-evidence.yml
```
