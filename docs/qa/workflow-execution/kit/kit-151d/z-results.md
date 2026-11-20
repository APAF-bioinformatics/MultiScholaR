# KIT-151 d, Q2 results: per-feature z-score and scale_views, scored once

Declaration `z-declaration.md` (5940475a, pushed alone). Lifecycle 44785248037 (fc5e8042). Counted run 44786925438
(35e50294): 30 jobs (2 N x 15 blocks of 2 seeds; the declaration's "60 jobs" is an arithmetic slip in its text, the
job list fixed in the lane is what ran), 1,110 of 1,110 cells ok, 0 failed, 0 at the 10,000-round limit (most 2,366).
Job wall time 12,295 s: **13.7 core-hours** (ceiling 25). Processors: AMD EPYC 7763, 9V45, 9V74; Intel Xeon 6973P-C,
8370C, 8573C. Tables: `z-rows.tsv`, `z-summary.tsv`, `z-pairs.tsv`, `z-scale-invariance.tsv`, `z-tie-k151b.tsv`;
scorer output `z-score.md` (unchanged). No correction pass.

Checks: with `scale_views` on, all four scale patterns gave exactly the same |r| (difference 0 in every seed at both
N), as predicted. Tie to k151b: 156 of 180 |r| equal within 1e-6, 42 of 42 on the same CPU model; the 8 seed-cells that
differ all ran here on the Intel Xeon 6973P-C (unpinned evidence runs; see Q3). Both arms of a paired comparison
always ran in the same job, so the pairing holds.

Selection fixed at the decided default throughout (observed in half the samples, 1,000 by variance minus the trend).
F1 and F3 are unaffected by every arm (all differences within +/-0.001); everything below is F2 (the weak factor of
the two small data sets, D3 and D4), which is also the worst of three in nearly every seed.

## z-score (mean paired difference in F2 |r| [95% interval], seeds better/worse)

| comparison | N = 20 | N = 60 |
|---|---|---|
| Z1 all vs none, missing on | **+0.301 [+0.200, +0.406]** better (26/4) | -0.000 no practical difference |
| Z2 all vs none, missing off | +0.057 [-0.002, +0.131] inconclusive | -0.001 no practical difference |
| Z3 interaction (Z1 minus Z2) | **+0.244 [+0.125, +0.359]** better (25/5) | +0.001 no practical difference |
| Z4 small data sets only vs none, missing on | **+0.317 [+0.214, +0.422]** better (29/1) | -0.000 no practical difference |
| Z5 large data sets only vs none, missing on | +0.044 [-0.023, +0.123] inconclusive | (not run) |
| Z6 all vs small only, missing on | -0.016 [-0.050, +0.004] inconclusive | +0.000 no practical difference |
| Z7 small only vs none, missing off | +0.061 [+0.002, +0.135] better (20/10) | -0.001 no practical difference |
| Z9 all vs none, `scale_views` off, missing on | +0.240 [+0.143, +0.342] better | (not run) |
| Z10 z all: `scale_views` on vs off | -0.000 no practical difference | (not run) |

Medians of F2 at N = 20, missing on: no z-score 0.569 (12 of 30 seeds >= 0.9); z small 0.967 (28); z all 0.960 (26);
z large 0.821 (15). Missing off: 0.965 to 0.973 in every arm. At N = 60 F2 is found in 30 of 30 seeds in every arm.

Reading: the gain is there only with missing values (Z3), comes from z-scoring the small data sets that carry the weak
factor (Z4 about equal to Z1; Z5 not distinguishable from 0), does not depend on `scale_views` (Z9, Z10), and vanishes
at N = 60, where the factor is found without it. It carries over from the raw cap of k151b to the decided selection
(k151b +0.347; here +0.301). Not tested: why equal per-feature weight in a small, partly missing data set helps (the
features there are kept whole, so their noisy, left-censored low-abundance features keep their full weight without
it); that is an inference.

## scale_views (`scale_views` off -> on, no z-score; F2 |r|)

| scales (D1, D2, D3, D4) | N = 20 missing on | N = 20 missing off | N = 60 missing on |
|---|---|---|---|
| none (1, 1, 1, 1) | -0.014 [-0.084, +0.046] inconclusive | **+0.107 [+0.026, +0.204]** better | no practical difference |
| orig (1, 2, 0.25, 4), k151b | -0.062 [-0.145, +0.010] inconclusive | -0.038 [-0.101, +0.006] inconclusive | no practical difference |
| lup (4, 2, 0.5, 0.25) | +0.031 [-0.079, +0.143] inconclusive | +0.084 [-0.013, +0.195] inconclusive | no practical difference |
| ldown (0.25, 0.5, 2, 4) | +0.028 [-0.083, +0.140] inconclusive | -0.017 [-0.079, +0.045] inconclusive | no practical difference |

Interaction (effect at the pattern minus effect at none), missing off: orig **-0.145 [-0.249, -0.054]**, ldown
**-0.124 [-0.228, -0.033]**, lup -0.022 inconclusive; missing on: all inconclusive (-0.048, +0.045, +0.042).

Reading: with `scale_views` on the result does not depend on how the data sets happen to be scaled (identical models
across the four patterns). With it off, the result depends on the scales: where the small data sets are scaled up
(orig's D4 x 4, ldown), switching it off is relatively better; where nothing is scaled, switching it on is better
(missing off) or not distinguishable (missing on). This fits the predicted cause of the k151b loss (undoing an arbitrary
scale that happened to favour the weak factor's data sets) for missing off; with missing values the intervals are too
wide to tell. No pattern showed `scale_views` on worse than off with an interval excluding 0. At N = 60 there is no
difference anywhere.

## What is shown and what is not

Shown (one generator, one training seed, 30 paired seeds): at N = 20, a per-feature z-score on the small data sets
(or all) recovers a weak factor confined to small data sets much better when values are missing in an
abundance-dependent way, and makes no practical difference otherwise or at N = 60; `scale_views` makes the model
invariant to data-set scale, and the earlier loss is explained by the arbitrary scales favouring the weak factor when
it was off. Not shown: real data, other missing-value mechanisms, other factor layouts (for example a weak factor in a
large data set), the z-score's effect on variance-explained figures (z-scoring changes what "variance explained"
refers to), or why the z-score helps under missing values.

Methods sentences the results support: "Each data set was divided by its overall standard deviation (`scale_views`),
which makes the model independent of the arbitrary scale of each data set." And, if the user adopts it: "In
simulations with 20 samples and abundance-dependent missing values, scaling each selected feature to unit variance
improved recovery of a weak factor confined to the smaller data sets (mean |r| +0.30, 95% interval +0.20 to +0.41,
30 data sets); without missing values or with 60 samples it made no practical difference."
