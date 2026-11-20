# KIT-151 SCI-01 / SCI-18: re-score of the existing rows (no new compute)

Written 2026-11-20 on branch `run-mofa-evidence-k151b`. Input: the committed `followup-rows.tsv` and
`followup-composition.tsv` of run 44654405980 (byte-identical to the run's artifacts, sha256 659f2069... and
a732196d...). Script `.github/mofa-smoke/k151-rescore.py` (sha256 7189272b...; Python standard library; 20,000
seed-bootstrap resamples, `random.Random(151)`). Outputs: `sci01-rescore-arms.tsv` (per scenario, mode, arm, factor:
median, mean, 10th percentile, worst seed, iterations), `sci01-rescore-pairs.tsv` (every pair of arms, every factor:
mean paired difference with 95% bootstrap interval, median difference, seeds better / worse, exact two-sided sign
test, worst seed difference), `sci01-rescore-composition.tsv`. This is a re-reading, not a new pre-declared test; the
intervals are unadjusted and 10 seeds is few.

## What the pre-declared rules returned

- First run (44650673965; N = 20, convergence fast, maxiter 1,000, F2 weight 1, i.e. per-feature SNR 0.44): **no cap
  met the rule** in the main scenario (closest cap 1,000: F2 2/10, F3 9/10). At N = 60 cap 1,000 met it; with 15%
  missing values no cap did.
- Follow-up (44654405980): F2 was strengthened to weight 1.5 (SNR 1) in a new scenario, and three convergence modes
  were added. Rule (c) returned cap 1,000 under slow, decided by the N = 20 strong scenario alone (N = 60 passes
  every cap; N = 20 weak passes none). Cap 1,000's ninth F2 success there is |r| = 0.902; at 0.899 the rule would have
  returned 250.
- At the original strength (n20weak, slow) no cap recovers F2 (at most 1 of 10 seeds for any cap); medians fall
  steadily from 0.777 (cap 250) to 0.673 (cap 5,000).

## N = 20 strong F2, slow: matched |r| median / 10th percentile / worst seed

| arm | F2 | F3 | worst of F2, F3 per seed | iterations median; at the 1,000 limit |
|---|---|---|---|---|
| no selection, scale_views off | 0.352 / 0.085 / 0.081 | 0.951 / 0.595 / 0.583 | 0.352 / 0.085 / 0.081 | 391; 0 |
| no selection, scale_views on | 0.297 / 0.075 / 0.053 | 0.967 / 0.595 / 0.586 | 0.297 / 0.075 / 0.053 | 548; 0 |
| cap 250 | 0.957 / 0.917 / 0.900 | 0.989 / 0.979 / 0.976 | 0.957 / 0.917 / 0.900 | 966; 5 |
| cap 500 | 0.957 / 0.882 / 0.871 | 0.991 / 0.964 / 0.954 | 0.956 / 0.882 / 0.871 | 846; 2 |
| cap 1,000 | 0.947 / 0.897 / 0.854 | 0.994 / 0.963 / 0.940 | 0.943 / 0.897 / 0.854 | 834; 2 |
| cap 2,000 | 0.939 / 0.791 / 0.232 | 0.993 / 0.939 / 0.754 | 0.919 / 0.701 / 0.232 | 938; 3 |
| cap 5,000 | 0.923 / 0.710 / 0.114 | 0.989 / 0.587 / 0.559 | 0.863 / 0.515 / 0.114 | 1001; 8 |

## Paired differences, slow, worst of F2 and F3 per seed (later arm minus earlier arm)

Mean difference [95% seed-bootstrap interval]; seeds better / worse; exact sign-test p.

| comparison | N = 20 strong | N = 20 weak |
|---|---|---|
| 500 vs 250 | -0.013 [-0.032, +0.002]; 3/7; p 0.34 | -0.043 [-0.127, +0.031]; 3/7; p 0.34 |
| 1,000 vs 250 | -0.019 [-0.036, -0.003]; 4/6; p 0.75 | -0.102 [-0.208, +0.006]; 3/7; p 0.34 |
| 1,000 vs 500 | -0.006 [-0.023, +0.010]; 4/6; p 0.75 | -0.059 [-0.164, +0.032]; 5/5; p 1.00 |
| 2,000 vs 1,000 | -0.093 [-0.220, -0.009]; 1/9; p 0.02 | -0.018 [-0.102, +0.082]; 4/6; p 0.75 |
| 5,000 vs 1,000 | -0.178 [-0.335, -0.052]; 1/9; p 0.02 | -0.074 [-0.184, +0.044]; 3/7; p 0.34 |
| 1,000 vs no selection (scale_views on) | +0.567 [+0.414, +0.708]; 10/0; p 0.002 | +0.351 [+0.249, +0.453]; 10/0; p 0.002 |
| 250 vs no selection (scale_views on) | +0.586 [+0.421, +0.731]; 10/0; p 0.002 | +0.452 [+0.342, +0.562]; 10/0; p 0.002 |

Per factor (strong, slow): F2 1,000 vs 250 -0.017 [-0.036, -0.001], 5/5 seeds; F3 1,000 vs 250 -0.002 [-0.011,
+0.004], 8/2 seeds in 1,000's favour but by at most a few thousandths.

## Signal-carrying features kept (median over 10 seeds; identical for weak and strong, which share D1 and D2)

| cap | D1 loaded features kept (of 3,805) | share of D1's kept that are loaded | D1 F3-loaded kept | D2 loaded kept | share of D2's kept that are loaded |
|---|---|---|---|---|---|
| 250 | 6.0% | 91% | 7.0% | 46% | 37% |
| 500 | 11% | 86% | 12.5% | 64% | 26% |
| 1,000 | 20% | 77% | 22% | 84% | 17% |
| 2,000 | 34% | 65% | 36% | 100% | 10% |
| 5,000 | 60% | 46% | 61% | 100% | 10% |

## Reading

- Selecting versus not selecting: a large paired gain at N = 20 in this simulator (every seed better, both strengths).
  The simulator makes the most variable features the loaded ones, so this gain mixes balance and signal enrichment.
- 250, 500 and 1,000 are not separated. The only interval excluding 0 among them (strong, 1,000 vs 250, -0.019) is
  small, its sign test is 4/6, and at weak strength the intervals span 0. If anything the smaller caps do slightly
  better on F2, the opposite of what the "largest passing cap" rule chose.
- 2,000 and 5,000 are worse than 1,000 in the strong scenario, driven by a few collapsed seeds (worst seed 0.23 and
  0.11); medians differ little.
- Recovering F3 needed only a small part of D1's F3 features (7% at cap 250), because those features were the most
  variable ones. That will not hold when unrelated features are more variable (the follow-up's stress arm).
