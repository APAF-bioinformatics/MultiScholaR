# KIT-151 selection evidence, run b: results

Declaration `selection-b-declaration.md` (b0d4e9bf, pushed alone before the counted run). Lifecycle 44676343317
(24babf60; seeds 151091 and 151092, not counted). Counted run 44677751982 (cfab7960): 30 jobs, 1,020 of 1,020 cells,
0 failed, 0 missing; scored once by `k151b-score.py` (sha256 d500aa9f..., printed by all 30 jobs and the score job;
`k151b.R` cc313336..., printed by all 30 jobs; both as declared). No correction pass. Job wall time 12,379 s in total
(longest job 570 s): **about 13.8 core-hours** (ceiling 30; estimate 15). CPU models: AMD EPYC 7763 11 jobs, AMD EPYC
9V74 10, Intel Xeon Platinum 8573C 5, AMD EPYC 9V45 2, Intel Xeon Platinum 8370C 2 (`selection-b-jobs.tsv`).
**No cell reached the 10,000-round limit**; the largest count was 3,486 rounds (full, raw250); medians 300 to 1,650.
Realised missing share in `full`: about 18 to 21% per data set.

Files: `selection-b-rows.tsv` (seed x scenario x arm x planted factor), `selection-b-summary.tsv`,
`selection-b-pairs.tsv` (every declared comparison, every factor and the worst-of-three),
`selection-b-composition.tsv` / `-summary.tsv`, `selection-b-tie.tsv`, `selection-b-jobs.tsv`,
`selection-b-score.json` / `.md` (the score job's output, unchanged), `selection-b-lifecycle-rows.tsv`.

Readings below use the declared rule on the 95% seed-bootstrap interval of the mean paired difference (better /
worse / no practical difference within +/-0.02 / inconclusive); unadjusted, descriptive.

## Matched |r|, median (worst seed); worst-of-three per seed

| arm | prev F2 | prev F3 | trend F2 | trend F3 | full F2 | full F3 | full worst-of-three median / 10th pct |
|---|---|---|---|---|---|---|---|
| none (scale_views off) | 0.405 (0.081) | 0.963 (0.583) | 0.657 (0.112) | 0.989 (0.555) | 0.347 (0.166) | 0.989 (0.565) | 0.347 / 0.211 |
| scale_views, no selection | 0.378 (0.078) | 0.974 (0.586) | 0.322 (0.124) | 0.991 (0.550) | 0.313 (0.104) | 0.988 (0.562) | 0.313 / 0.193 |
| raw 250 | 0.958 (0.887) | 0.988 (0.945) | 0.971 (0.890) | 0.548 (0.205) | 0.931 (0.159) | 0.522 (0.137) | 0.444 / 0.238 |
| raw 1,000 | 0.950 (0.527) | 0.994 (0.940) | 0.969 (0.156) | 0.922 (0.821) | 0.486 (0.177) | 0.971 (0.922) | 0.486 / 0.226 |
| raw 2,000 | 0.939 (0.232) | 0.994 (0.754) | 0.952 (0.184) | 0.972 (0.919) | 0.480 (0.118) | 0.987 (0.936) | 0.480 / 0.201 |
| random 1,000 | 0.953 (0.309) | 0.938 (0.801) | 0.969 (0.240) | 0.967 (0.865) | 0.702 (0.124) | 0.953 (0.773) | 0.702 / 0.221 |
| trend (difference) 1,000 | 0.949 (0.838) | 0.994 (0.946) | 0.965 (0.335) | 0.987 (0.949) | 0.721 (0.066) | 0.986 (0.933) | 0.721 / 0.218 |
| trend (ratio) 1,000 | 0.949 (0.835) | 0.994 (0.939) | 0.962 (0.211) | 0.995 (0.941) | 0.618 (0.148) | 0.994 (0.949) | 0.618 / 0.306 |
| raw 1,000 + z-score | 0.960 (0.852) | 0.994 (0.949) | 0.969 (0.870) | 0.923 (0.681) | 0.959 (0.399) | 0.969 (0.942) | 0.951 / 0.919 |
| raw 250, obs >= 50% | | | | | 0.933 (0.303) | 0.939 (0.851) | 0.884 / 0.486 |
| raw 1,000, obs >= 50% | | | | | 0.916 (0.113) | 0.986 (0.966) | 0.916 / 0.161 |
| raw 1,000 + z, obs >= 50% | | | | | 0.961 (0.312) | 0.985 (0.954) | 0.961 / 0.867 |

F1 (the group factor, in every data set) is recovered by every selection arm except raw 250 in trend / full (median
0.92 to 0.94). Full per-arm tables: `selection-b-score.md`.

## Declared comparisons, mean paired difference [95% interval] (b minus a)

Worst-of-three per seed:

| comparison | prev | trend | full |
|---|---|---|---|
| raw 1,000 vs no selection | +0.486 [+0.397, +0.571] better | +0.418 [+0.299, +0.528] better | +0.204 [+0.082, +0.325] better |
| random 1,000 vs no selection | +0.456 [+0.365, +0.542] better | +0.391 [+0.280, +0.499] better | +0.265 [+0.130, +0.396] better |
| trend 1,000 vs no selection | +0.498 [+0.409, +0.582] better | +0.462 [+0.349, +0.569] better | +0.286 [+0.169, +0.403] better |
| trend ratio 1,000 vs no selection | +0.498 [+0.408, +0.581] better | +0.354 [+0.221, +0.481] better | +0.278 [+0.162, +0.395] better |
| raw 250 vs no selection | +0.509 [+0.418, +0.595] better | +0.080 [-0.042, +0.196] inconcl. | +0.056 [-0.040, +0.149] inconcl. |
| random vs raw 1,000 | -0.030 [-0.079, +0.015] inconcl. | -0.027 [-0.142, +0.089] inconcl. | +0.061 [-0.052, +0.177] inconcl. |
| trend vs raw 1,000 | +0.012 [-0.004, +0.040] inconcl. | +0.044 [-0.035, +0.122] inconcl. | +0.083 [-0.004, +0.179] inconcl. |
| trend vs random 1,000 | +0.042 [+0.013, +0.086] better | +0.071 [-0.007, +0.154] inconcl. | +0.022 [-0.094, +0.135] inconcl. |
| ratio vs difference trend | +0.000 [-0.002, +0.003] no diff. | -0.109 [-0.202, -0.029] worse | -0.008 [-0.110, +0.103] inconcl. |
| raw 250 vs raw 1,000 | +0.023 [+0.005, +0.054] better | -0.338 [-0.423, -0.239] worse | -0.148 [-0.262, -0.035] worse |
| raw 2,000 vs raw 1,000 | -0.061 [-0.134, +0.000] inconcl. | -0.124 [-0.234, -0.025] worse | -0.037 [-0.151, +0.074] inconcl. |
| z-score vs raw 1,000 | +0.024 [+0.008, +0.050] better | +0.040 [-0.014, +0.114] inconcl. | +0.347 [+0.250, +0.444] better |
| scale_views vs none | -0.024 [-0.080, +0.008] inconcl. | -0.149 [-0.232, -0.077] worse | -0.086 [-0.166, -0.017] worse |

By factor, where the worst-of-three hides the mechanism (trend and full):

| comparison | trend F2 | trend F3 | full F2 | full F3 |
|---|---|---|---|---|
| trend vs raw 1,000 | -0.011 [-0.082, +0.064] | +0.074 [+0.061, +0.088] (30/0 seeds) | +0.082 [-0.005, +0.176] | +0.016 [+0.011, +0.022] |
| trend ratio vs raw 1,000 | -0.120 [-0.234, -0.008] | +0.081 [+0.068, +0.095] (30/0) | +0.074 [-0.035, +0.186] | +0.023 [+0.018, +0.029] |
| random vs raw 1,000 | -0.067 [-0.181, +0.046] | +0.048 [+0.035, +0.062] | +0.069 [-0.041, +0.188] | -0.027 [-0.045, -0.012] |
| trend vs random 1,000 | +0.056 [-0.022, +0.141] | +0.026 [+0.018, +0.036] | +0.013 [-0.102, +0.126] | +0.043 [+0.029, +0.060] |
| raw 250 vs raw 1,000 | +0.051 [+0.002, +0.120] | -0.386 [-0.451, -0.324] (0/30) | +0.221 [+0.113, +0.331] | -0.489 [-0.554, -0.423] (0/30) |
| raw 2,000 vs raw 1,000 | -0.175 [-0.279, -0.078] | +0.057 [+0.045, +0.070] | -0.038 [-0.153, +0.075] | +0.015 [+0.009, +0.022] |
| z-score vs raw 1,000 | +0.047 [-0.001, +0.120] | -0.005 [-0.023, +0.007] | +0.350 [+0.252, +0.447] | +0.002 [-0.002, +0.005] |

Eligibility (full only), observed in >= 10 of 20 against >= 2, worst-of-three: raw 250 +0.358 [+0.267, +0.446]
better (F3 +0.450, 30/0 seeds); raw 1,000 +0.097 [+0.000, +0.199] better (F3 +0.019); raw 2,000 +0.007 inconclusive;
random -0.033, trend -0.066, trend ratio -0.034, raw 1,000 + z -0.016, all inconclusive (F3 differences within
+/-0.02 for these).

## What each selection keeps (D1, median over 30 seeds)

| arm | prev: kept that are loaded | trend: kept that are loaded; noise sd of kept | full: kept that are loaded; noise sd; observed fraction |
|---|---|---|---|
| raw 250 | 91% | 26%; 3.37 | 30%; 2.71; 0.35 |
| raw 1,000 | 77% | 27%; 2.83 | 36%; 2.28; 0.50 |
| raw 2,000 | 64% | 28%; 2.56 | 37%; 2.07; 0.60 |
| random 1,000 | 19% | 19%; 1.51 | 19%; 1.50; 0.90 |
| trend (difference) 1,000 | 77% | 42%; 2.35 | 39%; 2.27; 0.50 |
| trend (ratio) 1,000 | 77% | 80%; 1.23 | 57%; 2.01; 0.65 |
| raw 1,000, obs >= 50% | | | 46%; 1.87; 0.70 |
| trend ratio, obs >= 50% | | | 71%; 1.32; 0.90 |

(D1's loaded share is 19%; its median noise sd is 1.5.) With a mean-variance trend, raw variance fills the cap with
the noisiest low-abundance features; at 250 it keeps about 2% of D1's F3 features and F3 is lost in every seed.
Variance above the trend as a ratio keeps the most signal (80% loaded in trend), the difference form much less (42%),
as expected from the noise of the variance estimate growing with the trend. With missing values, the raw cap
also prefers features observed in half the samples or fewer.

## Instrumentation check (not a gate): tie to the follow-up

prev seeds 1 to 10, arms none, scale_views, raw 250 / 1,000 / 2,000, against the follow-up's n20strong slow cells
(150 |r| values). Of the 120 values whose follow-up cell stopped by the rule (before 1,000 rounds), 94 agree within
1e-6 (90 of all 150 identical to 1e-9). The differing converged cells are all in seeds 6 and 10, whose follow-up jobs
ran on an Intel Xeon 6973P-C and whose run-b jobs ran on AMD EPYC 9V45 and 7763: uncapped arms there differ in
iteration count and |r|, up to 0.79 (seed 10, no selection, F2: 0.086 in the follow-up, 0.874 here). Seeds run on
other CPU pairs agree to rounding level. The other differences are cells that hit the 1,000-round limit in the
follow-up and ran on here (for example seed 5 raw 1,000 F2 0.915 at 1,000 rounds, 0.951 at 1,911). The data and
selection are therefore tied to the earlier run; the remaining differences are the round limit and CPU-dependent
divergence of uncapped models near the detection limit (SCI-05).

## Not as declared

Nothing. All cells ran; no correction pass; the scorer's output is reported unchanged.

## Limits

One generator and one choice each of trend shape, missing-value mechanism, scales and eligibility share; N = 20;
one training seed; F2 recovery is bimodal across seeds (recovered or lost), so its intervals are wide; about 30
comparisons per scenario with unadjusted intervals. The z-score's large gain for F2 in `full` was not predicted and
its mechanism (inferred: equal weight for the noisy, partly missing low-abundance features of D3 and D4, which no cap
touches) was not tested.
