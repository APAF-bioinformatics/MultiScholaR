# KIT-151 cap evidence follow-up: declaration (fixed before the counted run)

Authority: user, 2026-11-20 ("run follow up"), relayed by the coordinator. Evidence only; no package code or default
changes. Written 2026-11-20 after the follow-up lifecycle run (44653902387, commit f8f4761b) and before the counted
run; pushed on its own before the counted run's push. The first run (`cap-evidence-declaration.md`,
`cap-evidence-results.md`) stands as read and is not re-scored; both are reported.

## Purpose

The first run could not separate three things: (1) imbalance between data sets, (2) early stopping under convergence
mode "fast" (uncapped arms stopped after a median of 21 to 26 iterations), (3) a small-data-set factor F2 that may be
undetectable at N = 20 whatever the treatment. This run separates them by adding convergence mode, smaller caps and a
stronger F2, and records which features the variance cap keeps (balance and signal enrichment are different
mechanisms).

## Scripts

- `.github/mofa-smoke/k151-followup.R`, sha256 `631959aa90c1248e6a5bfb68a61101f727827ff6e457036f49400def552d7647`
- `.github/mofa-smoke/k151-followup-score.py`, sha256 `20f4abbdc41a7523aa38aff2592741ba16dd19610998a820f5904935e5a97c9b`
- Workflow `.github/workflows/k151-followup.yml`; parameters `.github/mofa-smoke/k151-followup-params.json`.
- Both scripts are byte-identical to the lifecycle's (same sha256 printed in its job). The first run's scripts are
  unchanged. The scorer runs once in the counted run's score job; its output is read once.

## Fixed design

Everything as in the first run unless stated: generator (D1 20,000 features with F1 and F3; D2 2,000 with F1; D3 200
and D4 100 with F1 and F2; 10% loaded features per factor and data set; weight +/-1; noise sd 1.5; per-feature
offset; F1 = group + noise sd 0.5), data seeds 151001 to 151010 (same random draws, so weak data are the first run's
data), training seed 20261119, double precision, 5 factors, maxiter 1000, MOFA2 called directly, view preparation as
`integrationAlignViews`, greedy maximum matching, recovered = |r| >= 0.9.

Added dimensions:

1. Convergence mode: fast, medium, slow. In the installed mofapy2 (core/BayesNet.py:466-479) a run stops when
   100 x |delta ELBO / first ELBO| is below 0.0005 (fast), 0.00005 (medium) or 0.000005 (slow) three times.
   Recorded per cell: iterations, stopped by convergence (iterations < 1000) or by maxiter (iterations = 1000), the
   last finite ELBO.
2. Caps 250 and 500 (with scale_views) beside 1,000, 2,000, 5,000. Cap 10,000 is not repeated.
3. F2 strength. weak = the first run (F2 weight +/-1). strong = F2 weight +/-1.5 in D3 and D4, nothing else changed
   (same draws; only F2's loadings scaled). Principle: the per-feature signal variance of an F2-loaded feature
   (1.5^2 = 2.25) equals the noise variance (1.5^2 = 2.25), a signal-to-noise ratio of 1 per feature, against 0.44
   for weak. Chosen from that principle, not from the first run's results.
4. N = 20 (main) and N = 60 (sensitivity); no missing values.

Arms (each under the modes listed below): none; scale_views; cap250; cap500; cap1000; cap2000; cap5000 (caps with
scale_views); weight_views (with scale_views); zscore (the workbook's: per-feature z-score, center_groups FALSE,
scale_views FALSE).

Inflated-noise sensitivity (a caution, no verdict): data variant "noisy" = the same data plus extra noise on 10% of
D1's features (2,000 features drawn, from seed data seed + 200000, among those loaded on no factor) so that their
noise sd is 2 x 1.5 = 3.0. Arms noisy_cap1000 and noisy_cap2000 (with scale_views), all three modes, N = 20, weak and
strong.

### The cross (reduced for budget, fixed now)

| scenario | N | F2 | modes | arms | cells per seed |
|---|---|---|---|---|---|
| n20weak | 20 | weak | fast, medium, slow | 9 main + 2 noisy | 33 |
| n20strong | 20 | strong | fast, medium, slow | 9 main + 2 noisy | 33 |
| n60weak | 60 | weak | fast, slow | 9 main | 18 |
| n60strong | 60 | strong | fast, slow | 9 main | 18 |

10 seeds each: 1,020 trainings, 40 jobs (one per scenario x seed). N = 60 has no "medium" mode: the lifecycle's
uncapped slow cell at N = 60 took 58 s for 531 iterations, and the full cross was estimated above the 30 core-hour
ceiling. No cell is dropped after seeing results.

## Metrics

As in the first run (per cell and planted factor: matched |r|, variance explained per data set by the matched
factor, features kept per data set, iterations, learned factors) plus stopped-by and the final ELBO. For every
capped arm (cap250 to cap5000, noisy_cap1000, noisy_cap2000), per seed and data set: how many features loaded on F1,
F2, F3, on none, and (noisy) with inflated noise were kept and dropped (`composition.tsv`).

## Pre-declared readings (not changed after reading)

Marks: a factor is recovered in at least 9 of 10 seeds; "losing F3" = F3 |r| < 0.9 in more than 1 of 10 seeds (as
the first declaration). A failed or missing cell counts as not recovered.

- (a) Early stopping. On n20weak (the first run's setting) under "slow": if arms none and scale_views each reach
  F1 >= 9 and F3 >= 9 of 10, the first run's poor uncapped result was early stopping. If not, it was not (or not
  only). Also reported for the other scenarios, without verdict.
- (b) Convergence default. Recommend the fastest mode (fast, then medium) for which, for every capped arm (cap250,
  cap500, cap1000, cap2000, cap5000), every scenario where both that mode and slow ran, and every planted factor,
  seeds recovered under the mode >= seeds recovered under slow minus 1. If neither passes, slow. ("medium" is
  therefore judged on N = 20 only.)
- (c) Cap. Under the mode from (b) (for N = 60, where that mode did not run, slow): per scenario, the caps among 250,
  500, 1,000, 2,000, 5,000 that recover F2 (>= 9) without losing F3 (>= 9). A scenario where any cap meets this
  "discriminates". The recommended cap is the largest cap that meets it in every discriminating scenario; if no cap
  is common to them, "no single cap" with the per-scenario results; if no scenario discriminates, "no cap meets the
  rule anywhere". If n20weak does not discriminate but n20strong or an N = 60 scenario does, the conclusion is that
  F2 at weak strength is not detectable at N = 20 in this design (stated as a limitation of the claim), and the cap
  is chosen from the discriminating scenarios as above.
- (d) weight_views. Same two marks per scenario under the same mode; reading "meets" only if it meets both in every
  discriminating scenario; it is offered as an option only on "meets".
- Inflated-noise arm: reported descriptively (recovery beside the clean cap1000 / cap2000 at the same cells, and the
  kept composition), as a caution on variance-based selection when the most variable features are not the
  informative ones.
- Instrumentation check, not a gate: the fast-mode weak cells repeat the first run's cells (same data, same
  settings); the scorer reports how many |r| values match the first run exactly and the largest difference.

## Budget

- Lifecycle (done): run 44653902387, one job (n60strong, seed 1) with three cells: none@slow (531 iterations, 58 s,
  stopped by convergence), noisy_cap1000@fast (composition written), probe_fail@fast (failed on purpose, written as
  failed; the job was marked failed and the score job still ran and counted it). Job 2 min 20 s. It showed those two
  cells' results for seed 151001 at N = 60 strong; nothing was changed after it.
- One counted run: the cross above.
- One bounded correction pass for an instrumentation fault only, rerunning only affected jobs.
- Ceiling 30 core-hours (job wall time x 4 cores). Estimate: about 10 to 11 min per N = 20 job and 8 to 12 min per
  N = 60 job including about 2.3 min of installation: about 20 to 30 core-hours; the upper end is near the ceiling.
  An overrun is reported as such.

## Limits stated in advance

One generator (Gaussian, linear, equal effect sizes except F2's strong setting, all data sets on one scale); N = 20
and 60; one training seed; N = 60 without "medium"; one inflated-noise pattern.
