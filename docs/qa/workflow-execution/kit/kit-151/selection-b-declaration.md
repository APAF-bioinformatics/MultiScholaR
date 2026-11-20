# KIT-151 selection evidence, run b: declaration (fixed before the counted run)

Authority: user, 2026-11-20, relayed by the coordinator (audit `tickets/kit/AUDIT-integration-series-20261120-scientific.md`,
SCI-01, SCI-02, SCI-03, SCI-13, "smallest next evidence"). Evidence only; no package code or default changes; the
user decides. Written after the lifecycle run (below) and before the counted run; committed and pushed on its own
before the counted run's push. The first run and the follow-up stand as read; `sci01-rescore.md` re-reads the
follow-up's rows without new compute.

## Purpose

The earlier simulator made the most variable features the signal-carrying ones (equal noise everywhere), so a
raw-variance cap both balanced the data sets and enriched them for signal, and the two could not be separated.
Real proteomics, metabolomics and RNA-seq data have a mean-variance trend and abundance-dependent missing values,
where the most variable raw features are often the noisiest. This run asks, on data with those properties, at N = 20
under convergence mode slow:

1. Does selecting features before MOFA2 help at all, against no selection?
2. Which statistic: raw variance, a random subset of the same size (balance without enrichment), or variance above
   the fitted mean-variance trend?
3. Does the number (250, 1,000, 2,000) matter?
4. Does a stricter eligibility rule (observed in at least half the samples) change anything against the current one
   (observed in at least 2 samples)? (SCI-03)
5. Does a per-feature z-score after a raw cap 1,000 change anything?

## Scripts

- `.github/mofa-smoke/k151b.R`, sha256 `cc3133360c47dd3446365526399a592441fb55aa8a8cd1fa911a4f567877492e`
- `.github/mofa-smoke/k151b-score.py`, sha256 `d500aa9fb8e6f10f5db5db408b0652d99f9d85d918096f3540decd09635173f5`
- Workflow `.github/workflows/k151b.yml`; parameters `.github/mofa-smoke/k151b-params.json`.
- Both scripts are byte-identical to the lifecycle's. The scorer runs once in the counted run's score job; its output
  is read once.

## Data (generator)

The follow-up's n20strong generator (D1 20,000 features with F1 and F3; D2 2,000 with F1; D3 200 and D4 100 with F1
and F2; 10% loaded features per factor and data set; weights +/-1, F2 +/-1.5; per-feature offset ~ N(20, 2) = the
feature's abundance; F1 = group + noise sd 0.5; N = 20, two groups of 10), with three switchable features:

- **trend**: noise sd of feature j = 1.5 x 2^(-(offset_j - 20) / 4). The median noise sd stays 1.5 (the earlier
  design) and doubles for every 4 units lower abundance (about 0.75 at offset 24, 3.0 at 16). Chosen as a moderate
  trend of the shape seen in log-intensity data (low-abundance features noisier), not fitted to any data set.
  Signal loadings do not depend on abundance.
- **missing**: each value x is missing with probability 1 - (1 - plogis((17.5 - x) / 0.5)) x (1 - 0.02): a soft
  detection limit near 17.5 (left-censoring of low values) plus 2% missing completely at random, drawn from its own
  random stream (data seed + 400,000) on the unscaled values, in every data set. The realised share per data set and
  seed is recorded in `job.json` (lifecycle: see below).
- **scales**: after the missing values, D1 x 1, D2 x 2, D3 x 0.25, D4 x 4 (a 16-fold range of spread across data
  sets), so `scale_views` has something to do.

Scenarios:

| scenario | trend | missing | scales | purpose |
|---|---|---|---|---|
| prev | off | off | off | the earlier design: seeds 1 to 10 are the follow-up's n20strong data exactly (tie) |
| trend | on | off | on | mean-variance trend and scales, no missing values |
| full | on | on | on | trend, abundance-dependent missing values and scales |

Data seeds 151001 to 151030 (30 paired seeds; the same seed gives the same factors and loadings in every scenario).
Training seed 20261119, double precision, 5 factors, convergence mode slow, round limit (maxiter) 10,000,
`center_groups` TRUE, MOFA2 called directly as in the earlier runs (`weight_views` FALSE).

## Arms

Preparation as `integrationAlignViews`: eligibility, then (z-score arm only) the zero-spread filter, then selection
of at most `cap` features per data set (data sets at or below the cap are kept whole: D3 and D4 always, D2 for cap
2,000), then (z-score arm only) per-feature z-score. All arms have `scale_views` TRUE except `none`.

| arm | selection | cap | eligibility | z-score | scale_views |
|---|---|---|---|---|---|
| none | none | - | >= 2 observed | no | FALSE |
| scale_views | none | - | >= 2 observed | no | TRUE |
| raw250, raw1000, raw2000 | raw variance (largest) | 250, 1,000, 2,000 | >= 2 observed | no | TRUE |
| random1000 | random subset (own stream, data seed + 300,000 + data set index) | 1,000 | >= 2 observed | no | TRUE |
| trend1000 | variance minus the fitted trend (difference, as scran's "biological" component) | 1,000 | >= 2 observed | no | TRUE |
| trendratio1000 | log variance minus the fitted trend of log variance (ratio to the trend, as Seurat's vst standardised variance) | 1,000 | >= 2 observed | no | TRUE |
| raw1000_z | raw variance, then per-feature z-score | 1,000 | >= 2 observed | yes | TRUE |
| *_obs50 (seven arms) | as raw250, raw1000, raw2000, random1000, trend1000, trendratio1000, raw1000_z | same | >= 10 of 20 observed (half) | as base | TRUE |

Trend: `lowess(mean, statistic, f = 0.3)` per data set over the eligible features (mean and variance of the
observed values), interpolated at each feature's mean. The task named one trend-adjusted arm; both common forms
(difference and ratio) are run because they rank differently when the noise of the variance estimate grows with the
trend (the difference favours noisy low-abundance features; the ratio does not). This is stated here, before the run,
as a second definition, not a choice made after seeing results.

The `*_obs50` arms run only in `full` (without missing values they are identical to the base arms).
Cells per seed: prev 9, trend 9, full 16 = 34; 30 seeds = 1,020 trainings.

## Metrics (per seed, scenario and arm)

- For each planted factor F1 (group, all data sets), F2 (D3 and D4 only), F3 (D1 only): |r| with its best-matched
  learned factor (greedy maximum matching, as before), and "worst" = the smallest of the three in that seed.
- Variance explained in each data set by the matched factor, of the **selected** features (MOFA2's r2 per factor).
- Iterations; stopped by the convergence rule or by the round limit (iterations >= 10,000); the last finite ELBO and
  its relative change over the last five recorded rounds.
- Composition per data set: features kept and dropped by what they carry (F1, F2, F3, none, all), the median true
  noise sd and the median observed fraction of the kept features. From this: share of signal-carrying features kept,
  share of kept features that carry signal.

## Analysis (paired, continuous; no pass/fail threshold on single seeds)

Primary measure: the per-seed worst-of-three |r|. Secondary: each factor separately, variance explained, composition.
For each declared comparison below and each scenario where both arms ran: the per-seed difference (b minus a), its
mean with a 95% seed-bootstrap interval (20,000 resamples, `random.Random(151)`), median difference, seeds better /
worse, worst seed difference. Reading of one comparison: interval inside [-0.02, +0.02] = "no practical difference"
(0.02 |r| is below what the matched |r| moves between neighbouring caps in the follow-up); else lower bound > 0 =
"better"; upper bound < 0 = "worse"; otherwise "inconclusive". The intervals are not adjusted for the number of
comparisons (about 30 comparisons x 3 scenarios x 4 measures); the readings are descriptive evidence for the user's
decision, not a test that selects a default. The count of seeds with |r| >= 0.9 is reported descriptively only.
A failed or missing cell drops that seed from the comparisons it is in; the number dropped is reported.

Declared comparisons (a -> b):

1. Selecting at all: scale_views -> raw250, raw1000, raw2000, random1000, trend1000, trendratio1000.
2. Statistic: raw1000 -> random1000, trend1000, trendratio1000; random1000 -> trend1000, trendratio1000;
   trend1000 -> trendratio1000.
3. Number: raw1000 -> raw250, raw2000.
4. Eligibility (full only): each of the seven base capped arms -> its `_obs50` version.
5. z-score: raw1000 -> raw1000_z; raw1000_obs50 -> raw1000_z_obs50.
6. scale_views: none -> scale_views.

Instrumentation check (not a gate): the prev scenario, seeds 1 to 10, arms none, scale_views, raw250, raw1000,
raw2000 repeat the follow-up's n20strong slow cells on the same data (the follow-up stopped at 1,000 rounds, this run
at 10,000). Cells that stopped by the rule before 1,000 rounds there are expected to agree to rounding level on the
same CPU model; the scorer reports how many agree and the largest difference, with the CPU models of both runs.

## What the possible outcomes would mean for the default

- **Selection helps everywhere** (selection arms "better" than scale_views in prev, trend and full): supports
  selecting before MOFA2 at N = 20; the statistic is read from 2.
- **Raw variance helps in prev but not in trend / full, while trend-adjusted or random does**: raw variance is the
  wrong statistic on realistic data; the default should be trend-adjusted (or the cap's benefit is balance only, if
  random matches trend).
- **Random ~ raw ~ trend ("no practical difference")**: the benefit is balancing the data sets; which features are
  kept matters little for these factors; then the simplest defensible statistic is a matter for the user, with the
  composition tables showing what each keeps.
- **Trend-adjusted "better" than random**: selection on variance above the trend adds signal enrichment beyond
  balance on realistic data.
- **Every selection arm "worse" or "inconclusive" against scale_views in trend / full**: no evidence that selecting
  helps on realistic data; selection should not be a default on this evidence.
- **Number**: if raw250 / raw2000 against raw1000 read "no practical difference" (or inconclusive) the evidence
  supports a band, not a value; the run does not pick a value by any threshold.
- **Eligibility**: "better" for obs50 in full supports the SCI-03 repair (a minimum observed share); "no practical
  difference" means the current rule does no measurable harm in this design; "worse" argues against it.
- **z-score**: "worse" or "no practical difference" for raw1000_z leaves the module's current choice (no per-feature
  z-score) unchallenged; "better" argues for offering or defaulting to it.
- **Round limit**: the share of cells stopping at 10,000 is reported; if it is material (more than 5% of cells), the
  "slow" default needs a different stopping rule or limit, and the |r| of those cells is not final.

## Budget

- Lifecycle (done, not counted): run 44676343317 (commit 24babf60), data seeds 151091 and 151092 (outside the counted set), every one of the 34 cells plus a deliberately failing cell per job. Both jobs ran every cell to the end; the failing cell was written as failed, each job was marked failed (exit 1), and the score job still ran and counted 68 of 68 cells ok, 2 failed (the probes), 0 missing. Instrumentation seen: realised missing share in `full` 18 to 21% per data set (0 in prev and trend); iterations up to 1,811, no cell at the 10,000 limit; training 362 s and 344 s per job on AMD EPYC 7763; job wall time about 7.4 min including installation. Nothing in the scripts, the parameters of the generator or the arms was changed after it.
- One counted run: 30 jobs (one per data seed, all 34 cells), free GitHub-hosted runners (4 cores, 16 GB), job
  timeout 300 minutes.
- One correction pass at most, for an instrumentation fault only, rerunning only the affected jobs.
- Expected: about 7.5 min per job x 30 jobs x 4 cores = about 15 core-hours (longer if other seeds need more rounds). Ceiling: 30 core-hours (job wall time x 4 cores). An overrun is reported as such.

## Limits stated in advance

One generator (Gaussian, linear; equal effects except F2; one trend shape; one missing-value mechanism; one set of
scales); N = 20 only; one training seed; one round limit; the two trend fits with one span; one eligibility share
(half). The trend and missing-value parameters are judgments, not estimates from real data.
