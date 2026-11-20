# KIT-151 cap evidence: declaration (fixed before the counted run)

Ticket: KIT-151 (plan `tickets/kit/PLAN-integration-workbook-20261120.md`, section E, approved by the user
2026-11-20). Evidence only: no package code and no default changes here. Written 2026-11-20, after the lifecycle run
(run 44650224302 (commit 494a7361)) and before the counted run; pushed on its own before the counted run's push.

## Question

With MOFA2 1.22.0, mofapy2 0.7.3, double precision, when one data set has 10 to 200 times more features than the
others, which treatment lets the model find a factor that exists only in the small data sets (F2) without losing a
factor that exists only in the large one (F3)?

## Where and what runs

- Branch `run-mofa-evidence-k151` (from `run-mofa-smoke-k062` 610d8534), workflow `.github/workflows/k151-evidence.yml`.
  GitHub-hosted `ubuntu-22.04` runners; R 4.6.1, MOFA2's closure from the bundle's snapshot (versions checked
  against the bundle's), the bundle's Python with mofapy2 0.7.3. No MultiScholaR, no token.
- Evidence script `.github/mofa-smoke/k151-evidence.R`, sha256 `d8c4b5e71093e576e430165a19d94e42a6525eb5329a77b1bba0c80790d7150e`.
- Scoring script `.github/mofa-smoke/k151-score.py`, sha256 `cba9761b436253945d3560df9e3a489c2b4b79f46748366008841e89c62820af`. It runs once in the counted run's
  score job; its output is read once.
- One job per scenario x data seed; each job trains every arm in turn.
- MOFA2 is called directly (`create_mofa` -> default data, model and training options with the declared ones set ->
  `prepare_mofa` -> `run_mofa(save_data = TRUE, use_basilisk = FALSE)` -> `load_model`), as
  `integrationTrainMofa` does (R/utils_workflow_integration_model.R, origin/kit-142-integration-record b20d4fd7).
  View preparation follows `integrationAlignViews` there: features observed in fewer than 2 samples dropped; under
  the z-score, features without spread dropped; the cap keeps `min(cap, D)` features by the variance of each feature
  on its own scale (`var(na.rm = TRUE)`), ties to the earlier row, kept rows in their original order; the z-score is
  `(x - mean) / sd` per feature with `na.rm = TRUE`, applied after the cap.

## Generator (self-contained Gaussian latent model, in the evidence script)

- RNG: Mersenne-Twister, Inversion, Rejection. Data seed = 151000 + seed index (seed index 1 to 10: 151001 to
  151010). Missing-value mask seed = data seed + 100000.
- Samples: N = 20 (two groups of 10) in the main scenario; N = 60 (two groups of 30) in a sensitivity scenario.
- Planted factors (sample scores), standardised to mean 0 and sd 1 over the samples:
  - F1 = group indicator (-1 / +1) + Normal(0, 0.5^2) noise: the group factor, in all four data sets.
  - F2 ~ Normal(0, 1): only in D3 and D4 (the two small data sets).
  - F3 ~ Normal(0, 1): only in D1 (the large data set).
- Data sets: D1 20,000 features (F1, F3); D2 2,000 (F1); D3 200 (F1, F2); D4 100 (F1, F2). D1 is 10 to 200 times
  larger than the others.
- Loadings: for each data set and each factor active in it, 10% of the features (drawn independently per factor, so a
  feature can carry two factors) get weight +1 or -1 (equal per-feature effect size 1, random sign); all others 0.
- Value = loadings x scores + Normal(0, 1.5^2) noise per value + a per-feature offset Normal(20, 2^2) (an intensity
  level that centring must remove). All data sets on the same scale.
- Missing values: none in the main and N = 60 scenarios. Sensitivity scenario "missing": N = 20 with 15% of values
  set missing completely at random in every data set.

## Arms (every arm: 5 factors, convergence mode fast, maxiter 1000, training seed 20261119, `use_float32 = FALSE`)

| arm | features kept | scale_views | center_groups | weight_views |
|---|---|---|---|---|
| none | all | FALSE | TRUE | FALSE |
| scale_views | all | TRUE | TRUE | FALSE |
| cap1000, cap2000, cap5000, cap10000 | `min(cap, D)` most variable per data set | TRUE | TRUE | FALSE |
| weight_views | all | TRUE | TRUE | TRUE (MOFA2 training option, located by name in the option lists) |
| zscore (the collaborator workbook's) | all, per-feature z-score | FALSE | FALSE | FALSE |

## Metrics per arm and seed

- For each planted factor F1, F2, F3: |Pearson r| with the learned factors over the samples; matching is greedy
  maximum (the largest |r| in the planted x learned matrix first, then the largest among the remaining rows and
  columns; one learned factor per planted factor). A planted factor left without a learned factor gets |r| = 0.
- Variance explained (MOFA2's `get_variance_explained`, percent, per data set) by the matched factor, on the data the
  model was given (after cap, z-score and MOFA2's own centring and scaling).
- Features kept per data set; iterations; learned factors and whether MOFA2 dropped any (5 minus learned);
  MOFA2 warnings count; seconds.

## Pre-declared rule (from the plan; not changed after reading)

- A planted factor is recovered in a seed when its matched |r| >= 0.9.
- An arm "recovers F2" when F2 is recovered in at least 9 of 10 seeds.
- "Losing F3" means F3 |r| < 0.9 in more than 1 of 10 seeds; an arm "keeps F3" when F3 is recovered in at least 9
  of 10 seeds.
- A failed arm-seed counts as not recovered.
- Recommendation (decided on the main scenario, N = 20, no missing values only): the largest cap in
  {1000, 2000, 5000, 10000} that recovers F2 and keeps F3. 10,000 is the upper bound considered. If none does:
  "no cap meets the rule", reported with the closest cap (most F2 seeds, then most F3 seeds, then the larger cap).
- `weight_views` is judged by the same two marks (recovers F2, keeps F3) on the main scenario. It is reported as
  meeting or not meeting them; it is not a candidate for the cap value.
- F1 is reported for every arm but is not part of the rule.
- The N = 60 and 15%-missing scenarios are sensitivity: scored and reported the same way, never used to choose.
- An unfavourable or ambiguous result is reported as it is.

## Other checks recorded in the same run

- The installed mofapy2 source lines that mention `scale_views`, `weight_views`, `center_groups`, `scale_groups`
  (file and line, with context), and the installed MOFA2 R functions whose deparsed body mentions these or
  `variance_explained` / `save_data`, printed into `source-excerpts.txt` per job.
- `save_data = FALSE`: a tiny two-view model trained twice (save_data TRUE and FALSE, same seed), loaded with
  `load_model`, `get_variance_explained` called on each; recorded: whether it works with FALSE, its error if not,
  the HDF5 top-level groups, and the largest difference in r2 per factor between the two.

## Budget

- Lifecycle run (done before this declaration): 1 seed, 2 arms (`scale_views` and `probe_fail`, an arm that fails on
  purpose) on the main scenario: run 44650224302 (commit 494a7361), the job ran through to its end and was marked failed; `probe_fail` was written as a
  failed row with its error; the score job still ran and counted the failed rows (6 of 6 rows, 3 failed); the source
  excerpts and the save_data check were written. The lifecycle also showed one arm's metrics for seed 151001
  (`scale_views`: F1 0.85, F2 0.35, F3 0.59, 26 iterations, 5 factors). Nothing in the generator, arms, metrics or
  rule was changed after it; the scripts are byte-identical to the lifecycle's (same sha256).
- One counted run: scenarios main (N = 20), n60 (N = 60), missing (N = 20, 15% MCAR), 10 seeds each, all 8 arms: 30 jobs.
- One bounded correction pass, for an instrumentation fault only (a job or arm that fails for a reason other than
  the model), rerunning only the affected jobs with the same scripts unless the fault is in a script.
- Ceiling 40 core-hours on hosted runners (job wall time x the runner's cores). Expected from the lifecycle:
  the lifecycle job took 1 min 42 s on a 4-core runner (install about 1.5 min, one arm 2.7 s training); 30 jobs of 3 to 6 min on 4 cores are about 6 to 12 core-hours.

## Limits stated in advance

One generator, Gaussian, linear, equal effect sizes, all data sets on one scale (so `scale_views` has little to do
here by design); N = 20 and 60 only; one training seed; one noise level and loaded share. The result says what
happens in this design, not in every real study.
