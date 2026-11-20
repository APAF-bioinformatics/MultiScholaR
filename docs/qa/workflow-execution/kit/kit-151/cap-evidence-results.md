# KIT-151 cap evidence: results

Declaration: `cap-evidence-declaration.md` (commit 6953941f, pushed alone before the counted run). Lifecycle run
44650224302 (494a7361); counted run 44650673965 (08e87cc6), scored once by `k151-score.py` (sha256 cba9761b...,
unchanged; evidence script sha256 d8c4b5e7..., unchanged). 30 jobs, 720 of 720 rows, 0 failed arms, one identical
source excerpt in all jobs. No correction pass was needed. Hosted time: 3,540 job-seconds on 4-core runners, about
3.9 core-hours (ceiling 40). MOFA2 1.22.0, mofapy2 0.7.3, R 4.6.1, Python 3.12.16, numpy 1.26.4; double precision.

Files: `results-rows.tsv` (one row per scenario x seed x arm x planted factor), `results-summary.tsv` (per scenario x
arm x planted factor: seeds at |r| >= 0.9, median and minimum |r|, median variance explained per data set, median
features kept, median iterations, dropped factors), `score.json` / `score.md` (the rule applied), `jobs.tsv` (CPU and
versions per job), `source-excerpts.txt`, `save-data-check.json`, `lifecycle-rows.tsv`.

## Main scenario (N = 20, no missing values; decides)

Seeds (of 10) with |r| >= 0.9, median |r|; variance explained (median %, per data set D1/D2/D3/D4) by the matched
factor. D1 20,000 features (F1, F3), D2 2,000 (F1), D3 200 and D4 100 (F1, F2).

| arm | kept D1/D2/D3/D4 | F1 | F2 | F3 | R2 F1 | R2 F2 | R2 F3 | iterations |
|---|---|---|---|---|---|---|---|---|
| none | 20000/2000/200/100 | 2, 0.853 | 0, 0.382 | 1, 0.734 | 7.2/4.9/5.0/5.3 | 1.2/1.3/1.9/1.0 | 5.5/2.5/2.7/2.4 | 26 |
| scale_views | 20000/2000/200/100 | 2, 0.855 | 0, 0.381 | 1, 0.736 | 7.2/4.9/5.0/5.3 | 1.2/1.2/1.8/1.1 | 5.5/2.6/2.7/2.4 | 26 |
| cap1000 | 1000/1000/200/100 | 10, 0.997 | 2, 0.623 | 9, 0.990 | 18.9/9.5/6.3/6.6 | 1.1/1.9/4.7/6.0 | 18.5/1.4/1.7/1.4 | 81 |
| cap2000 | 2000/2000/200/100 | 10, 0.976 | 0, 0.631 | 6, 0.940 | 17.4/5.8/5.8/5.8 | 1.5/1.6/3.7/3.8 | 14.6/1.6/2.1/1.4 | 48.5 |
| cap5000 | 5000/2000/200/100 | 5, 0.902 | 0, 0.497 | 3, 0.826 | 13.7/5.1/5.1/5.4 | 1.4/1.6/3.9/2.9 | 10.6/2.2/2.5/2.0 | 36 |
| cap10000 | 10000/2000/200/100 | 2, 0.874 | 0, 0.478 | 2, 0.782 | 10.6/5.1/5.1/5.5 | 1.4/1.5/2.7/2.3 | 8.2/2.5/2.7/2.3 | 31 |
| weight_views | 20000/2000/200/100 | 3, 0.892 | 0, 0.655 | 0, 0.499 | 5.3/4.9/5.8/8.3 | 0.9/1.1/4.6/6.2 | 2.1/1.3/1.3/5.4 | 108.5 |
| zscore (workbook) | 20000/2000/200/100 | 2, 0.855 | 0, 0.367 | 2, 0.714 | 6.9/5.3/5.4/5.7 | 1.4/1.5/1.6/1.2 | 5.7/3.3/3.4/3.3 | 21 |

No arm-seed dropped a factor (5 learned factors everywhere).

**Rule (pre-declared): no cap meets the rule.** No cap recovers F2 in 9 of 10 seeds. Closest: cap 1,000 (F2 in 2 of
10 seeds, median |r| 0.62; F3 kept in 9 of 10). `weight_views` meets neither mark (F2 0 of 10, F3 0 of 10). The
workbook's z-score arm meets neither (F2 0, F3 2).

## Sensitivity (reported, not used to choose)

| arm | N = 60: F1 / F2 / F3 seeds (median F2 |r|) | N = 20, 15% missing: F1 / F2 / F3 seeds (median F2 |r|) |
|---|---|---|
| none | 8 / 0 / 4 (0.267) | 2 / 0 / 2 (0.342) |
| scale_views | 8 / 0 / 5 (0.267) | 2 / 0 / 2 (0.328) |
| cap1000 | 10 / 10 / 10 (0.941) | 10 / 0 / 9 (0.533) |
| cap2000 | 10 / 8 / 10 (0.936) | 8 / 0 / 5 (0.476) |
| cap5000 | 10 / 7 / 10 (0.930) | 5 / 0 / 2 (0.467) |
| cap10000 | 9 / 5 / 9 (0.621) | 2 / 0 / 2 (0.366) |
| weight_views | 10 / 8 / 8 (0.925) | 1 / 0 / 0 (0.508) |
| zscore | 7 / 0 / 4 (0.250) | 2 / 0 / 2 (0.413) |

At N = 60 the rule would pick cap 1,000 (the only cap meeting both marks).

## What the installed source does (`source-excerpts.txt`, file:line on the runner)

- `scale_views` (mofapy2/build_model/utils.py:86-88): after centring, `data[m] /= np.nanstd(data[m])` - one divisor
  per data set, the standard deviation of all its values. It does not change how many features a data set has, so
  it does not balance a large data set against small ones. MOFA2 default FALSE (`get_default_data_options`).
- `center_groups` (utils.py:80-84): subtracts each feature's mean within each group (`np.nanmean`, axis 0), only for
  Gaussian data. MOFA2 default TRUE.
- `weight_views` (MOFA2 `get_default_training_options`: default FALSE; passed by `.run_mofa_reticulate` line 33 to
  `set_train_options`, mofapy2/run/entry_point.py:876, 1022-1025, message "Weighting the ELBO ... based on the number
  of features per view"): each view m gets weight `total_D / (M * D_m)`, renormalised to sum to M; the weights enter
  the factor (Z) update (mofapy2/core/nodes/Z_nodes.py:114-118) and the ELBO (mofapy2/core/BayesNet.py:513-522).
  Here D = 20,000/2,000/200/100 gives weights of about 0.013/0.13/1.29/2.57 (computed from that formula): the large
  data set counts about 200 times less than the smallest in the factor update.

## save_data = FALSE (`save-data-check.json`)

Yes, `get_variance_explained` works on a model trained with `save_data = FALSE`. The HDF5 still holds the
`variance_explained` group (it lacks only `data` and `intercepts`); `load_model` caches it (MOFA2::load_model lines
250-251 compute it only when the cache is empty) and `get_variance_explained` reads the cache (lines 8-9). On the tiny
check model the r2 per factor with FALSE equals that with TRUE exactly (largest difference 0), in all 30 jobs.

## Observations (not part of the rule)

- The uncapped arms stopped after a median 21 to 26 iterations under `convergence_mode = "fast"` (the plan's fixed
  setting); capped arms ran 31 to 81, `weight_views` 108. Part of the poor recovery without a cap may be early
  stopping rather than the imbalance itself. This was not tested here.
- Smaller caps help monotonically in this design (1,000 > 2,000 > 5,000 > 10,000 for F1, F2 and F3), the opposite of
  "the largest cap that still works": with N = 20 no cap reaches the F2 mark, and with N = 60 only 1,000 does.
- The cap of 1,000 also cut D2 (2,000 -> 1,000).

## Not as declared

Nothing. All 30 jobs and 8 arms ran; no correction pass. The lifecycle run showed one arm's result for seed 151001
before the declaration was pushed (stated in the declaration); nothing was changed after it.

## Limits

One generator: Gaussian, linear, equal effect sizes (1, noise sd 1.5), 10% loaded features per factor, all data sets
on one scale (so `scale_views` had almost nothing to do by design); N = 20 and 60 only; one training seed; "fast"
convergence only. The result describes this design, not every real study.
