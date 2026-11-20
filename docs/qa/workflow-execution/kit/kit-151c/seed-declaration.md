# KIT-151 c, task 1 (SCI-06): does the training seed give different models? Declaration

Fixed before the counted run. Branch `run-mofa-evidence-k151c`; lane `.github/workflows/k151c-evidence.yml`, script
`.github/mofa-smoke/k151c.R` (task `seed`), scorer `.github/mofa-smoke/k151c-score.R` (task `seed`).
Script hashes: see the end of this file. Lifecycle run: 44677064790 (all three tasks; deliberately failing cells reported as failed; two scorer defects fixed in 1 commit after it, evidence script unchanged).

## Purpose

The audit (SCI-06) infers from the mofapy2 source that training starts from a PCA of the joined data and that the
seed only perturbs scikit-learn's randomized SVD, so several seeds would return the same model. This run measures
it. No default, wording or package code changes here.

## Bindings

- MOFA2 1.22.0, mofapy2 0.7.3, numpy 1.26.4, scipy 1.12.0, scikit-learn 1.4.0, R 4.6.1: installed as the smoke lane
  installs them (bundle versions, checked against `kit-134/p2/linux/bundle-sources.tsv` and `requirements.txt`).
- Data: the k151 generator (copied unchanged from `k151-followup.R`), so the data equal the k151 runs' data:
  - hard: `n20weak@7` (N = 20, weak F2 in the two small data sets, data seed 151007; the case that diverged);
  - easy: `n60strong@1` (N = 60, F2 weight 1.5, data seed 151001).
- Arms: `scale_views` (uncapped) and `cap1000` (the follow-up rule's cap); modes `fast` and `slow`.
- Training seeds 1 to 10; everything else as k151: 5 factors, double precision, maxiter 1,000.
- One job (one processor) trains all 80 models.
- The job prints and saves (`source-excerpts.txt`) the installed mofapy2 lines that build the start and use the seed
  (`build_model.py`, `init_model.py`, `entry_point.py`, `BayesNet.py`) and scikit-learn 1.4.0's PCA solver choice.
  It also records, per model, the start values of the factors (hook on `initModel.initZ`) and the PCA solver
  actually used (hook on `PCA._fit`).

## Measures

Per cell (data set x arm x mode), over the 45 seed pairs: factors matched greedily by |r|; the smallest matched
|r| of the factor values (Z) and of the weights (W, all data sets stacked); the largest absolute difference of Z and
W after matching and sign alignment; bit-identity; iterations; ELBO range; the largest difference between start
values. Per factor of seed 1: median and smallest matched |r| against the other nine seeds, its variance explained
and the closest planted factor. Planted recovery (|r| to F1 to F3) per seed.

## Tiers and readings (declared)

- A: every pair bit-identical. B: largest abs factor difference after matching <= 1e-8 (machine precision).
  C: every matched factor |r| >= 0.99 (same solution, small numerical differences). D: otherwise (different
  solutions).
- All cells A or B: seeds give the same model; multi-seed training cannot measure stability; resampling would be
  needed.
- Any cell C or D: say by how much (pairs below 0.99, smallest |r|, which factors, whether planted recovery crosses
  0.9 between seeds) and in which cells. A seed difference is not a stability measure under resampling; that is not
  tested here.
- Start values: largest abs difference across seeds <= 1e-8 confirms that the seed changes the start only at
  rounding level.

## Seen in the lifecycle (not counted, stated before the counted run)

Seeds 1 and 2 on `n20weak@7`, `cap1000`, `fast`: same iterations (76), ELBO differs by 0.09, smallest matched |r|
0.984, start values differ by up to 1.95 (on mofapy2's -1 to 1 scale); the PCA solver used was `randomized` on
20 x 2,300. The counted run measures this over 10 seeds, two data sets, two arms and two modes.

## Budget

One job, about 80 small trainings; ceiling 60 minutes wall (4 vCPU: at most 4 core-hours; expected about 1.5).
One counted run; one correction pass at most, only for an instrumentation failure, reported as such.

## Script hashes (sha256, fixed for the counted run)

```
3d97c5f55a56c9415a3d1257bbe736350ccdce5a31dace25ef54f52da3d97afa  .github/mofa-smoke/k151c.R
84083aa17d79227309f92da2b9837ebea715d42f6f8ba369ff66ad2b38cd0b12  .github/mofa-smoke/k151c-score.R
405d84a1b78f7cdf2324ef7d85dc309539e88a1264a1ad2a29ca448c0bea2d17  .github/mofa-smoke/install-ruv-library.R
88084c4062a4b70b7a8a100ceec0198ca7a64c455ab37105ccca4e3ae37ccfad  .github/workflows/k151c-evidence.yml
```
