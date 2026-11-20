# KIT-151 c, task 2 (SCI-05): does pinning the maths library make processors agree? Declaration

Fixed before the counted run. Branch `run-mofa-evidence-k151c`; lane `.github/workflows/k151c-evidence.yml`, script
`.github/mofa-smoke/k151c.R` (task `pin`), scorer `.github/mofa-smoke/k151c-score.R` (task `pin`).
Script hashes: see the end of this file. Lifecycle run: 44677064790 (all three tasks; deliberately failing cells reported as failed; two scorer defects fixed in 1 commit after it, evidence script unchanged).

## Purpose

The k151 follow-up found that data seed 151007 (N = 20, weak, convergence `fast`) trained to different models on
two processors, with different iteration counts and |r| differences up to 0.29. From the job tables: the first run
(`jobs.tsv`, k151-main-7) ran on an Intel Xeon 6973P-C; the follow-up (`followup-jobs.tsv`, k151f-n20weak-7) on an
AMD EPYC 7763. This run asks whether fixing the BLAS/LAPACK code path and thread count before Python starts makes
processors give the identical model. No default, wording or package code changes here.

## Bindings

- Software as the k151 runs (MOFA2 1.22.0, mofapy2 0.7.3, numpy 1.26.4, scipy 1.12.0, scikit-learn 1.4.0, R 4.6.1;
  the bundle's Python). numpy and scipy wheels each ship their own OpenBLAS (DYNAMIC_ARCH); the job records which
  BLAS libraries are loaded, their version, thread count and the kernel type in use (`threadpoolctl`), numpy's
  `show_config` and numpy's CPU features (`python-env.json`).
- Case: the k151 generator, `n20weak@7` (data seed 151007), training seed 20261119, mode `fast`, the first run's
  arms: none, scale_views, cap1000, cap2000, cap5000, cap10000, weight_views, zscore (8 models per condition).
- Conditions, each a separate R process with the variables set before R (and so Python) starts:
  - U, unpinned: no variable set (OpenBLAS default threads, kernel chosen per CPU).
  - P, pinned: `OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_CORETYPE=Haswell`. Haswell
    (AVX2 and FMA) is the kernel OpenBLAS already uses on the AMD EPYC 7763 and runs on every hosted-runner CPU seen
    so far; MKL is set only in case it is loaded (the job shows whether it is).
  - PN, secondary: P plus `NPY_DISABLE_CPU_FEATURES` = numpy's dispatched AVX-512 targets that the CPU has (empty
    on a CPU without AVX-512), to separate numpy's own SIMD paths from BLAS if P is not enough.
- 12 identical jobs (pin-01 to pin-12). The processor cannot be chosen on hosted runners; 12 jobs landed on at least
  4 models in each earlier k151 run. If fewer than two models are seen, or no CPU without AVX-512 together with one
  with it, the report says so.

## Measures

Per condition and arm: classes of bit-identical models (factors, weights and iteration count) and the CPU models in
each; for every job pair, same or different CPU model: identical or not, iterations, the largest abs difference of
factors and weights as stored and after matching, the smallest matched |r|. History check (instrumentation): the
unpinned iterations on a CPU model that ran a committed k151 job against that job's committed row.

## Readings (declared)

- Instrumentation first: every cell trained; under P and PN every loaded OpenBLAS reports 1 thread and kernel
  Haswell; under PN numpy reports no AVX-512 feature on. A failure here is an instrumentation failure, not a result.
- U: does the divergence reproduce (more than one class across CPU models, different iterations)?
- P and PN, each: "identical across the CPU models seen" only if every arm has one class across all jobs
  (the audit's acceptance: identical output only); "agree to rounding" if iterations agree and the largest raw
  difference is <= 1e-10 but not bitwise; otherwise "processors still give different models".
- Same CPU model, same condition: expected identical; a difference is reported as such (it would mean run-to-run
  variation, for example from threading).

## Budget

12 jobs x 3 conditions x 8 small trainings; ceiling 60 minutes wall per job (expected about 6 minutes each, about
5 core-hours in all). One counted run; one correction pass at most, only for an instrumentation failure.

## Script hashes (sha256, fixed for the counted run)

```
3d97c5f55a56c9415a3d1257bbe736350ccdce5a31dace25ef54f52da3d97afa  .github/mofa-smoke/k151c.R
84083aa17d79227309f92da2b9837ebea715d42f6f8ba369ff66ad2b38cd0b12  .github/mofa-smoke/k151c-score.R
405d84a1b78f7cdf2324ef7d85dc309539e88a1264a1ad2a29ca448c0bea2d17  .github/mofa-smoke/install-ruv-library.R
88084c4062a4b70b7a8a100ceec0198ca7a64c455ab37105ccca4e3ae37ccfad  .github/workflows/k151c-evidence.yml
```
