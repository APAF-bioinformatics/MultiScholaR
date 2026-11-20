# KIT-151 c results (seeds, pinning, group-corrected inputs), scored once

Branch `run-mofa-evidence-k151c`. Lifecycle 44677064790 (probe cells failed as intended; two scorer defects fixed in
331d1d81 before any declaration). Each counted run follows its declaration commit:

| task | declaration | counted run | jobs | core-hours (4 vCPU x wall) | failed cells |
|---|---|---|---|---|---|
| 1 seed (SCI-06) | e9340fae | 44678126809 | 1 (AMD EPYC 9V74) | 1.5 | 0 / 80 |
| 2 pin (SCI-05) | 7a5d11ae | 44678183859 | 12 x 3 conditions | 1.9 | 0 / 288 |
| 3 ruv (SCI-04) | 6bb289f5 | 44680786131 | 8 | 2.8 | 0 / 200 |

Lifecycle 0.65 core-hours; total about 6.9. Tables: `seed-*.tsv`, `pin-*.tsv`, `ruv-*.tsv`; scorer output
`*-score.md`; installed mofapy2 / scikit-learn source lines `seed-source-excerpts.txt`.

## Task 1: seeds do give different models, in the weak factors

- Source (mofapy2 0.7.3, as installed): `initModel.__init__` calls `np.random.seed(seed)`; `initZ(qmean="pca")`
  fits `PCA(n_components=K, whiten=True)` on the concatenated data (mean-imputed), then rescales each component to
  -1..1; W starts at 0 (`initW(qmean=0)`, `initSW(qmean_S1=0)`); the entry point's global seed line is commented out.
  scikit-learn 1.4.0 chose `randomized` in every model (20 or 60 samples x 2,300 or 22,300 features), so the seed
  enters only through the randomized SVD.
- That is not a rounding-level effect: start values differ by up to 2.0 (sign flips) and some start components do not
  match between seeds at all (smallest matched |r| 0.0001 to 0.22 per cell). Trailing PCA components are nearly
  tied, and the randomized sketch picks different ones.
- Every cell is tier C or D (none A or B). Strong factors agree across seeds (F1, F3 and, where detectable, F2:
  matched |r| >= 0.99); the weaker factors change (smallest matched |r| down to 0.002).
- Seed changes recovery of a real planted factor: n60strong `scale_views` F2 >= 0.9 in 4 of 10 seeds (fast) and
  8 of 10 (slow). Iterations vary with the seed (e.g. n20weak scale_views slow 516 to 756; n60strong cap1000 slow
  796 to 1,001, one at the round limit).
- Only cell near-identical: n20weak cap1000 slow (tier C, smallest |r| 0.999). ELBO spread is small relative to the
  ELBO (2e-7 to 1.6e-4).

## Task 2: pinning made the processors seen give bit-identical models; the reported divergence did not recur

- CPU models: AMD EPYC 7763 (6 jobs, Zen kernel), EPYC 9V74 (4; 3 with AVX-512 -> Cooperlake kernel, 1 without ->
  Zen), EPYC 9V45 (2; Cooperlake, R's system OpenBLAS 0.3.20 on Barcelona). No Intel CPU was assigned, so the
  historical pair (Intel Xeon 6973P-C vs EPYC 7763) is not tested.
- Three OpenBLAS copies are loaded (R's system 0.3.20, scipy's 0.3.21.dev, numpy's 0.3.23.dev); under P all three
  report Haswell and 1 thread; no MKL. Under PN numpy's six AVX-512 dispatch targets were off.
- U (unpinned): 3 classes per arm by kernel (Zen; Cooperlake on 9V74; Cooperlake on 9V45), same iterations
  everywhere, largest difference 5.8e-13: agree to rounding, not bit-identical. Same-CPU-model pairs differ only where
  the 9V74 instance lacked AVX-512 (different kernel).
- P and PN: one class per arm across all 12 jobs and 3 models: bit-identical.
- History check: unpinned EPYC 7763 jobs reproduce the committed follow-up row (k151f-n20weak-7) exactly in
  iterations and |r| for all 7 shared arms. The first run's Intel row (weight_views 166 iterations) could not be
  checked.

## Task 3: RUV-III with k >= 2 removed the shared within-group factor; the group factor grew

Declared baseline rule not met: the oracle arm recovers S in 14 of 20 (N = 20) and 11 of 20 (N = 60) seeds, under
the declared 15, so the declared reading is "the design cannot show removal". The measured effects:

| N | arm | S recovered | median abs r S | S kept in data | group factor VE % | R^2 vs group | -log10 p |
|---|---|---|---|---|---|---|---|
| 20 | raw / oracle | 11 / 14 | 0.93 / 0.96 | 1 / 0.95 | 14.2 / 19.3 | 0.88 / 0.95 | 9.1 / 12.3 |
| 20 | ruv k1 / k2 / k3 | 12 / 0 / 0 | 0.93 / 0.017 / 0.018 | 0.85 / 0.016 / 0.015 | 19.3 / 23.6 / 25.2 | 0.99 / 0.995 / 0.996 | 19.4 / 21.6 / 22.3 |
| 60 | raw / oracle | 12 / 11 | 0.93 / 0.91 | 1 / 0.98 | 12.1 / 15.9 | 0.86 / 0.84 | 25 / 24 |
| 60 | ruv k1 / k2 / k3 | 14 / 0 / 0 | 0.95 / 0.009 / 0.010 | 0.92 / 0.009 / 0.009 | 16.0 / 18.7 / 19.2 | 0.96 / 0.996 / 0.996 | 43 / 71 / 72 |

- k = 1 takes the stronger unwanted factor (U kept 5 to 7 %) and keeps S; k = 2 and 3 also remove S (about 1 % left
  in the data) in every seed, at both N. Group signal in the data is kept (>= 99.6 %).
- Group factor against oracle, paired: k = 2 raises variance explained in 20 of 20 seeds (x1.22 at N = 20, x1.18 at
  N = 60), R^2 +0.05 / +0.16, -log10 p +8.8 / +49. k = 1 leaves variance explained unchanged (x1.00) but raises R^2
  and -log10 p (+2.5 at N = 20, +12 at N = 60).
- RUVIII_C_Varying agrees with ruv::RUVIII to 9e-10 on complete data.
- Limits: one simulated design in which S is weaker than each omic's unwanted factor; raw and oracle recover S in
  only about 11 to 14 of 20 seeds (the group and S factors have equal strength and partly mix in MOFA2).
