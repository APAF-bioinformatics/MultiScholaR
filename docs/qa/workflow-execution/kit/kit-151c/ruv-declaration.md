# KIT-151 c, task 3 (SCI-04, INT-004's required test): group-corrected inputs. Declaration

Fixed before the counted run. Branch `run-mofa-evidence-k151c`; lane `.github/workflows/k151c-evidence.yml`, script
`.github/mofa-smoke/k151c.R` (task `ruv`), RUV installer `.github/mofa-smoke/install-ruv-library.R`, scorer
`.github/mofa-smoke/k151c-score.R` (task `ruv`). Script hashes: see the end of this file. Lifecycle run:
44677064790 (all three tasks; deliberately failing cells reported as failed; two scorer defects fixed in 1 commit after it, evidence script unchanged).

## Purpose

Ledger M2: the default RUV design in all four omics uses the groups as replicate sets and every feature as a
negative control. M10 / INT-004: MOFA2 reads the corrected matrix. INT-004 requires "a planted-factor test on
corrected and uncorrected inputs". The audit (SCI-04) infers that this correction removes each omic's largest
within-group axes, including shared within-group biology, and inflates the group factor. This run measures it on
one simulated design. No default, wording or package code changes here.

## RUV-III as the package calls it

- Proteomics, lipidomics and metabolomics call `RUVIIIC::RUVIII_C_Varying(k, Y, M, toCorrect = colnames(Y),
  potentialControls = <controls>)` with Y samples by features and M the 0/1 replicate matrix of the groups
  (`getRuvIIIReplicateMatrixHelper`; R/func_prot_s4_norm_methods.R:248). Transcriptomics calls
  `ruv::RUVIII(Y, M, ctl, k)` (R/func_trx_ruv.R:178).
- This run calls `RUVIIIC::RUVIII_C_Varying` directly with every feature in `toCorrect` and `potentialControls`, M
  from the groups, without installing MultiScholaR. RUVIIIC 1.0.19 from the bundle's source
  (github.com/cran/RUVIIIC at c534151f), its dependencies and ruv 0.9.7.1 at the bundle's versions (checked
  against `kit-134/p2/linux/bundle-sources.tsv`). With complete data RUV-III-C reduces to RUV-III; the job checks
  `RUVIII_C_Varying` against `ruv::RUVIII` on its first data set for each k and records the largest difference.
- Not reproduced: proteomics' `restoreRuvFeatureLevels` (adds back each feature's level; MOFA2 centres every
  feature, so it does not change the model), and the package's choice of k (here k is set to 1, 2 and 3).

## Simulation (declared)

- N = 20 and N = 60, two groups of N/2; data seeds 153001 to 153020 (20 per N).
- Two omics: A (1,000 features), B (500 features). Per-feature offset N(20, 2); noise sd 1.
- Latent factors, each loaded with weight +1 or -1 on a random share of features (independently drawn sets):
  - G, the group (+1 / -1): 20 % of features in A and in B;
  - S, the planted shared within-group factor: N(0, 1), centred within each group and scaled to sd 1 (so exactly
    independent of the group in the sample), 20 % of features in A and in B;
  - UA and UB, unwanted factors, one per omic (technical), drawn and centred like S, 40 % of that omic's features.
- So S is weaker than each omic's unwanted factor; RUV-III with every feature as a control is expected to take the
  unwanted factor first (k = 1) and S next (k = 2). The result depends on this ordering, which is a design choice.

## Arms

- raw: uncorrected.
- oracle: each feature's component on its omic's true unwanted factor removed by least squares (the "known batch
  regressed" reference).
- ruv_k1, ruv_k2, ruv_k3: RUV-III as above with k = 1, 2, 3 per omic.

MOFA2 as k151 with `scale_views` on, 5 factors, double precision, convergence `slow` (the follow-up rule's mode),
maxiter 1,000, training seed 20261119. 2 N x 20 seeds x 5 arms = 200 models.

## Measures

Planted factors G, S, UA, UB matched greedily to MOFA factors by |r|. Per model: |r| for S (and the best |r| of S
with any factor); the group-aligned factor (matched to G): variance explained in A and B, R^2 against the group and
the F-test p-value; data-level retention in the corrected matrix (mean |slope| on S, G and U over the features
loaded on it, corrected over raw).

## Readings (declared)

- S recovered: matched |r| >= 0.9 (the k151 mark). Baseline: the oracle arm must recover S in >= 15 of 20 seeds at a
  given N, otherwise the design cannot show removal at that N (stated).
- Per k: "removes the within-group factor" if S is recovered in <= 5 of 20 seeds; "keeps" if within 2 seeds of the
  oracle count; otherwise "partly removes".
- Group factor, paired by data seed against oracle (and, reported, against raw): "inflates the group factor's
  variance explained" if it is higher in >= 15 of 20 seeds, with the median ratio and the median change in R^2 and
  -log10 p; otherwise "no inflation shown".
- These readings describe this simulated design only (one strength ordering, Gaussian, linear, two groups).

## Budget

8 jobs (N x 4 blocks of 5 seeds), each 25 trainings plus the RUV install; ceiling 60 minutes wall per job (expected
about 10 minutes each, about 5 core-hours in all). One counted run; one correction pass at most, only for an
instrumentation failure.

## Script hashes (sha256, fixed for the counted run)

```
3d97c5f55a56c9415a3d1257bbe736350ccdce5a31dace25ef54f52da3d97afa  .github/mofa-smoke/k151c.R
84083aa17d79227309f92da2b9837ebea715d42f6f8ba369ff66ad2b38cd0b12  .github/mofa-smoke/k151c-score.R
405d84a1b78f7cdf2324ef7d85dc309539e88a1264a1ad2a29ca448c0bea2d17  .github/mofa-smoke/install-ruv-library.R
88084c4062a4b70b7a8a100ceec0198ca7a64c455ab37105ccca4e3ae37ccfad  .github/workflows/k151c-evidence.yml
```
