# KIT-151 d, Q1: which RUV output should feed integration? Declaration

Fixed before the counted run. Authority: methods ledger INT-13 (user, 2026-11-21); audit SCI-04 (INT-004's required
planted-factor test). Branch `run-mofa-evidence-k151d`; lane `.github/workflows/k151d-evidence.yml`, script
`.github/mofa-smoke/k151d.R` (task `ruv`), installers `install-ruv-library.R` (RUVIIIC 1.0.19 from cran/RUVIIIC@c534151f,
ruv 0.9.7.1, the bundle's versions) and `install-sva-library.R` (sva 3.60.0, the version of the package lockfile at
e8a2ac67), scorer `.github/mofa-smoke/k151d-score.R` (task `ruv`). Script hashes: end of this file. Lifecycle:
run 44785248037 (commit fc5e8042): all three tasks, 8 jobs; every deliberately failing cell was written as failed and its job marked failed, every other cell trained to the end (0 at the round limit), the score job ran on all of them. One scorer wording defect was fixed after it (35c41ca3, pin history line); the evidence script and the lane are unchanged. Q1 lifecycle (4 jobs: strength C1 and C2 x N 20 and 60; ideal on data sets 91 to 100, every arm on 91 in both scenarios): the ideal arm recovered S in 10 of 10 data sets for both candidates at both N (median abs r 0.999), so the rule chose C1. On data set 91 every arm ran; the calibrated count gave k = 2 in both omics in every job; RUVIII_C_Varying against ruv::RUVIII largest difference 9.8e-10; sva 3.60.0 installed (Matrix 1.7-6 against the lockfile's 1.7-5, already in MOFA2's closure at the bundle's version; reported, not fatal). Noted: the share of U left in the ideal arm is not 0 (about 0.04 to 0.09): it is the noise floor of the coefficient on U after the U term is removed, so the U share is read against that floor. Evidence only: no package code, default or wording changes; the user decides.

## Purpose

Ledger M2: the default RUV design in every omic takes the groups as replicate sets and every feature as a negative
control; M10: the integration input is the RUV-corrected matrix; M5, M15, M30: k by the estimated count. k151c task 3
(run 44680786131) used the group-based design with k = 1, 2, 3 only, and its baseline was not met (the ideal
correction found the planted shared within-group factor S in 14 and 11 of 20 data sets; the failures were S and the
group factor rotating into each other, both of equal strength). This run strengthens S so that the ideal correction
finds it reliably, adds the package's calibrated default with its own k, more data sets, the round limit 10,000 and a
scenario where the unwanted factor is partly confounded with the groups.

## The calibrated default as the package defines it (e8a2ac67), and how it is reproduced

- Design: `runOmicsRuvFactors(ruv_design = "groups_all_features", ruv_k_rule = "num_sv_be")`
  (R/func_omics_ruv_factors.R:240-385, lipidomics and metabolomics; classic proteomics' `computeProtRuvModelFactors`,
  R/func_prot_norm_optimization_helpers.R:1274-1359, takes k the same way): the groups as replicate sets, every feature
  a control, k per omic from the shared resolver.
- k: `estimateRuvFactorCount` (R/func_general_ruv_factor_count.R:89-150): `sva::num.sv(dat, model.matrix(~ group),
  method = "be", B = 200, seed = 20261118)` on the features observed in every sample, under RNG kind
  Mersenne-Twister / Inversion / Rejection. Limits (`computeOmicsRuvAssayFactors`, R/func_omics_ruv_factors.R:137-178):
  k = 0 gives no correction; k at or above the number of controls, or above samples - groups - 1, is refused and the
  assay is analysed without RUV (automatic rule).
- Correction: `ruvIII_C_Varying` with the same replicate sets, k and every feature as a control, which calls
  `RUVIIIC::RUVIII_C_Varying(k, Y, M, toCorrect, potentialControls)` (R/func_prot_s4_norm_methods.R:248,
  R/func_metab_s4_ruv_methods.R:691); then each feature's level restored (`restoreOmicsRuvAssayLevels`).
- Reproduced here outside the package with the same underlying functions and arguments: `sva::num.sv` with the
  arguments above (sva 3.60.0 required, else the job stops) and `RUVIIIC::RUVIII_C_Varying` as in k151c (checked
  against `ruv::RUVIII` on each job's first data set; transcriptomics calls `ruv::RUVIII`). Not reproduced: the level
  restoration (a constant per feature; MOFA2 centres every feature, so it does not change the model). The package was
  not installed because both steps are calls of these two library functions with fixed arguments; the data here have
  no missing values, so "features observed in every sample" is every feature.

## Simulation

- Two omics: A (1,000 features), B (500 features); per-feature offset N(20, 2); noise sd 1; two groups of N / 2;
  N = 20 and N = 60; data seeds 154001 to 154030 (30 per N; same index, same draws for both scenarios).
- Latent factors, each loaded with weight +/- w on a random share of features (independent draws):
  G the group (+1 / -1), 20% of features; S the planted shared within-group factor, N(0, 1) centred within each group
  and scaled to sd 1 (exactly unrelated to the group in the sample), 20% of features in A and in B; UA and UB, one
  unwanted factor per omic, 40% of that omic's features.
- Strength (chosen by the lifecycle under the rule fixed before it, see below): C1, weights G 1, S 1.5, U 2 (per-feature variance from the factor on the loaded features: G 1, S 2.25, U 4; noise 1). S is stronger than the
  group effect (so the two do not rotate into each other in MOFA2) and weaker than each omic's unwanted factor (as in
  k151c: RUV with every feature a control is expected to take the unwanted factor first, S next).
- Scenarios: `main`, UA and UB within-group noise like S (unrelated to the group); `conf`, each unwanted factor built
  as 0.6 x (group / sd) + 0.8 x (its within-group noise), so its correlation with the group is exactly 0.6 in the
  sample (KIT-124's high association), the same draws otherwise.

## Arms

- `raw`: uncorrected.
- `ideal`: each omic's planted unwanted term (its loadings times its factor) subtracted exactly (the data as if no
  unwanted variation existed; in `conf` this also removes the part of U that follows the group, which no data-driven
  method can tell apart from the group effect).
- `ruv_k1`, `ruv_k2`, `ruv_k3`: RUV-III, groups as replicate sets, every feature a control, k fixed per omic.
- `calibrated`: the package's default as above, k per omic by the estimated count (recorded).

MOFA2 as the integration default (INT-7, INT-12, INT-15): `scale_views` on, no z-score, 5 factors, double precision,
`center_groups` TRUE, convergence mode slow, round limit 10,000, training seed 20261119. Feature selection does not
bind here (A has exactly 1,000 features, B 500, no missing values). 2 N x 2 scenarios x 30 data sets x 6 arms =
720 trainings.

## Measures (per training)

- Recovery: |r| of S, G, UA, UB with their greedily matched learned factors; S recovered = |r| >= 0.9.
- Share of planted signal left in the corrected data: per omic, the joint least-squares coefficients of each feature
  on (G, S, U); mean |coefficient| over the features loaded on that factor, corrected over uncorrected; reported for S,
  G and U (averaged over A and B).
- The group factor (matched to G): variance explained in A and B (MOFA2's r2, of all features), R^2 and the F-test
  p-value of the factor against the group; |r| with G (recovery of the group effect).
- Calibrated arm: k per omic and whether it was refused; iterations; stopped by the rule or at the limit.

## Analysis and readings (declared)

- Per N, scenario and arm: S recovered (count of 30), median and mean |r| of S with a 95% seed-bootstrap interval,
  the medians of the other measures.
- Paired by data set against `ideal` and against `raw`: mean difference with a 95% seed-bootstrap interval (20,000
  resamples, `set.seed(151)`). |r| of S and of G: interval inside [-0.02, +0.02] "no practical difference", lower
  bound > 0 "better", upper < 0 "worse", else "inconclusive". Group variance explained, R^2, -log10 p and shares left:
  "higher" / "lower" when the interval excludes 0, else "inconclusive" (no margin; magnitudes reported).
- Baseline: the ideal arm must recover S in at least 24 of 30 data sets at a given N and scenario; otherwise "the
  design cannot show removal" there and only the measured effects are reported.
- Per arm, where the baseline holds: "removes the within-group factor" if S is recovered in at most 6 of 30;
  "keeps" if within 3 of the ideal count; otherwise "partly removes".
- These readings describe this simulated design only (one strength ordering, Gaussian, linear, two groups).

## What the possible outcomes would mean

- Calibrated (and k >= 2) "removes" S while k = 1 "keeps" it, with the estimated k mostly 2: the default correction
  removes shared within-group biology whenever it is among the largest within-group axes; the integration input route
  (corrected or uncorrected or user-chosen) is then a trade-off the user decides, and the Methods must say that the
  corrected input can lose within-group structure.
- Calibrated "keeps" S (k mostly 1) in this design: the estimated count stops before S; not a general guarantee, since
  it depends on S being weaker than U and on the count's threshold.
- Raw recovers S as well as ideal and the group factor is unharmed: correction is not needed for integration in the
  unrelated scenario.
- In `conf`: if raw's group factor is distorted (|r| with G lower, variance explained higher, U's group-aligned part
  absorbed) and the RUV arms do not repair it (they cannot see the group-aligned part of U), correcting does not help
  where it would be needed most; if RUV arms move the group factor towards ideal, correction helps there.
- Group variance explained, R^2 and p higher than ideal under RUV: the corrected input inflates the group factor
  (SCI-04's second effect), magnitude reported.

## Budget

12 jobs (N x scenario x 3 blocks of 10 data sets), each the installs plus 60 trainings; job timeout 150 minutes.
Ceiling 25 core-hours (job wall time x 4 cores); expected from the lifecycle (about 2.5 s per training at N = 20, 9 s at N = 60, installs about 4 minutes), about 7 to 9 core-hours. One counted run; one correction pass at most,
for an instrumentation failure only, rerunning only the affected jobs.

## Script hashes (sha256, fixed for the counted run)

```
8fe56e9ad458fce63443eafe532e9cd47439ca4de73f020f9ae1885024eb6f7e  .github/mofa-smoke/k151d.R
6261544a4ab84361fcfcb6fc7dbc3b0ee19040a61db5520abc3c47cbbc73ba4e  .github/mofa-smoke/k151d-score.R
04c250d48e264f120042f4e3030388fc308480e62c0d9f2ad13a79c577b0eed4  .github/mofa-smoke/install-sva-library.R
405d84a1b78f7cdf2324ef7d85dc309539e88a1264a1ad2a29ca448c0bea2d17  .github/mofa-smoke/install-ruv-library.R
f4caa249375b225ce1d1cf012b6f3c52d4518ee554e465dfc68159274a1051bd  .github/workflows/k151d-evidence.yml
```
