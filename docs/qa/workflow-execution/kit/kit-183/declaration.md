# KIT-183 simulation: declaration (fixed before the counted run)

Ticket `tickets/kit/KIT-183.md` revision 1.4, sections 6.2 to 6.5; owner decisions VAL-51, VAL-53. This file is
committed alone, before the counted run S; a change after the lifecycle run L1 is a new dated section below, made
before S. The code that implements it: `simulate-k183.R` (simulation), `score-k183.R` (scorer), lane
`.github/workflows/k183-evidence.yml` on branch `run-mofa-evidence-k183` (the lane's `.github/mofa-smoke/k183.R` is a
byte copy of `simulate-k183.R`, checked by the job).

## Vehicle and cells

- MOFA2 evidence lane: R 4.6.1, the KIT-134 snapshot library, MultiScholaR at the package commit named in the lane's
  `k183-params.json` (the implementation's head), MOFA2 1.22.0 and mofapy2 0.7.3 with the bundle's Python, reference
  threads (`OPENBLAS_NUM_THREADS = OMP_NUM_THREADS = MKL_NUM_THREADS = 1`), four single-thread workers per job.
- Cells: n = 12 biological samples (deciding), n = 20 (deciding, for the sample-size rule), n = 6 (reported only); two
  groups of n/2; one run per biological sample. 800 replicates per cell. In the n = 12 cell the transcriptomics data
  set lacks the last biological sample (kept as missing, `keep_missing`), so its tests run at n_d = 11.
- Jobs (S; re-split after L1, VAL-64): 12 jobs of 200 replicates (each cell's replicates 1-200, 201-400, 401-600,
  601-800), timeout 75 minutes each, no replicate started in the last 6 minutes (recorded `not_started`); determinism
  job: n = 12 replicates 1-10 again, timeout 20 minutes.

## Generator (per replicate)

Seed `set.seed(20261183 + r)` for n = 12, `+ 100000 + r` for n = 20, `+ 200000 + r` for n = 6; one stream, drawn in
the order: group scores, Z2, then per data set (transcriptomics, proteomics, metabolomics) abundances, loadings,
noise, censoring, blocks.

| | transcriptomics | proteomics | metabolomics |
|---|---|---|---|
| features | 4,000 | 1,800 | 250 (on 204 identifiers) |
| abundance mu_f | N(5, 2.5^2) | N(20, 2^2) | N(16, 2^2) |
| noise sd | 0.25 + 1.0 exp(-(mu - 0)/2) | 0.15 + 0.6 exp(-(mu - 14)/2) | 0.15 + 0.6 exp(-(mu - 14)/2) |
| Z1 loadings | 10 % of loading units, half +, half -, abs ~ U(0.5, 1.5) sigma_f | same | same |
| Z2 loadings | 5 % of the other units, random sign, abs ~ U(0.3, 1.0) sigma_f | same | none |
| censoring | none | P(missing) = 1/(1 + exp(1.5 (y - 17))) | P(missing) = 1/(1 + exp(1.5 (y - 13))) |
| null blocks per rho and size | 3 | 2 | 1 |
| block sizes | 10, 30, 100 | 10, 30, 100 | 5, 10 |
| set sizes (C0) | 10, 30, 100 | 10, 30, 100 | 5, 10, 20 |
| C4 sizes | 10, 30, 100 | 10, 30 | 5, 10, 20 |
| P sizes (each direction) | 10, 30, 100 | 10, 30 | 5 |
| C2 sizes | 10, 30 | 10, 15 | 5, 10 |

- Z1 = -1 (group A) or +1 (group B) plus N(0, 0.25^2); Z2 ~ N(0, 1) residualised on Z1 in the sample, scaled to sd 1.
- Values y = mu_f + sigma_f (l1 Z1 + l2 Z2 + e), e ~ N(0, 1); loadings in units of sigma_f.
- Metabolomics identifiers: 31 of 204 identifiers carry two or three features (alternately 2 and 3; 77 features);
  the first 16 are independent features (each its own Z1 loading, active with probability 0.10), the last 15 adducts
  of one compound (one loading unit shared, noise correlated at 0.9). Identifiers `C1xxxx` in a `kegg` column;
  transcriptomics identifiers `G00001...` (gene), proteomics `P00001...` (UniProt-shaped).
- Censoring is drawn on the values before the block latents; a feature is complete when it has no missing value.
- Null blocks: from the complete features with zero loading on Z1 and Z2 (single-feature identifiers), in random order,
  disjoint, in the order rho 0.1, 0.3, 0.6, then size, then block; each member's noise becomes
  sqrt(rho) u_b + sqrt(1 - rho) e (variance unchanged; u_b ~ N(0, 1)); a block that the pool cannot fill is not formed
  and counted. The mask is not redrawn (the marginal distribution of the values is unchanged).
- C2 pool: the lowest-abundance fifth (by mu_f) of the complete zero-loading features in no block.

## The shipped path

Per replicate the three data sets are written with `writeIntegrationInput()`, imported into a new integration store and
trained with `runIntegrationModel()` at every shipped default except `training_seeds = 1` (F3). Factors (each oriented
by the shipped rule): F_main = largest |cor| of the oriented scores with Z1; F_weak = among the others, largest |cor|
with Z2; F_none = among the remaining, smallest max(|cor Z1|, |cor Z2|); F_ind = one N(0, 1) value per biological
sample, `set.seed(20261187 + r)` (no model; its direction sign is +1 since it has no weights). Per data set and factor:
`integrationEnrichmentRank()` (Pearson, 5 samples), `integrationEnrichmentMappedList()`,
`integrationEnrichmentTestValues()`, then `integrationEnrichmentPathwayTest()` for each test on one synthetic set table
holding every class's sets; for metabolomics on a learned factor, camera and gsea through `runIntegrationEnrichment()`
with a synthetic KEGG table (the whole path).

## Tests and comparators

- Offered and fallback: `camera` (the test of record), `gsea` (seed 20261121), `camera_pr`; sensitivity arm
  `camera_floor` (camera with `allow.neg.cor = FALSE`); C3 arm `camera_rcollapse` (camera on the identifiers collapsed
  by the largest absolute correlation among complete features, the pre-ranked test's rule), C3 sets only.
- Comparators (VAL-53 (5); reported only): `pcgse_cor` (Fisher z of each U_cam feature's correlation with the oriented
  factor, set mean difference with pooled sd, VIF = 1 + (m - 1) mean marginal correlation of the set's features, t on
  n_d - 2 df, two-sided); `pcgse_perm` (the factor's scores permuted across S_d 1,000 times, `set.seed(20261185 + r)`,
  each feature's Fisher z recomputed, the uncorrected t's |t| against its permutation distribution,
  p = (1 + #{|t*| >= |t|}) / 1,001); `mofa_default` (MOFA2 1.22.0's default, read in its `R/enrichment.R`: the uncorrected
  t of the mean difference of the absolute weights of the model's features, pooled sd, df m1 + m2 - 2, sets with at
  least 10 model features; learned factors only).
- Settings of every call: sizes 5 to 500 members, BH within the call (one call per data set, factor and test).

## Classes (members are identifiers of U_cam; drawn per replicate, data set and factor)

C0 (10 per size), C1 (each formed block; split by rho), C1-out (block members the model did not select, at least 5;
transcriptomics and proteomics), C2 (5 per size), C3 (metabolomics: identifiers carrying 2 or 3 zero-loading features,
5 per size 5 and 10), C4 (zero-loading complete features in no block that the model selected, 5 per size), C4b
(block members the model selected, at least 5), P-up and P-down (half planted with positive or negative Z1 loading,
complete, single-feature identifiers; half complete zero-loading features in no block; 5 per size and direction).

## Estimands, criteria and outcomes (as the ticket's 6.3 to 6.5)

- s_r = share of the replicate's sets of the class with unadjusted p <= 0.05; T = mean of s_r over the replicates with
  a set; interval T +/- t(R - 1) sd(s_r)/sqrt(R); Monte Carlo error sd(s_r)/sqrt(R).
- E1: camera, upper limit <= 0.075 on F_main and F_ind for C0, C1 (each rho) and C2 in every data set (30 checks) at
  n = 12 and n = 20. E2: gsea and camera_pr, C0 and C2, F_main and F_ind, every data set, n = 12. E3: direction
  accuracy >= 0.99 for camera and gsea at n = 12 among planted sets called at BH <= 0.05 on F_main, decided when at
  least 50 were called; the expected direction is planted sign x sign(orientation sign x cor(F_main, Z1)), the
  orientation sign computed by the simulation's own code from the weights (not the kernel's).
- Outcomes O1 to O5 exactly as ticket 6.5; a cell with more than 5 % of its replicates not completed is not decided.
- Reported: the C1 rates of gsea and camera_pr; C3; C4, C4b on every factor; every class on F_weak and F_none; the
  diagnostic ratio (mean within-set variance of camera's z over camera's pooled variance); power per test, size,
  direction and data set on F_main; the sensitivity arm; n = 6; family-wise (any BH <= 0.05 among the C0 and C2 sets,
  BH over them, Clopper-Pearson interval); censoring losses (share of zero-loading and of loaded features lost); factor
  strength (r2_percent of each factor kind and of the largest factor, median and range); the comparators.
- The numbers written into the texts (5.5): R_cam, I_cam, S_cam, R_gsea_0.3, I_gsea_0.3, R_cpr_0.3, R_cam_6, as the
  scorer's `placeholders.tsv` gives them.

## Attempts

One counted run S; one instrumentation correction S' only if S's instrumentation fails (not for an unfavourable
result). An unfavourable valid result is recorded and its outcome applied.

## Dated section 2026-11-24 (after the lifecycle runs, before S)

- Lifecycle: runs 45148992728 and 45151737782 (instrument failures of the script's last aggregation, located and
  corrected) and the check 45219821547 (15 of 15 replicates completed; 57 to 62 s per replicate under four workers).
  The stop rule of ticket 11 fired at 55 to 70 s; the owner accepted the re-split above (VAL-64).
- Classes with fewer than 90 % of their sets formed in the check, redefined (table above): proteomics C2 sizes 10, 15
  (size 30 formed 0 % at n = 20); proteomics C4 and P sizes 10, 30 (size 100: C4 0 to 60 %, P 80 %); metabolomics P
  size 5 (sizes 10 and 20: 0 to 80 %). No deciding check changes: C2 stays deciding in every data set.
- The simulation script and scorer are those of the package commit named in the lane's `k183-params.json` for S
  (`simulate-k183.R`, `score-k183.R`); the script's generator, training step, classes and arms are lists passed as data
  (KIT-184's request), which changes how it is organised, not what it computes.

