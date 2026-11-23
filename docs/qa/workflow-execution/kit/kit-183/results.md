# KIT-183 simulation: results (counted run S, 45225854502, package fae77e5f)

Declaration `declaration.md` (committed alone at fae77e5f, before the run). Every table of the scorer is in `s/`.
2,400 of 2,400 replicates completed (800 per cell); determinism job (n = 12, replicates 1 to 10) identical to the
counted replicates in every table.

## Outcome: O1

| check | result |
|---|---|
| E1, camera, n = 12 (30 checks: C0, C1 at rho 0.1, 0.3, 0.6, C2; F_main, F_ind; three data sets) | holds; largest upper limit 0.0737 (proteomics C2 on F_ind: T = 0.0657, Monte Carlo SE 0.0041) |
| E1, camera, n = 20 | holds; largest upper limit 0.0711 (proteomics C1 rho 0.1 on F_ind: T = 0.0641) |
| E2, gsea, n = 12 (C0, C2) | holds; largest upper limit 0.0599 |
| E2, camera_pr, n = 12 (C0, C2) | holds; largest upper limit 0.0637 (matters only under O2) |
| E3, direction, n = 12 | camera 33,877 of 33,877 called planted sets in the expected direction; gsea 34,891 of 34,891 |
| Operating characteristic (measured SEs), P(cell passes E1) if every check's true rate were 0.050 / 0.055 | n = 12: 0.998 / 0.959; n = 20: 0.990 / 0.879 |

Consequence (ticket 6.5, O1): `camera` is the test of record and the default for every data set; `camera_pr` is not
offered; the few-samples sentence below 11 biological samples. (The code was written to O1; nothing changes.)

## The numbers written into the texts (5.5)

| name | value | from |
|---|---|---|
| R_cam | 6.6 % | camera, C2, proteomics, F_ind, n = 12 |
| I_cam | 5.8 to 7.4 % | the same |
| S_cam | 12.8 % | camera, C4b (rho 0.6), proteomics, F_none, n = 12 |
| R_gsea_0.3 | 43.0 % | gsea, C1 rho 0.3, transcriptomics, F_main, n = 12 |
| I_gsea_0.3 | 41.9 to 44.1 % | the same |
| R_cpr_0.3 | 29.9 % | camera_pr, C1 rho 0.3, transcriptomics, F_main, n = 12 |
| R_cam_6 | 5.0 % | camera, C1 rho 0.1, proteomics, F_ind, n = 6 |
| S_gsea | 60.1 % | gsea, C1-out rho 0.6, transcriptomics, F_none, n = 12, for `enrich_factor_strength` |

## Reported, not decided (n = 12; `s/error-rates.tsv` holds every rate with its interval)

- Overstatement on correlated null blocks (C1, F_main, largest over data sets; rho 0.1 / 0.3 / 0.6): gsea 0.237 /
  0.430 / 0.577; camera_pr 0.129 / 0.299 / 0.424; camera 0.010 / 0.006 / 0.005. Sensitivity arm (camera with
  `allow.neg.cor = FALSE`): 0.009 / 0.006 / 0.005, C0 at most 0.023 (camera 0.044).
- On weaker and noise factors (largest over data sets and rho): camera C1-out 0.033 (F_weak), 0.071 (F_none);
  C4 (model-selected, zero-loading features): 0.094 to 0.101 already on F_ind; 0.085 on F_weak and 0.109 on F_none;
  not attributable to circularity. C4b at rho 0.6: 0.042 on F_ind, 0.128 on F_none (circularity). gsea C1-out 0.588,
  0.601; C4b 0.545, 0.576.
- C3 (metabolomics identifiers carrying 2 or 3 null features): F_ind camera (largest-mean collapse) 0.046, camera with
  the |r| collapse 0.068, gsea (|r| collapse) 0.065; F_main 0.008, 0.013, 0.036.
- Comparators (VAL-53 (5), reported only): PCGSE correlation-adjusted C1 at most 0.003 on F_main, C0 at most 0.021,
  power 0 at every set size (inferred: the marginal correlation counts the factor's own pull as noise; the t on n_d - 2
  degrees of freedom and BH within the call also act); PCGSE permutation C1 at most 0.011 on F_main and 0.0525 on F_ind
  (metabolomics rho 0.1), C4b on F_none 0.127, power 0.39 to 1.00 by size; MOFA2's default (uncorrected, absolute
  weights) 64.5 to 69.8 % on the model's features of correlated null blocks (transcriptomics; proteomics 46 to 61 %,
  metabolomics 3 to 18 %; its universe is the model's features, so these are the C4b sets), C0 at most 0.073.
- Power (P-up, BH <= 0.05 and the expected direction, F_main): camera proteomics 0.61 (size 10), 0.99 (30);
  transcriptomics 0.60, 0.99, 1.00 (100); metabolomics size 5 0.02. gsea 0.63, 1.00; 0.72, 1.00, 1.00; 0.01.
- Family-wise (any BH <= 0.05 among C0 and C2): camera 0.019 to 0.028; gsea 0.026 to 0.071.
- Censoring losses (share of features lost to the complete-case filter, zero-loading / loaded): proteomics 0.33 /
  0.33 (n = 12), metabolomics 0.35 / 0.37; transcriptomics none.
- n = 6 (reported only): the largest camera rate over C0, C1, C2 on F_main and F_ind is 5.0 % (R_cam_6).
- Factor strength, MOFA2's messages at n = 6 and 12, counts, the diagnostic ratio: `s/factor-strength.tsv`,
  `s/mofa-warnings.tsv`, `s/counts.tsv`, `s/diagnostic.tsv`.

## Against the expectation stated before the run (ticket 6.4)

The ticket expected camera "near or slightly below 0.05 on every class on F_ind". It was not, for two classes.
Camera on F_ind (the covariate drawn independently of the data), n = 12:

| class | proteomics | transcriptomics | metabolomics |
|---|---|---|---|
| C4 (zero-loading, in no block, selected by the model) | 0.0939, upper 0.1018 | 0.1009, upper 0.1069 | 0.0448 |
| C4b rho 0.1 | 0.0687, upper 0.0765 | 0.0831, upper 0.0918 | 0.0431 |

At n = 20 proteomics C4 is 0.1035 (upper 0.1138); at n = 6 transcriptomics C4 is 0.099. Diagnostic ratio (within-set
variance of camera's z over its pooled variance), F_ind, n = 12: C4 1.19 (proteomics) and 1.41 (transcriptomics); C0
1.00; C2 0.96 to 0.97. On the same C4 sets: gsea 0.041 and 0.061; PCGSE permutation 0.049 and 0.053; camera with
`allow.neg.cor = FALSE` 0.055 and 0.054.

Mechanism (inferred, not tested): the model selects features whose observed variance is high for their abundance; the
moderated variance shrinks those variances towards the abundance trend, so their moderated t is too large given the
data; with `allow.neg.cor = TRUE` a set of independent features often gets a variance inflation factor below 1. The
floor halves C4's rate, which fits.

Over the 30 deciding checks at n = 12 camera's T ranges from 0.004 (metabolomics C1 rho 0.6, F_main) to 0.066
(proteomics C2, F_ind; 95 % interval 0.058 to 0.074, which excludes 0.05). At n = 20 two intervals exclude 0.05:
proteomics C1 rho 0.1 on F_ind (0.064, 0.057 to 0.071) and proteomics C2 on F_ind (0.061, 0.052 to 0.071). On F_main
every null class is at most 0.044. Plug-in pass probabilities (inferred): if the true rates equal the measured ones, a
repeat of S passes E1 with probability about 0.62 at n = 12 and about 0.68 at n = 20.

## Observed against the declaration

- Proteomics C2 was redefined to sizes 10 and 15 in every cell, n = 12 included, where size 30 had formed in L1, so
  that one class definition holds across cells; 15 is the largest size the n = 20 pool formed. The L1 artifact held
  class-level counts of 5 replicates and nothing by size. Metabolomics C2 size 10 at n = 20 formed 62 % of its sets;
  the class is decided on all 800 replicates, and formation depends on the censoring draw, not on the covariate tested.

- Sets formed (`s/sets-formed.tsv`): every class and size at or above 90 % except metabolomics C2 size 10 at n = 20
  (62 %; the L1 check's five replicates had not shown it). T is computed over the replicates with a set, as declared;
  metabolomics C2 at n = 20 is decided on its size-5 and size-10 sets (800 replicates).
- Time: median 51 s per replicate at n = 12 (training 18 s), 48 s at n = 20, 45 s at n = 6.

## Limits

One generator (6.7): real data may have stronger blocks, other missingness, other trends, other rates of
multi-feature identifiers. In the deciding classes the error rate stayed within the programme's tolerance on the
main-effect factor and on a covariate independent of the data, but above 5 % with intervals excluding 5 % in two
proteomics classes (6.6 % at n = 12; 6.4 % and 6.1 % at n = 20); sets of the model's most variable features were called
in up to 10.1 % (10.4 % at n = 20) on the independent covariate, a reported class; no error rate is claimed for weaker
or noise factors. Null blocks were formed among complete features and their values redrawn without censoring, so
correlated null sets in the simulation are not censoring survivors, while every other complete feature is; real
low-abundance correlated sets are censored too.
