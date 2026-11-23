# KIT-183 simulation: results (counted run S, 45225854502, package fae77e5f)

Declaration `declaration.md` (committed alone at fae77e5f, before the run). Every table of the scorer is in `s/`.
2,400 of 2,400 replicates completed (800 per cell); determinism job (n = 12, replicates 1 to 10) identical to the
counted replicates in every table.

## Outcome: O1

| check | result |
|---|---|
| E1, camera, n = 12 (30 checks: C0, C1 at rho 0.1, 0.3, 0.6, C2; F_main, F_ind; three data sets) | holds; largest upper limit 0.0737 (proteomics C2 on F_ind: T = 0.0657, Monte Carlo SE 0.0041) |
| E1, camera, n = 20 | holds; largest upper limit 0.0710 (proteomics C1 rho 0.1 on F_ind: T = 0.0641) |
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

## Reported, not decided (n = 12; `s/error-rates.tsv` holds every rate with its interval)

- Overstatement on correlated null blocks (C1, F_main, largest over data sets; rho 0.1 / 0.3 / 0.6): gsea 0.237 /
  0.430 / 0.577; camera_pr 0.129 / 0.299 / 0.424; camera 0.010 / 0.006 / 0.005. Sensitivity arm (camera with
  `allow.neg.cor = FALSE`): 0.009 / 0.006 / 0.005, C0 at most 0.023 (camera 0.044).
- Circularity on weaker and noise factors (largest over data sets and rho): camera C1-out 0.033 (F_weak), 0.071
  (F_none); C4 0.085, 0.109; C4b 0.094, 0.128. gsea C1-out 0.588, 0.601; C4b 0.545, 0.576.
- C3 (metabolomics identifiers carrying 2 or 3 null features): F_ind camera (largest-mean collapse) 0.046, camera with
  the |r| collapse 0.068, gsea (|r| collapse) 0.065; F_main 0.008, 0.013, 0.036.
- Comparators (VAL-53 (5), reported only): PCGSE correlation-adjusted C1 at most 0.003 on F_main, C0 at most 0.021,
  power 0 at every set size (the marginal correlation counts the factor's own pull as noise); PCGSE permutation C1 at
  most 0.011 on F_main, 0.048 on F_ind (rho 0.3, transcriptomics), C4b on F_none 0.127, power 0.39 to 1.00 by size;
  MOFA2's default (uncorrected, absolute weights) C1 0.645 / 0.648 / 0.698 on F_main, C0 at most 0.073.
- Power (P-up, BH <= 0.05 and the expected direction, F_main): camera proteomics 0.61 (size 10), 0.99 (30);
  transcriptomics 0.60, 0.99, 1.00 (100); metabolomics size 5 0.02. gsea 0.63, 1.00; 0.72, 1.00, 1.00; 0.01.
- Family-wise (any BH <= 0.05 among C0 and C2): camera 0.019 to 0.028; gsea 0.026 to 0.071.
- Censoring losses (share of features lost to the complete-case filter, zero-loading / loaded): proteomics 0.33 /
  0.33 (n = 12), metabolomics 0.35 / 0.37; transcriptomics none.
- n = 6 (reported only): the largest camera rate over C0, C1, C2 on F_main and F_ind is 5.0 % (R_cam_6).
- Factor strength, MOFA2's messages at n = 6 and 12, counts, the diagnostic ratio: `s/factor-strength.tsv`,
  `s/mofa-warnings.tsv`, `s/counts.tsv`, `s/diagnostic.tsv`.

## Observed against the declaration

- Sets formed (`s/sets-formed.tsv`): every class and size at or above 90 % except metabolomics C2 size 10 at n = 20
  (62 %; the L1 check's five replicates had not shown it). T is computed over the replicates with a set, as declared;
  metabolomics C2 at n = 20 is decided on its size-5 and size-10 sets (800 replicates).
- Time: median 51 s per replicate at n = 12 (training 18 s), 48 s at n = 20, 45 s at n = 6.

## Limits

One generator (6.7): real data may have stronger blocks, other missingness, other trends, other rates of
multi-feature identifiers. Calibration is shown for a covariate independent of the data and on the main-effect factor;
on weaker or noise factors the excess is reported as measured and no error rate is claimed.
