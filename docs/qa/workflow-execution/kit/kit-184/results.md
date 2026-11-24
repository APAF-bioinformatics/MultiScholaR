# KIT-184 simulation: results (counted run S, 45353248796, package 81ee3052)

**Lifecycle disclosure: L1 (45333818847) was scored on real truth at five replicates per cell before this declaration;
what was seen and what was changed afterwards are in `declaration.md` sections 1, 2 and 4.**

Declaration `declaration.md` (committed alone at 81ee3052, before the run; the owner's decisions on it are methods
ledger VAL-150). Every table of the scorer is in `s/`. The run's jobs read the declared scripts (`simulate-k184.R`
58333b1e..., `score-k184.R` 282b844e..., `k184-evidence.yml` 15378427...). 2,000 of 2,000 replicates completed (cells
A and B 800 each, replicates 6 to 805; cell D 400, replicates 6 to 405); none failed, none was left unstarted. All 17
evidence jobs succeeded at the first attempt, so there was no S' and no re-score. Nothing below comes from L1 or L1b.

## Outcome: R1, with R7; determinism not met as declared (to be ruled on, below)

| Check | Result |
|---|---|
| **K-E1**, arm 3 (`camera`, removed factors in the design), cell A: 33 checks, upper 95 % limit at most 0.075 | **holds, 33 of 33.** Largest upper limit 0.0639 (proteomics C5 on F_ind, every size formed: T = 0.0545, Monte Carlo SE 0.0048). T over the checks: 0.003 to 0.055; on F_main at most 0.045 |
| **K-E1**, cell B | **holds, 33 of 33.** Largest upper limit 0.0593 (metabolomics C1 rho 0.1 on F_ind: T = 0.0488, SE 0.0054). T: 0.004 to 0.049 |
| Proteomics C5 on F_ind in both compositions (VAL-150) | holds in both, in both cells. Cell A: sizes 10 and 30 T = 0.0517 (SE 0.0043, upper 0.0601); every size formed T = 0.0545 (SE 0.0048, upper 0.0639). Cell B: 0.0432 (0.0047, 0.0525); 0.0464 (0.0052, 0.0566) |
| **K-E3**, direction accuracy of arm 3 (at least 0.99) | **holds.** Cell A 30,606 of 30,606 called planted sets in the expected direction; cell B 26,160 of 26,160 |
| **K-E4** (declared; does not hold the merge), arm 3 on C5 on F_main, cell A | holds in every data set: transcriptomics T = 0.0411 (SE 0.0036, upper 0.0481); proteomics 0.0086 (0.0018, 0.0121) and 0.0095 (0.0021, 0.0137) in the two compositions; metabolomics 0.0199 (0.0032, 0.0261) |
| **K-E4**, cell B (the unwanted factor tied to the groups at 0.6) | **fails in every data set**: transcriptomics T = 0.2956 (SE 0.0090, interval 0.278 to 0.313); proteomics 0.1123 (0.0074, 0.098 to 0.127) and 0.1328 (0.0083, 0.116 to 0.149); metabolomics 0.0807 (0.0071, 0.067 to 0.095) |
| **K-E5**, arm 2 (the form without the adjustment), the same 33 checks | **holds, 33 of 33 in both cells.** Largest upper limit 0.0606 in cell A (proteomics C2 on F_ind: T = 0.0539, SE 0.0034) and 0.0544 in cell B (metabolomics C1 rho 0.1 on F_ind: T = 0.0444, SE 0.0051) |
| K-E2 (`camera_pr`; decides nowhere under KIT-183's O1) | reported: C0 and C2 at most 0.057 (upper 0.068) in cell A and 0.042 in cell B; C1 up to 0.578 |
| Cells decided | A 800 of 800, B 800 of 800, D 400 of 400 |
| R6 | not met: no cell undecided; the written factors `exact` in 2,000 of 2,000 replicates of every data set; the identity check at most 6.2e-14 (none above 1e-8) |
| Operating characteristic (measured SEs): P(a cell passes K-E1) if every check's true rate were 0.050 / 0.055 / 0.060 | cell A 0.994 / 0.908 / 0.398; cell B 0.992 / 0.885 / 0.346 (34 rows per cell: the proteomics C5 check twice) |

**What the outcome rows of ticket 6.5 say.**

- **R1 holds** (K-E1 in A and B, K-E3, every cell decided): the test with the removed factors in its design ships as
  designed, and the measured numbers go into the texts.
- **R7 applies** (K-E4 fails in a data set): reported to the owner before merge, with arms 1, 2 and 3o on C5 and what
  the removed factors caught (next section). By the ticket it does not hold the merge (VAL-34), and the Methods prints
  [R5] and [I5] whatever K-E4 gave.
- **R8 does not apply** (K-E5 holds): the form without the adjustment stays offered.
- **R4 does not apply**: on F_ret the interval of T1 minus T3 excludes 0 in 29 of the 30 deciding rows (C0, C1, C2;
  cells A and B), arm 1 above arm 3 by 0.003 to 0.014. So at the set level too the test of corrected values without
  the factors called measurably more.
- **R5** ("arm 2 shows no loss of power against arm 3, or no excess on C5 in cell B": reported as measured, the Methods
  sentence on that form worded from the numbers). Arm 2 does lose power (below). On C5 on the main factor in cell B it
  is above 5 % in two data sets (proteomics 0.112, metabolomics 0.101) and below in transcriptomics (0.020); in two of
  the three data sets it is not above arm 3. I report the numbers and do not decide which half of R5's sentence that
  is; what it means for the two sentences that compare the forms is in `wording.md`.

## K-E4 in cell B: what the ticket asks to be reported with it (R7)

Sets of features that share unwanted variation and have no link to the simulated biology (C5), on the model's main
factor, `camera`, T (Monte Carlo SE). Proteomics: sizes 10 and 30 / every size formed.

| Cell | Data set | Arm 3 (removed factors in the design) | Arm 2 (no adjustment) | Arm 1 (corrected values, no factors) | Arm 3o (the planted factors in the design) |
|---|---|---|---|---|---|
| A | transcriptomics | 0.041 (0.004) | 0.000 | 0.055 (0.004) | 0.013 (0.001) |
| A | proteomics | 0.009 (0.002) / 0.010 (0.002) | 0.000 / 0.000 | 0.011 / 0.012 | 0.016 / 0.022 |
| A | metabolomics | 0.020 (0.003) | 0.000 | 0.024 (0.003) | 0.014 (0.003) |
| B | transcriptomics | **0.296 (0.009)** | 0.020 (0.003) | 0.345 (0.009) | 0.021 (0.001) |
| B | proteomics | **0.112 (0.007) / 0.133 (0.008)** | 0.073 (0.006) / 0.112 (0.007) | 0.176 (0.009) / 0.206 (0.010) | 0.021 (0.002) / 0.028 (0.003) |
| B | metabolomics | **0.081 (0.007)** | 0.101 (0.009) | 0.119 (0.009) | 0.016 (0.003) |
| D | transcriptomics | 0.050 (0.005) | 0.000 | 0.092 (0.007) | 0.015 (0.002) |
| D | proteomics | 0.009 / 0.011 | 0.000 / 0.000 | 0.019 / 0.021 | 0.013 / 0.015 |
| D | metabolomics | 0.020 (0.004) | 0.000 | 0.032 (0.005) | 0.010 (0.003) |

- With the planted factors in the design (arm 3o) these sets are called in 1 to 3 % everywhere. With the factors the
  correction actually removed (arm 3) they are called in 8 to 30 % when the unwanted factor follows the groups (cell
  B), and in 1 to 5 % when it does not (cells A and D).
- In cell B the adjusted test calls them more often than the form without the adjustment in transcriptomics (0.296
  against 0.020; paired difference arm 2 minus arm 3 -0.275, interval -0.292 to -0.259) and in proteomics, and less
  often in metabolomics (0.081 against 0.101).
- By set size the rate rises: transcriptomics, cell B, arm 3: 0.118 (size 10), 0.282 (30), 0.486 (100); proteomics
  0.091, 0.133, and 0.221 at size 100 over the 438 replicates in which that size formed
  (`s/error-rates-by-size.tsv`).
- What the removed factors caught (the smallest canonical correlation between the written and the planted unwanted
  factors, median over replicates; smallest replicate): transcriptomics 0.991 (0.008), 0.985 (0.202), 0.996 (0.799) in
  cells A, B, D; proteomics 0.917 (0.001), **0.846** (0.010), 0.984 (0.456); metabolomics 0.977 (0.022), 0.957 (0.003),
  0.985 (0.647). The share of the main factor that lies in the removed factors (R2, median): about 0.03 to 0.05 in
  cell A, 0.25 to 0.29 in cell B, 0.06 to 0.11 in cell D (`s/caught.tsv`).
- Inferred, not shown: the removed factors describe the planted ones well on the whole (medians above), so what is
  left is the part of the unwanted factor that runs with the groups, which a correction that takes the groups as its
  replicate sets cannot tell from the group difference and leaves in the corrected values the model was trained on.
  The simulation was not built to separate that from other causes.

## The numbers written into the texts (ticket 5.5 and 6.4; `s/numbers.tsv`)

| Name | Value | From |
|---|---|---|
| [R3] | 5.4 % | arm 3, `camera`, cell A, proteomics, F_ind, C5, every size formed (the largest T over K-E1's checks in cells A and B) |
| [I3] | 4.5 to 6.4 % | the same |
| [R5] | 29.6 % | arm 3, cell B, transcriptomics, F_main, C5 |
| [I5] | 27.8 to 31.3 % | the same |
| [R6] | 11.2 % | arm 2, cell B, proteomics, F_main, C5, every size formed (the largest T of arm 2 on C5 on F_main) |
| [I6] | 9.8 to 12.6 % | the same |
| [R5m] | 13.3 % | arm 3 in the data set, cell and composition of [R6] |
| [PF] | 13 to 32 % | arm 3, `camera_floor`, sets of 10, main factor, transcriptomics and proteomics, cells A and B |
| [PC] | 30 to 48 % | arm 3, `camera`, the same |

## Determinism: 7 of 12 tables identical; the declaration asks for 12

The determinism job (cell A replicates 6 to 15 again) against the same replicates in job s01 (`s/determinism.tsv`):
`errors` (52,683 rows), `power`, `formed`, `feature`, `floor`, `factors` and `status` are identical; `counts`,
`identity`, `exact`, `caught` and `orientation` are not.

- **Size.** In those five tables every count, label, status and sign is the same; only columns of real numbers
  differ, in their last digits: `factor_r2` by at most 9.9e-14 (9.4e-12 of its value); `raw_cor_main_z1` by at most
  7.0e-15; `smallest_canonical_correlation` by at most 6.4e-13; `ratio` by at most 7.5e-16 and
  `largest_relative_residual` by at most 5.2e-13 (quantities that are themselves rounding error); the identity check's
  two differences by at most 5.2e-14. Compared from the two jobs' worker tables, column by column.
- **Cause.** The two jobs ran on different processors: the determinism job on an AMD EPYC 7763, job s01 on an Intel
  Xeon Platinum 8573C (each job's `runner.txt`). The scorer compares the tables as text, so a difference in the 13th
  digit of a real number is a difference. Within one job on one processor (L1b) all twelve were identical.
- **What reads those tables.** Every rate of K-E1, K-E4 and K-E5 is computed from `errors`, and K-E3's counts from
  `power`: both identical. K-E3 also reads the sign of two columns of `orientation` (identical signs; the correlation
  differs in its 15th digit). Outcome R6 reads `exact` (statuses identical) and `identity` (both far below 1e-8).
- **What the declaration says.** "Twelve tables must be identical ... A difference is reported with its size and
  cause; a difference in a table that a criterion reads is a defect, handled as R3 within S'; one that needs more is
  a stop to the owner." By its letter `orientation` is a table a criterion reads. I have not re-read the rule, not
  re-run anything and not changed the comparison: this is reported for the coordinator's and the auditor's ruling.
  A re-run of the determinism job would land on whichever processor the runner has, so S' cannot make the tables
  identical on purpose.

## Reported, not decided (`s/error-rates.tsv` holds every rate with its interval)

- **Arm 1 (corrected values, no factors) against arm 3.** On the covariate that lies in the space of the corrected
  values (F_ret), feature level: arm 1 rejects 6.4 to 7.4 % of clean features in cells A and B and 12.3 to 14.2 % in
  cell D (four factors removed); arm 3 3.3 to 4.2 %; arm 3o 3.7 to 4.4 %. Set level on F_ret: T1 minus T3 between
  +0.003 and +0.014 in the deciding classes, interval excluding 0 in 29 of 30. On the independent covariate (F_ind)
  arm 1 rejects 3.7 to 4.5 % of features and arm 3 3.3 to 4.1 %.
- **Arm 2 (no adjustment) against arm 3.** Feature level 2.3 to 3.2 % on F_ind and F_ret (conservative). Power on the
  main factor, sets of 10 / 30: cell A transcriptomics 0.40 / 0.90 against arm 3's 0.45 / 0.94, proteomics 0.40 / 0.89
  against 0.48 / 0.95; cell B transcriptomics 0.17 / 0.63 against 0.34 / 0.84, proteomics 0.17 / 0.64 against 0.30 /
  0.80. In cell D (two more factors removed than planted) arm 3 has less power than arm 2: 0.30 / 0.79 against 0.39 /
  0.89 (transcriptomics), 0.34 / 0.84 against 0.38 / 0.89 (proteomics).
- **Power of the other tests** (main factor, sets of 10, cells A / B; transcriptomics, proteomics): arm 3
  `camera_floor` 0.26 / 0.16, 0.32 / 0.13 ([PF]); arm 3 `camera` 0.45 / 0.34, 0.48 / 0.30 ([PC]); arm 2 permutation
  0.30 / 0.07, 0.24 / 0.05; arm 2 pre-ranked 0.46 / 0.41, 0.34 / 0.33. Metabolomics (sets of 5): at most 0.02 under
  every arm and test. Direction: every called planted set in the expected direction under every arm and test
  (`s/ke3-direction.tsv`).
- **The oracle arm (3o)** on K-E1's classes: at most 0.059 in cell A and 0.064 in cell B (proteomics C5 on F_ind,
  upper limits 0.066 and 0.071); in cells A and B its power is that of arm 3 within 0.03 (in cell D it is higher, 0.43
  against 0.30 for transcriptomics sets of 10).
- **`camera_floor`** (reported, VAL-97 (3)) on C0, C1, C2: at most 0.045 (arm 3) and 0.046 (arm 2) in cell A, 0.047
  and 0.042 in cell B; p_floor at least p_camera in every one of 7,334,788 set comparisons. **Permutation** (arm 2):
  at most 0.053 in cell A and 0.063 in cell B (proteomics C1 rho 0.6 on F_ind, upper 0.072). **Pre-ranked** (arm 2):
  C0 and C2 at most 0.063 and 0.058; correlated null blocks up to 0.631 and 0.637.
- **The default model** (trained on the values before RUV), C5 on its main factor: without the adjustment (D0) 0.000
  to 0.002 in cell A and 0.175 to 0.193 in cell B; with it (D1) 0.014 to 0.039 in cell A and 0.080 to 0.254 in cell B.
  On its factor closest to the planted unwanted factors (F_unw): D0 0.071 to 0.379 (cell A) and 0.109 to 0.311 (B);
  D1 0.044 to 0.274 and 0.064 to 0.357.
- **Cell D** (four factors removed, two planted; reported only): arm 3 on K-E1's classes at most 0.058 (proteomics C5
  on F_ind, upper 0.069); C5 on F_main 0.009 to 0.050.
- **Classes defined by the model's selection** (not calibration evidence, declaration section 3): arm 3, C4 on F_ind
  0.033 to 0.094 over data sets and cells, C4b up to 0.092, C1-out at most 0.031; C3 (metabolomics identifiers with
  two or three features) 0.044 to 0.048 on F_ind.
- **Formation in S** (`s/pools.tsv`): proteomics C5 at size 100 formed in 446, 438 and 213 replicates of cells A, B, D
  (56, 55, 53 %); metabolomics C5 in 744, 748 and 364 (93, 94, 91 %), so the metabolomics C5 check of K-E1 rests on
  those replicates; transcriptomics C4 under the default model in 20 % of replicates at size 10 and almost never
  above (not reported, as declared); one metabolomics block of size 10 did not form in one replicate per cell; every
  other deciding class and size formed in every replicate.
- **Exactness**: the written factors `exact` in every replicate; largest s_(k+1)/s_1 5.4e-15; largest relative
  residual 3.2e-11. Identity check: largest difference of `correlation` 6.2e-14, of `pvalue` (proteomics,
  metabolomics) 2.1e-13.
- **Versions**: ruv 0.9.7.1, RUVIIIC 1.0.19, limma 3.68.4, MOFA2 1.22.0, one row each (`s/versions.tsv`).
- **Seconds per replicate** (median; largest): cell A 102; 122. Cell B 103; 135. Cell D 119; 176.

## Against the expectation stated before the run (ticket 6.4, last paragraph)

- As expected: arm 3 at or below 5.5 % on F_ind in every class, C5 included, and at or below 4.5 % on F_main in the
  clean classes; arm 1 above arm 3 on F_ret at the feature level (about 7 % against about 4 %, about 13 % in cell D;
  the ticket expected about 8, 5 and 13.5); arm 2 at or below 5.4 % on K-E5's checks and with less power than arm 3
  in cells A and B; the default model's unwanted-looking factor calling C5 sets (D0 on F_unw).
- **Contradicted:** "for transcriptomics the written factors should be close to the planted ones, arm 3 near nominal
  on C5 on F_main in cell B and arm 2 well above". Measured: arm 3 0.296, arm 2 0.020. The written factors are close
  to the planted ones (median 0.985) and the adjusted test still calls these sets.
- **Contradicted:** "D1 on F_unw showing few or no calls". Measured: 0.044 to 0.357.
- Not as simple as expected: arm 1 against arm 3 on F_ind ("indistinguishable at both levels"): the interval of the
  set-level difference excludes 0 in 12 of 38 rows, the largest -0.015 (proteomics C5, arm 1 lower).

## Not checked

- Anything outside this one generator, its two unwanted factors and RUV-III with the groups as replicate sets (ticket
  6.7); real data.
- Why the adjusted test calls C5 sets on the main factor in cell B (the explanation above is inferred).
- Whether the five tables that differ between processors would be identical on the same processor in S (they were in
  L1b).
- The scorer's tables were not recomputed by a second implementation.
