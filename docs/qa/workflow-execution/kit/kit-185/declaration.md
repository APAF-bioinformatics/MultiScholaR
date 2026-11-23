# KIT-185 declaration (before the counted run S)

Committed alone, after the blind lifecycle (`lifecycle.md`) and the owner's D1 to D5 (ledger VAL-67) and VAL-71, before
any counted job. Nothing below changes after results are seen; a later change is a dated section labelled as made after
the results. The design is KIT-185 revision 1.2 (`tickets/kit/KIT-185.md`, sections 5 to 7 and 10); this file pins it.

## 1. Bound code

- Package: `v0.9-staging` at ffdfdfb76a9eacdc681f2098f75ac574e2cf8af4 (the branch base). Every counted request's
  `package_commit` is this declaration's commit on `kit-185-ruviiic-missing`, which adds only files under
  `docs/qa/workflow-execution/kit/kit-185/`.
- Scripts (SHA-256):
  - `sim-k185.R` 3ee4aac433a694ca633975479835bf3a1b0daca368d611886f4b6c572c5895e7
  - `ruviiic-ref.R` 019ee0a6a039f3bf1b842dbff6543e964dda08ce393908a933ade1ea40e39130
  - `score-k185.R` 5a708c6e3e34feecec60b36836eaa4d090b32ded35ef315a79eea4a0593f71eb
- Runner: R 4.6.1, the snapshot library of the bound commit, `extra_packages: ["sva"]`, reference BLAS,
  `OMP_NUM_THREADS=1`, one RUV-III-C OpenMP thread, 4 workers (`parallel::mclapply`), one stream per replicate.

## 2. Generator, cells, masks, corrections, estimands (as revision 1.2)

- Generator: 5.2 (levels N(20, 2²), P0lo 14, P0hi 26; proteomics residual variance 4 x 0.25² / chi²_4, lipidomics and
  metabolomics 0.25 x / median(x), x = 10 / chi²_10; 20 % of features +-1 log2 in group 2, half up; k = 2 planted factors
  of loading variance 1.5 x the median residual variance, orthogonal to each other and to the groups, the lipid/metab
  first factor the run-order drift exp(-2u); P0r rho = 0.3; P5 75 % / 25 % unit / injection; D5 cells: centred batch
  offsets of variance 3 x the median residual variance (plus the drift in L1b, M1b), planted k 4, 2, 15). Draw order:
  levels, variances, effects, factors, loadings, (offsets), residuals, U, V, (W). Seeds:
  `set.seed(185000000 + 10000 c + r, kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")`,
  c: P0 1, P0r 2, P1 3, P2 4, P3 5, P4 6, P5 7, L1 8, L2 9, M1 10, M2 11, P0lo 12, P0hi 13, P4b 14, L1b 15, M1b 16.
- Masks: `none`; value-dependent 0.1, 0.5, 2, 10, 25, 50 %; at random 0.5, 2, 10 %; D5 cells `none`, VD 2, 10 %,
  batch-block 2, 10, 25 %; run-aligned VD 2, 10 % (P0, P1, L2, M2; reported).
- Filter: the route's own (proteomics `updateMissingValueParameters(2, 2)` then `removeRowsWithMissingValuesPercent`;
  lipidomics / metabolomics the intensity-filter S4 method, percentile 10, proportion 0.5, on 2^y).
- **G2 decided by L1: no lipidomics or metabolomics DA failed on the shipped filter, so the presence rule is not taken.**
- Corrections S (shipped S4 method and the route's post-processing), C (`ruviiic-ref.R`, `"complete"`), O (truth
  subtracted), N (none); Z (`"zero_filled"`) every 20th replicate in P4, L1, M1, P4b, L1b, M1b. Part A: k planted, every
  feature a potential control, replicate sets as the route resolves them, normalisation none.
- Estimands (a) to (f), contrasts Delta_tot, Delta_C, Delta_attr, oracle reference of the DA measures (5.5 with B01).
  Part B, the sensitivity panel, P0lo, P0hi and the run-aligned masks are reported, not classified.

## 3. Reading rules (the owner's D1, D3, D5 (a))

- Bands (material / serious): B 0.05 / 0.25; C* 0.02 / 0.10; FDP +0.010 / +0.025; power loss 0.02 / 0.10 (A1 and the null
  rate reported with 0.05 / 0.25 and +0.010 / +0.025).
- "No material effect": one-sided 95 % upper limit of the mean paired Delta_tot below the material edge. "Material" /
  "serious": Bonferroni one-sided lower limit (t on R - 1 df at 1 - 0.05 / 456) at or above the edge. Otherwise "not
  determined". The same for Delta_C. Refused above 5 %: "refused by the route" (Clopper-Pearson interval printed); did
  not complete above 5 %: "not determined"; below 5 % the refused share is printed beside each classification.
- Deciding family: 456 = 11 cells x 9 masks x 4 + 3 D5 cells x 5 masks x 4. The divisor stays 456 while P2 is held
  (section 4), so P2's later classifications use the same limits.
- Attribution, onset table, safe rate and outcomes O1 to O4 as 5.7. Attribution uses the replicates with at least
  max(20, 3 k) complete features (the others counted as "not assessable"); "not assessable" when fewer than two remain.
  The safe rate of a cell and mechanism stops at the first rate whose deciding classifications are not all "no material
  effect" ("refused by the route" included).

## 4. P2 held for KIT-187 (VAL-71)

The lifecycle found that P2 (2 x 3 runs) receives no missing value: `updateMissingValueParameters` rounds (3-2)/3 to
33.333 % and the filter removes a group whose share is greater, so one missing value of three fails the group and only
fully observed proteins pass. E1-neg failed on P2 at every L1 replicate (C equals S to 3e-14). This failure is the finding
that led to VAL-71. Decision (VAL-71): the rounding is a package defect fixed by KIT-187. P2 is not run now. It runs after
KIT-187 on the fixed filter under the rules of this declaration, with its own declaration addendum committed before its
run and E1-neg checked then. Until then every conclusion that involves P2 (its classifications, its onset and safe rate,
and the overall outcome, which P2 can only make worse) is provisional, and KIT-161's PXD000279 reading waits for it.

## 5. Lifecycle outcomes carried into S

- G1 passed in every part-A cell; G3: no re-split (lifecycle.md), jobs and timeouts as section 10 of the ticket.
- Mechanics changes after L1 (no computed value changes): E3's parameters, the part-B comparator's assay split and
  control order (run L1fix), G1 not printed for sensitivity variants, rows written to 7 significant digits (size), E8 on
  the N rows of the test on F_filt (the other N strata depend on the variant's own S and C).
- E6 for lipidomics (L2) mirrors the declared DA's input rule (stated as a limit).

## 6. Jobs (ceilings from section 10; P2's two jobs not sent)

Part A: P0 1-200, 201-400; P0r same; P1 1-100 ... 301-400; P3 and P4 80 each (5 jobs); P5 1-134, 135-267, 268-400; L1,
L2, M2 1-400; M1 1-134, 135-267, 268-400; P0lo, P0hi 1-200; P4b 1-134, 135-267, 268-400; L1b 1-400; M1b 1-200, 201-400.
Part B: P0 1-50, 51-100; L2 and M2 1-100 (one job); P1 and P4b 1-50, 51-100; L1b and M1b 1-100 (one job). Sensitivity:
P0 1-100, 101-200; M2 1-200. Run-aligned: P0, L2, M2 1-400 (one job); P1 1-200, 201-400. Determinism D: P0 and M2
replicates 1-5 (rows renamed `d-` for E5). Scoring SC once over the committed rows.

## 7. Dated section (2026-11-24), made after the results: scoring defect fixed and rows re-scored once

The first scoring run (SC, run 45281426216, script SHA-256 5a708c6e…71eb) reported E5 and E8 as failed. Checked outside R
on the committed text rows, the compared rows are identical (E5: P0 and M2 replicates 1-5, 2,360 rows each; E8: k1 and k4
N rows, 1,600 each; run-aligned `none` rows, 35,760). The failure was in `score-k185.R`'s comparison: `identical()` on
data frames read separately, whose inferred column types differ (an all-NA column is logical in one table, numeric in
another). Fix: the compared columns are converted to text before `identical()`; no estimand, contrast, interval, band or
classification code changed. New SHA-256 of `score-k185.R`: see the commit that adds this section. All rows re-scored once.
