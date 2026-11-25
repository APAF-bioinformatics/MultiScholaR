# KIT-211 declaration (committed before the lifecycle run L1 and before any counted job)

Ticket `tickets/kit/KIT-211.md` revision 1.4 (git-ignored; read in the main checkout), its register
`audits/AUDIT-KIT-211-plan-20261124.md` (Re-check 1 and Re-check 2 of 2026-11-26), methods ledger VAL-142 (2), VAL-171,
VAL-315, VAL-204. This is a measurement: no package code changes, no direction is chosen here. Everything below is fixed
before any run. After L1 nothing declared changes except in a dated VAL-54 section that changes no estimand, margin, rule,
seed or range. Rules, margins, families and stop rules below are the ticket's, quoted or restated without change; where the
ticket was under-specified the resolution is listed in "Resolutions" and nothing else is tuned.

Order followed (the coordinator's instruction of this task; the ticket's A3 says "after L1" for the declaration, the
coordinator's instruction is to commit it before L1, which is stricter): scripts and e4 tables -> this declaration with L1's
predictions -> L1 (one run, expected 30 job-minutes = 2.00 core-hours, ceiling 55 job-minutes = 3.67 core-hours; core-hours =
job-minutes x 4 / 60) -> stop for the coordinator's go before any counted job.

## Bindings

| Item | Value |
|---|---|
| Bound commit (branch point of `kit-211-threshold-filter-test`; a `v0.9-staging` head; holds KIT-200 and KIT-201) | 894b9cc821a70b03baab1e1fa47b3d47440c7ed5 (its package tree equals 4d588ba3ed5c6fee76a5a1450c028b409ba419ba's: `git diff --stat 4d588ba3 894b9cc8 -- R/ src/ inst/ tests/ DESCRIPTION NAMESPACE` is empty; only `docs/methods-decision-ledger.md` differs) |
| G0 | `l1/g0.txt`: the three protein filter files are unchanged since 19997323 (KIT-187's state): no stop; `R/func_lipid_qc.R` (KIT-202) scopes E9 only; the DA files (KIT-200, KIT-201, KIT-220, KIT-249) scope E1 |
| Every request's `package_commit` | a commit of this branch whose package tree equals the bound commit (`git diff --stat 894b9cc8 <package_commit> -- R/ src/ inst/ tests/ DESCRIPTION NAMESPACE` empty, recorded in `runs.tsv`); `tests` `[]`, `base_commit` "none", `blas` "reference", `extra_packages` `[]`, one `scripts` entry with `timeout_minutes`, `out` its directory under `kit-211/`, `out_pattern` `"\\.(tsv|gz|json|txt|md)$"` |
| KIT-197's scripts (extended, reused) | `sim-k197.R` SHA-256 71001de9b7b1c547e1e27fd92fb1b57a067121e336fa09e134f3317f1fa7c599; `score-k197.R` f546947985617c5a6c0e79edd2d2f4a062db2ed0bea922b9927af6a1a9e9e426; `l1-k197.R` 50d55e99c5eb99a4534a4ea93c39aca9fbe4738ca6a783e086b69cd78294da5e (merged ef63c918) |
| `sim-k211.R` SHA-256 | 51c2d4672fc7e758b565c4e74e886eac5f805751871e0a47feaa382c518011b4 |
| `score-k211.R` SHA-256 | aefc66baa7dc42a86f5f69247911e51b8701200b1e7d95851bec676ae1110d54 |
| `l1-k211.R` SHA-256 | 6079031f48a3c465b576c6d625d67fe93f1ddefcaffabd9c1d4ef796a0634562 |
| `e0-k211.py` / `e4hand-k211.py` SHA-256 | a1c20199b2111bf3626ed1db7cca34291958ef5d94c508fd9ed16eeccf3a830e / b82ff3112f1a19198758bb59a596ea9dc33adb91515ad1a88169bafb7bb2e734 |
| `e4/hand.tsv`, `e4/expect.tsv` SHA-256 | e1942e9efffd328228601889ed48af5ae68474092af81434fa02b598dd978f70 / 1f888116e00e867269e318bc5032500c3f82c482ac367eadbeffc4b1529e50cd; SHA-256 of `sha256sum` over the 73 files of `e4/` in sorted path order: 8154cf67bcb15223bcb5c518da8da8c62336a67cd81e8a574c361fc205b80012 |
| E0 (local text check, `e0-k211.py`, `l1/e0.txt`) | pass: `sim-k211.R`'s verbatim region equals `sim-k197.R` lines 26 to 279 plus exactly two inserted lines (`if (identical(cellname, "Av")) s2 <- s2 * 2^((20 - mu) / 2)` after the variance draw; `if (identical(cellname, "Ag")) beta <- 3 * beta` directly after `beta[where] <- signs`); the three verbatim blocks of `score-k211.R` and the two of `l1-k211.R` equal their KIT-197 lines |
| Software | R 4.6.1, the bound commit's snapshot library, reference BLAS (the request default) |

## Design (ticket 5, 6.1, 6.2, 18.3)

- Generator and seed rule KIT-197's: `set.seed(185000000 + 10000 idx + r, "Mersenne-Twister", "Inversion", "Rejection")`;
  cells A = 2 x 3 (idx 201), B = 4 x 3 (202), C = 5 x 3 (203), D = 3 x 4 (204), 3,000 proteins, 600 planted at +-1 log2 in
  `treated`, strength 0; masks vd5, vd10, vd25, mcar5, mcar10, mcar25 on the same values. Reported cell **Av**: cell A's layout, idx
  **211**, vd masks only, the marked line `if (identical(cellname, "Av")) s2 <- s2 * 2^((20 - mu) / 2)` (a stated assumption: the
  replicate SD halves for every 4 log2 of level), `OLD_RULE$Av = c(m = 3, g = 2)`, C0 run once per replicate over its three
  masks. Reported cell **Ag** (VAL-315): cell A's layout, idx **212**, vd masks only, 8-fold changes by the marked line
  `if (identical(cellname, "Ag")) beta <- 3 * beta` directly after `beta[where] <- signs` (it draws nothing; `y` and the returned
  `beta` both carry +-3), `OLD_RULE$Ag = c(m = 3, g = 2)`.
- Ranges: counted r = 201 to 400 (R = 200 per cell; two jobs per cell, r 201 to 300 and 301 to 400); L1 r = 9101 and 9102
  (E6 9101 to 9104); E1 reads KIT-197's r = 101 to 105 (S only, instrument check). Ag: r = 201 to 400, L1 9101 and 9102.
- Arms per replicate and mask: **S** (the package filter at (2, 2), percentile 1; DA on the kept rows of the masked matrix),
  **A** (the filter is S's; every finite value at or below the filter's own threshold tau, the type-7 1st percentile of every
  finite value of the masked table, is set to NA in the kept rows at the DA input), **B** (the filter counts a value as missing
  only when it is not finite, (m, g) = (2, 2); DA on `masked[H_B, ]`), **B'** (B's count, then a protein is removed when the share of
  ALL samples holding a finite value strictly below tau is at least 0.5; DA on `masked[H_B', ]`), references **N** (S's DA call on
  every row of `masked`, no filter) and **C0** (S's DA call on every row of the complete values, once per replicate). O = the
  pre-KIT-187 rule in-script (`perGroupKept`), fixed with S's definition of present. Ag: S, A, B, B' at (2, 2) and B1, B'1 = B, B' at
  (2, 1); no N, no C0. DA through `differentialAbundanceAnalysis(~ 0 + group, treat_lfc_cutoff = 0, eBayes_trend = TRUE,
  eBayes_robust = TRUE)` (`daTable`, verbatim). Tested = finite p and logFC in that arm; called = tested and q < 0.05; t_r,X =
  the largest p called (0 with no call).
- Subgroup: ell = 1 when the masked matrix holds a finite value at or below tau in `control` or `treated`; strata one_treated,
  one_control, both, other_only, none; s = +1 (NA in `treated` only), -1 (NA in `control` only), 0 otherwise. Delta_X = H_X \ O;
  L_X = proteins of Delta_X with ell = 1. ell0 = the complete-data counterpart (a complete value at or below that mask's tau in
  `control` or `treated`).
- Rows stored (sufficient statistics, nothing per protein): per (cell, mask, r, arm) the sums of n, n_tested, n_p05, n_call,
  n_bh (Benjamini-Hochberg at 0.05 on the arm's raw p), sum_err, sum_s_err, sum_sgn_err over the subsets added (not in O), kept_both
  (in O), whole, whole0, common (added and in S's kept set), onoff (Ag only), keyed by class, stratum, ell, band, ell0 as the
  subset uses them; per (cell, mask, r, arm) a metrics row (kept, added, calls, t_r, F1 violations, E3 flag, the E10 deviation,
  planted counts, calls shared with S, proteins gained and lost against S, seconds on arm S only).

## Estimands (ticket 6.3, as `score-k211.R` computes them)

Ratios of sums; "null" = planted change 0; "tested" = finite p and logFC in that arm; alpha_X(P) = null p < 0.05 in P tested in X /
null in P tested in X; rho^l_X = alpha_X(L_X) / alpha_X(O in H_X); the unit of resampling is the cell's replicate (all masks, arms and
references move together).

- **Primary: Phi_X = rho^l_X / rho^l_N**, pooled over the ten deciding settings (A vd5/10/25, B vd5/10/25, C vd5/10/25, D vd25), for
  X = S, A, B, B'. The reported counterparts: rho^l_X against O (R-a), rho^l_C0 = alpha_C0(ell0 = 1) / alpha_C0(ell0 = 0) (R-i), the share
  of S's excess a candidate removes (rho^l_S - rho^l_X) / (rho^l_S - rho^l_N), Phi^H_X = R-c_X / R-c_N (R-c = alpha_X(H_X, ell = 1) /
  alpha_X(H_X, ell = 0)), rho^l on L_S intersect H_X (R-b, a check), KIT-197's rho_O,X, rho^{ell0}_X, per stratum x ell values, the mcar
  pool, paired differences rho^l_S - rho^l_X (R-a).
- R-d: pow_X = sum T_X / sum 600; FDR_X = sum V_X / sum (V_X + T_X); S2's mean-of-ratios FDP; kept, tested, called and untested-in-contrast
  counts per data set; paired differences against S (dpow_X = pow_S - pow_X; dfdr_X = FDR_X - FDR_S).
- R-e: m6_l and m6_NA on Delta_X (null, one-sided strata) and m6_l on H_X; M5 and attenuation A_X over planted proteins of H_X and Delta_X.
- R-f: share of S's calls also called by X, X's new calls, proteins gained and lost against S. R-g: KIT-197's E_cal and E_comp per arm.
- R-h (reported): for S, A, B, B', N, C0: alpha among null proteins tested in X whose complete-data row mean lies within 1 log2 above tau
  (resolution 2 below), its ratio to C0's in the same band, power (called over all planted in the band) and attenuation in the band.
- R-v: every estimand for Av's three settings with S, A, B, B', N, C0, not pooled, no label.
- Ag (reported only; never labelled, pooled or in the outcome table), per arm and mask (vd5, vd10, vd25): (i) R-d; (ii) power on the
  300 planted -3 log2 and, separately, the 300 planted +3 log2 proteins, called over all planted of that direction, by q < 0.05 and by
  Benjamini-Hochberg at 0.05 on the arm's raw p; (iii) the paired change from 2 to 1 group required within B and within B' in (i) and (ii);
  (iv) the planted proteins on/off (complete values below the mask's lod, `solveLod(sim$y, maskRate(mask))`, in every sample of `control`
  or of `treated`): their number per data set, the share each arm removes, each arm's power on them by both rules; (v) attenuation over
  planted proteins of H_X by direction; (vi) calls gained or lost against S. Each with a two-sided 95 % percentile interval and a
  jackknife SE from Ag's own index matrix. S's and A's filter at 1 group is not run.

## Uncertainty (6.5)

Bootstrap B = 10,000; `set.seed(211000001, "Mersenne-Twister", "Inversion", "Rejection")` once; then for cells A, B, C, D, Av, then Ag in that
order `sample.int(R, B x R, replace = TRUE)` filled by row into a B x R index matrix (a cell absent from the rows draws nothing, so the
deciding cells' indices do not depend on Ag); one index set per cell, shared by its masks, statistics and arms; two-sided percentile limits
`quantile(type = 7)` for the reported estimands. The two ranges of a cell are joined into one cell of R = 200 before resampling. A range still
incomplete after its one continuation leaves its cell scored on the replicates that completed, R stated beside every figure and the incomplete range
listed. A resample in which a statistic is not finite (no tested unchanged protein in a subgroup) holds NA there: the limits are taken over the
finite resamples and the count of the others is printed (`boot_nonfinite`). Jackknife: delete one replicate of one cell with all its masks and arms;
SE = sqrt((n - 1)/n sum (theta_(i) - theta-bar)^2); jackknife-normal limits beside the bootstrap ones; a label that differs is printed as such.

## The primary decision, margins and classes (6.4, 6.5)

- **One primary endpoint:** Phi_X at R = 200 for X in {A, B, B'}; S the positive control. Margin **M = 1.20**.
- **PC** (outside the family): the one-sided 97.5 % lower limit of Phi_S is above M = 1.20. If not, outcome N0; each candidate's primary state is still
  printed, marked "does not discriminate from S".
- **removes**: the upper limit of Phi_X is at or below 1.20; **excess remains**: the lower limit of Phi_X is above 1 and the upper limit above 1.20;
  **not determined**: otherwise; precedence: lower limit above 1 and upper at or below 1.20 is "removes", printed "non-zero, within the margin". One-sided
  limits at 1 - 0.025/3 (Bonferroni over the three candidates): the 0.025/3 and 1 - 0.025/3 quantiles.
- **Costs power** (dpow_X = pow_S - pow_X paired; one-sided at 0.025/3): **costs power** if the lower limit is above **0.005**; **no material cost** if the
  upper limit is at or below 0.005 (a gain included); otherwise not determined. A's cost is judged on the estimates (R-h beside), not on the class alone.
- **False-share guard** (dfdr_X = FDR_X - FDR_S paired; one-sided at 0.025/3): **raised** when the lower limit is above **0.0025**; **not raised** when the
  upper limit is at or below 0.0025; otherwise not determined.
- **Class per candidate:** **E** (removes, no material cost, not raised); **C** (removes and not E); **R** (excess remains); **U** (primary not determined).
  A candidate withheld by F3 or F5 counts as U.
- Per-setting labels of Phi (24 settings x 3 candidates) at one-sided 0.05/72, reported, not deciding; a deciding setting whose label differs from the
  pooled one is named. R-a to R-v: two-sided 95 % percentile intervals and jackknife SEs, no labels besides these; Av and Ag never labelled.

## Falsifiers and stop rules (6.7)

G0 passed (above). **F1**: in any replicate-mask O not in H_S, H_A not equal H_S, H_S not in H_B or H_B' not in H_B (Ag also H_B not in H_B1 and
H_B' not in H_B'1): stop, no outcome read. **F2**: A's DA input differs from the definition (E3), or an E10 identity fails, in any counted replicate-mask:
stop (a mechanics failure before any estimand is an instrument defect under VAL-54; after scoring a result). **F3**: the pooled interval of m6_l,A excludes 0:
A's labels are withheld; the owner is told. **F4**: KIT-197's F4 on S (m6_NA,S), judged on the pool: stop. **F5**: the pooled interval of m6_l,B computed on
H_B (not Delta_B) excludes 0: B's labels are withheld; the same for B' on H_B'. F3, F4 and F5 are judged on the pool of the 24 settings (ratios of
sums; the cell's replicate resampled) with two-sided bootstrap limits at 0.05/4 each (quantiles 0.00625 and 0.99375; jackknife-normal beside); the 24
per-setting Bonferroni intervals (quantiles 0.05/48 and 1 - 0.05/48) are reported, not deciding. A counted job that does not complete is "did not
complete" for its cell and range (never a denominator), continued once inside S'; an unfavourable or inconclusive valid result is not re-run (VAL-54); no
replicates are added under this declaration; R is never reduced.

## The total outcome table (6.8; the first row that holds applies)

| Outcome | Condition |
|---|---|
| Stop | G0, F1, F2 (after scoring) or F4 |
| N0 | PC fails: S's excess beyond the unfiltered analysis is not shown to exceed the margin at the bound commit |
| K2 | two or three candidates are E (the owner chooses; declared tie-break recommendation: higher whole-list power when the lower limit of the two-sided 95 % interval of the paired difference is above 0.005, otherwise the next criterion; then fewer changed calls (R-f); then the smaller code change; not applied between B and B') |
| K1 | exactly one candidate is E (confirm it) |
| K3 | none E; at least one C |
| K4 | all three R |
| K5 | otherwise (some U, none E or C, not all R) |

Expected before the run (inferred; the result can contradict it): PC holds; A: Phi_A <= 1 and rho^l_A above 1 (indicative 1.4); B: Phi_B near 1 and rho^l_B
near rho^l_N (indicative 1.7); B' close to B with fewer proteins kept; A conservative and B nominal near tau (R-h); A loses a little power and attenuates more near
tau, B gains a little. **Expected outcome: K2** (B and B' E, the owner choosing between them; A E or C, its power cost being near the margin).

## The decision table (6.9; filled by `score-k211.R` as `decision-table.md`; the owner edits the wording, VAL-41)

One row per arm (S, A, B, B'; N and C0 reference rows). Columns: Subgroup (Phi_X with limits and class; rho^l_X against O; rho^l_N and rho^l_C0 beside; the share
of S's excess removed; per-setting values, the mcar pool, the strata; Phi^H_X beside); Near the threshold (R-h per arm); Whole list (kept, tested, called per data
set; power and its paired change; false share and its paired change); Fold change (m6_l, attenuation over the whole list and the added set); Intensity-dependent
variance (Av's figures, reported only); What changes in users' results (share of S's calls kept; new calls; proteins gained or lost; not-tested contrasts); Routes
affected; What users are told and see; Flag column (SD8); Groups required (VAL-204 (2), (7); VAL-315): Ag's figures, B and B' at 2 groups against 1, paired, with
KIT-217's cell Le beside (BH power 0.712 against 0.819 at vd25), not comparable in level (Le is lipid, 500 features in two assays, median SD 0.5, 2 x 6, the 10th
percentile, R = 100; Ag is proteomics, 3,000 features, median SD 0.27, 2 x 3, the 1st percentile, R = 200); in 2 x 3 a protein kept at 1 group with no value in a compared
group is still not tested there (KIT-204); not run: S's and A's filter at 1 group; Code change it needs; Measured where (one generator whose variance does not depend on
intensity; Av apart).

## Checks (7.3), jobs and budget (11, 18.3)

E0 to E10, D and G0 as the ticket lists them. E1: S on KIT-197's r = 101 to 105 (120 replicate-masks): the filter-dependent columns (membership counts per subset, class,
stratum and ell; kept and added counts; tau) equal KIT-197's committed rows exactly; the DA-dependent columns are written as one flag (G0 lists DA-file changes). E2: B's kept set
equals the package helper with the threshold switched off (sentinel protein of -1e6, percentile 0), Ag also at (2, 1). E3: A's DA input holds NA exactly at the masked NA cells and
the finite cells at or below tau of the kept rows. E4: the scorer on constructed tables against hand values. E5: the scorer twice on the L1 rows, byte-identical. E6: killed-and-continued
run. E7: F1's identities. E8: determinism mechanics (`K211_DETERMINISM=1`, cell A r = 9101, 9102). E9: B''s share rule equals `lipidIntensityFilteringHelper` on the B-kept table
(inside `quiet()`, `sample_columns` = the run columns). E10: `logFC` of B equals S's on every protein of H_S, B''s equals S's on H_S intersect H_B', N's equals B's on every protein of H_B
(Ag: B1's equals B's on H_B, B'1's equals B''s on H_B'), maximum absolute deviation 0.

| Job | Request | `timeout_minutes` | Expected (min -> core-h) | Ceiling (min -> core-h) |
|---|---|---|---|---|
| L1 | blind lifecycle (cells A to D, Av, Ag; r 9101, 9102 as two ranges; E1 on 120 S-only replicate-masks; E2, E9 helper runs; E6's two runs; E4; E5) | 45 | 30 -> 2.00 | 45 + 10 = 55 -> 3.67 |
| S-A1 ... S-D2 | eight jobs, one per cell and range (r 201 to 300 and 301 to 400; 600 replicate-masks each); S-A1 also runs D (the determinism check) | 120 | 8 x 74.5 = 596 -> 39.73 | 8 x 130 = 1,040 -> 69.33 |
| S-Av | cell Av, r 201 to 400 | 120 | 74.5 -> 4.97 | 130 -> 8.67 |
| S-Ag | cell Ag, r 201 to 400 (planned at the ceiling rate 92.3 min) | 139 | 83.7 -> 5.58 | 149 -> 9.93 |
| SC | scoring | 30 | 15 -> 1.00 | 40 -> 2.67 |
| S' reserve | one continuation, VAL-54 instrument rerun, or a repeat of L1 after a script change | | 0 | 130 -> 8.67 |
| SC' reserve | one re-scoring after a defect in `score-k211.R` | | 0 | 40 -> 2.67 |
| **Total** | | | **799.2 min -> 53.28** | **1,584 min -> 105.60** |

G2: each counted job's plan is recomputed from L1's seconds per replicate-mask; a job planned above 150 minutes is split at the same R with each timeout 1.5 times its plan; if the
total ceiling then exceeds the accepted one, stop: the coordinator decides. Ag cannot hold the deciding reading: if G2's re-plan of S-Ag alone takes the total over the accepted ceiling,
S-Ag is not sent and Ag is reported as not run; a re-run of S-Ag that does not fit S' is not made, and Ag is reported as "did not complete". Attempts (VAL-54): L1 once (repeated only
after a script change, inside S'); each counted job once; scoring once. An instrument defect may be fixed and re-sent under VAL-54's four conditions (no expected value, threshold, declared
comparison or stop rule changes; each run recorded in `runs.tsv` and `outcome.md`; the ticket's total ceiling unchanged; at most three further runs per item).

## Resolutions of what the ticket left open (none tunes an estimand, margin, rule, seed or range)

1. **F3's set.** The ticket names F5's set (H_B, not Delta_B) and says F3 tests "m6_l,A"; the register says F3's exact null holds on both H_A and Delta_A. F3 is judged on **H_A**
   (parallel to F5); m6_l on Delta_A is reported.
2. **R-h's band.** "Within 1 log2 above tau": the complete-data row mean `rowMeans(sim$y)` greater than tau and at most tau + 1, tau being the mask's filter threshold (arm-independent); the
   planted proteins "in the band" for power are counted over all 3,000 proteins of the data set (filter losses count against power).
3. **L1's two ranges.** The ticket says separate `rows/<cell>-<range>/` directories; L1 writes `rows/r9101` and `rows/r9102`, each holding every cell (one simulation process per range rather than
   twelve). The scorer joins every directory's cells exactly as it will join the counted jobs' (the same code path); a cell present in two directories is joined.
4. **The constructed tables of E4** are committed files written by `e4hand-k211.py` (which also computes the hand values and the designed outcomes from the definitions, with exact fractions) rather than
   built in R by `l1-k211.R`: the hand values stay independent of any R code. `K211_B_BOOT` (default 10,000) is read by `score-k211.R` and set to 1,000 only for the degenerate constructed scenarios (every
   replicate alike, so every limit equals its estimate); the base scenario that carries the hand values runs at 10,000; the scorer records B in `checks.tsv` (`bootstrap_B`, pass only at 10,000) and
   `outcome.tsv`, and no counted scoring sets it.
5. **Rows stored.** KIT-197's `ref` subset and residual-variance columns are dropped (no estimand here uses them); `band`, `ell0`, `n_bh`, the `common`, `whole0` and (Ag) `onoff` subsets and per-arm metrics are
   added. E1 therefore compares S's `added`, `kept_both` and `whole` rows with KIT-197's after collapsing `band`.
6. **Incomplete ranges.** `K211_PARTIAL=1` assembles the tasks that completed when a range is incomplete after its one continuation; `range.tsv` records expected and completed replicates per cell and range and
   the scorer lists incomplete ranges (`ranges.tsv`) and states each cell's R (`replicates.tsv`).
7. **G2's seconds** come from L1's lifecycle run, which carries E2's package filter calls and E9's helper; they overstate a counted replicate-mask and so plan conservatively.
8. **E9 at this commit.** `lipidIntensityFilteringHelper` takes the threshold and the proportion explicitly and the KIT-217 rule is absent (no classed warning `omics_intensity_filter_groups_guessed` in `R/`); the
   warnings muffled inside `quiet()` are counted and expected to be 0.
9. **F2 in counted rows.** Any E3 failure or non-exact E10 deviation stored in the metrics stops the scorer (outcome Stop, no estimand read).
10. **Ag's direction denominators.** 300 planted decreases and 300 increases (`n_down_all` per data set from the generator).

## L1: predictions (written before the run; blind: L1 prints and uploads only the files below)

L1 writes `l1/checks.tsv` (id, units, passed, pass, details), `l1/timing.tsv`, `l1/hashes.tsv`, `l1/console.txt`, `l1/g0.txt` and `l1/e0.txt` (the last two made locally). Expected: **every check passes at the first
attempt; run once.**

| Check | Expectation |
|---|---|
| E0 | pass (`l1/e0.txt`, local) |
| E1 | 120 of 120 units pass (cells A to D, six masks, r 101 to 105: S's filter-dependent columns equal KIT-197's committed rows); the DA-dependent flag is expected FALSE (KIT-201 moves every p-value) |
| E2 | 60 of 60 replicate-masks (two ranges x 30: cells A to D six masks, Av three, Ag three; Ag's (2, 1) variant included) |
| E3, E7 | 60 of 60 |
| E9 | 60 of 60; classed warnings counted: 0 |
| E10 | 60 of 60; maximum absolute deviation 0 |
| E8 | 2 of 2 identical |
| E4 | all sub-checks pass (56 expected: per scenario exit and outcome (11 x 2); classes, withheld, positive control and F4 flags for the ten scenarios that exit 0 (10 x 3); the hand estimates, the hand jackknife SEs, the join's R and the empty-resample statistic) |
| E5 | two runs exit 0, every output file identical |
| E6 | killed with at least one and fewer than four tasks complete; continued with the missing tasks; rows and metrics (seconds excluded) identical to an uninterrupted run; each replicate once with six masks |
| G2 | the planned minutes of a deciding-cell job and of the Ag job are each at most 150 |
| Cost | 30 job-minutes (2.00 core-hours) expected, ceiling 55 (3.67); the ticket's totals 53.28 expected and 105.60 at the ceiling |

A failed check is classified by VAL-54: an instrument defect (the log shows a mechanics failure and no result was judged) may be fixed and re-sent under the four conditions, at most three times per item; a result is never
re-run and no expectation here is changed after it is seen; when it is not plainly mechanics, the work stops and is reported.
