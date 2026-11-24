# KIT-217 measurement M: declaration (before any run)

Committed with its two scripts, alone, before the lifecycle ML and before any counted run (plan section 3 rule 1). The
design is ticket KIT-217 revision 1.3, section 4.1 (`tickets/kit/KIT-217.md`), with the owner's SD-10 (methods ledger
VAL-184, 2026-11-25). Nothing below changes after a result is seen; a later change is a dated section labelled as made
after the results, and any change to a script after ML goes through ML again. Status: **revision 2 (repairs D-01 to
D-03 of the declaration audit); waiting for the auditor's diff-check before ML is sent.**

**Revision 2 (2026-11-25, before any run).** The declaration audit (`tickets/kit/audits/AUDIT-KIT-217-plan-20261125.md`,
"Declaration audit") found three gaps in what the lifecycle checks and in the completeness rules; revision 1 is commit
fc644ffb. Changed: (D-01) the decision script's self-test now plants the audit's scenarios (a) to (f) with their
expected labels, including every outcome row through `score()`; (D-02) the lifecycle writes its rows through the
counted writing code to `$RUNNER_TEMP`, re-runs L2's first replicate as a determinism file, and scores them blind
(E7, added to G-M); the determinism range is read from the det file's name; (D-03) a unit whose DA failed for an assay
is not complete, the 5 % rule is applied per mask over the masks a statistic uses, paired statistics print the units
they used and dropped (the guard and Q2 read "not determined" above 5 % dropped), and a non-finite bootstrap value
makes a classification "not determined". These are instrument and completeness rules only: no estimand, margin, level,
cell, mask, arm, seed or outcome consequence changed.

## 1. Bound code and scripts

- Package: `v0.9-staging` at e6791275d2121a1ce2d00cf0be79b4b52395c8a6. Every request's `package_commit` is the commit
  on branch `kit-217-measure` that carries this file (it adds only files under
  `docs/qa/workflow-execution/kit/kit-217/`; nothing under `R/` changes). KIT-202 is not in it and need not be
  (KIT-185's tables have text ids; E1 shows the package's S set equals the oracle).
- Scripts (SHA-256 of the committed files):
  - `m/sim-k217.R` e80b7d5d597299f43c76be8663a6ee47a4ba675549d5d9e7bc35855aa44db1b3 (generator and measurement)
  - `m/score-k217.R` c3b97ebd9bb7697ef4a753ee23886ab73abda3c269bf88fc72165178a2386d89 (decision script)
  - (revision 1, fc644ffb: 18db7551…8ef6 and d6f59789…528f)
- Generator source: `docs/qa/workflow-execution/kit/kit-185/sim-k185.R` at 735c2ea5 (the KIT-185 merge; unchanged at
  e6791275, SHA-256 5ffc96049efbcba0c7c465c78201465aa7a6cc13490b20afa6e705524aa234ce), copied by line range under
  "verbatim" headers: 43-52 (`warning_count`, `quiet`), 71-76 (`GROUPS`, `cell`), 90-91, 101-102 (`maskMech`,
  `maskRate`), 105-108 (`scale1`), 109-199 (`generateReplicate`), 202-205 (`solveLod`), 206-234 (`maskMatrix`),
  257-273 (`ID_COL`, `assayTable`, `omicSlot`, `omicObject`), 303-328 (`filterArm`), 429-465 (`CONTRAST`, `qMethod`,
  `daTable`). Marked additions: (i) the cells; (ii) the mask `vdg` inside `maskMatrix`; (iii) `beta <- 3 * beta` in
  Le and Me after the effect draw; and the Lv variance line after the variance draw. `strength = 0`.
- Runner: `run-tests` script requests, R 4.6.1, the snapshot library of the bound commit, `blas: "reference"`,
  `ubuntu-22.04`, `tests: []`, `base_commit: "none"`, 4 workers (`parallel::mclapply`), one random stream per replicate.

## 2. Cells, replicates, seeds

| Cell | Layout (KIT-185) | idx | Role | Masks | R counted |
|---|---|---|---|---|---|
| L2 | lipid, 500 features in two assays of 250, 2 x 6 | 3201 | deciding | none, vd5, vd10, vd25, mcar10, vdg10, vdg25 | 200 |
| L1 | lipid, 500 in two assays of 250, 2 x 32 | 3202 | deciding | as L2 | 200 |
| M2 | metabolite, 150, 4 x 3 | 3203 | deciding | as L2 | 200 |
| M1 | metabolite, 150, 4 x 45 | 3204 | deciding | as L2 | 200 |
| Lv | L2's, s2 x 2^((20 - mu)/2) | 3205 | reported | vd10, vd25 | 100 |
| Le | L2's, effects +-3 log2 | 3206 | reported | vd10, vd25, vdg25 | 100 |
| Me | M2's, effects +-3 log2 | 3207 | reported | vd10, vd25, vdg25 | 100 |

- Seeds: KIT-185's line `set.seed(185000000 + 10000 idx + r, kind = "Mersenne-Twister", normal.kind = "Inversion",
  sample.kind = "Rejection")`, that is 217000000 + 10000 i + r: counted r = 1 to 200 (1 to 100 in Lv, Le, Me),
  lifecycle r = 401 to 403. Range 217,010,001 to 217,070,403, distinct from KIT-185, KIT-187/197/211, KIT-188 and
  KIT-210 (ticket G4).
- Unchanged features: 400 per L data set, 120 per M data set. Scored contrast `treated - control`; in 4-group cells
  the other two groups enter the fit only.
- Masks (KIT-185's `maskMatrix`): `vd<rate>` value-dependent with one lod for every sample; `mcar10` at random;
  `vdg<rate>` as `vd` with the lod raised by 0.5 log2 for the `control` samples before solving for the rate (same
  uniforms `U`); `none` an instrument check (E2, F1), never deciding. Exchangeable stratum (Q1, deciding): vd5, vd10,
  vd25, mcar10. Group-linked stratum (Q2, reported): vdg10, vdg25.

## 3. Arms (per replicate and mask, on the same masked matrix)

Every in-script rule uses each assay's own threshold tau = the type-7 10th percentile of every finite raw (2^y) value
of the assay (S's threshold, L-1), pi = 0.5, m = 2; b = finite values strictly below tau, f = finite values, per group
f_i (finite) and a_i (finite and not below tau), need_i = min(n_i, m).

| Arm | Keep a feature when | Role |
|---|---|---|
| NF | always (no filter) | reference |
| S | the package's `lipidIntensityFiltering` / `metaboliteIntensityFiltering`, percentile 10, proportion 0.5 (`filterArm`) | reference; its cells decide O0 |
| PC | NF tests it and its \|observed log2 fold change\| is at or above the median among NF's tested features | positive control (gate) |
| B | f >= 1 and b/n < pi | reported; carries S's figures (E3) |
| D | (b + n - f)/n < pi | primary family |
| Dk | f >= k and b/n < pi; k = `trxFilterMinSampleSize`: L2, Lv, Le 6; L1 25.4; M2, Me 3; M1 34.5 | primary family |
| Atau_g2 | at least 2 groups with a_i >= need_i | primary family |
| Aobs_any | at least 1 group with f_i >= need_i, and b/n < pi | primary family |
| Aobs_g2 | at least 2 groups with f_i >= need_i, and b/n < pi | primary family |
| Atau_any, Atau_all, Aobs_all | as named (`all` only in 4-group cells; with two groups g2 = all) | reported |

- DA for every arm: KIT-185's `daTable` (the package's `runLipidsDA` / `runMetabolitesDA`, `~ 0 + group`, q 0.05,
  TREAT 0, trend and robust TRUE) on the arm's kept rows of the masked log2 matrix, rows without any value included
  (KIT-196's removal in force).
- Observed fold change (PC, F2, mean fold change): the in-script difference of the observed group means
  (treated minus control). For `~ 0 + group` it equals limma's logFC by algebra; E6 checks it against the DA's logFC.

## 4. Estimands (per data set = cell, mask, replicate; "tested by X" = kept by X and given a finite raw p)

- alpha_X,c (Q1): sum of unchanged features tested by X with raw p < 0.05 over sum of unchanged features tested by X,
  pooled over the cell's four exchangeable masks and its replicates (ratio of sums).
- alpha^G_X,c and e^G_X,c = alpha^G_X,c - alpha^G_NF,c (Q2): the same over vdg10 and vdg25, paired with NF.
- Power and false share at the matched rule (VAL-11): Benjamini-Hochberg on each arm's raw p per assay at 0.05;
  pow_X = called planted / all planted of the data set (filter losses count against power), FDP_X = false calls /
  calls (0 with no call), FDR_X = mean FDP over data sets. Pooled over the deciding settings (power a ratio of sums,
  FDR a mean); the shipped call rule (q < 0.05, `qMethod`) printed beside, deciding nothing.
- Printed per arm and cell, never deciding: kept per data set and of them with no value, zero residual df (f >= 1 and
  f = the number of groups with a value) and no value in `control` or `treated`; realised missing share before and
  after; E-l (alpha among unchanged tested features holding a finite value below tau in `control` or `treated`, and its
  ratio to those without); mean fold change of unchanged kept features; planted features removed by direction and, under
  vd and vdg masks, those "on/off" (complete values below the mask's effective lod in every sample of `control` or of
  `treated`); alpha_X - alpha_S (paired); the fixed-set ratio on K_X (unchanged features X tests that hold a missing or
  below-tau value in the scored groups): their rate under X over their rate under NF; changed BH calls against S; the
  null rate per setting (16 per candidate).

## 5. Intervals and the "not determined" rule

- Resampling unit: the cell's replicate; every arm, mask and statistic of a replicate moves together.
- Percentile bootstrap, B = 100,000 resamples, type-7 percentiles. One index set per cell, shared by every arm, mask
  and statistic: `set.seed(s_c, kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")`,
  `idx = matrix(sample.int(n, n * 100000, replace = TRUE), nrow = 100000)`, row b = resample b; s_c = 217990001 (L2),
  217990002 (L1), 217990003 (M2), 217990004 (M1), 217990005 (Lv), 217990006 (Le), 217990007 (Me). Statistics pooled
  over cells resample each cell with its own index set (stratified).
- Beside it, the delete-one-replicate jackknife: estimate +- t(df) x SE; in one cell df = n - 1; pooled over cells the
  stratified form, var = sum_h (n_h - 1)/n_h sum_j (theta_hj - mean_h)^2, df = sum_h (n_h - 1).
- Agreement rule: every classification below (Q1 per cell, NF's reference check, the PC gate, the guard, Q2, F2's
  interval) is made by each method; where they differ the classification is "not determined". Where a rule says "both
  methods", both must hold. Clopper-Pearson over features is not used (VAL-11).
- Paired statistics (X against S, X against NF, X against D-k or B) use only the (replicate, mask) units both arms
  completed.

## 6. Q1 under SD-10 (a), the positive control, the guard, the class

- Family: D, Dk, Atau_g2, Aobs_any, Aobs_g2 (Bonferroni over five). S and NF are classified at the same levels.
- Per (X, c): **meets** when the one-sided upper limit of alpha_X,c at 0.025/5 is below 0.075; **fails** when the
  one-sided lower limit at 0.025/(5 x 4) is at or above 0.075; otherwise **not determined**.
- Reference check: NF "meets" in a cell when both methods' one-sided 97.5 % upper limits of alpha_NF,c are below 0.075.
  Where it does not, every X (S and PC included) is classified in that cell on the paired excess alpha_X,c -
  alpha_NF,c with edge +0.025 at the same levels, and the output says so.
- Verdict per candidate over the reading cells (intersection-union, no adjustment over cells): calibrated = meets in
  all; miscalibrated = fails in at least one (cells named); otherwise not determined (cells named).
- Positive-control gate: in every deciding cell, both methods' one-sided 97.5 % lower limits of alpha_PC,c are above
  0.075, or, where NF does not meet, of the paired excess over NF above +0.025. Passing shows that a 2-fold excess is
  seen, not one of the band's size.
- Guard (paired against S, one-sided at 0.025/5, pooled over the reading cells' exchangeable settings): costs power if
  the lower limit of pow_S - pow_X is above 0.005; no material cost if its upper limit is at or below 0.005; false share
  raised if the lower limit of FDR_X - FDR_S is above 0.0025; not raised if its upper limit is at or below 0.0025;
  otherwise not determined.
- Class: E (calibrated, no material power cost, false share not raised); C (calibrated, not E); R (miscalibrated);
  U (not determined). SD-10 (iii) is not adopted (VAL-184 does not include it): alpha_X - alpha_S is printed only.

## 7. Q2 (SD-10 (ii-b): reported, with a fixed consequence)

X in the family plus S and B: **robust** in a cell when the one-sided upper limit of e^G_X,c at 0.025/5 is below
+0.025; **not robust** when the one-sided lower limit at 0.025/(5 x 4) is at or above +0.025; otherwise not
determined. Verdict: robust in all four cells / not robust in at least one / not determined. Printed: alpha^G of every
arm and NF; e^G_Aobs_g2 - e^G_Dk with its two-sided 95 % interval. **Fixed consequence:** "If A-obs-g2 is 'not robust'
in Q2 in any cell and D-k is 'robust' in every cell, the decision table names it and the owner is told before choosing;
no class changes."

## 8. Falsifiers, determinism, "did not complete"

- **E1 in counted rows:** in any counted (replicate, mask) the package's S set differs from the oracle: Stop.
- **F1:** in any counted (replicate, mask), B, D, Dk or an Aobs arm keeps a feature S removes, Aobs_g2 keeps one
  Aobs_any removes, or (mask none) any of B, D, Dk, Aobs_* differs from S: Stop.
- **F2:** for B, Dk, Aobs_any, Aobs_g2 and Aobs_all (4-group cells only), per deciding setting (cell x exchangeable
  mask), the mean fold change of unchanged kept features with its two-sided interval at level 0.01/80 (99 %
  Bonferroni over 16 settings x 5 arms): fires when an interval of either method excludes 0 **and** the absolute
  estimate exceeds 0.02 log2. That arm gets no label and the owner is told; the rest is read.
- **Determinism:** L2's replicates 1 to 5 are run again in a separate job (`K217_JOB=det`); their rows (seconds are
  kept out of the rows) must be byte-identical to the counted L2 job's. A difference: Stop.
- **Did not complete:** a unit (cell, mask, replicate, arm) that errors is "did not complete", listed with its error
  class, never a denominator. Above 5 % of a cell's units for an arm (over all the cell's masks), that arm's cell is
  "no result" and counts as not determined, the cell named. Revision 2 (D-03): (i) a unit whose DA failed for one
  assay (`da_failed > 0`) is not complete (KIT-188 r1.9's rule); (ii) the share is computed over the masks the
  statistic uses (the exchangeable masks for Q1, the guard and F2; vdg10 and vdg25 for Q2), and any single mask above
  5 % also makes that arm's cell "no result"; (iii) every paired statistic (excess, guard, Q2, O2's power tie-break,
  large effects) prints the (replicate, mask) units it used and, per arm, the units it completed that were dropped
  because the other arm did not; the guard and Q2 read "not determined" when more than 5 % of their units are dropped;
  (iv) every interval prints its count of non-finite bootstrap values, and any such value makes the classification
  "not determined" (an F2 setting not determined counts as firing: the arm gets no label).

## 9. Outcomes (the first row that holds applies) and their consequence

| Outcome | Condition | Consequence |
|---|---|---|
| Stop | E1 in counted rows, F1 or determinism fails | the owner is told; nothing else is read |
| N0 | the positive-control gate fails in some cell | no candidate is labelled; the owner is told |
| O0 | S does not meet in some deciding cell | the owner is told first with NF's and every candidate's per-cell figures; O1 to O3 are read on the cells where S meets; the Methods names the others. With no such cell: no recommendation, the owner decides on the per-cell figures |
| O1 | Aobs_g2 is class E | recommend A-obs-g2 |
| O2 | Aobs_g2 not E; Aobs_any or Dk is E | the owner chooses; recommendation: the E candidate with higher whole-list power (the paired difference's two-sided 95 % interval excludes 0 by both methods), else fewer changed BH calls against S per data set, else the smaller code change (D-k: no new user control) |
| O3 | no per-group or count candidate is E | recommend B (it changes no tested result and removes only rows with no value); the Methods states that a per-group rule was not shown calibrated on these shapes |

In every outcome: A-tau's and D's classes are reported; a class-E A-tau does not displace A-obs. F2's named arms,
Q2's verdicts with the fixed consequence, the large-effect consequence and Lv are printed. **Large effects (Le, Me):**
if, in any of their masks, the recommended candidate removes a larger share of planted features than D-k, or than B,
with both methods' paired two-sided 95 % lower limits of the difference above 0.05, the decision table names the loss
and the owner is told before choosing; no class changes. An unfavourable valid result is recorded, never re-run (VAL-54).

## 10. Lifecycle ML (blind; the plan auditor confirms this file and the scripts before it is sent)

One request, two scripts in one job:
1. `sim-k217.R L2,L1,M2,M1,Lv,Le,Me 401 403 docs/qa/workflow-execution/kit/kit-217/m/ml` with `K217_LIFECYCLE=1`,
   `K217_FAULT=1` (arm D raises in L2, replicate 401, mask vd10), `timeout_minutes: 120`. Writes only completion per
   (cell, mask, arm), seconds per replicate, the SHA-256 of each replicate's rows, the job-minute plan, and pass/fail of:
   **E1** (oracle = package S, every unit), **E2** (mask none: B, D, Dk, Aobs_* = S), **E3** (B's tested rows and raw
   p = S's to 1e-12; reported), **E4** (the fault is "did not complete", every other row present), **E5** (the
   counted jobs' plan from these seconds), **E6** (the rules reproduce 2.4's ten-feature table in every column, T3 and
   T4; `trxFilterMinSampleSize` gives the declared k; the in-script fold change equals the DA's logFC to 1e-10).
   Revision 2 (D-02): it also writes every cell's rows through the counted jobs' writing function to
   `$RUNNER_TEMP/k217-ml-rows` (never under `out`) and runs L2's replicate 401 a second time as
   `det-rows-L2-401-401.tsv.gz` (the fault fires identically, so the rows must match).
2. `score-k217.R selftest docs/qa/workflow-execution/kit/kit-217/m/ml`: (a) to (f) of the self-test, each with its
   expected label, as listed in `score-k217.R` (the declaration audit's D-01: a planted calibrated arm, a 2-fold arm
   failing every cell, NF meeting, the PC gate; a failure in one cell only; an arm at the band reading not determined;
   NF failing in one cell, with the excess basis, meets and fails there; the guard's no material cost, costs power and
   raised; and `score()` on planted directories giving Stop, N0, O0 then O1, O1, O2 for A-obs-any, O2 for D-k and O3
   for B); any mismatch fails G-M.
3. `score-k217.R e7 docs/qa/workflow-execution/kit/kit-217/m/ml`: runs `score()` on `$RUNNER_TEMP/k217-ml-rows` into
   `$RUNNER_TEMP/k217-ml-score` with its console discarded, computing every table even where an earlier outcome row
   holds. **E7** passes when it completes, `stops.tsv`, `outcome.tsv`, `q1-cells.tsv`, `classes.tsv`, `q2.tsv`,
   `f2.tsv`, `beside.tsv` and `per-setting.tsv` exist, and `stops.tsv` shows 0 determinism failures. Only E7's pass or
   fail is written (`e7.tsv`); nothing under `$RUNNER_TEMP` is kept. (The audit wrote this step as a shell line; the
   run-tests runner passes arguments without a shell, so the script resolves `$RUNNER_TEMP` itself.)

Gate G-M: E1, E2, E4, E6, E7 and the selftest pass; otherwise an instrument defect under VAL-54 (at most three further
instrument runs, each recorded). Timing rule: a counted job whose plan exceeds 150 minutes is split at the same R; R is
never reduced; if the total ceiling would be exceeded, stop: the coordinator decides.

## 11. Counted runs MC and scoring MS (after G-M)

| Job | Script arguments | Env | timeout_minutes |
|---|---|---|---|
| L2 | `L2 1 200 <out>` | `K217_JOB=m` | 150 |
| L1 | `L1 1 200 <out>` | `K217_JOB=m` | 150 |
| M2 | `M2 1 200 <out>` | `K217_JOB=m` | 150 |
| M1 | `M1 1 200 <out>` | `K217_JOB=m` | 150 |
| reported | `Lv,Le,Me 1 100 <out>` | `K217_JOB=m` | 150 |
| determinism | `L2 1 5 <out>` | `K217_JOB=det` | 60 |
| MS | `score-k217.R score <rows dir> <out>` on the committed rows | none | 120 |

The rows are committed under `kit-217/m/rows/` before MS; MS binds that commit.

## 12. Budget and stop rules

Core-hours = job minutes x 4 / 60. Planning from KIT-185's measured seconds (`kit-185/jobs/*/timing-*.tsv`: per L2
replicate 9.4 s of DA for about 56 calls, 0.4 s per filter call; M2 0.05 s and M1 0.07 s per DA call; about 11
minutes of set-up per job), to be replaced by ML's seconds:

| Step | Expected (min) | Ceiling (min, ticket section 7) |
|---|---|---|
| ML (lifecycle, selftest with seven planted `score()` runs, E7) | 40 | 197 |
| MC: L2 25, L1 27, M2 16, M1 19, reported 17, determinism 12 | 116 | 967 |
| MS | 25 | 88 |
| Instrument reruns (VAL-54) | 0 | 450 |
| **M** | **181 = 12.1 core-hours** | **1,702 = 113.47 core-hours** |

The ticket's expected 33.47 core-hours (502 min) was scaled from r1.0's table; this plan's 12.1 is lower because
KIT-185's per-call seconds are smaller than r1.0 assumed. The ceiling is the ticket's and does not change. Stop rules:
a counted run that would take the total past 1,702 job-minutes is not sent (the coordinator decides); a script or
mechanics failure with no result is an instrument defect (VAL-54, at most three further runs per item); anything else
is a result.

## 13. Limits (as ticket 4.4)

One generator (normal levels; variance independent of level except Lv; logistic or uniform missingness; no unwanted
factor; no normalisation, ITSD or RUV between filter and test; equal groups); effects +-1 log2 in deciding cells (+-3
only in Le and Me, so presence/absence changes are rare in the deciding cells by construction); percentile 10 and pi
0.5 only, m = 2; group-linked missingness in one stylised form (0.5 log2 lod shift in one group at two rates); no pool
groups, no design without a group column; the DA at the bound commit; no real data. Under (a) a rule at about 1.3 times
nominal is not separated from nominal, and one at 1.5 times is "not determined" rather than "fails" (VAL-184).
