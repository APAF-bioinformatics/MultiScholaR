# KIT-275: the counted declaration (32 counted jobs; tiers D, K, F+S, M, ID, T, E)

Committed alone, before any counted request, after the auditor's reading of it. Tier V's six jobs and L1-V are held (they need KIT-242 on staging) and are declared separately. This file fixes what the counted run does; it states **no prediction of any outcome**, only of the checks.

## 1. Binding

- **Bound commit** `c4e5be4091488cd1f4b40013c3f37e730074d685` (branch `kit-275`; every counted request names it as `package_commit`, or a later commit of the branch for which items (a) to (c) below are re-verified and recorded in `runs.tsv`). It has staging `ec27f69e0` (which holds KIT-255 with its I01 repair, merge aa119308f) as an ancestor.
  - (a) `git diff --stat f3c1358b c4e5be4091488cd1f4b40013c3f37e730074d685 -- R/ src/ DESCRIPTION NAMESPACE` is **empty** (verified when this file was written).
  - (b) The scripts are those audited at `d6b95a05`: `git diff --stat d6b95a05 c4e5be4091488cd1f4b40013c3f37e730074d685` over `sim-k275.R`, `score-k275.R`, `fixture-k275.R`, `e0-diff.py` and `selftest/` is **empty** (verified). SHA-256, first 16 hex digits: `sim-k275.R` 2cc87ea569eef7fd, `score-k275.R` a43a2c10c9a41ec2, `fixture-k275.R` 407949a9a6c0f90f; each job's `provenance-<job>.json` records the full hashes, which must equal these.
  - (c) KIT-242 (`origin/kit-242`) and KIT-210 (`origin/kit-210-ruviiic-remedy`) are **not ancestors** of the bound commit (`git merge-base --is-ancestor`: both absent; verified). G0's reading stands (`g0.txt`); if either lands on staging before a request, the request is not sent without a re-reading of G0 and the coordinator's word.
- G1 (run 45755206721) and L1 (run 45756994368) passed; the proving send (run 45785308371) showed the continuation paths. The measurement's R, seeds (generator 275000000 + 100000 c + 10000 s + r; planting 276000000 + 100000 c + r; bootstrap 277000000 + 100 c + t, ID none, G1 t = 10), arms, tiers, designs, estimands, margins, levels, falsifiers and the outcome table are those of ticket sections 4 to 7 (quoted in section 7 below); nothing was changed after L1.

## 2. The ceiling and the margin rule

- **Ceiling 406.90 core-hours (VAL-382)**: used 12.43 (smokes 4.35, G1 0.73, L1 7.35) + the proving send 3.27 + remaining 391.20. Core-hours = job-minutes x 4 / 60. Timeouts follow the **1.5 x rule**: `timeout_minutes` = ceil(1.5 x the job's L1-measured compute), compute = data sets x seconds per data set / 4 / 60 (`lifecycle.md` section 18, `budget-replan-r15.py`); ceiling per job = timeout + 10 minutes of set-up; every value is set explicitly (the default of 120 would kill the longer jobs). The 32 jobs' ceilings sum to 4555 job-minutes (= 303.67 core-hours); expected about 3,229 job-minutes.
- Each counted job is its own request (scripts of one request run in one job of at most 300 minutes).
- **The running-ceiling row**: before each request the implementer writes in `runs.tsv` the core-hours already used + this request's ceiling + the ceilings of every request still planned (the 391.20 less the requests already sent, plus the continuation and reserve rows below); the request is sent only if the sum is at most 406.90, otherwise stop and report to the coordinator.
- Budget rows outside the 32 jobs: the interim controls read 30 job-minutes (20 + 10 set-up); final scoring 40; the reserve 150 (inside it: one continuation per job); the correction pass 230 (largest job 190 + rescoring 40); tier V's six jobs (822) and L1-V (41), held.

## 3. The jobs (one request per job; no lifecycle, smoke, fault or guard flag; no `K275_ONLY_UNITS`)

Every job: `Rscript sim-k275.R <tiers> <design> <first> <last> <out>` from the package root, environment `OMP_NUM_THREADS=1`, `K275_WORKERS=4`, `K275_JOB=<job>`, `K275_DECLARATION=docs/qa/workflow-execution/kit/kit-275/declaration.md`; `out` = `k275-<job>`; `timeout_minutes` as below; `extra_packages` `["sva"]`; `base_commit` none; `fetch_refs` `[]`. **Upload pattern** (`out_pattern`): `\.(tsv|gz|json|txt)$` (covers `rows/*.tsv.gz`, `rows/part-*.tsv`, `checks-*.tsv`, `progress/*.tsv`, `timing-*.tsv`, `provenance-*.json`).

| # | Job | Command arguments | `timeout_minutes` | Ceiling (job-min) |
|---|---|---|---|---|
| 1 | `D-D3-1-50` | `sim-k275.R D D3 1 50 k275-D-D3-1-50` | 141 | 151 |
| 2 | `D-D3-51-100` | `sim-k275.R D D3 51 100 k275-D-D3-51-100` | 141 | 151 |
| 3 | `D-D3-101-150` | `sim-k275.R D D3 101 150 k275-D-D3-101-150` | 141 | 151 |
| 4 | `D-D3-151-200` | `sim-k275.R D D3 151 200 k275-D-D3-151-200` | 141 | 151 |
| 5 | `D-D4-1-50` | `sim-k275.R D D4 1 50 k275-D-D4-1-50` | 151 | 161 |
| 6 | `D-D4-51-100` | `sim-k275.R D D4 51 100 k275-D-D4-51-100` | 151 | 161 |
| 7 | `D-D4-101-150` | `sim-k275.R D D4 101 150 k275-D-D4-101-150` | 151 | 161 |
| 8 | `D-D4-151-200` | `sim-k275.R D D4 151 200 k275-D-D4-151-200` | 151 | 161 |
| 9 | `D-D6-1-50` | `sim-k275.R D D6 1 50 k275-D-D6-1-50` | 168 | 178 |
| 10 | `D-D6-51-100` | `sim-k275.R D D6 51 100 k275-D-D6-51-100` | 168 | 178 |
| 11 | `D-D6-101-150` | `sim-k275.R D D6 101 150 k275-D-D6-101-150` | 168 | 178 |
| 12 | `D-D6-151-200` | `sim-k275.R D D6 151 200 k275-D-D6-151-200` | 168 | 178 |
| 13 | `K-D3-1-50` | `sim-k275.R K D3 1 50 k275-K-D3-1-50` | 146 | 156 |
| 14 | `K-D3-51-100` | `sim-k275.R K D3 51 100 k275-K-D3-51-100` | 146 | 156 |
| 15 | `K-D4-1-50` | `sim-k275.R K D4 1 50 k275-K-D4-1-50` | 162 | 172 |
| 16 | `K-D4-51-100` | `sim-k275.R K D4 51 100 k275-K-D4-51-100` | 162 | 172 |
| 17 | `K-D6-1-50` | `sim-k275.R K D6 1 50 k275-K-D6-1-50` | 180 | 190 |
| 18 | `K-D6-51-100` | `sim-k275.R K D6 51 100 k275-K-D6-51-100` | 180 | 190 |
| 19 | `FS-D3-1-50` | `sim-k275.R F,S D3 1 50 k275-FS-D3-1-50` | 91 | 101 |
| 20 | `FS-D3-51-100` | `sim-k275.R F,S D3 51 100 k275-FS-D3-51-100` | 91 | 101 |
| 21 | `FS-D4-1-50` | `sim-k275.R F,S D4 1 50 k275-FS-D4-1-50` | 97 | 107 |
| 22 | `FS-D4-51-100` | `sim-k275.R F,S D4 51 100 k275-FS-D4-51-100` | 97 | 107 |
| 23 | `FS-D6-1-50` | `sim-k275.R F,S D6 1 50 k275-FS-D6-1-50` | 107 | 117 |
| 24 | `FS-D6-51-100` | `sim-k275.R F,S D6 51 100 k275-FS-D6-51-100` | 107 | 117 |
| 25 | `M-D6-1-100` | `sim-k275.R M D6 1 100 k275-M-D6-1-100` | 148 | 158 |
| 26 | `ID-D4-1-20` | `sim-k275.R ID D4 1 20 k275-ID-D4-1-20` | 22 | 32 |
| 27 | `T-D3-1-100` | `sim-k275.R T D3 1 100 k275-T-D3-1-100` | 114 | 124 |
| 28 | `T-D4-1-100` | `sim-k275.R T D4 1 100 k275-T-D4-1-100` | 124 | 134 |
| 29 | `T-D6-1-100` | `sim-k275.R T D6 1 100 k275-T-D6-1-100` | 138 | 148 |
| 30 | `E-D3-1-100` | `sim-k275.R E D3 1 100 k275-E-D3-1-100` | 76 | 86 |
| 31 | `E-D4-1-100` | `sim-k275.R E D4 1 100 k275-E-D4-1-100` | 97 | 107 |
| 32 | `E-D6-1-100` | `sim-k275.R E D6 1 100 k275-E-D6-1-100` | 110 | 120 |

## 4. The order

1. **Tier D's twelve jobs** (1 to 12), then **every tier D continuation, if any** (section 5). Nothing else is sent before tier D is complete.
2. **The interim controls read**, once, only after all tier D continuations: `score-k275.R <tier D rows dir> <out>` with `K275_CONTROLS_ONLY=1` only (it writes `controls-judged.tsv` with PC1 and NC1 and prints only their pass or fail; it computes, writes and prints no outcome row, no tier D cell and no "outcome" line). Its own budget row: 30 job-minutes (limit 20 + 10 set-up). **F1 or F2 is read as: exit 1 AND `controls-judged.tsv` present AND a FAILS label (PC1 or NC1); then every other tier's request is cancelled (not sent), the owner is told, and nothing else is read.** A pass is exit 0 with all six labels `passes`.
3. **Then the remaining tiers**: K (13 to 18), F+S (19 to 24), M (25), ID (26), T (27 to 29), E (30 to 32), each with its continuations as needed. Tier V waits for L1-V.
4. **Final scoring once**, all tiers in one directory (each job directory under its own name, which is its source; tier E's comparators come from tier K; crossNO needs D and E), `K275_G1=docs/qa/workflow-execution/kit/kit-275/g1/g1.tsv`, `K275_BOOT_B` at its default 10,000, `K275_SCORE_R` at its default 200, 100, 20. Timeout 30, ceiling 40 job-minutes.
5. **Procedural rule: no tier D row is opened or summarised before the final scoring**, beyond the continuation line counts of section 5 and the interim controls read of step 2. The primary outcome is not seen early.

## 5. The continuation rule (verbatim) and the missing-unit procedure

**Verbatim (the coordinator's ruling, VAL-382):** a timed-out job's missing data sets are run ONCE; an interrupted data set is discarded whole; no partial re-run of a cell; a data set complete in two sources stops scoring; one continuation per job, inside the reserve and under the running-ceiling rule.

- **Missing units, found without reading any value.** A unit `<tier>-<design>-<r>` is missing when `rows/part-<tier>-<design>-<r>.tsv` is absent from the job's directory (a job that finished has its parts assembled into `rows/<tier>-<design>-<first>-<last>.tsv.gz`) or holds fewer data lines (line count minus the header) than its tier's rows: D 16, K 17, F 5, S 4, V 12, M 7, T 6, E 6, ID 5. Only line counts are used.
- **The continuation**: the original's tiers, design, first and last, with `K275_ONLY_UNITS=<the missing tags, comma-separated>` (tags outside the job stop the run), the job name gains `-cont`; its own out directory (`k275-<job>-cont`) and scoring name; `timeout_minutes` = ceil(1.5 x its units' L1 compute); its ceiling row is taken from the reserve and written in `runs.tsv` first. **One continuation per job.** A continuation that itself times out is not repeated: its complete units count, the rest is reported to the coordinator.
- **Scoring input for a timed-out job**: the part files of the killed job directory stay in the records and are scored as their own source (an incomplete data set is dropped whole and listed as "interrupted"). For F4, a job without `checks-<job>.tsv` contributes its `progress/checks-*.tsv`. The scorer stops on a data set complete in two sources, on a unit row, on a false `unit_completed` check, on an instrument row and on duplicate or out-of-range rows.
- **Pre-scoring line-count check of the checks** (counts only): per complete tier D data set 6 E1 and 1 E2 check rows, i.e. **3,600 E1 and 600 E2 over the 600 tier D data sets** (200 per design); a shortfall is a missing or interrupted unit, handled by the procedure above before scoring.

## 6. What is handed back, and when

- **After tier D and the interim read** (steps 1 and 2 of section 4): to the coordinator, the auditor's reading follows: the tier D completeness (units complete, continuations used, line counts of the checks), the runs.tsv rows, and the controls read's exit status, the label of each of the six PC1 and NC1 cells and `controls-judged.tsv`'s presence; no other tier D content.
- **After the last tier's completeness and BEFORE final scoring**: the completeness table of every job (units, continuations, line counts), the runs.tsv ledger and the running-ceiling figures; final scoring is not run until the coordinator says.
- After final scoring: the outputs of ticket section 8.1, `outcome.md`, `results.md`, and the ledger row `VAL-<merge>` are written (not part of this declaration).

## 7. The check predictions (checks only; no outcome is predicted)

For a job that ran to the end, `checks-<job>.tsv` and its `progress` files hold exactly: **E3** one row per omic the job's units use (lipid; metab for ID), all TRUE, listing `correlation_exclusion` at `within_assay`; for the **tier D jobs** E1 6 rows and E2 1 row per data set (so 300 E1 and 50 E2 per tier D job of 50 data sets), all TRUE; **E12** rows for LG journeys of data sets 1 to 20 only (K: 3 per data set; T: 1 per data set; none elsewhere), all TRUE; no `unit_completed` row and no `instrument` row; `timing-units-<job>.tsv` holds one row per unit with an empty error column; `provenance-<job>.json` has `final: true`, `smoke: false`, `lifecycle: false` and the audited hashes; the assembled `rows/<tier>-<design>-<first>-<last>.tsv.gz` hold D 16, K 17, F 5, S 4, M 7, T 6, E 6, ID 5 rows per data set. A job killed by its timeout has `final: false`, part files, and no timing or final checks file. A check row that differs is recorded and handed back; none is repaired and no counted result is re-run (VAL-54: an unfavourable valid result is a result; a mechanics failure with no result may be re-run only inside the reserve and the four conditions).

## 8. The gates, margins, falsifiers and the outcome table (quoted from KIT-275, sections 7.1 to 7.3, as fixed)

The ticket text below is quoted verbatim from `tickets/kit/KIT-275.md` (sections 7.1, 7.2, 7.3), unchanged since before G1 (the copy read here carries the field `ticket_revision: 1.2`; the coordinator's later revision 1.3 changed the sections 4 and 6.1 emulation reference, R-2, which is applied in `score-k275.R`'s `REF_ARM`; the auditor should confirm that sections 7.1 to 7.3 are unchanged in r1.3). Where the quoted text names data sets or arms it is read with this declaration's job list.

### 7.1 Margins, levels and the multiplicity

Margins (declared now; they are fixed unless the owner changes them before the declaration commit): **M_abs = 0.05** (the filter "fails nothing" if its failure probability is at most 0.05); **M_act = 0.50** ("acts" if it fails the sample in more than half of the data sets); **M_ρ = 0.95** (a standard error at least 5 % too small, a t at least 5.3 % too large, is "understated"); **M_match = 0.03** (the derived ρ is "matched" if the interval lies within the derived value ± 0.03; the three derived values are 0.07 and 0.06 apart); **M_any = 0.10** (any-call rate raised by at least 10 percentage points, the margin VAL-185 and KIT-190 used for the any-call; **a non-inferiority margin used here as a harm threshold, so VAL-344 OQ-2 prints the "any rise" reading beside it**); **M_α = 0.0125** (the per-feature null rate raised by a quarter of the nominal 0.05). Intervals: exact Clopper-Pearson for a single proportion (the decision of VAL-243 OQ-S1 (c)); paired percentile bootstrap over data sets (B = 10,000, type-7 percentiles, the shared index matrices of 5.3) for differences and ratios, with the delete-one jackknife-normal limit beside, a label that differs printed as such. A label is **"not determined"** exactly when none of its bounds holds; no data sets are added, nothing is re-run (VAL-54), and the owner gets the estimates and intervals. A journey that did not complete (`crashed`, `refused`) is never a denominator (except `refused_after_exclusion` in the emulation arms, a result counted in their denominators); a cell with fewer than 0.90 R valid journeys is labelled "not determined (coverage)"; `dropped`, `uncorrected` and `no_da` are results (RUV not applied) and stay in the denominators of the arms where RUV is requested. There is **no tie-break**: no row of the outcome table selects among arms.

| Id | Estimand | Cells | Direction, margin | Level (one-sided unless stated) | Role |
|---|---|---|---|---|---|
| **PC1** | P_pre(0.5), V1g, over the data sets of tier D | D3, D4, D6 | lower limit above 0.90 | 0.01/3 per design | positive control; fails the ticket (F1) |
| **NC1** | P(fail(0.5)) of sample a in V0 **on GK2 after RUV (the drift removed)**, tier D | D3, D4, D6 | upper limit at most M_abs | 0.01/3 per design | negative control; F2 |
| **PC2** | ρ_unmod of the hand-made pseudo-replicate (G1) | X3, X4, X6 | within 0.03 of sqrt((2n-3)/(2n-1)) | two-sided 0.01/3 per design | positive control of the reading; F3 |
| **NC2** | ρ_unmod of genuine data against the generator's truth (G1) | X3, X4, X6 | within 0.03 of 1 | two-sided 0.01/3 per design | negative control of the reading; F3 |
| **A** | P_post(0.5), V1g | 12 cells: D3, D4, D6 x LM2, GK2, LA, GA | **absorbed** if the upper limit is at most M_abs; **acts** if the lower limit is above M_act; **mixed** otherwise | 0.025/12 | primary |
| **B** | ρ_mod, V1g against VD (k - 1) | 6 cells: D3, D4, D6 x (LM2 against LM1; GK2 against GK1) | **understated** if the upper limit is below M_ρ; **not understated** if the lower limit is at least M_ρ; else not determined | 0.025/6 | primary (on the legacy rule-(d) route it measures the increment only, 2.6) |
| **B′** | ρ_unmod **of LX** (legacy test on GK2's corrected matrix, V1g, against the same test on VD's GK1 corrected matrix) against the derived value | 3 cells: D3, D4, D6 | **as derived** if the two-sided interval lies in derived ± M_match; else "differs from derived" with the value | two-sided 0.05/3 | reported label (not in 7.3). LG2's ρ is reported beside (no label); LM2's is **not** compared with the derived value (rule (d), 2.4) |
| **π** (A7) | the imprint retained, V1g | 12 cells: D3, D4, D6 x LM2, GK2, LG2, LX | **imprint retained** if the two-sided 95 % lower limit is above 0.5; **removed** if the upper limit is below 0.5; else not determined | two-sided 95 % per cell | reported label (not in 7.3 and not in the pooled level); printed under 7.3 beside B and C |
| **C** | Δany = P(any \| V1g) - P(any \| V0) | 6 cells: D3, D4, D6 x LM2, GK2 | both readings (VAL-344 OQ-2): **raised (material)** if the lower limit is above M_any; **not raised** if the upper limit is at most M_any; else not determined; and **any rise shown** if the lower limit is above 0; P(any \| V0) and P(any \| VD) printed beside | 0.025/12 (C and C′ together) | primary |
| **C′** | Δα̂ against V0 | the same 6 cells | **co-primary with C** (VAL-344 OQ-2): raised (material) if the lower limit is above M_α, not raised if the upper limit is at most M_α, else not determined; and any rise shown if the lower limit is above 0; V0's and VD's rates beside | 0.025/12 (C and C′ together) | primary |

The two readings of C and of C′ are read off the **same** interval at the same level, so they add no multiplicity. Everything else (the k grid, the other thresholds, V1p, V1s, V1t, V2, tier V, LM2@20, LG2's A, GSC, A2 to A6, B3, the shipped arms' B and C, the emulation arms' H1 to H3, X's checks) is **reported with two-sided 95 % intervals and no label** (the labelled cells are the ones above only). The shipped arms LA and GA are labelled in family A (so a user's actual route is read) but the total outcome (7.3) classifies on the fixed-k arms LM2 and GK2, where k is not a data-dependent choice.

### 7.2 Falsifiers and the pooled level

- **F1** PC1 fails in any design (the planted anomaly does not behave: wrong assay, not planted, mis-scaled).
- **F2** NC1 fails in any design (the genuine sample is failed by the filter after RUV in more than the margin: a data or instrument defect). **The same sample's pre-RUV failure rate in V0 (the NO arm) is a reported result, not a falsifier** (6.1 A4; indicatively about 8 % at D3 because of the generator's own drift, AUDIT-KIT-275-plan P04), and so is its post-RUV failure rate on the other RUV arms (A1 on V0).
- **F3** PC2 or NC2 fails in any design (the standard-error reading is wrong).
- **F4** an identity check fails in any counted journey (E1, E2, E5 of 8.3). After scoring has begun a failure is a result and stops the ticket; before any estimand is produced, a mechanics failure is an instrument defect under VAL-54 (the four conditions).
- **Level.** Twelve tests (PC1, NC1, PC2, NC2, each in three designs) each at 0.01/3 (one-sided for PC1 and NC1; two-sided for PC2 and NC2 on the 20 fixture data sets), so when the instrument is right the probability that any falsifier fires is at most 12 x 0.01/3 = 0.04 by the union bound; in fact far less, because the right values sit far from the limits: PC1 expects 1 against a limit near 0.90 (it fires only with 9 or more non-failures in 200, impossible at the true value 1 and very unlikely at 0.97), NC1 (on GK2, the drift removed) expects 0 against 0.05, and PC2 and NC2 have a sampling SD of about 0.002 against a tolerance of 0.03. 7.4 states how many deviations each tolerates.
- **Stop rules.** G0; F1, F2, F3, F4: stop, the owner is told and nothing else is read. A counted job that does not complete is "did not complete" for its seeds (never a denominator), continued once inside the reserve. An unfavourable or inconclusive valid result is not re-run (VAL-54). No data sets are added under this declaration.

### 7.3 The total outcome table (the first row that holds applies; fixed-k cells only)

| Outcome | Condition | What goes to the owner |
|---|---|---|
| **Stop** | G0, F1, F2, F3 or F4 | the falsifier and its numbers; nothing else is read |
| **R1 absorbed on both routes** | all six cells D3, D4, D6 x LM2, GK2 of family A are "absorbed" | the filter after RUV cannot fail a single gross anomaly at k = T + 1 on either route; the consequence readings by route (B, C, C′, π); P_pre and the pre-RUV false fails beside, and the emulation arms' H1 to H3 (what option (b) of 7.5 would buy and cost) |
| **R2 absorbed on the declared route only** | the three GK2 cells "absorbed", not all three LM2 cells "absorbed" | the legacy route's behaviour by design with κ_c, PR and LG2 (2.4's mechanism observed or not); the declared route's consequence |
| **R3 absorbed on the legacy route only** | the three LM2 cells "absorbed", not all three GK2 cells | same, reversed; the owner is told this is not what the derivation expects |
| **R4 acts on both routes** | all six cells "acts" | the derivation is refuted for these designs: the filter does act after RUV; P_post by threshold; the sentence adopted under VAL-344 OQ-A would be wrong (NEWS, 11) |
| **R5 mixed** | anything else (including any cell "mixed") | the estimates by cell, k, κ_c, PR and π; no label is stretched |

The consequence labels (B, B′, π, C, C′) are printed per route and design under every row except Stop: "legacy standard errors understated and as derived (LX)", "understated, differs from derived", "not understated", "not determined"; "imprint retained / removed / not determined"; "false-positive rate raised (material) / any rise shown / not raised / not determined". They are read unconditionally (a legacy cell whose filter acts in many data sets can still be "not understated"). The reported (unlabelled) results are a part of the owner's table whatever the row.


