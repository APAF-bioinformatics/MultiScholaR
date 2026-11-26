# KIT-262 M3: lifecycle L3 declaration (DRAFT for the audit; to be committed alone before the lifecycle is sent; nothing is sent from this draft)

**Why.** M3 of KIT-262 (ticket section 6): the present rule P against candidate C1 on VAL-235's eight cells and the two 24-sample stress cells. L3 is its blind rehearsal (VAL-11; ticket 9.4): it proves that the instrument observes the intended work through completion and failure, and it measures the seconds the counted plan is made from. **L3 reads no estimate.** Owner decisions in force: VAL-302 (O-1 (ii): M3 tests P and C1 only, not C2; O-2 (a): "no worse" on the cell average; O-3: the stress cells and the skip rule, yes; O-4 (a): no 2 x 2 error rates), VAL-381 (stress cell unchanged: M2's cell 79 with KIT-172's floor order), VAL-390 (ticket ceiling 250.32 core-hours; M3's plan is made from L3's own timings and returns to the owner if it exceeds the ceiling). No package code, no NEWS, no version change.

**Binding.** `package_commit` of the L3 request is the commit that adds this file to branch `kit-262`. It is `origin/v0.9-staging` a3c69901ea9b50b8de875d79ca89fb872912c092 (KIT-190, KIT-242 Part B with VAL-373 and VAL-290, KIT-260) plus `docs/qa/workflow-execution/kit/kit-262/` only; `git diff` over `R/`, `src/`, `inst/`, `DESCRIPTION`, `NAMESPACE` from that staging commit to the declaring commit is empty (`m3-binding-a7.md`, which also holds the A7 statement). Script SHA-256, recorded at the declaring commit (the values below are the draft's): see section 8.

## 1. The scripts, in one request, in this order; every timeout explicit

| # | Script and arguments (from the package root) | Environment | Plan (expected wall, at spread 1.0 / 1.93) | `timeout_minutes` |
|---|---|---|---|---|
| S1 | `reestimate-m3-k262.R selftest out/kit-262/l3/s1` | | 1 min | 5 |
| S2 | `score-m3-k262.R selftest out/kit-262/l3/s2 docs/qa/workflow-execution/kit/kit-262/e5hand-m3-expected.tsv` | | about 5 min | 20 |
| S3 | `m3-k262.R out/kit-262/l3/s3` (every cell, both replicates, every arm forced, the E9 units, E8 on one replicate per cell, the scorer rehearsal) | `K262_LIFECYCLE=1`, `K262_M3_CELLS=N-r,N-s,E-r,E-s,E-p,E-g,N-w,E-w,S-N,S-E`, `K262_M3_REPS=9001,9002`, `K262_M3_E9REPS=1,2`, `K262_M3_FORCE_ARMS=1`, `K262_M3_E8=1`, `K262_M3_E8_MAX=9001`, `K262_M3_WORKERS=4` | 14.0 / 27.0 min (3,363 CPU-s at 4 workers) | 55 |
| S4 | `m3-k262.R out/kit-262/l3/s4` (the failure check: four planted faults and one clean unit) | `K262_LIFECYCLE=1`, `K262_M3_FAULT=1`, `K262_M3_CELLS=N-r`, `K262_M3_REPS=9101:9105`, `K262_M3_E8=0`, `K262_M3_WORKERS=4` | 2 / 3 min | 10 |
| S5 | `m3-k262.R out/kit-262/l3/s5` (stage 1 of the re-estimation, rehearsed on 5 data sets per stress cell) | `K262_LIFECYCLE=1`, `K262_M3_STAGE1=1`, `K262_M3_REPS=9011:9015`, `K262_M3_E8=0`, `K262_M3_WORKERS=4` | 6.5 / 12.6 min (1,568 CPU-s) | 30 |

Every timeout is at least twice the plan at the larger spread (S3: 2 x 27.0 = 54, rounded up to 55, five minutes above the ticket's 50; S5: 2 x 12.6 = 25.2, rounded up to 30). The sum of the timeouts (120 script-minutes) is a limit for a hung script, not the plan; see section 6. Each script writes under `out/kit-262/l3/sN`: `checks.tsv`, `timing.tsv`, `hashes.tsv` and its console; **the row files stay under `tempdir()`** (S3 to S5). No share, rate, k distribution, FDP, power, calibration, interval or test statistic leaves L3. The console of S3 to S5 shows per unit only its tag, its status and its seconds; S5 shows the stage-1 estimate's cell, endpoint and pair count and never its value (its file is removed). `out_pattern` `\.(tsv|txt)$`.

## 2. What each script does

- **S1, S2.** The self-tests of the re-estimation (`reestimate-m3-k262.R`: planted inputs give byte-identical estimates; the plan's R values equal the hand values of `e5hand-m3-k262.py`: effect SD 0.265 gives 250, 0.35 gives 450, 0.40 gives 550; null share 0.30 gives 300; the floors 200 and 250; R is never reduced) and of the scorer (`score-m3-k262.R selftest`: constructed cases for every path, the expected values read from `e5hand-m3-expected.tsv`, which Python worked independently).
- **S3.** For each cell and replicate (the eight cells at r = 9001, 9002; the stress cells likewise) one unit: the P journey (`kit172Journey` with KIT-172's legacy RUV, Automatic), the rules read from its trace (`rules-k262.R`: P and C1 only, no hook), C1's Manual journey at C1's percentage and k where its choice differs; every arm forced (`K262_M3_FORCE_ARMS=1`: where C1 re-uses P's row, the Manual journey at P's own choice is run as E11). The 16 E9 units (the eight cells at r = 1, 2, on KIT-190's seeds, P only) compare P's percentage and sweep k with KIT-190's committed D rows and P's k used with its D+B rows. The stress cells also run the direct M2-form sweep before between-sample normalisation (E10). E8: the unit through the metabolomics class at r = 9001. After the units, the lifecycle runs the blind scorer twice on its own rows (E6), checks every estimand in it is NA, and runs the full scorer once on the rows (so the format the simulation writes is the format the scorer reads); only pass or fail leaves.
- **S4.** Units 1 to 4 carry planted faults: 1 a journey that stops (a data-set error row), 2 a partial match inside a journey, 3 a partial match in unit code outside a journey (both instrument rows with the class `partial_match_guard`, a trace file each), 4 a worker that dies with no row and no trace (an instrument row written by the parent). Unit 5 is clean.
- **S5.** The stage-1 mechanics: 10 stress units, the capped rows to `reestimate-m3-k262.R estimate`, the schema and pair counts checked.

## 3. Predictions (mechanics and controls only; every one derived from the code at the binding under the shipped settings, or marked OBSERVE)

| Item | Prediction | Basis |
|---|---|---|
| S1 | every case `ok` | `e5hand-m3-k262.py` and the formula of KIT-190 section 8 at k_cells = 2 |
| S2 | every case passes (the checks of `checks.tsv`) | constructed inputs; Python hand values |
| S3 COMPLETION | 36 rows of 36 planned (20 main, 16 E9); INSTRUMENT_ROWS 0 | the unit appends its line before it returns; the parent writes a row for any unit without one |
| S3 E1 | 36 of 36: the rules' usable set, composites and winner equal the shipped trace and the step's `ruv_applied`, k used = min(k sweep, cap) | `k262Rules` is M2's (E1 passed in 7,920 + 3,520 data sets); `omicsRuvKBound` `utils_workflow_omics_declarations.R:1052-1087` |
| S3 E7 (reduced) | 36 of 36: P within C1 row by row; C1's k used at most the cap; RUV under P implies RUV under C1 | ticket 3 (nesting), M2's E7; C2 and the curves are not used in M3 |
| S3 E9 | 16 of 16 compared and equal: P's percentage and sweep k at the D row's, k used at the D+B row's | the ticket's F3 expectation: same generator and seeds (E0), the sweep's R3 to R6 and ranking not edited (A7), rule (d) shipped by KIT-190; **OBSERVE**: this is the first comparison of the head with KIT-190's D tree (a1c12d28). A difference is not an instrument defect: it is what F3 is for. If L3 shows one, no counted job is declared; the numbers go to the coordinator |
| S3 E11 | every run passes (a Manual journey at P's own percentage and k gives P's DA table field for field); the number run is the number of re-used units (**OBSERVE**) | the bounded arms of KIT-190 used this route (audit `:152-156`); `omicsRuvCorrectByCause`: with nothing left out the result is the earlier rule's |
| S3 F2 | every Manual arm passes (percentage, k and controls as declared and as rule (d)'s oracle at that percentage) | KIT-190's `oracle_ok` (falsifier `oracle` 0 in its bounded arms) |
| S3 E8 | 10 of 10 identical rows between lipidomics and metabolomics (all fields but the class and the seconds) | KIT-190's metabolomics identity check; M1 and M2's E8 |
| S3 E10 | 4 of 4 stress units ran the direct sweep | the comparison itself is the scorer's (McNemar counts) and is not read in L3 |
| S3 seeds | SEEDS 36 of 36; SEED_RANGE 4 of 4 (stress seeds 262299001 to 262309002, above M2's 262119895 and below the re-measurement's 263100000); FLOOR_ORDER_STRESS 1 of 1 | `m3Seed`; `e0-m3-k262.py` |
| S3 scorer rehearsal | E6_BLIND_SCORER, BLIND_ESTIMANDS_NA, FULL_SCORER_ON_ROWS each 1 of 1 | `score-m3-k262.R`'s self-test |
| S3 error rows | 0 (ERROR_ROWS_REPORTED 36 of 36) | KIT-190: 0 journey errors in the eight cells; **OBSERVE** for the stress cells (n = 24, no standards, MCAR 10 %): a stopped stress journey is a finding to hand back, not something to re-run |
| S4 | FAULT_ROWS 5 of 5: unit 1 `error`, units 2 and 3 `instrument` with class `partial_match_guard`, unit 4 `instrument` with class `no_row`, the other unit `ok`; four trace files for units 1 to 3 and no trace for unit 4 | the guard is a local calling handler: the stop reaches the enclosing `tryCatch` (KIT-275 CC2) |
| S5 | STAGE1_ESTIMATE_SCHEMA 1 of 1 (columns cell, endpoint, n_pairs, value; both cells present; n_pairs 5 each if all units are ok) | `reestimate-m3-k262.R estimate` |
| shapes (OBSERVE) | the step's value carries `optimisation`, `ruv_applied`, `ruv_not_applied_reason`, `best_k`, `k_sweep`, `k_cap` per assay `LCMS_Pos`; the normalised state exists as `<prefix>_normalized`; the DA table has `feature`, `raw_pvalue`, `fdr_qvalue` | read in the code at the binding (`utils_workflow_omics_declarations.R:884-950`; `helper-kit172.R:141-146`); first run of these accessors: a NA or empty read shows as an E1 or E9 failure with a mechanics failure in the log, an instrument defect under VAL-54 |

**Not predicted:** the share of units in which C1 re-uses, rescues or re-ranks P; any FDP, any-call rate, power, calibration or k; the number of E11 runs; every second.

**KIT-242's behaviours on the route (cited, A7):** an assay kept without RUV with the reason `no_usable_controls` or `no_factor_supported` (`utils_workflow_omics_declarations.R:959-1018`) is P's "no RUV" in M3 (`ruv_applied` FALSE, the reason kept in `p_reason`, the DA on the uncorrected assay); the sweep status `skipped_no_controls_selected` (`func_lipid_norm_ruv_optimization_helpers.R:190-209`) is a rejected percentage: its composite is NA (`:371-372`), `t_n_rejected` counts every status but `ok`, and under both P and C1 it is unusable. The scripts count both per data set; whether either occurs is OBSERVE.

## 4. What L3 measures for the plan, and the re-estimation of R

1. **Seconds.** `timing.tsv` of S3 gives seconds per data set by cell (all arms together; `seconds_p` for the P journey alone) with four workers side by side; S5 gives the stress cells' seconds in the stage-1 form. The counted jobs' plan is recomputed from these (G2): per cell, mean seconds per data set x units x 1.93 / 4 workers, jobs cut at 150 minutes of plan, `timeout_minutes` 2 x plan (`m3-continuation.md` section 3). The ticket's basis (39.2 s Automatic, 17.8 s Manual at n = 12, doubled at n = 24) is replaced by L3's.
2. **The blinded re-estimation (VAL-185's form; ticket 6.1).** Stage 1: `m3-k262.R` with `K262_M3_STAGE1=1`, cells S-N and S-E, replicates 1 to 50, P and C1 paired; rows private; `out` holds `reestimate.tsv` and `provenance.json` only. `reestimate-m3-k262.R plan` then gives R per stress cell: the smallest multiple of 50 at or above the ticket's minimum (S-N 250, S-E 200), at most 1,000, with P_c(R) at least 0.80^(1/2) at one-sided 0.05 (effect: Phi(0.05 sqrt(R) / SD - z); null: Phi(0.10 sqrt(R) / sqrt(discordant share) - z)). N-r and N-s are not re-estimated (C2 is not tested; k_cells = 2).
3. **The route fixed now if a re-estimated R exceeds 400** (the ticket's "R doubled"): the plan summary says `above_400 TRUE`; nothing is sent; the formula's R, the projected core-hours (spent + sum of R x stage 1's median seconds per unit x 1.93 / 3,600) and the ceiling line go to the coordinator under G2; the ceiling is restated before any M3 job. **R is never reduced**: `reestimate-m3-k262.R` has no loop that lowers the target (KIT-190's has one; it is removed and E0 checks it is absent), and a cheaper stress design is not an option here. Planning numbers (the audit's, R about 3,354 x SD^2 at this target): 235 at SD 0.265, 411 at 0.35, 536 at 0.40, so R exceeds 400 above an SD of about 0.35.

## 5. Pass rule and what L3 decides

L3 passes when every row of every `checks.tsv` is `pass` for the hard checks (E1, E7, E9 as OBSERVE above, F2, COMPLETION, INSTRUMENT_ROWS, SEEDS, SEED_RANGE, FLOOR_ORDER_STRESS, FAULT_ROWS, the scorer rehearsal, STAGE1_ESTIMATE_SCHEMA), S1 and S2 pass, and the 16 E9 units agree with KIT-190's D rows. E8, E10 and E11 are reported (a failed E11 means the ticket's remedy, every candidate row run, is costed and put to the coordinator; a failed E8 means lipidomics decides). It decides only whether the instrument may be used and what the counted plan is; it decides nothing about C1 or P. A failure whose log shows a mechanics failure and no result judged is an instrument defect (VAL-54: fixed, the same request sent again, at most three further runs for the item, recorded in `runs.tsv`); a failure that is a comparison of the product with its committed past (E9, F2, E7) is a result and goes to the coordinator.

## 6. Budget and the running-ceiling line (core-hours = job-minutes x 4 / 60; set-up 10 job-minutes)

| Row | Ticket line | Expected (spread 1.0) | At spread 1.93 |
|---|---|---|---|
| L3 (S1 to S5) | expected 2.00, ceiling 4.00 (60 job-minutes) | 38.5 job-minutes = 2.57 | 58.7 job-minutes = 3.91 |
| Stage 1 (100 stress units, rows private) | not a ticket line (the ticket's stress ceiling "R doubled" holds headroom, not this) | 60.5 = 4.03 | 107.5 = 7.16 |
| M3 counted jobs (O-1 (ii), O-4 (a), the stress cells at their minima) | expected 43.82; ceiling 128.03 is the widest option set (C2, the 2 x 2 pair, the fallback sweep) | 657.25 job-min = 43.82 | at R doubled for the stress cells and spread 1.93: 1,730.6 job-min = 115.37 |
| Scoring | 2.67 line, 0.80 used, 1.87 left | 1.00 | 1.87 |

Arithmetic (stdlib Python, `round(m * 4 / 60, 2)`): L3 expected 1 + 5 + 14.01 + 2 + 6.53 + 10 = 38.55 min; at 1.93, 1 + 5 + 27.05 + 3 + 12.61 + 10 = 58.66 min (S3 3,363 CPU-s: eight cells 16 units x (39.2 + 17.8 + 17.8) s, plus E8 on 8, stress 4 x (114 + 7.2 + 35.6) s plus E8 on 2, E9 16 x 39.2 s; S5 10 x 156.8 s). `m3-budget-k262.py` reproduces every figure here. Stage 1: 100 x (114 + 7.2) = 12,120 CPU-s / 240 = 50.5 min + 10 = 60.5 min; at 1.93, 107.5 min. M3 expected 145,740 CPU-s (1,550 x 57 + 450 x 121.2 + 2,850) / 240 + 5 x 10 = 657.25 min; at 1.93 and R doubled (900 stress units) 386,540 CPU-s / 240 + 12 x 10 = 1,730.6 min.

**The ticket's L3 line (ceiling 4.00) holds at the larger spread (3.91, 0.09 to spare); L3's expected 2.57 is above the ticket's 2.00 because of the forced arms (E11 on every re-used unit), the 16 E9 units, E8 and the stage-1 rehearsal, which the ticket's 30 job-minutes did not itemise; the explicit timeouts (120 script-minutes with set-up 130 job-minutes = 8.67) are limits for a hung script.** A run that reaches a timeout is an instrument matter charged to the instrument reserve (7.43 left on paper) and counted in the line below.

**Running-ceiling line against the ticket's 250.32 (VAL-390):** used 87.24 (M1 0.88, L2 5.56, M2 jobs 46.77, continuation 4.56, M2 scoring 0.80, re-measurement 28.68, rounded; recomputed from `runs.tsv`); remaining 163.08. After L3 at its ceiling 91.24; after stage 1 at spread 1.93 98.40; after M3's counted jobs at their ceiling for this option set 213.77; after scoring 215.64; the two reserves on paper (correction 8.77, instrument 7.43) bring it to 231.84, **18.48 below 250.32**. M3's plan is made from L3's own timings and returns to the owner if it exceeds the ceiling.

**Attempts (VAL-54).** L3 once; repeated only after a script change for an instrument defect, inside the instrument reserve, at most three further runs. Stage 1, the counted jobs and the continuation are not part of this declaration.

## 7. What is taken unchanged from M2 and what is new

`m3-e0-statement.md` (E0 for M3: each block, its source, its marks, the check that proves it).

## 8. Script SHA-256 (draft values, recomputed at the declaring commit)

See `m3-script-hashes.txt` at the declaring commit.
