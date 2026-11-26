# KIT-262 M3: E0 statement (what each script takes unchanged, what is new, with the marks)

Checked locally by `e0-m3-k262.py` (standard library; 35 checks, all PASS at this draft; run from the repository root, output in `m3-local-checks.txt`). Rule of the check: each line of a copied function in the destination is a line of the source function (same text after stripping, in the source's order) or carries a mark `# [X]`; each source line the destination lacks is counted against the marked lines. "Verbatim" below means the check found 0 marked lines and 0 replaced lines. The sources are on the branch (it holds the binding's tree): KIT-190's `effect.R`, `decide.R`, `bounds.R`, `reestimate.R`; `tests/testthat/helper-kit172.R`; KIT-262's own `m2-k262.R`, `rules-k262.R`.

## `m3-k262.R` (the simulation)

| Block | Source | Taken | Marks |
|---|---|---|---|
| `e2PickCells` | `kit-190/effect.R` | verbatim (8 of 8 lines) | none |
| `e2Values` | `kit-190/effect.R` (VAL-235's generator, `helper-kit172.R:34-62` with KIT-190's two added arguments) | 37 of 37 lines verbatim | **[L]** one added line: `level = ...` in the returned list (the abundance thirds of the calibration need the true level); no draw moves |
| `E2_CELLS`, `ALIGNED_CELLS`, `e2Data` | `kit-190/effect.R` | the two constants verbatim; `e2Data` 19 of 20 lines verbatim, the seed rule `190000000 + 10000 x idx + r` unchanged | **[L]** `level = sim$level` added to the returned list (one list line replaced by two) |
| `k190Oracle`, `k190Candidate`, `CUTOFF` | `kit-190/effect.R` | verbatim (26, 8 lines; the constant) | none |
| `m2Design`, `m2Object`, `M2_RUV_PARAMS`, `m2Prepare`, `m2Sweep` | `m2-k262.R` | verbatim (6, 16, 2, 7, 3 lines); used only by E10 and the stress design | none |
| `m3Object` | `helper-kit172.R` `kit172Object` | 13 of 15 lines | **[D]** the design is a parameter (signature line; the source's `design <- kit172Design()` removed, stated in a marked comment line) |
| `m3StateMatrix` | `helper-kit172.R` `kit172StateMatrix` | 5 of 7 lines | **[D]** the runs are a parameter (signature; the source's `runs <- kit172Design()$Run` removed) |
| `m3StressValues` | `m2Values` (`m2-k262.R`, itself a declared generalisation of `e2Values`) and `e2Values` | every line not marked is a line of one of them (comments aside); the draws before the mask are M2's, in M2's order | **[N]** the floor order is **named and fixed**: `M3_STRESS_FLOOR_ORDER <- "after"`, `floor_order = M3_STRESS_FLOOR_ORDER` in the signature, `stopifnot(identical(floor_order, "after"))`; the mask first, then `e2PickCells` on the observed cells, then `values[ones] <- 1` (KIT-172's order, VAL-381's ruling; the check asserts the order of the three lines); **[H]** `effect_share` (0 for S-N, 0.1 for S-E: e2Values' own argument form); **[I]** MCAR only (the abundance-dependent branch removed; the stress cells are MCAR 10 %); **[L]** the level returned; **[P]** the design fixed at n = 24, g = 2, 300 features, rate 0.10 |
| `kit172Quiet` in the helper environment | `helper-kit172.R:105` | replaced by `m3Quiet` | **[M]** the helper's `suppressWarnings(suppressMessages(.))` would hide a partial-match warning from the guard; the replacement is a calling handler that stops on a partial match and muffles the rest. No other helper function is changed |
| `rules-k262.R` | M2 | **sourced unchanged** (SHA-256 9547d350...; unchanged since M1): `k262Rules` with `curves = NULL` (no hook), `k262ParseReason`, `k262Cause` unused | none; M3 adds no rule |
| `kit172Journey`, `kit172LegacyRuv`, `kit172Object`, `kit172StateMatrix` (eight cells), `kit172DaTable`, `kit172Design` | `tests/testthat/helper-kit172.R`, read through `testthat::source_test_helpers` | **called, not copied** (KIT-190's route) | none |

New in M3, not copied from anywhere (audited by reading and by the self-test, not by E0): the unit (`m3Unit`: the P journey, the rules from its trace, C1's Manual journey where its choice differs, E1, E7, E9, E10, E11, E8, F2), the endpoint extraction (`m3Endpoints`: the step's `ruv_applied` and reason, the DA calls, the calibration counts by abundance third), the committed-D reader (`m3CommittedD`), the row format (70 columns) and its one-line append, the local partial-match guard (`m3GuardMuffle`), the unit wrapper and its trace files, the plan and continuation settings (`K262_M3_*`), the planted faults, the lifecycle's scorer rehearsal. M3's own differences from M2's script: **no hook** (`k262WithCurveCapture` is not used; C2 is not an arm), the unit is a data set with journeys instead of a sweep, rows are appended by the unit, and the guard is local.

## `score-m3-k262.R` (the scorer)

| Block | Source | Taken | Marks |
|---|---|---|---|
| `tangoZ`, `bergerBoos`, `skewness`, `nullTest` | `kit-190/decide.R` | verbatim (7, 16, 5, 7 lines) | none |
| `effectTest` | `kit-190/decide.R` | 37 of 38 lines | **[S]** the sign-flip seed is this ticket's (`SIGNFLIP_SEED` 262000002) in the signature; Johnson's comparison stays `t_johnson < -Z` with Z = qnorm(1 - 0.05), the declared alpha (O-1 (ii): 0.05, no split) |
| `tangoZb` | `kit-190/bounds.R` `tangoZ` (the version that allows negative margins) | the body verbatim; renamed | **[S]** name only. The upper bound is `bounds.R`'s form (inverted Tango, at 1 - alpha); the lower bound is its mirror, used only to label a not-NW cell |
| `ALPHA`, `M_NULL`, `M_EFFECT` | `decide.R` | the same values (0.05, 0.10, 0.05) | none |

New: the reader and its stops (instrument rows, a unit in two sources, a duplicate, a row outside the plan, an unknown cell, a seed that is not the formula, a header that is not the row format; a short last line is dropped), the blind scorer, the verdicts and their labels (`NW`, `shown worse`, `inconclusive`, `not determined`), 6.3's figures (changed data sets with Newcombe's and the paired-t intervals; calibration by third with the percentile bootstrap; power; coverage), 6.4's table, E10's McNemar counts, the outcome (Stop, Stop (E11), K1, K4), the `missing` mode, the self-test (50 cases; expected values read from `e5hand-m3-expected.tsv`, which `e5hand-m3-k262.py` worked in Python). Newcombe's method 10 and the correlation term are implemented from the paper's formulae in both R and Python; the auditor should check that formula against the source (INFERRED from memory of the paper, not read here).

## `reestimate-m3-k262.R`

A copy of `kit-190/reestimate.R`: `estimate` 25 of 27 lines and `writeEstimate`, `SCHEMA` verbatim. **[R]** marks: the cells (`S-N`, `S-E`), the arms (`P` the reference, `C1` the candidate), the floors (250, 200), `LIMIT` 1,000, `ROUTE_ABOVE` 400, `SPREAD` 1.93 (KIT-190 used 1.25; M2's counted jobs ran up to 1.93 times the lifecycle's sweep seconds), the target 0.80^(1/2). **The ceiling loop of `plan` is removed** (the check asserts that `t - 0.0005` and `while (core` are absent): R is never reduced, and above 400 the plan only says so.

## What M3 takes from M2's records, not its code

The grid-independent machinery of M2 that M3 does not use: the abundance-dependent mask, `solveLod`, the C2 hook, the per-percentage rows file, the 198-cell grid. The seeds of M3 do not overlap M2's (262100000 + 100 x cell + r), the re-measurement's (263100000 + 100 x cell + r) or the lifecycles': the stress cells use `262200000 + 10000 x idx + r` with idx 9 and 10 (E9's disjointness clause is `SEED_RANGE` in the checks); the eight cells reuse KIT-190's seeds `190000000 + 10000 x idx + r` by design (E9 compares them with KIT-190's D rows).
