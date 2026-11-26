# KIT-281 lifecycle record: the scripts, the self-test, the job plan and the budget (scripts only; nothing has run)

Status: **the scripts of the calibration gate are written and self-checked locally; nothing has been run or sent** (no R on this machine, no smoke, no lifecycle, no counted run,
no `gh`). The ticket's ceiling (now 74.67 core-hours with the S09 send; not yet accepted by the owner) awaits the owner; the scripts get an independent audit first. Part C (the default change) is not here.
Branch `kit-281` from `origin/v0.9-staging` a5de1c2e; no change under `R/`, `tests/`, `inst/`, `man/`.

## 1. Files

| File | What |
|---|---|
| `sim-k281.R` | the simulation: KIT-211's script with exactly the ticket's changes (`script-diff.md`) |
| `score-k281.R` | the scorer: KIT-211's reused functions, adapted input, the new decision section |
| `l1-k281.R` | the blind lifecycle L1: E1', E2 to E10, E4 (the self-test), E5, E6, G0, G2, G3 (the one count, rule G3, its tables) |
| `e0-k281.py` | E0, local: verbatim blocks and named functions against KIT-211's and KIT-197's, located by markers |
| `check-request-k281.py` | refuses a counted request whose R differs from rule G3's (and one that does not cite the table) |
| `selftest/` | `make-selftest.py` (constructed tables and every prediction, Python stdlib), `cases/`, `cases.tsv`, `expected.tsv`, `g3-cases.tsv`, `requests/`, `requests.tsv`, `check-selftest.py`, `run-selftest.sh` |
| `make-script-diff.py`, `script-diff.md` | the list of every changed line range against KIT-211's files |
| `lifecycle.md` | this file |

SHA-256: `sim-k281.R` 16d67e46e65ccbbd606c25cfbe4f7f59f0f0f4aaefec8aaada42ee2e545d6b27; `score-k281.R` 4e543a78318f6f65e7d9fae50bb8c72301bfd029a87e76ba4bc3a49f08730898; `l1-k281.R` 9a899c6112bb91d463ca935ed2b6012444135d568c609ed1e0d3aec200979251 (at this commit; the declaration records the final ones).

## 2. Requirement table (ticket requirement to script lines)

| Ticket | Requirement | Where |
|---|---|---|
| 4, 6 | arms S, B, B1, N | sim-k281.R:296 |
| 4, 6 | twelve settings: A to D at vd5, vd10, vd25 (no mcar, no Av, no Ag) | sim-k281.R:290; score-k281.R:32; score-k281.R:36 |
| 6 | `g_added` subset (H_B1 without H_B, arm B1 only); `kept_b`/`in_b` | sim-k281.R:299; sim-k281.R:368; sim-k281.R:397; sim-k281.R:427 |
| 6 | `kept$B1` unconditional; `kept$A`, `kept$Bp` removed | sim-k281.R:463 |
| 6, 7 | E2 at 2 and 1 group, unconditional; E9 removed; four checks in the summary | sim-k281.R:517; sim-k281.R:687 |
| 7 | E3-type identity: B1's DA input is `masked[H_B1, ]` (not in the ticket's list of changed lines: see section 9) | sim-k281.R:494 |
| 6 (a) | provenance written before the first task | sim-k281.R:647 |
| 6 (b), lesson 3 | one line per completed task (cell, r, seconds) in `progress.tsv` | sim-k281.R:639; sim-k281.R:659 |
| 6 (c), K281-P05, lesson 2 | the partial-match guard around the script's own table handling only (never package, limma or qvalue calls); `options(warn = 2)` not used | sim-k281.R:571; used at sim-k281.R:683, sim-k281.R:700, sim-k281.R:711, sim-k281.R:728 |
| 3 (commit rule) | the bound tree carries no code of KIT-280: `threshold_values_in_filter` absent from the replay marker table | sim-k281.R:589; sim-k281.R:652 |
| 8 (G3) | a counted job is refused unless its range lies inside the G3 table's r = 401 to 400 + R | sim-k281.R:615 |
| scorer | a smoke flag has a producer (`K281_SMOKE=1`) and the scorer refuses it | sim-k281.R:611; score-k281.R:159 |
| 5, 6 | new atoms `n0_gadd`, `p05_gadd`, `*_l0`, `*_l1`, `bhV_whole`, `bhT_whole` | score-k281.R:297-304 |
| 5, 6 | new stats `alpha_gadd`, `alpha_gadd_l0/l1`, `alpha_whole`, `fdr_bh` | score-k281.R:324 |
| 5, 6 | `rhog.B1`, `rhog_l0.B1`, `rhog_l1.B1`, `dfdrbh.B1`, `dpowbh.B1`, `rhoratio.B1`, `Phi.S`, `Phi.B`, `Phi.B1` | score-k281.R:362-386 |
| 5 | the deciding limit = the bootstrap 0.975 quantile (type 7) of `rhog.B1` over the twelve settings' pool; state `pass` iff upper <= 1.10 as written; jackknife-normal beside, deciding nothing | score-k281.R:480; score-k281.R:621 |
| 5 | non-finite resamples dropped and counted (`boot_nonfinite`) | resampleStats (identical to KIT-211's, E0) and the `boot_nonfinite` column of `decision.tsv` |
| 5 | guards G1a (FDP_BH(B1) <= 0.05), G1b (paired <= 0.0025), G2 (power cost <= 0.005) | score-k281.R:622-624 |
| 5, 9 | F5': pooled two-sided 95 % interval of `B1.m6_l_H` excludes 0 | score-k281.R:491 |
| 5 | positive control reported (Phi_S lower limit against 1.20), decides nothing | score-k281.R:632 |
| 5 | per design, per setting (labels at one-sided 0.05/12), by flag, Phi: reported | score-k281.R:637; score-k281.R:475; `designs.tsv`, `settings.tsv`, `pooled.tsv`, `decision-table.md` at score-k281.R:699 |
| 10 | outcome rule: Stop, Withheld, Calibrated, Not calibrated (the first row that holds) | score-k281.R:647 |
| 10 | Stop: F1, F2, instrument statuses (a stop writes `outcome.tsv` and exits 1) | score-k281.R:464; score-k281.R:713 |
| coordinator S02 | a NON-FINITE ESTIMATE of rho_g or of a guard is an instrument Stop ("could not be formed"); resample-level non-finite values stay dropped and counted | score-k281.R:628 |
| coordinator S03 | the count of rows whose recorded q-method is not qvalue is reported (per cell, mask, arm; checks; outcome notes), never a stop | score-k281.R:559 |
| 7, 8 | G3 as code: the one count (`g_added`, arm B1, null, `n_tested`, mean over L1's two replicates), N-hat = 200 x the sum, R by the rule, one R for all cells; `g3-count.tsv`, `g3-rule.tsv` | l1-k281.R:274; l1-k281.R:282; tables written after it |
| 8 | the rule function against Python's Fraction table (24 values, both edges of every branch) | l1-k281.R:281; `selftest/g3-cases.tsv` |
| 8 | a check refusing a counted request whose R differs from the rule's (and one that does not cite the table's hash); the scorer's declared range | `check-request-k281.py`; score-k281.R:127; sim-k281.R:615 |
| 7 | E1': S and B filter-dependent columns at r = 201 to 205, cells A to D, twelve vd settings = KIT-211's committed rows exactly (120 units) | l1-k281.R:182 |
| 7 | E2, E3, E7, E10 (F1 identities, logFC identity); E8 determinism; E5 scorer twice; E6 kill and continue (three masks); G0; G2 (<= 75 minutes) | l1-k281.R:118; l1-k281.R:132; l1-k281.R:207; l1-k281.R:247; l1-k281.R:100; l1-k281.R:268 |
| 6 | the scorer self-test runs inside L1 (a defect: fix, re-hash, re-declare before any counted job) | l1-k281.R:194 |

## 3. The script lessons (KIT-188 script audit 2), each with its lines

| Lesson | Implemented at |
|---|---|
| 1. every table read as text (`colClasses = "character"`, `na.strings` empty) and converted by name; no value turns into NA | scorer: score-k281.R:88, score-k281.R:93, score-k281.R:104 (logical), score-k281.R:120 (labels); the same in `l1-k281.R` (`readText`, `asNum`); the simulation's `readPart` was already text-only (unchanged) |
| 2. `options(warnPartialMatchDollar = TRUE)`; the warning is an error only around the scripts' own table handling; text matched case-insensitively (`tolower`) | scorer (no package code: global): score-k281.R:27, handler score-k281.R:713, guard-stop = instrument status via `stopRecord` score-k281.R:464, probe score-k281.R:509; simulation (scoped, K281-P05): sim-k281.R:571; `l1-k281.R`: l1-k281.R:28 |
| no `$` inside `suppressWarnings` | none of the three scripts calls `suppressWarnings` (the scorer converts with a regex check, then `as.numeric`) |
| guard stop in the simulation is an instrument status that stops the scorer | the simulation writes `instrument-stop.txt` before aborting (sim-k281.R:577); the scorer refuses any job directory holding it or `errors.txt` (score-k281.R:157) |
| `qmethod == "qvalue"` in every metrics row: COORDINATOR'S RULING S03: counted and REPORTED, not a stop (a BH fallback of `qMethod` is deterministic and the decision reads no q-value) | score-k281.R:558 (`qmethod-fallback.tsv` per cell, mask and arm; the `qmethod` check; `n_qmethod_fallback` and `notes` in `outcome.tsv`); the guard stays out of package code (K281-P05), the way a silent switch could arise |
| stop on rows outside the declared replicate range (RANGE-1) | score-k281.R:127 (`declaredRange`: the G3 table, or `K281_RANGE` for mechanics runs); score-k281.R:206 |
| stop on any input carrying a smoke flag (E16's `smokeflag`) | score-k281.R:159 |
| progress files and provenance before any unit runs | simulation: sim-k281.R:647, sim-k281.R:639; scorer: score-k281.R:522, `progress-score.tsv` (score-k281.R:461); L1: l1-k281.R:71 |
| bound to a tree WITHOUT KIT-280 (BIND-1 in the ticket's form): `threshold_values_in_filter` absent from `workflowReplayRuleMarkers` | sim-k281.R:652 (every simulation process, hence every L1 child and every counted job) |
| positional arguments landing in the wrong field (KIT-275's scorer) | S08: every call of `addCheck`, `stopRecord`, `progress` (scorer), `addCheck`, `run`, `failureNote` (L1), `decisionRow`, `summariseStats` and every `data.frame(...)` in the three scripts' new code names its arguments (a mechanical pass: 16 scorer calls and 20 L1 calls were converted); the two-argument condition raiser `k211Abort(class, message)` and KIT-211's verbatim functions are the declared exceptions; the `fields` case plants a distinct value in every field of every table (section 4) |

## 4. The scorer self-test (E16 manner): 2007 checks (64 cases, 1936 predicted values, 7 request cases)

Constructed tables, no simulated value, no emulation of R's random-number generator. Most tables are degenerate (every replicate of a setting alike) so every
bootstrap and jackknife limit equals its estimate; the predictions are derived in `make-selftest.py` (Python 3 stdlib, `python3 -I`) from the planted counts with the
ticket's operations in their order (IEEE doubles are correctly rounded, so a ratio of totals equals the ratio of the per-unit counts and the boundary cases decide exactly as written).
**S01 (script audit): the `q_*` cases are NOT degenerate.** They use two replicates (`K281_RANGE=401,402`); cells B, C, D are alike, cell A differs between its replicates, so per resample cell A's weights are (2,0), (1,1) or (0,2) and every pooled statistic takes exactly three values, each an exact ratio of integer sums. The 0.975 quantile (type 7, position 975.025 of 1,000) is the MAXIMUM of the three unless fewer than 26 resamples take its pattern (probability below 1e-50), the 0.025 quantile the MINIMUM, the estimate the (1,1) value: predicted exactly without knowing a draw, so a scorer that decided on the estimate, the lower limit, the median or a mis-set alpha fails these cases. The jackknife SE is the delete-one formula by hand. `partial` asserts only "`boot_nonfinite` > 0". `K211_B_BOOT` = 1,000 for all cases except `counted_ok`
(B = 10,000, the counted path), recorded in `checks.tsv` (`bootstrap_B`). The generator asserts the outcome of every named case by hand (`HAND`) before writing a prediction.
Runner: `selftest/run-selftest.sh <selftest dir> <out dir>` (scores each case with `score-k281.R` under its own environment, four at a time); checker: `selftest/check-selftest.py`.

| Case | Covers |
|---|---|
| `cal105` | rho_g 1.05 (planted 105 of 2,000 against 50 of 1,000), guards pass: Calibrated; Phi_S 1.0 (positive control fails): reported, outcome unchanged; also the checks table, replicates, counted FALSE |
| `bnd_below` | deciding limit just below 1.10 (5,499 of 100,000 over 500 of 10,000 = 1.0997999999999999): pass, Calibrated |
| `bnd_at` | the limit at 1.10 as written: 550 of 10,000 over 500 of 10,000 = 1.0999999999999999 in doubles: pass, Calibrated (no tolerance) |
| `bnd_above` | just above: 5,501 of 100,000 over 500 of 10,000 = 1.1002: fail, Not calibrated |
| `rho112` | rho_g 1.12: Not calibrated |
| `split` | the split by the flag: flagged 90 of 1,000 (1.8), unflagged 4,950 of 99,000 (1.0), pooled 1.008: Calibrated |
| `g1a_fail` | G1a fails: FDP_BH(B1) = 26/500 = 0.052 (paired difference 0): Not calibrated |
| `g1a_pass_edge` | G1a at its limit: 25/500 = 0.05 exactly: pass |
| `g1b_fail` | G1b fails: 99/2,000 against 23/500, difference 0.0035: Not calibrated |
| `g1b_pass` | G1b passes: difference 0.002 |
| `g2_fail` | G2 fails: power cost 0.006: Not calibrated |
| `g2_below` | G2 just below: 0.0049: pass |
| `g2_at` | G2 at 0.005 as written: 0.2 - 0.195 = 0.0050000000000000044 in doubles: fail (no tolerance) |
| `g2_above` | G2 above: 0.0051: fail |
| `pc_pass` | Phi_S = 3: the positive control passes; the outcome equals cal105's (it decides nothing) |
| `pc_fail_ok` | Phi_S = 1: the positive control fails and stops nothing: Calibrated |
| `f5_fires_pos` | F5' fires (bias +0.2, degenerate interval above 0): Withheld although every other row says Calibrated |
| `f5_fires_neg` | F5' fires (bias -0.2): Withheld |
| `f5_zero` | F5' does not fire (bias 0): Calibrated |
| `f5_and_fail` | F5' fires and rho_g fails: Withheld (the first row of the outcome table that holds) |
| `noadded` | B1 adds nothing in any setting: the ESTIMATE of rho_g is not finite: Stop (instrument), coordinator ruling S02 (nothing claimed; the default stays 2) |
| `partial` | only replicate 401 of A at vd5 holds added proteins (the pooled estimate is finite): non-finite resamples dropped and counted (> 0), limits over the finite ones (2.0), jackknife limit not finite, labels "exceeds margin" / "not determined" (UNLABELLED) per setting |
| `fields` | distinct planted values in every setting: every column of settings.tsv, designs.tsv, pooled.tsv and labels-settings.tsv lands in its own field (the positional-argument defect); 377 values |
| `miss_rep` | a replicate absent for one cell: Stop |
| `smoke` | the provenance carries smoke = true: Stop |
| `smoke_absent` | the provenance carries no smoke flag: Stop |
| `range_out` | rows outside the declared replicate range: Stop |
| `guardstop` | a job directory holding instrument-stop.txt (a guard stop in the simulation): Stop (instrument) |
| `probe` | the scorer's deliberate `$` prefix match (K281_GUARD_PROBE=1): the guard stops it: Stop (instrument), never an ordinary failure |
| `f1` | an F1 identity violated in one replicate-mask: Stop (falsifier) |
| `f2_e10` | an E10 logFC identity not exact (1e-09): Stop (falsifier F2) |
| `f2_e3` | E3 false for B1's input: Stop (falsifier F2) |
| `incomplete_range` | a range whose complete count is below its expected count: Stop |
| `cell_T` | a cell labelled T: Stop (a label, not a logical) |
| `cell_F` | a cell labelled F: Stop |
| `cell_NA` | a cell labelled NA: Stop (not read as a missing value) |
| `arm_NA` | an arm labelled NA: Stop |
| `mask_mcar` | a random-missingness mask (not a deciding setting): Stop |
| `text_num` | a numeric column holding text: Stop (instrument) |
| `na_num` | a numeric column holding NA: Stop (instrument) |
| `logical_na` | a logical column holding NA: Stop (not coerced) |
| `tf` | metrics columns spelled T and F: accepted as logical, the result equals cal105's (no value lost to a logical guess) |
| `join` | two jobs per cell (r 401-402 and 403-404) joined into one cell of replicates: result equals the single-job result |
| `dup` | a replicate joined twice (two jobs overlap): Stop |
| `f2_det` | a determinism difference: Stop (falsifier F2) |
| `counted_ok` | counted mode (range from the G3 table, B = 10,000, determinism present): Calibrated, counted TRUE, range_source g3-table |
| `counted_nodet` | counted mode without a determinism table: Stop |
| `counted_g3five` | counted mode, G3 table R = 5 against rows R = 4 (R differs from the rule's): Stop |
| `counted_g3none` | counted mode with no G3 table: Stop |
| `q_rho` | S01 (two replicates, cell A differs): the pooled rho_g takes three values (0.9, 1.075, 1.25); estimate 1.075 passes but the bootstrap 0.975 quantile (the maximum, 1.25) decides: Not calibrated; also one SETTING (A vd5: lo 0.6, estimate 1.3, hi 2.0, label "not determined") and the jackknife SE by hand |
| `q_pc` | S01: Phi_S estimate 1.3 (> 1.20) but the lower limit is the minimum 1.0: the positive control fails (reported; outcome Calibrated) |
| `q_g1b` | S01: G1b estimate 0.00226 (pass) but the maximum 0.00898 decides: fail (G1a likewise): Not calibrated |
| `q_g2` | S01: G2 estimate 0.00417 (pass) but the maximum 0.00833 decides: fail: Not calibrated |
| `q_f5_fire` | S01: F5' interval from its minimum (all three values above 0): fires: Withheld |
| `q_f5_fire_neg` | S01: mirror (all three values below 0): fires: Withheld |
| `q_f5_straddle` | S01: F5' estimate positive but the interval straddles 0 (min -0.05, max 0.075): does not fire (an estimate-based scorer would fire) |
| `q_f5_touch` | S07: F5' lower limit exactly 0: does not fire |
| `g1b_at` | S07: G1b at 0.0025 as written (0.05 - 0.0475 = 0.0025000000000000022 in doubles): fail, no tolerance |
| `f5_na` | S07: no flagged one-sided protein, F5' statistic 0/0 not finite: does not fire (the NA path), outcome Calibrated |
| `errors_file` | S07: a job directory holding errors.txt: Stop (instrument) |
| `killed` | S07: a killed job (start provenance only, no provenance-sim.json): Stop (instrument) |
| `g1a_nonfinite` | B1 has no BH call at all: the estimate of guard G1a is not finite: Stop (instrument) |
| `comp_zero` | the comparator B has no p < 0.05 event: rho_g is not finite (division by 0): Stop (instrument) |
| `qmethod_bh` | S03 ruling: one metrics row with q-method BH is REPORTED (count 1 in checks, outcome notes, qmethod-fallback.tsv per cell, mask, arm), not a stop: Calibrated |
| `req_ok.json` (check-request) | exit 0 |
| `req_r276.json` (check-request) | exit 1 |
| `req_cellshort.json` (check-request) | exit 1 |
| `req_overlap.json` (check-request) | exit 1 |
| `req_toolong.json` (check-request) | exit 1 |
| `req_nocite.json` (check-request) | exit 1 |
| `req_range.json` (check-request) | exit 1 |

Interpretation (to be confirmed by the audit): "cells with labels T, F, NA" = a design cell, arm or mask labelled `T`, `F` or `NA` in the row files (the E16 `tf` defect: a column of
such labels read as logical or as missing); "UNLABELLED cells" = settings whose rho_g is not finite, so the setting label is "not determined" (`noadded`, `partial`).

## 5. Job plan, timeouts and the budget

Time basis (KIT-211's `runs.tsv`, request created to run last updated, an upper bound of the job): its eight cell jobs of 600 replicate-masks each took
33.2, 44.4, 41.5, 42.5, 52.9, 42.4, 40.1, 41.8 minutes (sum 338.8, **mean 42.35, maximum 52.9**); the five-arm DA work (S, A, B, B', N plus C0 once) against KIT-281's four DA arms
(S, B, B1, N) is a factor 0.774, taken as 0.8; set-up taken as 10 minutes. One counted job = one cell, up to 200 replicates, three masks (600 replicate-masks at 200 replicates), four workers:

`T(n) = 10 + (42.35 - 10) x 0.8 x n / 200 = 10 + 25.88 x n / 200` minutes, so T(200) = **35.88, taken as 36**; the slowest KIT-211 job projected the same way is (52.9 - 10) x 0.8 + 10 = 44.3.
Timeout rule (KIT-188's section 10, whole minutes): `timeout_minutes = min(75, 5 x ceiling(1.5 x T / 5))`. No job is above 75 minutes; the slowest projected job (44.3) is inside every timeout below.

| Job | Replicates | T expected (min) | `timeout_minutes` | Derivation |
|---|---|---|---|---|
| L1 (one request) | r 9201, 9202 (+ E1' r 201-205, E6 r 9201-9204) | 20 | 55 | ticket section 12 (KIT-211's L1 script took 832 s for a larger lifecycle; KIT-281 adds the self-test and E1'); 55 >= 1.5 x 20 and G2 checks the counted plan |
| counted, 200 replicates of one cell | r 401-600 | 36 | 55 | T(200) = 35.9; 1.5 x 35.9 = 53.8, rounded up to 55; 44.3 (slowest KIT-211 job projected) < 55 |
| counted, 171 replicates (R = 571, third job) | r 801-971 | 32 | 50 | T(171) = 32.1; 1.5 x 32.1 = 48.2 -> 50 |
| counted, 77 replicates (R = 277, second job) | r 601-677 | 20 | 30 | T(77) = 20.0; 30.0 -> 30 |
| scoring (one request) | all cells | 15 | 40 | KIT-211's scorer ran 63 s on 10 cells at B = 10,000; here 17 resampling units, R up to 571; ticket section 12: 15 expected, 40 ceiling |
| reserve: one continuation (55) and one re-scoring (40) | | 0 | 55 and 40 | the instrument-defect rule (VAL-54), ticket section 12 |

**S05:** the counted jobs' `timeout_minutes` are the ones L1's G2 derives from its own seconds (`g2-plan.tsv`: planned script minutes = n x 3 masks x mean seconds / 4 / 60, planned job = script + 10 set-up, timeout = min(75, 5 x ceiling(1.5 x job / 5))), passing when every planned job (200, 171, 77 replicates) is at most 75; the table above is the plan from KIT-211's times, and the requests use g2-plan.tsv's timeouts, with the running ceiling recomputed from them before the first counted request. L1's seconds include the lifecycle checks, so the plan errs long.

Jobs per cell: R <= 200: one (r 401 to 400 + R); R > 200: ceiling(R / 200) jobs of at most 200 replicates, consecutive ranges starting at 401 (R = 277: 401-600 and 601-677; R = 571: 401-600, 601-800, 801-971).
One cell per job (`K211_CELLS`), `K211_DETERMINISM=1` on the first job of cell A. Every job: `out_pattern` `"\\.(tsv|gz|json|txt|md)$"`, `blas` `reference`, `tests` `[]`, `base_commit` `"none"`, `extra_packages` `[]`.

Budget (core-hours = job-minutes x 4 / 60; arithmetic by `python3 -I`):

| Line | Jobs | Expected min | Ceiling min | Expected core-h | Ceiling core-h |
|---|---|---|---|---|---|
| parse-and-self-test send (S09, before L1; `smoke/`) | 1 | 20 | 30 | 1.33 | 2.00 |
| L1 | 1 | 20 | 55 | 1.33 | 3.67 |
| counted at R = 200 (r 401-600) | 4 x 36 | 144 | (4 x 55 = 220) | **9.60** | (14.67) |
| counted at R = 571 (the cap of rule G3) | 12 x 36 | 432 | 12 x 55 = 660 | **28.80** | **44.00** |
| scoring | 1 | 15 | 40 | 1.00 | 2.67 |
| reserve (a continuation 55, a re-scoring 40) | | 0 | 95 | 0 | 6.33 |
| **Gate at R = 200** (with the send) | | **199** | | **13.27** | |
| **Gate ceiling (R = 571)** (with the send) | | | **880** | | **58.67** |
| Part C (only on Calibrated; not built here): tests head and base 2 x 38.2, check 40 | 3 | 116.3 | 240 | 7.75 | 16.00 |
| **Ticket total** | | **315.3** | **1,120** | **21.02** | **74.67** |

Without the S09 send these are the ticket's figures (section 12: 19.69 expected, 72.67 ceiling) and reproduce: 179 min = 11.93 and 850 min = 56.67 for the gate; 116.3 + 179 = 295.3 min = 19.69 and 240 + 850 = 1,090 min = 72.67. **With the send (about 20 job-minutes = 1.33 core-hours expected, `timeout_minutes` 30 = 2.00 at the ceiling) the totals are 315.3 min = 21.02 expected and 1,120 min = 74.67 ceiling: the ticket's ceiling becomes 74.67, not yet accepted by the owner (no request of any kind before it is).** (Part C 7.75 / 16.00; the ticket's 7.76 is a rounding slip.)
At R = 277 (two jobs per cell, 200 + 77): expected 4 x (35.9 + 20.0) = 223.6 min = 14.91 core-hours for the counted jobs; at the ticket's 55 minutes per job 8 x 55 = 440 min = 29.33. Every R between 200 and
571 stays inside the R = 571 ceiling line (12 jobs at 55 minutes is the largest job count).

### The running-ceiling rule for `runs.tsv`

Before ANY request is sent: `used + this request's ceiling + everything still planned (each at its ceiling, the counted jobs at the cap R = 571 until the G3 table fixes R) <= 74.67 core-hours` (72.67 before the S09 send is added), where `used` is the sum of
`core_hours` already in `runs.tsv` (a request's ceiling is its `timeout_minutes` x 4 / 60 per job). The figures are written into the request's declaration and the new `runs.tsv` row is added with the result.

| Moment | used | this request's ceiling | still planned | sum (<= 74.67) |
|---|---|---|---|---|
| before the S09 send | 0.00 | 2.00 | L1 3.67 + counted 44.00 + scoring 2.67 + reserve 6.33 + Part C 16.00 = 72.67 | **74.67** |
| before L1 | <= 2.00 (1.33 expected) | 3.67 | counted 44.00 + scoring 2.67 + reserve 6.33 + Part C 16.00 = 69.00 | <= 74.67 |
| before the first counted job at R = 200 (S09 + L1 used 2.67 expected, at worst 5.67) | <= 5.67 | 3.67 | the other 3 jobs 11.00 + scoring 2.67 + reserve 6.33 + Part C 16.00 = 36.00 | <= 45.34 |
| before the first counted job at R = 571 | <= 5.67 | 3.67 | 11 jobs 40.33 + scoring 2.67 + reserve 6.33 + Part C 16.00 = 65.33 | <= 74.67 |
| before scoring | the used total | 2.67 | reserve 6.33 + Part C 16.00 | used + 25.00 <= 74.67, i.e. used <= 49.67 |

A request whose sum would exceed 74.67 is not sent; the coordinator is told which line to reduce (the reserve and Part C are the planned lines that can be left unspent).

## 6. Declarations (DRAFTS, not sent: the declaration is committed with the final script hashes before L1, and the owner's acceptance of the ceiling precedes it)

### 6.0 The parse-and-self-test send (S09; before L1): drafted in `smoke/` (`smoke-k281.R`, `declaration.md`, `request-draft.json`), NOT sent

- It must show three parses; `SELFTEST: 2007 checks, 0 differ`; the probe case stopped as instrument; `cal105` completing; and one smoke-flagged simulation task (cell A, r = 9301, a newly declared seed) with E2, E3, E7, E10 3/3 each, the marker assertion passed, provenance smoke true, and the scorer refusing that directory. Cost 20 job-minutes = 1.33 core-hours expected, `timeout_minutes` 30 = 2.00 at the ceiling; the ticket's ceiling becomes 74.67.

### 6.0a Statements for every declaration (S10, S02, S03)

- **Boundary rule (S10):** a limit passes iff its computed double is <= the double of the written limit; no tolerance (`bnd_at` passes at 1.0999999999999999, `g2_at` and `g1b_at` fail at 0.0050000000000000044 and 0.0025000000000000022).
- **Non-finite estimate (coordinator's ruling S02):** a non-finite ESTIMATE of rho_g or of a guard is an instrument Stop ("the deciding quantity or a guard could not be formed (not finite)"): nothing is claimed and the default stays 2; it is not Not calibrated. Resample-level non-finite values stay dropped and counted (`boot_nonfinite`). The ticket text is being updated by its writer (section 10, Stop row).
- **q-method (coordinator's ruling S03):** the scorer reports the count of metrics rows whose recorded q-method is not qvalue (per cell, mask and arm; checks table; outcome notes); it is not a stop, because the fallback is deterministic and the decision reads no q-value.

### 6.1 L1 (blind lifecycle, once)

- **Purpose:** prove the instrument before any counted job: the scripts run to completion and failure, the kept sets equal the package filter at 2 and at 1 group, S and B equal KIT-211's committed rows, the
  identities of F1 and E10 hold, the scorer behaves as predicted on constructed tables, and rule G3's one count and table are produced. It decides nothing about calibration.
- **Binding:** `package_commit` = a commit of branch `kit-281` whose `R/`, `src/`, `inst/`, `tests/`, `DESCRIPTION` and `NAMESPACE` tree equals the branch point's a5de1c2e (`git diff --stat a5de1c2e <package_commit> -- R/ src/ inst/ tests/ DESCRIPTION NAMESPACE` empty, recorded in `runs.tsv`); G0:
  `git diff --stat 894b9cc8 <package_commit> -- R/func_prot_qc_filtering_helpers.R R/func_prot_s4_qc_filter_methods.R R/func_general_parameter_helpers.R` empty (`l1/g0.txt`); `base_commit` "none", `tests` [], `blas` "reference", `extra_packages` [],
  one `scripts` entry `Rscript docs/qa/workflow-execution/kit/kit-281/l1-k281.R l1-out`, `out` `l1-out`, `out_pattern` `"\\.(tsv|gz|json|txt|md)$"`, `timeout_minutes` 55. If KIT-280 has merged when the commit is chosen, bind to the last staging commit before its merge.
- **Ceiling:** 55 job-minutes = 3.67 core-hours (expected 20 = 1.33); running-ceiling sum <= 74.67 (section 5).
- **Predictions (all fixed now):** G0 1/1; E1' 120/120; E2, E3, E7, E10 24/24 each (12 replicate-masks per range); E10 maximum absolute deviation 0; E8 2/2; E4 2007 checks, 0 differ; E5 2/2 (two scorer runs identical); E6 1/1 (killed with 1 to 3 of 4 tasks done, continued, identical to uninterrupted, each replicate once with three masks); G2 3/3 (the planned jobs of 200, 171 and 77 replicates, set-up plus script, each <= 75 min; `g2-plan.tsv` gives the timeouts); G3rule 24/24; G3 1/1; `g3-count.tsv` 12 rows, `g3-rule.tsv` 4 rows with one R in 200 to 571.
  The E1' DA flag may read FALSE (KIT-201/KIT-221 move p-values; flagged only, as KIT-211's E1).
- **Output (public):** `checks.tsv`, `timing.tsv`, `hashes.tsv`, `provenance-l1.json`, `progress-l1.tsv`, `g3-count.tsv`, `g3-rule.tsv`; nothing derived from a p-value, a call or a fold change leaves the blind rows except the count (K281-P08).
- **Decision it feeds:** R for the counted run (rule G3, mechanical); any check failing is an instrument defect (VAL-54: fix, re-hash, re-declare, at most three re-sends, no value, margin or seed changes) or, if it is a result, a stop.

### 6.2 The counted run (once; one request per job)

- **Purpose:** the calibration gate: pooled over the twelve settings, rho_g with its bootstrap 0.975 upper limit against 1.10, the guards G1a, G1b, G2, F5', reported per design, per setting, by the flag, Phi and the positive control.
- **Binding:** as 6.1, the same `package_commit` (or another commit whose tree is equal), `requested_by` citing the first 16 hex digits of the SHA-256 of `g3-rule.tsv` and the declaration's hash; each request one `scripts` entry:
  `Rscript docs/qa/workflow-execution/kit/kit-281/sim-k281.R k281-out <first r> <last r>` with `env` `{"K211_CELLS": "<one cell>"}` (and `"K211_DETERMINISM": "1"` on the first job of cell A), `out` `k281-out`, `timeout_minutes` 55 (50 for 171 replicates, 30 for 77).
- **R and seeds:** read from `g3-rule.tsv` (R per cell, `r = 401 to 400 + R`); `python3 -I check-request-k281.py --g3 g3-rule.tsv --count g3-count.tsv <every request file>` exits 0 before the first request is pushed; the simulation refuses a counted range outside the table and the scorer stops on any other set of replicates.
- **Ceiling:** jobs x 55 min x 4 / 60 core-hours (R = 200: 4 jobs, 14.67 core-h at the ceiling, 9.60 expected; R = 571: 12 jobs, 44.00 / 28.80); running-ceiling sum in the declaration.
- **Scoring (one request after every counted job is complete):** `Rscript docs/qa/workflow-execution/kit/kit-281/score-k281.R <rows root> <out>` on the downloaded job directories, no environment variable, `timeout_minutes` 40; `outcome.tsv` must read `counted` TRUE.
- **Decision rule (VAL-365 Q5, fixed):** Calibrated iff rho_g upper limit <= 1.10 and G1a, G1b, G2 pass and F5' does not fire; Withheld if F5' fires; Not calibrated otherwise; Stop on G0, F1, F2 or an instrument status. The positive control decides and stops nothing. One correction budget (VAL-54); an unfavourable valid result is a result and is never re-run.

## 7. Everything unrun (nothing below has been executed by anyone)

- **Every R script:** `sim-k281.R`, `score-k281.R`, `l1-k281.R` have never been parsed or run by R (no R on this machine). Local checks only: bracket and quote balance (Python), E0 (Python), the self-test generator and its prediction asserts (Python), `check-request-k281.py` on its fixtures (Python).
- The scorer self-test (`run-selftest.sh`, 2007 checks), L1 (E1' to G3), the counted jobs, the scoring, Part C, the G0 check on a runner, the replay-marker assertion on the bound tree, `provenance-sim.json` reading, `jsonlite`/`digest` under the guard, the 10,000-resample timing of the scorer, the self-test's `counted_ok` timing.
- No `gh` call, no push of a request, no smoke send, and the S09 parse-and-self-test send is drafted (`smoke/`) but NOT sent: no run of any kind before the coordinator says so.

## 8. Local checks done on this branch (no R)

- Bracket and quote balance of the three R scripts: balanced (a Python scan that ignores comments and handles strings; it does not parse R).
- E0 (`python3 -I docs/qa/workflow-execution/kit/kit-281/e0-k281.py`): KIT-211's three scripts carry the declared SHA-256; the simulation's verbatim region (256 lines) and eleven named functions are byte-identical; the scorer's three verbatim blocks and eleven named functions are identical to KIT-211's (and the blocks to KIT-197's); L1's two verbatim blocks equal KIT-197's.
- `make-selftest.py`: all predictions generated; the outcome of each named case (the HAND table) asserted against the ticket's rules by hand inside the generator; the G3 table asserted at the ticket's examples (96,786 -> 571; 200,000 -> 277; 276,220 -> 200).
- `check-request-k281.py` on its seven fixtures: 1 accepted, 6 refused as predicted.
- G0 locally: `git diff --stat 894b9cc8 HEAD -- <the three filter files>` is empty at the branch point.

## 9. Departures, interpretations and open points for the script audit

1. **E3 (not in the ticket's changed-line list).** Section 6 lists :497-498 as unchanged and section 7 requires an "E3-type identity that B1's DA input is `masked[H_B1, ]`"; with arm A gone, KIT-211's E3 would be a constant TRUE (and the E3 row of the summary a check not run). Implemented at sim-k281.R:494 and the metrics `e3_pass` on arm B1; marked ADAPTED in `script-diff.md`.
1a. **S06 (audit):** the E3 adaptation is accepted; the audit gives the exact ticket text for the changed-line list (E3 at :497-498 and :510, the added refusal, marker assertion, start provenance and `K281_SMOKE`; `e4hand-k281.py`/`selftest-k281.py` read `selftest/`; the scorer's "E1' flag" reads "F2 (E3, E10, determinism); E1' is L1's"). The ticket writer applies it.
1b. **Applied after the script audit (S01 to S10):** see the commit message and the coordinator's hand-back: S01 and S07 cases in `selftest/`, S02 and S03 rulings in the scorer, S04 and S05 in `l1-k281.R`, S08 by a mechanical naming pass, S09 in `smoke/`.
2. **"E1' flag" in the scorer's checks table.** The scorer has no E1' input (E1' compares S and B with KIT-211's committed rows and is an L1 check); the scorer's F2 row (E3, E10, determinism) stands in the checks table. Say if the ticket meant something else.
3. **`e4hand-k281.py` not produced.** KIT-211's e4 tables and hand values belong to a decision structure that no longer exists; the hand values are `selftest/make-selftest.py` (one generator, one checker) and L1's E4 runs them. `declaration.md` is not written here (drafts in section 6); the ticket lists it among the files to be committed before L1.
4. **Setting labels.** The ticket says "labels at one-sided 0.05/12 for the settings" without naming the labels; implemented as KIT-211's `costState` pattern: "within margin" (upper <= 1.10), "exceeds margin" (lower > 1.10), else "not determined" (also when not finite); reported, deciding nothing.
5. **A non-finite deciding quantity: RULED (S02, coordinator).** A non-finite estimate of rho_g or of a guard is an instrument Stop; resample-level non-finite values are dropped and counted. Cases `noadded`, `g1a_nonfinite`, `comp_zero`. (My first reading, state `fail` / Not calibrated, is withdrawn.)
6. **The determinism table is required in a counted scoring** (A3 "determinism identical", F2 "determinism fails"); a mechanics scoring (`K281_RANGE`) does not need it. A counted scoring is one whose range is the G3 table's and whose B is 10,000 (`outcome.tsv` `counted`).
7. **Counted jobs are recognised by their range** (it meets r = 401 to 971): the simulation then needs the G3 table and refuses ranges outside it. L1's ranges (9201 and above, 201-205) are not counted ranges.
8. **The L1 count** is the sum of `n_tested` over `g_added` null rows of arm B1: "tested" means a finite p and logFC (ticket), so the count depends on the finiteness of p-values, never on their values.
9. **Per-design and per-setting tables carry every statistic** (`settings.tsv`, `designs.tsv`, `pooled.tsv`: all statistics with two-sided 95 % bootstrap limits and jackknife SEs), more than the ticket's five `THIN` statistics, which are the ones used for the setting labels.
10. **Residual risks I could not check without R:** every R semantic above (partial-match behaviour of `jsonlite::write_json`/`digest` under the global option in the scorer; `globalCallingHandlers` in L1; `quantile` type 7 returning identical values exactly, which the degenerate self-test relies on; the cost of `weights %*% parts` at R = 571).

