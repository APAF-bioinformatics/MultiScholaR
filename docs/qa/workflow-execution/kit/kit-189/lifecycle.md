# KIT-189 lifecycle L1 and gates G1, G2

Ticket revision 1.2. Package code: `v0.9-staging` at 3ff60190448150fa13e06d624378cad1481d3820 (the branch's base); nothing under `R/`, `inst/`, `src/`, `tests/`, `man/`, `DESCRIPTION` or `NAMESPACE` differs on the branch.

## Runs

| Run | Commit | What | Result |
|---|---|---|---|
| 45275836889 (`k189-l1-426cebb0`) | 426cebb0 | `sim-k189.R l1 all 9001 9005` (5 lists in each of the 24 cells, 4 workers), then `e7-k189.R`, then `score-k189.R` blind on the L1 rows | 120 of 120 replicates written; 240 route calls, all `ok` (0 `empty`, 0 `crashed`); E0, E1, E2, E3, E4, E6 pass (`l1/checks-l1.tsv`); E7's driver stopped at a parse error in its report line before it ran anything; blind scoring completed on 1,680 rows (24 files), 108 classifications, outputs' SHA-256 in `l1/score-blind.tsv` (not read) |
| 45281701515 (`k189-l1e7-9f78d703`) | 9f78d703 | E7 instrument run 1 of 3 (VAL-54): one closing parenthesis moved in `e7-k189.R`'s report line; nothing declared changed; the limit read from L1's committed `l1/lifecycle-timing.tsv` | E7 pass (`l1/checks-e7.tsv`) |

Blindness: the L1 rows and the scoring tables were written under `k189-blind/` on the runner and never uploaded; the uploaded files are completion, timing, mechanics counts, check results and hashes (`l1/`). E0 reproduces a published count (KIT-173), as the ticket allows.

## Checks (`checks.tsv`)

| Id | Result | Observed |
|---|---|---|
| E0 | pass | KIT-173 null draw 5 under `"none"`: 156 tested, 120 with q < 0.05, 0 with p.adjust < 0.05; pi0_hat 0.0716, from stage 1 of enrichit's chain (the default smoother); no earlier stage failed |
| E1 | pass | the three archives' SHA-256 equal 2.1's; the five file hashes equal 2.1's; `deparse(formals)` and `deparse(body)` identical for `enrichit:::ora_gson`, `enrichit:::calculate_qvalue`, `qvalue::qvalue`, `qvalue::pi0est`, `qvalue::lfdr`, `clusterProfiler::enricher` (the detail cell is cut at 300 characters before the last name; the pass is the conjunction of all six) |
| E2 | pass | 240 calls: margins (n, M, N) as the result's ratios; p.adjust = BH over the kept terms; q / p.adjust one constant per call within 1e-12; the repeated chain reproduces every q-value exactly; Count > 0 on every row; p = phyper(Count − 1, M, N − M, n) within 1e-12; BHf calls within BH calls; BH calls within q-value calls |
| pairing | pass | both propagation values received the same list in all 120 replicates |
| E3 | pass | `@result` of the route equals a direct `clusterProfiler::enricher` call apart from the order inside `geneID` cells, for B1, B2, B3 under both values (160, 738, 101, 622, 323, 1,277 terms) |
| E4 | pass | BHf (0.01, 0.05, 2/3) and Qf (pi0 = 1, stage 1) equal the hand values |
| E6 | pass | planted case ER of T 2.5 and of a disjoint term 0.625; null lists ER = 1 exactly on every tested term of the 90 L1 null calls with a term |
| E7 | pass | killed by `timeout` at 178 s (exit 124) with 20 of 40 parts complete and no `.tmp` left; the continuation ran the other 20 from the driver's list; the union holds each replicate once (560 rows) and equals the uninterrupted run's rows without `part` and `seconds` |

**G1 (the mechanism): passed** (E0, E1, E2).

## Time (`l1/lifecycle-timing.tsv`, `l1/lifecycle-cells.tsv`)

Observed: R start-up and package load 3.1 s; background set-up B1 24.0 s, B2 11.4 s, B3 30.2 s; 120 replicates in 1,027 s on 4 workers. Median wall per replicate (two route calls, the analysis, the part written), each worker running beside three others: B1 29 to 31 s, B2 17 to 19 s, B3 50 to 55 s; per call under `"none"` 23.5 s (B1), 15 s (B2), 36 s (B3), with parent terms 6.2 s, 3.1 s, 16 s. The ticket's basis (KIT-173, sequential: 11.5 and 3.2 s per call on B1) and assumed 4-worker efficiency of 0.8 gave 14.7 s per replicate / 3.2 = 4.6 s of job wall per replicate on B1; measured, it is 30.2 / 4 = 7.6 s (the sum of the workers' walls over the job's wall is 3.9, but each call takes about twice its sequential time when four run at once). B3 is 1.7 times B1, not 1.5.

Inferred: per job, wall = replicates × median replicate wall / 4 + 1.5 minutes per script.

## G2 (time): the re-split exceeds the accepted ceiling — stop, the owner decides

Re-split by G2's rule (planned wall at most 170 minutes per job; timeout 1.5 times the plan rounded up to 5, at most 260; 10 minutes of set-up per job; R unchanged):

| Row | Jobs (replicates each) | Planned wall per job (min) | Timeout per job (min) | Expected core-hours | Ceiling core-hours |
|---|---|---|---|---|---|
| S null B1 | 6 (each size split 1,000 + 1,000) | 127.3 | 195 | 36.0 | 54.7 |
| S null B2 | 3 (one size each, 2,000) | 154.0 | 235 | 22.5 | 33.0 |
| S null B3 | 9 (each size split 667 + 667 + 666) | 147.1 | 225 | 94.3 | 141.0 |
| S planted B1 | 2 (cells 10, 11; cells 12, 13, 14) | 103.6, 105.2 | 160, 160 | 15.3 | 22.7 |
| S planted B2 | 1 (cells 15 to 19) | 129.6 | 210 | 9.3 | 14.7 |
| S planted B3 | 3 (cells 20 and 23; 21 and 24; 22) | 134.0, 134.0, 88.8 | 205, 205, 135 | 25.0 | 37.7 |
| **S** | **24 jobs** | | | **232.3** | **347.7** |

Totals by the ticket's rows: expected = L1 and its E7 run as used (3.77) + S 232.3 + D 1.33 + SC 1.33 = **238.7 core-hours** (accepted, D3: 151.22); ceiling = L1 6.67 + S 347.7 + D 2.67 + SC 2.67 + S′ (25 % of S's ceiling) 86.9 = **446.6 core-hours** (accepted, D3: 281.17).

The ceiling exceeds the accepted one, so by G2 the work stops here and the owner decides. Nothing has been declared or run beyond L1; R is not reduced (VAL-24). The question for the owner: accept about 239 core-hours expected and 447 at the ceiling at R as declared (D3 revised), or decide otherwise. A different job shape is not expected to lower the cost (inferred: the measured cost is per-call slowdown when four calls run at once, so fewer workers per job spread the same core-hours over more jobs; not measured). Any change of R or of cells is the owner's (VAL-24); none is proposed.
