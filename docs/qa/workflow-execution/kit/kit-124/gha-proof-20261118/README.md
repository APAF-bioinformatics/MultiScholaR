# KIT-124 on GitHub-hosted runners: the equivalence proof (2026-11-18)

Workflow on branch `sim-k124` (run-only; not merged). Runs: 44192523498 (proof), 44194745643 (lifecycle and failure check, all arms), 44197079835 (all arms), 44198837284 (run ranges).

Rule: `packaging/verify/compare-results.R` (relative 1e-6, absolute floor 1e-9, discrete values exact), against the laptop's confirm rows for the same seeds.

| Comparison | Result |
|---|---|
| Arms 1, 2, 7, 8; C0, C5 (and C6 in S1); runs 1 to 20; four simulators (1,120 rows) | Within the rule; no discrete difference; null rates, false discovery proportion and power identical; bias and rmse differ by at most 5.6e-8 relative |
| All arms, core cells and one extra-hard cell per simulator | Every arm within the rule except arm 3 (today's design): S1 C6 6 of 20 runs; S4 C0 5 of 20, C5 3 of 20, extra-hard S2 1 of 8. Calls differ by 1 to 4 per run; per-cell means close (S1 C6 mean false discovery proportion 0.0369 laptop, 0.0351 runner; S4 C0 0.4134, 0.4140). Cause not established (BLAS and compiler flags differ) |
| Dependency versions | R, sva, RUVSeq, ruv, RUVIIIC, limma, edgeR, limpa, fastlimpadpc and the package tree equal; 11 secondary packages at other versions (listed in the runs' library tables) |

Declaration of the new-arms run (launched 2026-11-18 after this proof): arms 1, 2, 22, 23, 24 (S4: 1, 2, 22, 23), cells C0 to C7, 1,000 runs, package tree 0b6fcb6d, simulators at 35e5a25a, one job per simulator and cell, 180 minutes per job, one attempt (a failed job may be rerun once as the same command: seeds are per run). Purpose: a measured figure for each selectable combination of control set and k rule. Read: arms 1 and 2 must equal the laptop's rows under the rule before arms 22 to 24 are scored. Arm 3 is not in this run.

Extra-hard cells on the runners wait for the user's decision on arm 3.

Extra-hard cells on the runners (user, 2026-11-18: "A AND B in parallel"): A, arm 3's platform dependence is accepted and stated, and the cells not finished on the laptop were launched (label `hard`: S1 S9 to S14, S2 S1 to S14, S3 S14; 500 runs, all arms, one job per cell, 300 minutes, one attempt). B, in parallel, one bounded test of the cause: the differing cells rerun on the runner with reference BLAS (one run, one correction pass).

## The bounded BLAS test (run 44201424399, 2026-11-18)

The cells where arm 3 differed, rerun on the runner with reference BLAS and LAPACK 3.10.0 (the laptop has reference 3.12.0), runs 1 to 20, all arms, against the laptop's rows under the fixed rule. Outputs in `refblas-compare/`.

| Cell | Arm 3 runs differing, OpenBLAS | Arm 3 runs differing, reference BLAS | Other arms outside the rule |
|---|---|---|---|
| S1 C6 | 6 | 6 (the same runs) | none |
| S4 C0 | 5 | 4 | none |
| S4 C5 | 3 | 6 | none |

In the differing runs the reference-BLAS values match neither the laptop's nor the OpenBLAS runner's (S1 C6 run 13, rmse: laptop 1.31, OpenBLAS runner 16.8, reference-BLAS runner 209). Reading: arm 3's corrected values in these runs are sensitive to rounding-level differences, so each platform gives its own values; BLAS is one input and not the sole cause; compiler flags and the BLAS version difference were not tested (one bounded test, as declared). No other arm moves under either BLAS.
