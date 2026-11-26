# KIT-275 smoke send 2: declaration (written before the send)

**Instrument check, blind, nothing judged.** The one second send the coordinator allowed (ruling after AUDIT-KIT-275-scripts-20261126). A miss is recorded and handed back without repair; a third send needs a new ruling.
Same blind shape as send 1: row files under `tempdir()`, no estimand uploaded. Package commit: the tip named in the request (scripts after the repairs S01 to S16 and rulings R-1 to R-4).

## What it must show (the audit's eight items)

1. All scripts parse; E0 (`e0-diff.py`) lists the new marks (K275-assaycols, K275-nogroup, the moved K275-REFUSALS) and `0 unmarked`.
2. The self-test with S16's cases: 34 cases, 2,745 predictions, 0 differ.
3. Unit `E-D3-9101` completes: V1g's PEL2 ends `uncorrected` or `dropped`, PEG2 ends `refused_after_exclusion` (class `ruv_k_above_limit`), E11 passes.
4. `K-D3-9101`'s VD NO is `no_da` with class `filter_removed_group:B`, accepted by `L1_status`.
5. The planted guard journey (`K275_FAULT_GUARD=1`) ends `instrument` (check `GUARD_planted` passes) and the blind scorer reports the guard stop (`instrument_stop` item false, exit 3).
6. E3 passes and lists `correlation_exclusion` at `within_assay`.
7. `timing-<job>.tsv` gives completed-journey means for every arm class at D3 and D6 (tier D at D6 has NO, LM, LA, GA, GK; LG, GSC, PEL2, PEG2 are at D3 only) and every `G2_coverage` check passes.
8. No other failed check.

## The request (one)

| # | Script | Arguments and environment | Limit (min) |
|---|---|---|---|
| 1 | `selftest/run-selftest.R` | `K275_SELFTEST_OUT=k275-smoke2-selftest` | 8 |
| 2 | `sim-k275.R` | `D,K,E,T D3 9101 9101 k275-smoke2-a`; `K275_LIFECYCLE=1 K275_SMOKE=1 K275_FAULT=1 K275_FAULT_GUARD=1 K275_REPEAT=T K275_JOB=smoke2-a` | 20 |
| 3 | `sim-k275.R` | `D D6 9101 9101 k275-smoke2-b` (12-sample timing); `K275_LIFECYCLE=1 K275_SMOKE=1 K275_JOB=smoke2-b` | 12 |

Ceiling: (8 + 20 + 12 + 10 set-up) = 50 job-minutes x 4 / 60 = **3.33 core-hours**. Running-ceiling rule: used 1.67 (send 1) + 3.33 + planned 321.20 = 326.20, at most 350.73: met (`runs.tsv`).
Declared decision: send L1 only after the coordinator's word; this send sets no timeout, ceiling or threshold.
