# KIT-275 smoke send 3: declaration (written before the send)

**Instrument check, blind, nothing judged.** The third send, allowed by the coordinator's ruling on send 2 (item 5 missed). A miss is recorded and handed back without repair; a fourth send needs a new ruling.

## The cause found by reading (not observed)
`sim-k275.R` (send 2) raised the guard's `stop()` from a GLOBAL calling handler (`globalCallingHandlers(warning = ...)`, header). R runs a global calling handler with every local handler disabled, so that `stop()` cannot reach any `tryCatch`: the planted `data.frame(abcdef = 1)$abc` in `journeyL` (the old :662-666) ended the child at top level, exit 1, before a row, a check, the result file or a trace was written. A guard stop raised from `quiet()`'s local handler was mapped correctly; one outside `quiet()` was not. The same would have happened to any real partial match outside `quiet()`. Reading is consistent with the 13.9 s, no row, no trace, exit 1 observed; it was not shown by a run.
Repair (only that path): `guardWarning`, one function, is now a LOCAL calling handler around each journey (inside `journeyL`'s `tryCatch`) and around each unit (the child's `tryCatch`), where its stop reaches the enclosing handler; the global handler stays as a last net together with an `options(error=)` writer that leaves a trace before the non-zero exit. A second planted variant is raised at unit level (`K275_FAULT_GUARD_UNIT`, in `runUnit` before any journey). The scorer stops on an instrument row (a unit row too), on a unit row of status crashed, and on a false `unit_completed` check.

## What it must show
1. All scripts parse; E0 reports 0 unmarked hunks.
2. The self-test: 37 cases, 2,751 predictions, 0 differ (the three new cases: instrument unit row, crashed unit row, false `unit_completed`).
3. Planted journey-level guard (`K275_FAULT_GUARD=1`, unit `D-D3-9101`): its journey row has status `instrument` and class `partial_match_guard: ...`, the check `GUARD_planted` passes, a trace file `progress/unit-trace-D-D3-9101-journey.txt` exists, the unit continues; the blind scorer writes `instrument_stop` false and exits 3.
4. Planted unit-level guard (`K275_FAULT_GUARD_UNIT=T`, unit `T-D3-9101`): the unit ends with a unit row `instrument` (class contains `partial_match_guard`), a trace file `progress/unit-trace-T-D3-9101.txt`, a failed `unit_completed` check, and exit non-zero only after those are written; the blind scorer stops on it.
5. The control unit (`S-D3-9101`, ordinary, tier S) completes with every check passing.
6. No other failed check.

## The request (one)
| # | Script | Arguments and environment | Limit (min) |
|---|---|---|---|
| 1 | `selftest/run-selftest.R` | `K275_SELFTEST_OUT=k275-smoke3-selftest` | 6 |
| 2 | `sim-k275.R` | `D,T,S D3 9101 9101 k275-smoke3-a`; `K275_LIFECYCLE=1 K275_SMOKE=1 K275_FAULT_GUARD=1 K275_FAULT_GUARD_UNIT=T K275_JOB=smoke3-a` | 12 |

Ceiling 6 + 12 + 10 set-up = 28 job-minutes x 4 / 60 = **1.87 core-hours** (expected about 18 minutes = 1.2). Running ceiling: used 1.67 + 1.53 + this 1.87 + planned 339.27 (the re-plan accepted by the coordinator) = 344.34, at most 350.73: met (`runs.tsv`). Decision: send L1 only after the coordinator's word.
