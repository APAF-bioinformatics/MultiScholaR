# CORE-010 outcome: accepted by user decision (2026-10-08), not by a protocol pass

Ticket: `tickets/openai/CORE-010.md`, protocol versions 6 and 7 (`core/protocol.md` §4a). Product commits:
- `0d5ea8ee`: CORE-017's per-batch whole-input checks re-applied;
- `77e8ac35`: worker cap bounded by free host memory, and adaptive P14 batches.

CORE-005 was reverted in `887cde28`. Candidate build `core_010` = `browser-smoke.vtJFQY` (source `9d9a8daa`, build `0f26956f…`).

## Lifecycle re-proof (harness `0085b897`)

CORE-017's lifecycle folders were renamed `*.core017`.
- `core-lifecycle-complete`: `lifecycle_complete_observed` (decision `a1128db7…`).
- `core-lifecycle-fail`: `lifecycle_failure_observed` (decision `8da986c2…`).

Both record anonymous phase peaks; the completed member also records the admission record.

## Stage `core-010-control` (16 GiB, PERF-108 build against `core_005`)

All four members were valid and verified, with monitor gaps of 73–86 ms. The decision (receipt `e3f3c09b…`) was **reject**, failing rule `guard_breached:peak`:
- per-pair peak ratios were 1.005 (7.98 against 7.94 GiB) and 0.961 (7.88 against 8.20 GiB);
- the median ratio was 0.983.

By the frozen rule CORE-005 is `reverted` (user 2026-10-08). Its change was memory-neutral by construction, so the breach is likely noise. The rule was applied as frozen.

B was fixed from the PERF-108 members' anonymous peaks (7.10 / 7.33 GiB): 0.90 × median, rounded down to 0.25 GiB, giving **6.0 GiB** (`budget.json`, committed `9d9a8daa` before any budget member).

## Stage `core-010` (placeholder leaves B + 2 GiB = 8 GiB free at each start)

Decision: **`budget_not_binding`** (receipts: stage `16a89a70…`, decision `6064cd58…`). All three members completed, verified against the oracle, and had zero swap.

| Member | Build | Free at start GiB | Batches | Worker cap GiB | Peak GiB | Anonymous peak GiB | Wall s | CPU s | Monitor gap s |
|---|---|---|---|---|---|---|---|---|---|
| control-1 | PERF-108 (0c38dce1) | 9.13 | 3 (fixed) | 13.38 (slice) | 7.56 | 6.51 | 95.4 | 322.5 | 0.104 |
| candidate-1 | core_010 (0f26956f) | 8.97 | 5 | 5.32 (host) | 6.02 | 5.15 | 94.3 | 302.8 | 0.080 |
| candidate-2 | core_010 (0f26956f) | 8.79 | 6 | 4.98 (host) | 5.64 | 4.77 | 93.6 | 307.3 | 0.076 |

**Why the control fit:**
- Free memory at start drifted 0.8–1.1 GiB above the 8 GiB target after the placeholder was sized. This is a harness limitation: the placeholder is sized once, and the host's other use changes afterwards.
- The control's anonymous peak this time (6.51 GiB) was below its peaks in `core-010-control` (7.10 / 7.33 GiB), so it stayed above the 2 GiB host floor.

The protocol does not allow choosing a new budget after this result.

## Decision (user, 2026-10-08)

**Accepted on this evidence.** This is not a protocol pass.
- Under reduced free memory the candidate chose 5 and 6 batches by itself.
- Its whole-run peak was 20–25% lower (7.56 → 6.02 / 5.64 GiB) and its anonymous peak 21–27% lower.
- Output was identical, with no measured time cost: wall time 94 against 95 s, CPU time about 5% lower.
- With ample free memory the cap is the slice cap and P14 keeps three batches, so normal runs are unchanged.

This build (`core_010`, including CORE-017's per-batch checks) is the next control.

## Open items

- Before the placeholder is reused, size it against a target and re-check free memory immediately before the controller launch, recording the drift.
- Recalibrate `diaRollupRowBytes` against these receipts. The candidates' estimates admitted 5–6 batches with real peaks well inside their caps (5.15 / 4.77 GiB anonymous against 5.32 / 4.98 GiB caps, whole slice), so the margin is real but thin.
- `multischolar.dia_rollup_batches` is removed in CORE-015.
- A noise tolerance for per-pair peak guards should be decided before the next stage that uses one.

## Corrections (2026-10-08, post-hoc audit `CORE-SERIES-POSTHOC-AUDIT-20261008c`)

- **P01, budget not reached.** Free memory at member start was 0.79–1.13 GiB above the declared B + floor (8.00 GiB), and protocol v6 had no upper-bound check. The declared budget was therefore never tested, and the control's completion is uninterpretable; it is not evidence that the old build fits 6.0 GiB. The candidates' evidence stands, since each adapted to the cap it observed. Before the budget stage is reused, the placeholder must re-check free memory just before launch and re-size or refuse beyond a stated drift. That fix, and recording batches and time ratios on every decision path (P09), are deferred to that reuse.
- **P03, margin.** From `9e540d8f` the host term of the worker cap also leaves the 1 GiB growth reserve free above the host floor. The measured runs predate this; with it, the same free memory gives a smaller cap and more batches.
- **P05, calibration.** `diaRollupRowBytes` comes from a two-point fit using the CORE-017 traces and PERF-108's six-batch run (another build, on a quiet host), not from the CORE-017 traces alone. The docstring is corrected.
- **P07, time.** "No measured time cost" means: in this stage, against one control member, under protocol v5, where time is reported and not judged. PERF-108 measured +46% for six batches on a quiet host.
