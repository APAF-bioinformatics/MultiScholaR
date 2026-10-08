# CORE-017 outcome: stage invalid; change reverted (user, 2026-10-08)

Ticket: `tickets/openai/CORE-017.md` revision 1.0, protocol version 3. Product change `7152942b`, reverted in `fea29bb6`. The phase-timing harness (`c325457c`) is kept.

## Lifecycle re-proof (protocol v3, harness `c325457c`, build `core_004`)

CORE-004's lifecycle stage folders were renamed to `runtime/core-lifecycle-complete.core004` and `runtime/core-lifecycle-fail.core004` (contents unchanged), because the coordinator refuses an existing stage folder.

- `core-lifecycle-complete`: `lifecycle_complete_observed`. Peak 2.70 GiB, 65.9 s; all seven P14 phases `complete`, each with a sampled peak except `query_engine_close` (0.02 s, no sample inside).
- `core-lifecycle-fail`: `lifecycle_failure_observed`. Status `worker_failed`; `relational_validation` and `working_memory_admission` `complete`, `relational_aggregate_order_summary_write` and `ordered_result_stage` `failed`.

## Stage `core-017` (`full_p14`, control `core_005`, candidate `core_017`)

Decision: **`invalid`**, failing rule `control-2:measurement_valid`. Every member completed, verified against all 1,200 oracle pages, ran with zero swap and stopped cleanly. Control-2's application monitor had a 131 ms maximum observation interval, above the 100 ms limit.

| Member | Build | Peak GiB | Elapsed s | Validation peak GiB | Admission peak GiB | Aggregate phase peak GiB | Validation / admission s | Max monitor gap s | Measurement valid |
|---|---|---|---|---|---|---|---|---|---|
| control-1 | 693afda9 | 7.50 | 144.5 | 4.45 | 4.21 | 7.47 | 2.6 / 1.4 | 0.092 | TRUE |
| candidate-1 | 7caf73f9 | 7.33 | 161.6 | 4.24 | 4.13 | 7.29 | 3.0 / 2.5 | 0.078 | TRUE |
| candidate-2 | 7caf73f9 | 7.65 | 159.3 | 4.27 | 4.18 | 7.58 | 3.1 / 2.6 | 0.085 | TRUE |
| control-2 | 693afda9 | 7.88 | 147.5 | 5.66 | 5.46 | 7.78 | 2.4 / 1.4 | 0.131 | FALSE |

Receipts: `runtime/core-017/stage.rds` (SHA-256 `cca41a96…bfdf7a`), `decision.rds` (`91c6b12f…a18f68`). Phase peaks are sampled `memory.current` maxima (protocol v3 §5).

## Cause of the invalid member

The host was not quiet. At the end of the stage the ACPI platform profile was `low-power` (power-saver), some cores ran at 1.4 GHz, a `CUO_MAIN_THREAD` process used about 110% CPU, Java and Discord were running, and the 1-minute load average was about 4.7. The protocol checked memory headroom only. All timings are therefore inflated, about 1.7 times PERF-108's 88 s for the same population, and the lifecycle member took 65.9 s against 34.6 s in CORE-004 on the same build.

## Reading (not a protocol decision)

- **Memory, the primary question:** the checks phase on the valid members fell from 4.45 to 4.24 / 4.27 GiB, about 5%. The in-process premise projected about 1 GB. Splitting the checks shrinks only their grouping state; most of the roughly 4.2 GiB process memory in that phase is elsewhere, so it is not the floor this ticket targeted.
- **Time:** the pair-1 elapsed ratio is 1.12, but under this contention it cannot be attributed in either direction.

## Decision (user, 2026-10-08)

Revert the product change and keep the phase-timing harness. No re-run: the memory effect is too small to justify the change whatever a re-run's validity. Protocol version 4 adds a quiet-host start precondition (performance profile, at most 5% busy CPU over 2 s). CORE-005's qualification moves to a quiet-host 16 GiB control pair run by CORE-010, which also sets CORE-010's smaller budget. `core_017` is recorded as `reverted`.
