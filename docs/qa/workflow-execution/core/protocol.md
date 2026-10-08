# CORE measurement protocol

**Version 3** (2026-10-08, user decision after CORE-010's premise check): adds CORE-017 (whole-input checks per batch) with phase timing and sampled phase peaks (§5), and replaces CORE-010's 16 GiB pair with a smaller-budget stage whose budget is fixed by a predeclared rule (§4a). CORE-009 closed without a stage. Version 2 had SHA-256 `f6bbb48a4e9799b2e953d8d1c8b6f87ae25e1ecf31c384131092f6117970baad`; no counted stage ran under it, so no control receipt is invalidated.

**Version 2** (2026-10-08, CORE-004 post-hoc audit, user decision): limits hold for each matched pair, host state is recorded at the start and end of every member, refusals are recorded as `setup_refused`, and the dry-run scope is stated. Version 1 (hash in `core-004-evidence.md`) governed only the two CORE-004 lifecycle members; no measured stage used it.

Frozen by CORE-004 before any counted CORE run. Later tickets cite this file by the SHA-256 recorded in `docs/qa/workflow-execution/core/bindings.json` (`harness.protocol_sha256`). A change to this file is a new protocol version: it needs a user decision and a new hash, and it invalidates control reuse across the change.

## 1. Measurement protocol (verbatim from the series plan v1.1)

- **Control definition.** The control for a measured ticket is the installed library of the most recent ticket whose recorded outcome is `accepted`; a reverted ticket's library is never a control. The first control is the PERF-108 `jEmE2F` library. A control pair is reused from its own receipts when input, library, protocol and harness bindings match; otherwise it is rerun and counted.
- **Vehicles.** `core:<ticket>:(control|candidate):[12]` members. Rollup tickets bind `full$p14` with `p14_oracle`; CORE-013 binds `full$primary` with `filter_oracles$primary` through the sample-filter action; CORE-014 binds `full$primary` adopted into a resident store for Save then Open; lifecycle and overflow members bind `screen$p14`.
- **Per-ticket primary metric.** Frozen here and in `protocol.md`:

| Ticket | Primary metric (must not worsen; expected direction) | Guards |
|---|---|---|
| CORE-006 | elapsed of the staged-result read phase (lower) | peak ≤ 1.00, total elapsed ≤ 1.05 |
| CORE-007 | combined peak (lower) | elapsed ≤ 1.05 |
| CORE-008 | combined peak (lower or equal) and validation-phase elapsed (lower) | total elapsed ≤ 1.05 |
| CORE-009 | closed without change (v3); no stage | — |
| CORE-017 | sampled peak of the checks phase, `phase_peak:checks` (lower) | combined peak ≤ 1.00, elapsed ≤ 1.05 |
| CORE-010 | (v3) completion with verified pages under the predeclared budget B (§4a) | none; elapsed against the control's 16 GiB median is a recorded trade-off for the user |
| CORE-011 | combined peak on the normal path (≤ 1.00); overflow member completes without OOM | elapsed ≤ 1.05; anon and file maxima reported separately |
| CORE-012 | decision only; prototype peak reported against the 0.85 adoption rule | science tolerance as predeclared |
| CORE-013 | `synchronous_wait_seconds` and `sync` fork count (lower) | elapsed ≤ 1.05, peak ≤ 1.00 |
| CORE-014 | project disk bytes after Save then Open (lower) and head publishes per Save (lower) | elapsed ≤ 1.05, peak ≤ 1.00 |

- **Peak rule.** Peak is the combined managed-process kernel peak including file cache, exactly as PERF-108 measured it. Spilled or dirty file pages count. Anonymous and file maxima are always reported alongside.
- **Acceptance.** Both candidate runs complete with valid monitoring (gap ≤ 100 ms), zero swap, output hashes identical to the oracle (CORE-012 excepted), the primary metric not worse than the control median, every guard satisfied, cleanup verified. A valid unfavourable result is recorded as the ticket outcome and the mechanism reverted in a follow-up commit; it is never retried or re-thresholded.
- A unit-test pass is never a performance claim. Builds are recorded with elapsed and peak but are not counted runs.

## 2. Preconditions for every counted member

- `MemAvailable` at least 12 GiB and no active `multischolar*` user unit at the start of each member. The coordinator checks both before every launch (`coreSetupRefusal`). On refusal the stage receipt records `status = "setup_refused"` with the reason and host reading; no member is launched and nothing is counted. The validate-only pass applies the same check and exits 76 on refusal. Host state at start is `host_at_start` in the stage receipt; host state after each member's teardown is `host_at_end` in its operation receipt.
- A validate-only dry run of the exact coordinator invocation has passed on the same library, bindings and harness revision. The coordinator's dry run checks bindings, protocol and harness hashes, inputs, libraries and host state. The reader-check `--validate-only` checks only selector admission, the input path and service paths: it exits before the guard block and the verifier source snapshot, which are evidenced by each counted verifier unit's `controller-request.nul`.
- The harness files listed in `bindings.json` (`harness.files`) hash to their recorded values. The coordinator refuses to start otherwise.
- `OMP_NUM_THREADS=1`, `OPENBLAS_NUM_THREADS=1` (set by the service launcher); DuckDB eight threads.
- Test results are compared with `docs/qa/workflow-execution/core/test-baselines.json`.

## 3. Execution path (user decision 2026-10-08)

Members run through the existing PERF coordinator, `tools/profiling/workflow_perf106_repair.R`, with a `core-*` stage. Each member is a managed `adapters` application in its own slice: 16 GiB hard cap, zero swap, eight threads. The coordinator samples the slice every 50 ms, stops it at a 15.5 GiB peak, a host `MemAvailable` below 2 GiB or the member ceiling, and records peak, anonymous and file maxima, swap, accounting and elapsed time on the shared uptime clock. The independent verifier is the `core:*` reader-check selector, run by `run_workflow_checkpoint_proof.sh` inside its own guarded unit (12 GiB startup refusal, 16 GiB cap, 15.5 GiB early stop, 2 GiB floor, 300 s).

```
# validate only (not counted)
Rscript tools/profiling/workflow_perf106_repair.R <stage> <library> <root> --validate-only
bash tools/profiling/run_workflow_checkpoint_proof.sh --service reader-check <library> <root-parent> <repo> <selector> --validate-only
# counted member(s) of one stage; the coordinator launches each verifier itself
Rscript tools/profiling/workflow_perf106_repair.R <stage> <library> <root>
```

`<root>` is `harness.root` in `bindings.json`, and `<root-parent>` is its parent directory. `<library>` is the stage's candidate library: for lifecycle stages, the `core_004` build.

## 4. Stages, member order and ceilings

| Stage | Members in launch order | Population | Ceiling per member | Verifier selector |
|---|---|---|---|---|
| `core-lifecycle-complete` | `member` (candidate, no fault) | `screen_p14` | 300 s | `core:lifecycle:complete` |
| `core-lifecycle-fail` | `member` (candidate, fault `after_aggregate`: the worker stops right after the first `ordered_result_stage`) | `screen_p14` | 300 s | `core:lifecycle:fail` |
| `core-006` … `core-013`, `core-017` | `control-1`, `candidate-1`, `candidate-2`, `control-2` | `full_p14` (`full_primary` for CORE-013) | 600 s | `core:<ticket>:<role>:<rep>` |
| `core-010` (v3) | `control-1`, `candidate-1`, `candidate-2`, all at budget B (§4a) | `full_p14` | 600 s | `core:010:<role>:<rep>` |
| `core-011-overflow` | defined by CORE-011 | `screen_p14` | 300 s | `core:011:overflow` |
| `core-014` | defined by CORE-014 (checkpoint Save/Open vehicle) | `full_primary` | 600 s | `core:014:<role>:<rep>` |

Each ticket's `control` and `candidate` builds are bound under `harness.tickets.<ticket>` before its stage runs. The control is the most recent build whose ticket outcome is `accepted` (`coreControl`); a control pair whose bindings still match is reused from its receipts rather than rerun.

## 4a. CORE-010 smaller budget (v3)

- **Budget rule, fixed before any CORE-010 run.** B = 0.90 × the median combined peak of CORE-010's control build on `full_p14`, taken from the CORE-017 stage receipts (the control pair if CORE-017 is rejected, the candidate pair if it is accepted), rounded down to a multiple of 0.25 GiB. B, its inputs and the receipt hashes are written to `core/core-010/budget.md` and committed before the first CORE-010 member.
- **Member policy at B.** The application slice `memory_max` is B and the coordinator's early stop is B − 0.5 GiB, mirroring 16 / 15.5 GiB. Swap stays zero, threads eight, monitoring 50 ms. The host start precondition (12 GiB available, no other managed unit) is unchanged.
- **Outcomes.** The control member is expected to stop at the guard; if it completes, B is not binding and the stage is `budget_not_binding` (recorded; no new budget is chosen). Pass (`fits_budget`) requires both candidates to complete with valid monitoring, zero swap, verified cleanup and pages identical to the oracle. A candidate stopped at the guard is a valid unfavourable result (`exceeds_budget`). Elapsed of each candidate is divided by the control build's 16 GiB median from the same CORE-017 receipts and recorded as the trade-off.

## 5. Decision rule as implemented (`tools/profiling/workflow_core.R`)

- **Eligibility of every member, in this order:** verifier accepted, status `complete`, measurement valid (gap ≤ 100 ms and no `monitor_lost_after_result`), lifecycle tests passed, zero swap, peak below 15.5 GiB, elapsed within the member ceiling. The first failing rule is named and the stage is `invalid`.
- **Verifier (P14):** logical 10,000-row page hashes identical to `p14_oracle` (screen: 48 pages; full: 1,200 pages and the frozen summary), 480,000 / 12,000,000 output rows, the build's exact-sum identity, query owner `fastduckdb_c`, profile `P14`, eight engine and writer threads, writer `fast_r`, an admitted worker cap, twelve non-empty figures, export ready, journal `complete`. The request is rejected if its journal is `conflict`, its completion is `guard_payload_invalid` or it carries `monitor_lost_after_result`.
- **Comparison:** candidate median divided by control median for each primary metric must be at most 1 (lower is better). Then every guard must hold **for each matched pair**: candidate-*k* divided by control-*k* (*k* = 1, 2) at most the guard value. A single breaching pair fails the guard; medians are not used for guards. Otherwise the result is `reject`, naming `primary_metric_worse:<metric>` or `guard_breached:<guard>`. Both the median ratios (`ratios`) and the per-pair ratios (`pair_ratios`) are recorded. A metric that is not recorded is `primary_metric_unavailable:<metric>`, never a pass. For CORE-010, an elapsed breach with a memory gain is `trade_off` and goes to the user. CORE-012 is `decision_only`.
- **Metric names:** `peak` = `peak_bytes`; `elapsed` = `elapsed_seconds`; `phase:<name>` = `phase_seconds$<name>`; `phase_peak:<name>` = `phase_peak_bytes$<name>`; `phase_peak:checks` = the larger of `phase_peak:relational_validation` and `phase_peak:working_memory_admission`; `execution:<name>` = the member's `execution$<name>`.
- **Phase measures (v3).** The worker writes `phase-events.tsv` (shared uptime clock; `start`, then `complete` or `failed`) through the existing `multischolar.dia_rollup_phase` hook; the coordinator writes its 50 ms samples to `memory-trace.tsv`. `phase_seconds` sums each phase's intervals; `phase_peak_bytes` is the largest sampled `memory.current` inside any of its intervals. Phase peaks are sampled, not kernel high-water marks, and a phase with no sample inside it has no peak (`primary_metric_unavailable`).

| Ticket | Metrics (all lower is better; medians) | Guard ratios (each matched pair) |
|---|---|---|
| 006 | `phase:staged_result_read_summary_write` | peak ≤ 1.00, elapsed ≤ 1.05 |
| 007 | `peak` | elapsed ≤ 1.05 |
| 008 | `peak`, `phase:relational_validation` | elapsed ≤ 1.05 |
| 009 | `peak` | elapsed ≤ 1.05 |
| 010 | v3: budget stage (§4a), not a ratio comparison | none |
| 017 | `phase_peak:checks` | peak ≤ 1.00, elapsed ≤ 1.05 |
| 011 | `peak` | peak ≤ 1.00, elapsed ≤ 1.05 |
| 012 | decision only | predeclared tolerance (CORE-012) |
| 013 | `execution:synchronous_wait_seconds`, `execution:sync_forks` | peak ≤ 1.00, elapsed ≤ 1.05 |
| 014 | `execution:project_disk_bytes`, `execution:head_publishes` | peak ≤ 1.00, elapsed ≤ 1.05 |

Since v3 the harness records phase seconds and phase peaks (CORE-017). CORE-013's and CORE-014's execution metrics are not recorded by the CORE-004 harness. The ticket that needs one must add its recording before its counted stage (a harness change: the harness files are re-hashed in `bindings.json`); until then its decision is `primary_metric_unavailable`.

- **Lifecycle stages:** `core-lifecycle-complete` is `lifecycle_complete_observed` only if every eligibility rule passes. `core-lifecycle-fail` is `lifecycle_failure_observed` only if the verifier confirms all of the following, and is never `complete`:
  - the recorded status is `worker_failed`, in the receipt, the journal and the completion;
  - the worker error is the forced-failure message;
  - the project head is unchanged and the source head is unchanged;
  - lifecycle tests passed, there was zero swap, and cleanup stopped both the worker and the monitor.

## 6. Install path (recorded, not counted)

```
bash tools/profiling/run_workflow_checkpoint_proof.sh --service install-corrected <previous accepted library> <proof root>/core
```

This archives the current `HEAD`, applies the working-tree diff of `R src configure DESCRIPTION NAMESPACE` (empty on a clean tree) and installs into `<run dir>/library`, alongside the pinned fastduckdb, duckplyr and nanoarrow. Each candidate ticket performs one install. The source commit, patch hash, build identity (`workflowExecutionBuildIdentity`) and exact-sum identity (`diaRollupNativeInfo()$manifest$binary_sha256`) are recorded in `bindings.json` under `builds` before any member uses it.
