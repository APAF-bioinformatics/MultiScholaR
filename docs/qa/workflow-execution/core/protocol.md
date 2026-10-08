# CORE measurement protocol

**Version 7.1** (2026-10-08, editorial, audit CORE-S3-P08 and P10; no rule changes): §1 is relabelled as amended, and §4a's cap statement is limited to builds with `workflowAdmissionCap`. Version 7 had SHA-256 `ef4c43b603c004a7ae1d163dd5558385282257a814272012b69f8f9b7698d10d`.

**Version 7** (2026-10-08, forced by the control-definition rule after stage `core-010-control` rejected CORE-005; user decision to apply the rule): CORE-005 is `reverted`, so CORE-010's control build is the PERF-108 build. B is taken from that build's members in `core-010-control` (role `control`) by the unchanged 0.90 × median anonymous peak rule, the budget stage's control member is the PERF-108 build, and its time ratios divide by the PERF-108 medians. The CORE-010 candidate is rebuilt without CORE-005. Version 6 had SHA-256 `f4838cdcb526f6bb1821c766156e532aab3374d8cd4a12363ae447f66239460a`.

**Version 6** (2026-10-08, user decision on CORE-010): the CORE-010 smaller budget is produced by a placeholder process outside the slice that holds memory so only B + the 2 GiB host floor is free when each member starts, because the product's 16 GiB policy is fixed and its worker cap now follows free host memory (§4a). B and phase peaks use anonymous memory, which file cache cannot inflate; phase peaks report `phase_peak_anon_bytes` beside `phase_peak_bytes` (§5). CORE-005 is qualified in stage `core-010-control`. Version 5 had SHA-256 `a64e2cda4d9612d51076a6668e92dc06cfe3034c93578fd62e11848f0b83cbee`; no counted member ran under it.

**Version 5** (2026-10-08, user decision: the measurement must work with the host's normal use, including other heavy programs and any power profile). Version 4's quiet-host refusal is withdrawn: the platform profile and busy CPU fraction are recorded per member (`host_at_start`, `host_at_end`), never a refusal. Time is reported, never judged: on a host in normal use neither wall-clock nor CPU time is robust (stage `core-017` under power-saver and a game took 445–488 CPU s against 284 s quiet), so every comparison records per-pair CPU and wall-clock ratios (`time_pair_ratios`) for the user, and memory alone decides acceptance. CORE-005 is qualified against the PERF-108 build run in the same session, not against PERF-108's quiet-host receipts (§4a). A member's measurement is valid with a monitor sampling gap up to 250 ms, because the peak is the kernel high-water mark however the monitor is scheduled (§5). Phase seconds and phase peaks remain wall-clock and sampled, so they are noisier on a busy host. Version 4 had SHA-256 `fe1dade154e61fee030da5fe194f8580afcd410c6bafc0bfe699d00b846fa300`; no counted member ran under it.

**Version 4** (2026-10-08, user decision after stage `core-017` was invalid under host CPU contention): adds a quiet-host start precondition (§2) and takes CORE-010's budget and CORE-005's qualification from a 16 GiB control pair that CORE-010 runs on a quiet host (§4a). Version 3 had SHA-256 `7ba6f7dca14b6c6330b8bed8b3b89531f9b2f3775adf7ce43ff1dea1273c68aa`; its only counted stage, `core-017`, is recorded as invalid and is not a control for later stages.

**Version 3** (2026-10-08, user decision after CORE-010's premise check): adds CORE-017 (whole-input checks per batch) with phase timing and sampled phase peaks (§5), and replaces CORE-010's 16 GiB pair with a smaller-budget stage whose budget is fixed by a predeclared rule (§4a). CORE-009 closed without a stage. Version 2 had SHA-256 `f6bbb48a4e9799b2e953d8d1c8b6f87ae25e1ecf31c384131092f6117970baad`; no counted stage ran under it, so no control receipt is invalidated.

**Version 2** (2026-10-08, CORE-004 post-hoc audit, user decision): limits hold for each matched pair, host state is recorded at the start and end of every member, refusals are recorded as `setup_refused`, and the dry-run scope is stated. Version 1 (hash in `core-004-evidence.md`) governed only the two CORE-004 lifecycle members; no measured stage used it.

Frozen by CORE-004 before any counted CORE run. Later tickets cite this file by the SHA-256 recorded in `docs/qa/workflow-execution/core/bindings.json` (`harness.protocol_sha256`). A change to this file is a new protocol version: it needs a user decision and a new hash, and it invalidates control reuse across the change.

## 1. Measurement protocol (from the series plan v1.1, as amended through v7; the v1.1 text is in git history at the CORE-004 freeze)

- **Control definition.** The control for a measured ticket is the installed library of the most recent ticket whose recorded outcome is `accepted`; a reverted ticket's library is never a control. The first control is the PERF-108 `jEmE2F` library. A control pair is reused from its own receipts when input, library, protocol and harness bindings match; otherwise it is rerun and counted.
- **Vehicles.** `core:<ticket>:(control|candidate):[12]` members. Rollup tickets bind `full$p14` with `p14_oracle`; CORE-013 binds `full$primary` with `filter_oracles$primary` through the sample-filter action; CORE-014 binds `full$primary` adopted into a resident store for Save then Open; lifecycle and overflow members bind `screen$p14`.
- **Per-ticket primary metric.** Frozen here and in `protocol.md`:

| Ticket | Primary metric (must not worsen; expected direction) | Guards |
|---|---|---|
| CORE-006 | elapsed of the staged-result read phase (lower) | peak ≤ 1.00; time reported (v5) |
| CORE-007 | combined peak (lower) | time reported (v5) |
| CORE-008 | combined peak (lower or equal) and validation-phase elapsed (lower) | time reported (v5) |
| CORE-009 | closed without change (v3); no stage | — |
| CORE-017 | sampled peak of the checks phase, `phase_peak:checks` (lower) | combined peak ≤ 1.00; time reported (v5) |
| CORE-010 | (v3) completion with verified pages under the predeclared budget B (§4a) | none; elapsed against the control's 16 GiB median is a recorded trade-off for the user |
| CORE-011 | combined peak on the normal path (≤ 1.00); overflow member completes without OOM | time reported (v5); anon and file maxima reported separately |
| CORE-012 | decision only; prototype peak reported against the 0.85 adoption rule | science tolerance as predeclared |
| CORE-013 | `synchronous_wait_seconds` and `sync` fork count (lower) | peak ≤ 1.00; time reported (v5) |
| CORE-014 | project disk bytes after Save then Open (lower) and head publishes per Save (lower) | peak ≤ 1.00; time reported (v5) |

- **Peak rule.** Peak is the combined managed-process kernel peak including file cache, exactly as PERF-108 measured it. Spilled or dirty file pages count. Anonymous and file maxima are always reported alongside.
- **Acceptance.** Both candidate runs complete with valid monitoring (gap ≤ 100 ms; 250 ms since v5, see §5), zero swap, output hashes identical to the oracle (CORE-012 excepted), the primary metric not worse than the control median, every guard satisfied, cleanup verified. A valid unfavourable result is recorded as the ticket outcome and the mechanism reverted in a follow-up commit; it is never retried or re-thresholded.
- A unit-test pass is never a performance claim. Builds are recorded with elapsed and peak but are not counted runs.

## 2. Preconditions for every counted member

- `MemAvailable` at least 12 GiB and no active `multischolar*` user unit at the start of each member. (v4 added a performance-profile and ≤ 5% busy-CPU refusal; v5 withdrew it. Both readings are recorded per member.) The coordinator checks both before every launch (`coreSetupRefusal`). On refusal the stage receipt records `status = "setup_refused"` with the reason and host reading; no member is launched and nothing is counted. The validate-only pass applies the same check and exits 76 on refusal. Host state at start is `host_at_start` in the stage receipt; host state after each member's teardown is `host_at_end` in its operation receipt.
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
| `core-010-control` (v6) | `control-1`, `candidate-1`, `candidate-2`, `control-2` | `full_p14` | 600 s | `core:010-control:<role>:<rep>` |
| `core-010` (v6) | `control-1`, `candidate-1`, `candidate-2`, each with the placeholder leaving B + 2 GiB free (§4a) | `full_p14` | 600 s | `core:010:<role>:<rep>` |
| `core-011-overflow` | defined by CORE-011 | `screen_p14` | 300 s | `core:011:overflow` |
| `core-014` | defined by CORE-014 (checkpoint Save/Open vehicle) | `full_primary` | 600 s | `core:014:<role>:<rep>` |

Each ticket's `control` and `candidate` builds are bound under `harness.tickets.<ticket>` before its stage runs. The control is the most recent build whose ticket outcome is `accepted` (`coreControl`); a control pair whose bindings still match is reused from its receipts rather than rerun.

## 4a. CORE-010 control stage and smaller budget (v6)

- **Control stage `core-010-control`.** At the normal 16 GiB slice on `full_p14`: the PERF-108 build (`control`) and `core_005` (`candidate`), alternating `control-1`, `candidate-1`, `candidate-2`, `control-2`. Its §5 decision (primary metric `peak`, per-pair peak ≤ 1.00, time reported) is CORE-005's outcome.
- **Budget rule.** B = 0.90 × the median whole-slice anonymous peak (`sampled_anon_bytes`) of CORE-010's control build in that stage (v7: the PERF-108 build, role `control`, because CORE-005 was reverted), rounded down to a multiple of 0.25 GiB (`coreBudgetBytes`). Both of that build's members must be eligible. B and the stage receipt hash are written to `core/core-010/budget.json` and the bindings, and committed, before the first budget member.
- **Budget stage `core-010`.** Members `control-1` (CORE-010's control build; v7: the PERF-108 build), `candidate-1`, `candidate-2` (the CORE-010 build), on `full_p14`, the normal policy and slice. Before each member the coordinator starts a placeholder unit outside the slice (no swap) holding `MemAvailable − (B + 2 GiB)` and records the free memory it leaves; it is refused (`setup_refused`) if less than B + 2 GiB is free. The placeholder stops after the member's teardown, or on any coordinator error. The product's own guards then apply: in builds with `workflowAdmissionCap` the worker cap at dispatch is at most free memory minus the 2 GiB floor (minus the growth reserve from `9e540d8f`); older builds take the slice cap. In every build the monitor stops the action if host free memory falls below the floor.
- **Outcomes (`coreBudgetDecision`).** `budget_not_binding` if the control completes (recorded; no new budget is chosen). `invalid` if the control ends in any status other than `host_headroom`, `memory_guard`, `worker_failed` or `worker_lost`, or a candidate fails eligibility for another reason. `exceeds_budget` if a candidate is stopped by one of those statuses (a valid unfavourable result). `fits_budget` if both candidates complete, verify and are eligible. The chosen batch counts and each candidate's CPU and wall-clock time divided by the control build's medians from `core-010-control` are reported to the user.

## 5. Decision rule as implemented (`tools/profiling/workflow_core.R`)

- **Eligibility of every member, in this order:** verifier accepted, status `complete`, measurement valid (v5, `coreMeasurementValid`: monitor present without a stop reason, sampling gap ≤ 250 ms, no `monitor_lost_after_result`, operation peak and elapsed at least the worker's, worker and monitor stopped), lifecycle tests passed, zero swap, peak below 15.5 GiB, elapsed within the member ceiling. The first failing rule is named and the stage is `invalid`.
- **Verifier (P14):** logical 10,000-row page hashes identical to `p14_oracle` (screen: 48 pages; full: 1,200 pages and the frozen summary), 480,000 / 12,000,000 output rows, the build's exact-sum identity, query owner `fastduckdb_c`, profile `P14`, eight engine and writer threads, writer `fast_r`, an admitted worker cap, twelve non-empty figures, export ready, journal `complete`. The request is rejected if its journal is `conflict`, its completion is `guard_payload_invalid` or it carries `monitor_lost_after_result`.
- **Comparison:** candidate median divided by control median for each primary metric must be at most 1 (lower is better). Then every guard must hold **for each matched pair**: candidate-*k* divided by control-*k* (*k* = 1, 2) at most the guard value. A single breaching pair fails the guard; medians are not used for guards. Otherwise the result is `reject`, naming `primary_metric_worse:<metric>` or `guard_breached:<guard>`. Both the median ratios (`ratios`) and the per-pair ratios (`pair_ratios`) are recorded. A metric that is not recorded is `primary_metric_unavailable:<metric>`, never a pass. Since v5 there is no time guard and no `trade_off` status; time ratios go to the user. CORE-012 is `decision_only`.
- **Metric names:** `peak` = `peak_bytes`; `elapsed` = `elapsed_seconds` and `cpu` = `accounting$cpu_seconds` of the managed slice, both reported per pair in `time_pair_ratios` and never judged since v5; `phase:<name>` = `phase_seconds$<name>`; `phase_peak:<name>` = `phase_peak_bytes$<name>` (cache-inclusive `memory.current`; since v6 `phase_peak_anon_bytes$<name>` reports anonymous memory beside it); `phase_peak:checks` = the larger of `phase_peak:relational_validation` and `phase_peak:working_memory_admission`; `execution:<name>` = the member's `execution$<name>`.
- **Phase measures (v3).** The worker writes `phase-events.tsv` (shared uptime clock; `start`, then `complete` or `failed`) through the existing `multischolar.dia_rollup_phase` hook; the coordinator writes its 50 ms samples to `memory-trace.tsv`. `phase_seconds` sums each phase's intervals; `phase_peak_bytes` is the largest sampled `memory.current` inside any of its intervals. Phase peaks are sampled, not kernel high-water marks, and a phase with no sample inside it has no peak (`primary_metric_unavailable`).

| Ticket | Metrics (all lower is better; medians) | Guard ratios (each matched pair) |
|---|---|---|
| 006 | `phase:staged_result_read_summary_write` | peak ≤ 1.00; time reported (v5) |
| 007 | `peak` | time reported (v5) |
| 008 | `peak`, `phase:relational_validation` | time reported (v5) |
| 009 | `peak` | time reported (v5) |
| 010 | v3: budget stage (§4a), not a ratio comparison | none |
| 010-control | `peak` (CORE-005's qualification, v6) | peak ≤ 1.00 |
| 017 | `phase_peak:checks` | peak ≤ 1.00; time reported (v5) |
| 011 | `peak` | peak ≤ 1.00; time reported (v5) |
| 012 | decision only | predeclared tolerance (CORE-012) |
| 013 | `execution:synchronous_wait_seconds`, `execution:sync_forks` | peak ≤ 1.00; time reported (v5) |
| 014 | `execution:project_disk_bytes`, `execution:head_publishes` | peak ≤ 1.00; time reported (v5) |

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
