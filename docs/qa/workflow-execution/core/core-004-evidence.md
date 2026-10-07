# CORE-004 evidence: harness binding, population binding and lifecycle proof

Ticket: `tickets/openai/CORE-004.md` revision 1.1. Harness commit `fcabd85c`. Protocol: `docs/qa/workflow-execution/core/protocol.md` (SHA-256 in `bindings.json` → `harness.protocol_sha256`). Receipts: `lifecycle-complete.json`, `lifecycle-fail.json`.

## Execution path (user decision 2026-10-08)

The ticket assumed that `core:*` selectors would run measured members inside the guarded runner. PERF-108's measured P14 members actually ran through the PERF coordinator (`tools/profiling/workflow_perf106_repair.R`): a managed `adapters` application per member, in its own 16 GiB slice, with the coordinator's 50 ms sampler. The `perf108:*` reader-check selectors were only the independent verifiers. The user chose to reuse that path.

- `core-*` coordinator stages launch the members.
- Each `core:*` reader-check selector is that member's verifier, run inside the bulk-rollup guard block: 12 GiB startup refusal, 16 GiB cap, 15.5 GiB early stop, 2 GiB floor, 300 s.
- The coordinator gained a validate-only mode, a per-member 12 GiB / no-active-`multischolar*`-unit start check (`setup_refused`), and a refusal when the harness files differ from `bindings.json`.

## What was built

| Item | Where |
|---|---|
| Cases, bindings reader, control selection, decision rule, independent verifier | `tools/profiling/workflow_core.R` (`coreCases`, `coreBindings`, `coreControl`, `coreDecision`, `coreLifecycleDecision`, `verifyCoreMember`) |
| `core-*` stages | `tools/profiling/workflow_perf106_repair.R` |
| `core:*` selectors: admission, guard block, verifier source snapshot | `run_workflow_checkpoint_proof.sh`, `run_workflow_checkpoint_proof.R` (16 GiB service memory), `workflow_checkpoint_dia_driver.R` (input validation and dispatch) |
| Forced failure right after the first aggregate (`ordered_result_stage`) | `tests/testthat/test-workflow-execution-adapters.R`, fault `after_aggregate` |
| Bindings | `bindings.json`: `builds` (control `perf108` = `jEmE2F`, lifecycle `core_004`), `outcomes`, P14 oracle page hashes (48 screen, 1,200 full) and full summary, `control_receipts`, `harness` (root, protocol hash, plan hash at freeze, 10 harness file hashes) |

## Install (recorded, not counted)

| Attempt | Member | Outcome |
|---|---|---|
| 1 | `core/browser-smoke.I5wqHM` | Failed before building: `install-corrected` ran `git apply` on the empty patch of a clean tree, and git refuses an empty patch. This runner defect was fixed in the harness commit by applying the patch only when it is non-empty. The member is kept. |
| 2 | `core/browser-smoke.zjxStr` | Exit 0. Source `b47d0bf8` (CORE-001 to CORE-003), empty patch. Build `30d4d7fe…09ec4`; exact-sum identity `a6ecc0e4…bd579`. It differs from the control's `77905b65…`, as expected, because the provider `.so` changed in CORE-001. |

The control library `jEmE2F` was re-probed while writing the bindings: its build and exact-sum identity match the CORE-001 record.

## Control receipts (PERF-108)

The closeout's `reused_rollup` stage receipt (`ordered-rollup/stage.rds`, SHA-256 matches) holds the two `full$p14` members on the `jEmE2F` build with every experimental option off: `p14-current-1` and `p14-current-2`. Both are complete, valid and verified, with peaks of 7.10 and 7.05 GiB and elapsed times of 87.5 and 87.0 s. They bind to `jEmE2F` and are recorded under `control_receipts` as a reference. CORE control pairs still run under the CORE protocol and harness.

## Dry runs (not counted)

- The first coordinator validate-only pass was refused, because the harness hash check caught a driver edit made after hashing. This proves that check works.
- After re-freezing and committing, all four validate-only passes exited 0, with about 13.3 GiB available and 0 active managed units:
  - coordinator `core-lifecycle-complete`;
  - coordinator `core-lifecycle-fail`;
  - reader-check `core:lifecycle:complete`;
  - reader-check `core:lifecycle:fail`.
- No stage directory was claimed by any dry run.
- Scope of these dry runs (audit CORE-004-P04): the coordinator's pass checks bindings, protocol and harness hashes, inputs, libraries and host state. The reader-check `--validate-only` checks only selector admission, the input path and service paths; it exits before the guard block and the source snapshot. The verifier units' guard limits (16 GiB, 300 s) and frozen source snapshot are evidenced by the counted runs: `controller-request.nul` of `browser-smoke.GiBsU9` and `browser-smoke.pCuyoq`.

## Counted lifecycle members (ledger: 2 small members, 5 min each; both first attempts, no correction used)

| Member | Decision | Status | Peak (anon / file) | Swap | Elapsed | Max monitor gap | Verifier |
|---|---|---|---|---|---|---|---|
| `core-lifecycle-complete` | `lifecycle_complete_observed` | `complete` | 2.74 GiB (2.66 / 0.05) | 0 | 34.6 s | 53 ms | accepted. All 48 screen logical pages are identical to `p14_oracle$screen_pages`; exact-sum identity `a6ecc0e4…`; `fastduckdb_c`; P14; eight engine threads; writer `fast_r`; 12 figures; journal `complete`; cleanup verified. |
| `core-lifecycle-fail` | `lifecycle_failure_observed` | `worker_failed` | 2.71 GiB | 0 | 27.0 s | 52 ms | accepted. The worker error is "CORE-004 forced failure after the aggregate stage". The journal and completion both record `worker_failed`, and the project head is unchanged (input entity `5b6a36c5…`). Worker and monitor were stopped; nothing was relabelled as success. |

Both stage receipts bind `bindings.json` by SHA-256 (`input_bindings_sha256` equals the committed file's hash). Every library and input path comes from `bindings.json`, with no `latest` lookup. Host at start: 13.31 GiB available, 0 active managed units (stage receipt `host_at_start`). **Deviation (audit CORE-004-P01):** the plan requires host state at the end as well. These two members predate end recording: their end state was observed by hand for units only (0 active `multischolar*` units after each run), and `MemAvailable` at the end was not recorded. From protocol v2 on, the coordinator records `host_at_end` in every operation receipt. Eight `failed` (not running) units from earlier sessions remain listed by systemd; they are not from this ticket and were left untouched.

**Exact P14 check carried from CORE-001:** the installed build with CORE-001's repaired exact-sum provider reproduces the PERF-108 screen oracle exactly. The CORE-001 deferral is closed.

## Regressions (source runner, control `jEmE2F` library first)

| File | Result | Baseline |
|---|---|---|
| `test-workflow-core-harness.R` (new) | 5 tests, 50 passed | — |
| `test-workflow-checkpoint-harness.R` | 780 passed, 2 failed, 1 error | identical to `test-baselines.json` |
| `test-perf108-*` (10 files) | 20 passed, 0 failed; the rest are installed-only skips | unchanged |
| `test-perf106-runner-audit.R` (reads the coordinator; added after the post-hoc audit, missed from the original regression set) | 2 tests, 24 passed | — |
| `test-workflow-core-harness.R` after the audit repairs | 7 tests, 59 passed | — |

**Regression found and fixed:** the first run of `test-workflow-checkpoint-harness.R` gave 759 passed, 2 failed, 2 errors. The new error came from `workflow_checkpoint_dia_driver.R`, where the recorded-QC diagnostic sequence calls `computeSavedWorkflowPrecursorRollup(store, input)` without persisting. CORE-002's new refusal broke it. CORE-002's caller trace searched only `R/` and missed this `tools/profiling` caller. The diagnostic legitimately uses the non-exact arm, so the call now passes `non_exact = TRUE`. The rerun matches the baseline exactly. The other callers outside `R/` are two mocks and a persist call, all reached only on the persist path, so they are unaffected.

## Known limits handed forward

- **Unrecorded metrics:** the harness does not record per-phase timings (CORE-006 and CORE-008) or the execution counters that CORE-013 and CORE-014 need. `coreDecision` returns `primary_metric_unavailable` for them, never a pass. Each of those tickets must add the recording, which is a harness change and requires re-hashing `bindings.json`, before its counted stage.
- **Undefined vehicles:** `core-011-overflow` and the CORE-014 checkpoint vehicle are admitted by the selector regexes, but their members are refused by `coreCases` until CORE-011 and CORE-014 define them.
- **Per-ticket bindings:** each measured ticket binds `harness.tickets.<ticket>` (control and candidate build keys) and records its outcome in `outcomes`.

## Post-hoc audit repairs (2026-10-08)

- **Protocol v2** (user decision on CORE-004-P02): limits hold for each matched pair; primary metrics stay median-to-median. v1 hash `cae16bdcbdfb…` governed only the two lifecycle members. The v2 hash is `harness.protocol_sha256`.
- **Refusals and end state:** the coordinator records `setup_refused` (with reason and host) through `coreSetupRefusal`, also in validate-only (exit 76), and records `host_at_end` after each member's teardown.
- **Generator:** `bindings.json` is generated by the tracked `bindings-generator.R`. Two consecutive runs produced byte-identical output.

## Handoff

CORE-005 to CORE-014: use the protocol and bindings above. The control for CORE-006 is CORE-005's build once that ticket is recorded `accepted`.
