# CORE-003 evidence: execution recovery and leak repairs

Ticket: `tickets/openai/CORE-003.md` revision 1.1. Source: `R/utils_workflow_execution.R`. Tests: `tests/testthat/test-workflow-execution.R` (five new cases that run from source; systemd, cgroup and application checks are mocked, while store, adapter, journal and control files are real). No counted runs, no installed build.

## Premise checks (before editing)

| Item | Premise | Verified |
|---|---|---|
| 1 | The recovery conflict aborts before `active_request` is cleared, so every later controller fails the same way. | Yes. Recording `conflict` alone would not be enough: on the next open the old code falls through to `controller_lost` and clears the request silently. Recovery therefore refuses explicitly on `conflict`. |
| 2 | A `cleanup_failed` abort inside `workflowManagedMain`'s `on.exit` skips `application-monitor-done.rds`. | Yes. The application monitor binds only to the supervisor and runs with `RuntimeMaxSec=86500`. The controller binds to the monitor, not the other way round. |
| 3 | The first `monitor_lost` return in `pollWorkflowAction` comes before the result file is read. | Yes for the first site. The second site is reached only when the result file is absent. |
| 4 | An unvalidated guard `reason` breaks teardown. | Yes, in two ways. A NULL reason removes `status` from the receipt and the journal, and the GUI observer's unguarded `if (result$status %in% ...)` then errors. A non-list payload makes `$` error on every poll, so the job never finishes. |
| 5 | The publication handler drops the condition message. | Yes. |

## Changes

1. When recovery finds that the durable head is neither the expected nor the proposed head, it writes journal `status = "conflict"` and `conflict_head`, then refuses with class `c("workflow_execution_conflict", "workflow_execution_stale")`. Existing `stale` handlers still catch it, and the message names the entry point. Later opens refuse on `conflict` without retrying recovery. `acknowledgeWorkflowExecutionConflict(root, request_id, head)` holds the mutation lease and checks three things: `head` must equal the durable head, `active_request` must equal `request_id`, and the journal must be in `conflict`. It then writes `acknowledged_head` to the journal and clears `active_request`. The status stays `conflict` (one new state). The proposed head is never applied.
2. `workflowManagedMain` exits through `closeWorkflowManagedApplication(app)`. It closes every controller, collects the first failure, always writes the GUI monitor-done file, sets `app$closed`, and then re-signals the failure.
3. A monitor loss with `worker-result.rds` present no longer ends the job. The existing worker-success and result checks still decide acceptance. The receipt gains `monitor_lost_after_result = TRUE` (present only when true) and `measurement_valid = FALSE`.
4. `readWorkflowGuardReason` is used at both guard sites. It returns the reason only if it is a single non-NA, non-empty string; an unreadable file, a non-list payload, NULL, NA or a vector gives `guard_payload_invalid`.
5. The publication handler stores `conditionMessage` (truncated to 1,024 characters, as at the existing launch-error sites) in `job$error`. The receipt already copies `job$error`.

## Tests (source runner, control `jEmE2F` library first)

| Case | Before the change | After |
|---|---|---|
| conflict holds controllers until acknowledged (refusal twice, wrong head refused, correct head opens with `active_request` clear) | error: abort classed `workflow_execution_stale` only (no `workflow_execution_conflict`); the test stopped before reading the journal, which stays `committing` by code | 12/12 |
| cleanup failure at exit still releases the application monitor (done file written, job stage `cleanup_failed`, terminal `controller_closed`, then the in-process `workflowMonitorMain` returns at once with no guard stop) | function absent. A mutant without the per-controller `tryCatch` fails the done-file expectation. Before that expectation was guarded, the mutant left the in-process monitor looping until the run had to be stopped, which reproduces the leak. | pass |
| worker result before monitor loss (control: no result, so still `monitor_lost`) | `monitor_lost`, `measurement_valid = TRUE` | `complete`, flag TRUE, `measurement_valid = FALSE`, flag in `completion.rds`. The worker unit's state is mocked as `available`, `Result = "success"`; see the open item below. |
| invalid guard payloads (NULL, vector, NA, non-list; valid `memory_guard` unchanged) | receipt and journal `status` NULL | `guard_payload_invalid` in receipt, journal status and `completion_status` |
| publication failure message | `job$error` NULL | message in `job$error`, the receipt and `completion.rds` |

`test-workflow-execution.R` from source: 15 tests, 46 expectations passed, 0 failed, 10 skipped (the existing installed-only lifecycle cases).

## Regressions (source)

| File | Result | Baseline |
|---|---|---|
| `test-workflow-execution-fixtures.R` | 55 passed, 0 failed | — |
| `test-workflow-execution-adapters.R` | 6 skipped (installed preflight or member only) | unchanged skip conditions |
| `test-workflow-saved-recovery.R` | 58 passed, 4 failed, 2 errors | identical to `test-baselines.json` (two SIGKILL/concurrency tests, `MultiScholaR.so` not found in the callr child) |

The installed-only lifecycle cases were not executed: the ticket allows no counted runs, and there is no installed CORE-003 build. By reading the code, the existing "monitor loss discards worker delivery" case kills the monitor while the worker sleeps, before any result exists, so it still ends `monitor_lost`. The CORE-004 lifecycle proof is the first installed execution of this code.

## Deviations

- The acknowledgement takes a third argument, `head`. The acceptance criterion requires that "acknowledging with a wrong head is refused", which needs the operator's asserted head; the plan's two-argument signature cannot express that.
- Rollback note correction: a previous build does **not** read a `conflict` journal as the old `stale` refusal. Its recovery treats any status other than `committing` or accepted as interrupted: it writes `controller_lost`, clears `active_request` and opens on the current durable head. No entity or head is lost or changed, because the proposed head was never applied, but the conflict marker is overwritten and the request is released without acknowledgement.

## Open item and known limitation

- **Item 3 in production (open, audit CORE-003-P01).** After a monitor loss, the result is accepted only if the worker unit still reports `Result = success`. In production the worker is `BindsTo` the monitor, so systemd stops it when the monitor goes. Whether a remain-after-exit unit stopped that way, or one already garbage-collected, keeps `Result=success` has not been observed. If it does not, the result is still discarded, now as `worker_failed`. CORE-004's lifecycle members never lose the monitor, so they cannot close this. It needs a one-second `systemd-run --user` probe, or an installed monitor-loss-after-result case, either of which needs explicit authorisation.
- **No durable head (audit CORE-003-P05).** If `dia-state-head.rds` has been removed outside the package while the journal names an `expected_head`, the conflict records `conflict_head = NULL` and cannot be acknowledged, because `head` must equal a durable head id. No package path removes the head file. Supporting this case would be new behaviour and needs a user decision.

## Handoff

CORE-004 (harness) and CORE-013: new journal status `conflict` with fields `conflict_head` and `acknowledged_head`; new receipt field `monitor_lost_after_result`; new terminal status `guard_payload_invalid`; error class `workflow_execution_conflict`. CORE-004 lists all of these (updated after audit CORE-003-P04). Readers of journals and receipts were searched in `tools/profiling` and `R/` (`journal.json`, `completion.rds`, receipt and journal name checks). The two found are unaffected: `readDIANativeJournal` reads the `dia_native_export.v1` schema, and the PERF-107 driver reads the GUI action receipt `gui-completion.rds`. (The search during the ticket covered only three field names; this wider search was done in the post-hoc audit.)
