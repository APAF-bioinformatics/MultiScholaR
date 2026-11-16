# KIT-131 outcome: the managed worker without systemd (record of 2026-11-17)

Status: **Linux part done; Windows runs not started; step 8 waits.** Branch `kit-131` from `performance-v2` at `a98c8a6b`,
not merged. Design: `w0-design.md` in this folder. Decisions D1 to D4 as in the ticket's revision 1.2; the user's answers
of 2026-11-17 (through the coordinator) are recorded under "Decisions".

## What was built

One protocol, two launchers. The systemd code path is the code it was; the plain launcher sits beside it, chosen once at
launch and carried in the frozen configuration (`config$launcher`: absent for systemd, `"process"` for the plain one).

| Piece | Where (`R/utils_workflow_execution.R` unless named) |
|---|---|
| The seam: `workflowLaunchService`, `workflowUnitState`, `workflowUnitStop`, `workflowServiceEmpty` branch on the launcher; `workflowApplicationMemory` for the two cgroup reads | step 1 (systemd side, a pure move), step 2 (process side) |
| The plain launcher: `workflowProcessStart`, `workflowProcessState`, `workflowProcessStop`, `workflowProcessEmpty`, `workflowProcessForget`, `workflowLaunchProcess` | step 2 |
| The choice: `workflowLauncherChoice`. The option `multischolar.workflow_launcher` (`"systemd"` or `"process"`) wins; without it Linux with a `systemctl` binary gets systemd, anything else the plain launcher | step 2 |
| Outside Linux: host memory (`ps::ps_system_memory()`) and the controller's identity with its start time (`ps::ps_create_time()`) | step 2 |
| The application under the plain launcher: `launchWorkflowProcessApplication` runs the GUI application in the calling process; `checkWorkflowApplication` checks the class, the owner and the launcher; no peak window, accounting or monitor waits | step 3 |
| `processx (>= 3.9.0)` and `ps` in `DESCRIPTION` Imports (D1); `NAMESPACE` unchanged (both called with `::`) | step 2 |
| The state owner's paths (D4): `winslash = "/"` at the 37 `normalizePath()` calls of `R/utils_workflow_*.R` that lacked it | step 6, Linux part |
| Lifecycle tests: `tests/testthat/test-workflow-execution-process.R` | steps 2, 3, 6 |
| CI: `.github/workflows/kit-131-launcher.yml` (ubuntu-latest and windows-latest) | steps 0 and 5, committed, **not yet run** |

No exported function changed (D2). `tools/profiling/*` and `tests/testthat/helper-workflow-checkpoint-fixtures.R` are
untouched (`git diff performance-v2...kit-131 -- tools/profiling tests/testthat/helper-workflow-checkpoint-fixtures.R`
is empty).

## Commits

| Commit | Step |
|---|---|
| `bc535cf5` | 1: the seam, systemd side only |
| `b7992a5c` | 2: the process side, the choice, D1; layer 1 tests written first |
| `2305a9d5` | 3 and 4: the application under the plain launcher; layer 2 tests written first |
| `57c8fc04` | 6, Linux part: D4 path repair; its test fails without the repair and passes with it |
| `24c346bf` | 0 and 5: the CI workflow (not run; see "Open") |

## Evidence (Linux, laptop, R 4.6.1, from source, one test file per process)

The lifecycle file forces the plain launcher through the configuration and the option. Final tree: **183 passed,
0 failed, 0 skipped.**

| # | Property | Result |
|---|---|---|
| L1.1 | exit status (`success`, `exit-code`), a command that cannot start is refused, one name is one unit | pass |
| L1.2 | stop ends the unit and a grandchild that left its process group; the unit is then empty | pass |
| L1.3 | a unit past its deadline is stopped by the state read and reported as `timeout` (2 s deadline, ended in under 6 s) | pass |
| L1.4 | a hard-killed controller takes its unit with it | pass |
| L1.5 | the contract's environment, the scratch directory as working directory, the log; the launcher's own command line asserted; the controller's random stream and environment unchanged by a start | pass |
| L1.6 | a unit never started: unavailable, empty, silent to stop | pass |
| L1.7 | the choice of launcher | pass |
| D4 | the containment checks where `normalizePath()` answers with backslashes (emulated): 8 of 8 with the repair; without it the store refuses its own directory (`saved-entity directory escapes project`) | pass |
| L2.1 | a GUI action on a plain worker: one commit, the result computed in another process and read back, the receipt without measurements, the GUI action closes without waiting for a monitor | pass |
| L2.2 | `fail` and `forbidden_commit` end as `worker_failed` with the worker's message; the head stays | pass |
| L2.3 | cancellation: the worker goes, head, history and action lock unchanged, scratch removed | pass |
| L2.4 | the per-step deadline through the protocol: `worker_failed`, `Result = timeout` | pass |
| L2.5 | closing the application stops a running worker and releases the project | pass |
| L2.6 | after a controller is killed its worker goes; a new controller finds the head unchanged, records `controller_lost` and runs an action | pass |
| L2.7 | the admitted cap is finite, positive and in the request | pass (inside L2.1) |
| L2.8 | the plain launch is the GUI application in the calling process; suites, a supervisor unit and a journey record are refused; without the option a system with `systemctl` takes the systemd launch | pass |

From source a child cannot load an installed build, so the test replaces each worker's command line with a script that
loads the source tree (after asserting the launcher's own command line). The installed-package mode of the same file is
what the CI workflow runs; it has not run yet.

### Linux with systemd unchanged

| File | Passed | Failed | Skipped |
|---|---|---|---|
| test-workflow-execution.R | 114 | 0 | 10 (same as before any edit) |
| test-perf-003.R | 137 | 0 | 1 |
| test-perf107-audit.R | 191 | 0 | 1 |
| test-perf-107.R | 21 | 0 | 0 |
| test-perf-102-audit-repairs.R | 74 | 0 | 0 |
| test-workflow-execution-adapters.R | 0 | 0 | 6 (installed lifecycle suite only) |
| test-workflow-execution-launch-record.R | 25 | 0 | 0 |
| test-workflow-omics-worker.R | 289 | 0 | 0 |
| test-workflow-omics-sessions.R | 150 | 0 | 0 |
| test-workflow-omics-application-reads.R | 230 | 0 | 0 |
| test-workflow-trx-declarations.R | 324 | 0 | 0 |
| test-prot-norm-worker-figures.R | 63 | 0 | 0 |
| test-prot-norm-core.R | 1854 | 0 | 0 |
| test-workflow-saved-da.R | 951 | 0 | 0 |
| test-omic-da-module.R | 1369 | 0 | 0 |

For the path repair (D4), on the step 6 tree: test-workflow-saved-identity.R 73/0/0, test-workflow-saved-resources.R
304/0/0, test-workflow-saved-recovery.R 57/0/2, test-workflow-saved-exports.R 276/0/4,
test-workflow-saved-projection.R 406/0/0, test-workflow-native-session.R 733/0/2.

Where Linux behaviour changes:

1. Linux **without** a `systemctl` binary now gets the plain launcher; it was refused before. Accepted by the user (design
   §4.1; the bundle's fallback order needs it).
2. An invalid value of the option is refused (`workflow_execution_invalid`).

Nothing else changes on Linux with systemd: every new branch is guarded by the launcher, and `winslash` is ignored outside
Windows. Claimed from the tests above; the installed suites and the journeys (step 8) have not been run.

## Claim limits (accepted by the user, 2026-11-17)

1. **No memory figure from the plain launcher.** It has no control group: no enforced worker limit (an exhausted
   allocation is the worker's own error, `worker_failed`), no slice limit, no swap control, no monitor guards, no
   accounting or peak window. Receipts carry `NA` and `measurement_valid = FALSE`. Memory figures are Linux with systemd
   only.
2. **After a controller crash on Linux only the worker process itself is killed** (`linux_pdeathsig`); what the worker
   started survives. On Windows the Job Object covers the tree (not yet run).
3. **The per-step deadline holds while the controller polls** (every 50 ms in the GUI). If the controller dies the
   lifetime binding applies instead.

## Decisions

| | Decision | By |
|---|---|---|
| D1 | `processx (>= 3.9.0)` and `ps` declared in Imports | best practice, rev 1.2 |
| D2 | no public API change; the plain launch returns when the application stops | best practice, rev 1.2 |
| D3 | no DIA lane on Windows (KIT-134 D4); the worker lane on Windows is lipidomics or metabolomics | user |
| D4 | the state owner's Windows paths repaired here, 5 hours and one correction pass | user; Linux part done in one pass |
| | Linux without `systemctl` keeps the plain launcher | user, 2026-11-17 |
| | the three claim limits above | user, 2026-11-17 |

## Open

1. **Windows (steps 0, 5, Windows half of 6): not run.** The workflow is committed (`24c346bf`); the push of branch
   `kit-131` that starts it was refused by the session's permission system, so no run exists. Budget declared: run 1 and
   at most one correction run per platform (the ticket allows 8 runs per platform).
2. **Step 7** (a lane with the worker, Linux then lipidomics or metabolomics on Windows): not started.
3. **Step 8** (R CMD check, install, the installed `protocol` and `resources` suites, the screen journeys): waits for the
   laptop (KIT-124 simulations hold 12 cores until about Wednesday).
4. Untested on Windows: the Job Object kill-on-close; `ps_kill_tree` reaching the R session behind `Rscript.exe`;
   `ps_system_memory`, `ps_create_time`; the `winslash` repair (emulated only); file locks and replacing open files;
   the GNU tool switch.
