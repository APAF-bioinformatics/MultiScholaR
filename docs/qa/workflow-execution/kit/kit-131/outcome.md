# KIT-131 outcome: the managed worker without systemd (record of 2026-11-17)

Status: **the launcher passes its lifecycle tests on Linux and on Windows (run 2). Step 8 has run on the laptop with
no KIT-131 defect found (see "Step 8 and step 7 on Linux"). Two results wait for the user's call: the `resources` suite
fails because of an old harness gap, and the worker lane differs from its oracle only in an internal store count.** Branch `kit-131` from `performance-v2` at `a98c8a6b`,
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
| CI: `.github/workflows/kit-131-launcher.yml` (ubuntu-latest and windows-latest) | steps 0 and 5; runs 1 and 2 below |
| The metadata read of a staged table: a plain stream closed on exit, not a memory map (`R/utils_workflow_saved_table_codec.R`) | correction after run 1 (Windows) |

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
| `24c346bf` | 0 and 5: the CI workflow |
| `3baf11e5` | 9: this record (first version) and NEWS |
| `2e870a4c` | the one correction after run 1 |

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
what the CI workflow runs (run 2: 185 on ubuntu, 179 on Windows; see below).

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

## Windows runs (GitHub Actions, `.github/workflows/kit-131-launcher.yml`, R 4.6.1)

Budget declared in `24c346bf`: run 1 and at most one correction run per platform (the ticket allows 8). Both used.

| Run | Commit | Link |
|---|---|---|
| 1 | `3baf11e5` | https://github.com/APAF-bioinformatics/MultiScholaR/actions/runs/43962809795 |
| 2 (the one correction) | `2e870a4c` | https://github.com/APAF-bioinformatics/MultiScholaR/actions/runs/43968604399 |

### Run 1

| | ubuntu-latest | windows-latest |
|---|---|---|
| `R CMD INSTALL` | 0 | 0 (the `configure` warning as in KIT-130) |
| Lifecycle, installed package | 183 passed, 0 failed | 98 passed, 1 failed, 1 skipped, 6 errors |
| Discovery (report only) | not run | no output: instrumentation failure |

On Windows, every test of layer 1, D4 and L2.8 passed except one assertion. That includes the Job Object taking a
unit down with a hard-killed controller (L1.4), the deadline (L1.3), and the environment and log (L1.5). The
systemd-only L2.8 test skipped.

| Failure | Class | Cause | Correction |
|---|---|---|---|
| L1.2: `identical(pids$unit, launched)` expected FALSE | instrumentation (a wrong assumption of the test) | on windows-latest the R session was the launched `Rscript.exe` itself, not its child | the platform assumption removed. In the same test the tree stop held: the unit, the grandchild and the launched process were gone, and the unit was empty |
| L2.1 to L2.6 errored at fixture setup (`manager$saveState`): `concurrent table request conflict` | **real defect of the state owner on Windows** | the staged table's commit rename failed with `cannot rename file ... reason 'Access is denied'`. `readSavedWorkflowTableMetadata` opened the payload with arrow's default memory map (`ParquetFileReader$create(path)`, `mmap = TRUE`). The mapping and its descriptor outlived `reader$.unsafe_delete()` until garbage collection, and Windows refuses to rename a directory around a mapped file. Linux allows it, so it never showed there | the metadata read opens a plain `ReadableFile` closed on exit, the pattern of `readProtImportParquetSchema`. Shown on Linux: after the old read `/proc/self/maps` and `/proc/self/fd` still list the payload; after the new one neither does. The new test fails without the repair (1 failure) and passes with it (3 of 3) |
| discovery printed nothing; its logs are empty | instrumentation | Windows cut the multi-line `Rscript -e` at its first newline | the step runs a script file |

Linux after the correction (laptop): test-workflow-execution-process.R 185/0/0, test-workflow-saved-projection.R
406/0/0, test-workflow-saved-resources.R 304/0/0, test-workflow-saved-identity.R 73/0/0.

### Run 2

| | ubuntu-latest | windows-latest |
|---|---|---|
| `R CMD INSTALL` | 0 | 0 |
| Lifecycle, installed package | **185 passed, 0 failed** | **179 passed, 0 failed, 2 skipped** |

**The plain launcher works on Windows against the installed package.** It passes every lifecycle property of the
ticket's contract:
- the worker per managed step, with the result read back and one commit (L2.1);
- failure and the forbidden commit (L2.2);
- cancellation (L2.3);
- the per-step deadline (L2.4);
- application shutdown (L2.5);
- controller death with the Job Object taking the worker, and recovery (L2.6);
- tree termination, exit status, deadline, environment, lifetime and choice (layer 1);
- D4.

The two skips are by design: the handle test reads `/proc`, which Windows has no equivalent of here (the repair
it guards is shown by L2.1 to L2.6 now passing), and the systemd launch test needs `systemctl`.

#### Discovery on Windows (report only; findings for KIT-052 and KIT-134, not corrected under this budget)

`MULTISCHOLAR_DISABLE_TOOLS=bash,dd,sha256sum,sync,tar,cp`, from source (`load_all`), one file per process.

| File | Passed | Failed | Skipped | Error |
|---|---|---|---|---|
| test-workflow-saved-identity.R | 69 | 0 | 1 | 0 |
| test-workflow-saved-resources.R | 0 | 0 | 0 | 0 |
| test-workflow-saved-recovery.R | 0 | 0 | 0 | 0 |
| test-workflow-saved-exports.R | 268 | 3 | 4 | 1 |
| test-gui-entry-lipidomics.R | not run | | | |

| Finding | Class | What |
|---|---|---|
| saved-identity: 1 skip | by design | `symlinks cannot redirect payloads outside the store` has `skip_on_os("windows")` |
| saved-resources, saved-recovery: no test ran | test gate (by design, now out of date) | both files open with a file-level `skip_if_not(sysname == "Linux", "durable saved-entity backend is Linux-qualified only")`. So on Windows the owner's resource and recovery contracts are untested, not failing. Lifting the gate and running the two files is KIT-052's to decide |
| saved-exports: `export table sharing preserves attributes and copy-on-modify history isolation` errors | test environment | `there is no package called 'lobstr'`: the test uses `lobstr`, which is not in `DESCRIPTION`, so CI does not install it. Not a product defect |
| saved-exports: `native artifact reuse skips copying and failed replacement preserves prior bytes` (2 failures) | Windows-specific test assumption | the test injects the copy failure by mocking `runDIANativeIO`, the GNU `cp` path. With `cp` disabled, `snapshotWorkflowFile` copies with `file.copy`, so the mock never fires and the copy succeeds. By reading, the test would fail the same way on Linux with `MULTISCHOLAR_DISABLE_TOOLS=cp`. The product's failed-replacement guarantee is therefore untested on the R fallback path. Finding for KIT-052 |
| saved-exports: `owned snapshots remain independent of later source writes` fails at `Sys.readlink(dangling)` | Windows-specific test assumption | `Sys.readlink` returns `""` on Windows. The product's refusal to snapshot onto the dangling link passed in the same test, just before. The product's symbolic-link guards (`Sys.readlink`) are inert on Windows, as recorded in the design (§5.1) |
| saved-exports: 4 skips | by design | withdrawn DIA-025 trials, kept for their installed candidates |
| lipidomics lane: no line | instrumentation | GitHub runs bash steps with `-e -o pipefail`. The failing `Rscript ... \| tee` of saved-exports ended the step (exit 1 at 22:25:32) before its own status check, so the lane never started. No evidence either way; the next run of the step needs `set +e` around the loop |

No discovery finding is a defect of the launcher or of the D4 repair. The one product defect found on Windows, the
memory-mapped metadata read, was corrected in run 2 and is shown by L2.1 to L2.6.

## Step 8 and step 7 on Linux (laptop, 2026-11-17, tree `ee888419`, R 4.6.1)

These checks are evidence of correctness only. They make no memory claim. Everything ran under `nice -n 19`, one R
process at a time, with thread variables set to 1. Tctl stayed between 64 and 66 °C. The evidence files are in
[step8/](step8/).

**`R CMD build`, then `R CMD check --no-manual --no-tests --no-build-vignettes`** (`_R_CHECK_FORCE_SUGGESTS_=false`,
proof library on `R_LIBS`): 3 WARNINGs and 3 NOTEs. These are the six baseline findings, and the log is line for line
the same as `kit-123/final/check.log` apart from the log path ([step8/check.log](step8/check.log)). No new finding.

**Install:** the checked tarball installed into a scratch library with exit status 0.

**Installed suites:**

| Suite | Launcher | Passed | Failed | Skipped | Error |
|---|---|---|---|---|---|
| `protocol` (`run_workflow_checkpoint_proof.sh --workflow-lifecycle`) | systemd (the frozen configuration has a slice and no `launcher`) | 223 | 0 | 7 | 0 |
| `resources` (same entry) | systemd | 296 | **1** | 1 | 0 |
| `test-workflow-execution-process.R`, installed mode, no option | systemd chosen (`workflowLauncherChoice` → `systemd`) | 185 | 0 | 0 | 0 |
| `test-workflow-execution-process.R`, installed mode, `multischolar.workflow_launcher = "process"` | plain | 185 | 0 | 0 | 0 |
| `protocol` and `resources` with the option `"process"` | plain | refused, `workflow_execution_unsupported` (by design, L2.8) | | | |

The `resources` failure is a gap in the harness. It is not a product defect, and KIT-131 did not cause it.
- **What failed.** The suite's last test, `affected saved-entity regressions preserve existing boundaries`, runs the
  five `test-workflow-saved-*.R` files against the installed package. It sources only
  `helper-workflow-checkpoint-fixtures.R`.
- **Why.** Two tests in `test-workflow-saved-recovery.R` call `skipIfNoInstalledRuntime()` (from
  `helper-installed-runtime.R`, added at KIT-034, `5250c88c`). Without that helper they stop with
  `could not find function "skipIfNoInstalledRuntime"`. The two tests are "SIGKILL at storage boundaries" and
  "concurrent writers".
- **The four other files** passed inside the suite: resources 304/0, identity 73/0, projection 406/0, exports 276/0
  with 4 skips.
- **The same file run properly.** Against the same installed build, with `helper-installed-runtime.R` sourced and
  `MULTISCHOLAR_INSTALLED_RUNTIME=1`, the file passes 108/0/0/0, both tests included.
- **Why this is not KIT-131's.** KIT-131 changes neither `test-workflow-execution.R` nor that helper
  (`git diff a98c8a6b` is empty for both). The suite fails this way on any tree since KIT-034.
- **What would fix it.** Source the helper in that test, a one-line change to a test file. That is not this step's to
  make: tests are not edited here.

**Dependents of the correction `2e870a4c`:** the test files that name `readSavedWorkflowTableMetadata` or one of
its five direct callers, minus those already passed after the correction. Run from source, one file per process:
15 files, **0 failures, 0 errors** ([step8/dependents.log](step8/dependents.log)).
- With tests run: perf-103-audit-repairs 122, perf-103 85, perf-105-audit-repairs 19, perf-105 62, perf-106 193,
  perf106-repair-v2 75, bulk-relational 45 (3 skipped), checkpoint-dia-state 174 (4 skipped), kit-harness 735,
  native-session 733 (2 skipped), partitions 438 (2 skipped), saved-qc-rollup 59 (1 skipped).
- All tests skipped, because they need the installed package: execution-adapters (6), memory-mechanism-member (13),
  rollup-qualification (7).
- The saved-entity files ran installed, inside the `resources` suite, as listed above.

**Step 7: the lipidomics entry lane with every managed step on a plain worker** (`lipid_canonical`, installed
package, option `"process"`).
- **The driver.** [step8/lane_worker.R](step8/lane_worker.R) is a scratch driver; no test file was edited. It reuses
  `helper-gui-entry.R` unchanged and gives the module's session a plain application, built as `workflowManagedMain`
  builds it. After each input it waits until the omic's runner is idle, as a user waits before the next click.
- **Attempts.** Attempt 1 stopped in the driver's own setup before the lane began (an invalid assignment target), an
  instrumentation failure. Attempt 2 completed.
  - 20 steps in the oracle's order.
  - 9 managed waits, 290 s on workers; lane 394 s.
- **Identical to the frozen oracle at all 20 steps:** head, history, state content digests, `tab_status`, every
  output, every published file and table hash, and every notification. The managed step's own "Applying lipidomics
  step ..." notice (19 of them) was set aside first; the in-process path has none.
- **The one difference.** The oracle counts the store's committed entities per masked path (KIT-101 F-13). The count
  of `dia-state-*` entities differs from `qc_intensity` on: 30 against 31 there, 107 against 144 at the end. Every
  other entity count is equal ([step8/lane-lipid-canonical-diff.txt](step8/lane-lipid-canonical-diff.txt)).
- **Reading, not tested.** This is how the managed path commits compared with the in-process one: one commit per
  action. It is not the launcher's doing, because the launcher only starts and stops processes. A systemd worker on
  this lane would settle it, but the lane runs under `testServer`, which gives no control group, so that comparison
  was not run.
- **Verdict.** The science is identical. Strict identity with the oracle is **not** met, on store internals only.
- **Not run.** Step 7 on Windows. The CI budget is spent, so it stays open.

Not run in this step: the four screen GUI journeys (design §6.3). This step's list did not include them.

## Claim limits (accepted by the user, 2026-11-17)

1. **No memory figure from the plain launcher.** It has no control group: no enforced worker limit (an exhausted
   allocation is the worker's own error, `worker_failed`), no slice limit, no swap control, no monitor guards, no
   accounting or peak window. Receipts carry `NA` and `measurement_valid = FALSE`. Memory figures are Linux with systemd
   only.
2. **After a controller crash on Linux only the worker process itself is killed** (`linux_pdeathsig`); what the worker
   started survives. On Windows the Job Object took the worker down with a killed controller (run 2, L1.4 and L2.6).
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

1. **Steps 0, 5 and the Windows half of 6: done.** The lifecycle tests pass on Windows (run 2), within the declared
   budget of two runs. Step 0's discovery is partial: the lipidomics lane did not run (instrumentation), and the
   state owner's resource and recovery files are gated to Linux. These are findings for KIT-052 and KIT-134.
2. **Step 7.** Linux is done with a scratch driver: identical except the store's `dia-state` entity count (above). The
   user decides whether that difference is accepted. Windows is not run.
3. **Step 8.** Done on the laptop, with no new check finding and no launcher failure under either launcher. The
   `resources` suite stays red until its test sources `helper-installed-runtime.R`. That gap is old (since KIT-034), and
   fixing it is a test edit the user decides on. The screen journeys were not run.
4. Still untested on Windows:
   - the state owner's resource and recovery contracts (gated to Linux);
   - the failed-replacement guarantee on the R fallback copy;
   - the lipidomics entry lane in process.
   Now shown on Windows (run 2):
   - the Job Object kill-on-close (L1.4, L2.6);
   - the tree stop (L1.2);
   - `ps_system_memory` (admission, L2.1);
   - `ps_create_time` (the controller identity, L2.6);
   - the `winslash` repair (D4 and L2.1 to L2.6);
   - file locks (the mutation lease and the action lock, L2.3, L2.5);
   - the GNU tool switch.
