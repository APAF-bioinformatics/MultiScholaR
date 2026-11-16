# KIT-131 W0: design of the second worker launcher (read-only step)

Status 2026-11-17: **W0 written, nothing implemented.** Read at `performance-v2` `5534aa97`. No R code, test or
harness file was changed; no test, journey or systemd unit was run. Every statement about Windows below is from
reading the code and the installed `processx` 3.9.0 and `ps` 1.9.3 help, not from a run. Line numbers are
`R/utils_workflow_execution.R` unless a file is named.

Bounded work item (AGENTS.md): outcome = this document; owner = `tickets/kit/KIT-131.md` (W0); files = this file
only; checks = none beyond reading; exit = the design is in front of the user.

## Summary

1. **Design choice.** One protocol, two launchers. The run, poll, accept and finish logic stays one copy. The
   systemd calls sit behind four unit functions (start, state, stop, empty) that branch on `config$launcher`. The
   plain launcher makes the launching R process the controller (the Shiny application) and starts each worker as a
   `processx` child of it. It has no slice, no monitor process and no cgroup reading.
2. **What is ported.** Worker in its own process, lifetime bound to the application, whole-tree termination, the
   per-step deadline, the environment, the log, the exit status, the result file, and the admission cap as a number
   handed to the step (it sizes the batches).
3. **What is not ported.** The 16 GiB slice limit, the worker's enforced memory limit and OOM classification, swap
   off, the three monitor guards, and every memory, CPU and I/O measurement. Consequence: no memory claim on
   Windows, and a worker there is bounded by the deadline and by the operating system only.
4. **Three findings that change the ticket's picture.**
   - Nothing in `R/` calls `launchWorkflowApplication`. Only `tools/profiling/*` does. A normal `run_app()` runs
     every heavy step in process on every platform, Linux included. The worker reaches users only when the Electron
     shell (KIT-050) starts R through `launchWorkflowApplication`.
   - The DIA precursor rollup, a step of both the DIA and the DIA_LIMPA route, supports only the native P14 path,
     which is gated to Linux x86_64 and needs the `fastduckdb` query runtime
     (`R/utils_workflow_dia_relational_rollup.R:11-13`, `:400-401`; `src/workflow_dia_ordered_reduce.c:13`, `:227-228`).
     By reading, the ticket's exit line "the DIA_LIMPA entry lane on Windows with the worker" cannot pass as the code
     stands, with or without a worker.
   - The on-disk state owner compares `normalizePath()` output with paths joined by `/` (39 calls without
     `winslash`, 9 containment checks). By reading, it refuses every payload path on Windows. The worker needs that
     owner (`:710-711`), so any worker lane on Windows needs this repaired first.
5. **Estimate.** 28 hours for the launcher, its tests on Linux and the process-level tests on Windows. A bounded
   5 hours more for the state owner's path handling on Windows, if the user puts it in this ticket. The rollup on
   Windows is not estimated here (KIT-052).

## 1. Inventory

### 1.1 Where systemd is used in `R/`

All of it is in `R/utils_workflow_execution.R`. Three commands are run: `systemctl --user`, `systemd-run --user`
and `busctl --user`, all through `workflowSystemd` (`:395-405`).

| Function | Lines | What it asks systemd for |
|---|---|---|
| `workflowSystemd` | 395-405 | Runs one command with `--user` and an argument vector. A missing binary or a failed required command raises `workflow_execution_unsupported`. |
| `workflowServiceProperties` | 411-423 | `systemctl show` of one unit: `ActiveState`, `SubState`, `Result`, `MainPID`, `ControlGroup`, `MemoryPeak`, `MemorySwapPeak`, `CPUUsageNSec`, `IOReadBytes`, `IOWriteBytes`, `MemoryMax`, `MemorySwapMax`. This is how the exit status is collected. |
| `createWorkflowApplicationSlice` | 2113-2124 | `busctl call ... StartTransientUnit` for the application slice: `MemoryAccounting`, `CPUAccounting`, `IOAccounting`, `MemoryMax` = 16 GiB (4 GiB diagnostic), `MemorySwapMax` = 0. Returns the slice's control group. |
| `workflowLaunchService` | 1167-1194 | `systemd-run` of one of three fixed entry points (`application`, `monitor`, `worker`) inside the slice. Properties in the next table. |
| `launchWorkflowApplication` | 1052-1155 | `show-environment` as the probe (1060); the slice (1121); the application monitor (1137, `gui` mode); the controller (1148); `stop` of the slice if the launch fails (1123). |
| `checkWorkflowApplication` | 991-1010 | `show` of the controller unit and, in `gui` mode, of the monitor unit: the controller's `MainPID` is this process, both sit under the slice, and the slice's `memory.max` and `memory.swap.max` equal the policy. |
| `workflowBrowserMembership` | 963-985 | `show` of the harness's proof browser unit, only when the harness wrote `application-browser-member.rds`. |
| `runWorkflowAction` | 1633-1638 | Starts the per-request monitor unit (modes other than `gui`). |
| `pollWorkflowAction` | 1672, 1704-1707, 1716, 1728 | `show` of the monitor; start of the worker unit with its cap and deadline; `show` of the worker. |
| `finishWorkflowExecutionJob` | 1910-1913, 1928, 1988 | `show` then `stop` of the worker; `stop` of the per-request monitor. |
| `cancelWorkflowAction` | 361 | `stop` of the worker. |
| `workflowMonitorMain` | 1414, 1449 | The monitor stops the slice (application scope) or the worker (request scope) when a guard fires. |

What `workflowLaunchService` sets on a unit (`:1173-1191`):

| Property | Value | Purpose |
|---|---|---|
| `--slice` | the application slice | The unit is charged to the 16 GiB aggregate limit. |
| `Type=exec` | | `systemd-run` returns once the binary has been executed, so a launch failure is seen at launch. |
| `KillMode=control-group`, `TimeoutStopSec=1` | | `stop` signals every process of the unit, and kills what is left after one second. |
| `MemoryAccounting`, `CPUAccounting`, `IOAccounting` | yes | Counters for the receipts. |
| `MemorySwapMax=0` | | No swap for the unit. |
| `RuntimeMaxSec` | worker: the per-step deadline (`workflowActionWorkerDeadline`, 2962-2964; policy 1800 s, diagnostic 90 s; no operation names its own today); per-request monitor: deadline + 60; controller: 86400 s or the journey's deadline; application monitor: 86500 s | systemd ends the unit at the limit, whatever the controller is doing. |
| `StandardOutput=append:<root>/<role>-<request file>.log`, `StandardError=inherit` | | One appended log per role. |
| `--setenv` | `R_LIBS`, `OMP_NUM_THREADS=1`, `OPENBLAS_NUM_THREADS=1`, `MKL_NUM_THREADS=1`, `DUCKPLYR_FALLBACK_AUTOUPLOAD=0`, `DUCKPLYR_FALLBACK_COLLECT=0`, `TMPDIR=<request dir>/<role>-scratch`, and for a worker `MULTISCHOLAR_WORKFLOW_WORKER=1` | The unit does not inherit the caller's environment; it gets the user manager's plus these. |
| `MemoryMax=<cap>`, `OOMPolicy=kill` | worker only | The admitted cap, enforced by the kernel; an OOM kill ends the whole unit with `Result=oom-kill`. |
| `BindsTo=`, `After=` | worker: controller and monitor; per-request monitor: controller; controller: application monitor; application monitor: the harness's supervisor unit | Lifetime binding: the unit stops when what it is bound to stops. |
| `--remain-after-exit` | worker only | The unit stays readable after the process ends, so `Result` can be collected. |
| `--wait` | controller in the `protocol` and `resources` suites and diagnostic `adapters` | The launcher blocks until the suite ends. |
| Working directory | not set | A user unit starts in the home directory. Nothing in the code relies on it. |
| Command | `<R bin>/Rscript --vanilla -e MultiScholaR:::workflow<Role>Main(commandArgs(TRUE)) <request path> <sha256>` | The entry reads a hashed, bounded request file (`readWorkflowExecutionRequest`, 1200-1209). |

Files of the control group that the code reads directly (they exist because the slice exists):

| Reader | Lines | File | Callers |
|---|---|---|---|
| `workflowMemoryCounter` | 429-434 | `memory.max`, `memory.swap.max`, `memory.current`, `memory.peak` | 939, 1006-1007, 1125-1126, 1344, 1431, 1686, 1972 |
| `workflowWorkerMemory` | 443-448 | the worker unit's `memory.current` | 1368, 1443 (monitor) |
| `openWorkflowPeakWindow`, `closeWorkflowPeakWindow` | 1016-1042 | `memory.peak`, reset by a write (Linux 6.12 or later) | 835, 932, 1547, 1966 |
| `workflowCgroupAccounting` | 2023-2035 | `cpu.stat`, `io.stat` | 847, 933, 1576, 1967 |
| `workflowServiceEmpty` | 2042-2049 | `cgroup.procs` of the unit | 1930, 1932 |
| `workflowMonitorMain` | 1346, 1433 | `memory.stat` (`anon`, `file`) | monitor loop |

Other Linux-only reads on the same path: `workflowHostAvailable` (`/proc/meminfo`, 453-461; callers 1058, 1345,
1432, 1560, 1686), `workflowExecutionProcess` (`/proc/<pid>/stat` and the boot id, 2054-2065), `/proc/self/cgroup`
(1062, 1489).

Outside this file, `assertDIACompatibilityGuard` (`R/utils_workflow_dia_sessions.R:1699-1718`) reads the control
group of the older proof service for the explicit S4 compatibility export. It is not on the worker path and is not
touched.

### 1.2 The processes of a managed application today

`launchWorkflowApplication(root, mode = "gui", port)` creates the slice, starts the application monitor and the
controller as units, saves the configuration to `<root>/application.rds`, and returns. The controller is an R
process that runs `workflowManagedMain` (1215-1251): it builds a `WorkflowExecutionApplication`, checks it, and
runs `shiny::runApp(run_app(workflow_execution = app))`. Each heavy step then starts one worker unit. In modes
other than `gui` each action also starts its own monitor unit; in `gui` mode the application monitor is shared.

### 1.3 Caller paths, from the GUI down

| Path | Code | What happens |
|---|---|---|
| Launch | `launchWorkflowApplication` 1052-1155 → `workflowLaunchService` → `workflowManagedMain` 1215 → `run_app(workflow_execution = app)` (`R/run_app.R:62-67`) → `workflowManagedServer` 2011-2017 | The session gets `session$userData$workflow_execution`. |
| Dispatch rule | `workflowSessionApplication` 509-512 | Returns the application, or NULL for an unmanaged session. |
| A module's Apply button | 19 `installWorkflowManagedStep` sites (17 in `R/mod_prot_qc_*.R`, `R/mod_prot_da_handlers_run_analysis.R:18`, `R/utils_workflow_omics_actions.R:54`), each behind `if (!is.null(workflowSessionApplication(session)))` | `step$apply` 692-747. |
| One heavy action per application | `beginWorkflowGUIAction` 819-879: `application$active` and the lock file `<root>/.heavy-action.lock` (836) | A second action is refused as busy. |
| Queue | `installWorkflowOmicsActions` (`R/utils_workflow_omics_actions.R:13-104`): a tab's actions run in order, a click during another action waits in `run$queued`; `protNormFiguresDispatcher` (`R/mod_prot_norm_workflow_helpers.R:181-197`) | The queue is in the controller and does not know the launcher. |
| Managed runner | `runWorkflowAction` 1533-1651 (request, token, journal, control record, monitor) → `pollWorkflowAction` 1657-1744 (admission cap 1685-1692, worker start 1704, state 1716, result 1737-1743) → `acceptWorkflowActionResult` 1764-1823 → `finishWorkflowExecutionJob` 1904-2005 | Polled every 50 ms by the step's observer (748-770). |
| Per-step deadline | `workflowActionWorkerDeadline` 2962-2964 → `RuntimeMaxSec` at 1707 | A unit that hits it ends with a `Result` other than `success`, read as `worker_failed` (1732-1735). |
| Cancellation | `step$cancel` 682-691 → `cancelWorkflowAction` 352-368 → `stop` of the worker (361) → `finishWorkflowExecutionJob`. Triggers: Revert (`R/utils_workflow_dia_state.R:1753`), project replacement (751-753), session end (808), a failed step (645) | The generation is advanced first, so a late result is stale. |
| Application shutdown | last session ends: `R/app_server.R:28-37` closes every controller; process exit: `closeWorkflowManagedApplication` 1281-1294; store close: `R/utils_workflow_saved_entities.R:175`; `closeWorkflowExecutionController` 374-387 cancels a running job | The slice and `BindsTo` take whatever is left. |
| Worker entry | `workflowWorkerMain` 1506-1523 → `computeWorkflowPeptideAction` 2420, `computeWorkflowOmicsAction` 2444, or the diagnostic fixture 1465 | Writes `worker-result.rds`; a failure also exits non-zero. |
| Callers of `launchWorkflowApplication` | `tools/profiling/workflow_kit_coordinator.R:205`, `:285`, `workflow_checkpoint_dia_driver.R:1253`, `workflow_perf106_repair.R:323`, `:334`, `run_workflow_checkpoint_proof.R:2556` | **No caller in `R/`.** |

The harness also binds these internals by name, so their names, arguments and Linux behaviour are frozen:
`workflowSystemd` (13 uses), `workflowMemoryCounter` (7), `launchWorkflowApplication` (6), `workflowHostAvailable`
(5), `workflowServiceProperties` (3), `workflowCgroupAccounting` (2), `workflowExecutionBuildIdentity` (4).

## 2. What depends on each property

"Application" means a step's result, its state, or the application staying usable. "Harness" means only the
receipts that `tools/profiling/*` reads.

| Property | Who depends on it | Evidence | Plain launcher |
|---|---|---|---|
| The worker is a separate process; only the controller moves the head | Application (results, and the memory design: the step's memory goes back to the system when the worker exits) | `assertWorkflowController` 106-110; `workflowWorkerMain` 1508; acceptance under the lease 1764-1823 | Ported. |
| Exit status and result file | Application | 1732-1743: `Result == "success"` and `worker-result.rds`; the worker's own error text 1751-1757 | Ported: exit code and the same file. |
| Whole-tree stop and the emptiness check | Application (Revert, session end, project replacement; ownership is kept while anything is left: `cleanup_failed`, 1930-1937) | 361, 1913, 2042-2049 | Ported. |
| Lifetime binding of the worker to the controller | Application, as hygiene. An orphan cannot move the head (it has no lease, and the token's generation is checked at 1799-1802), but it keeps computing and holds files | `BindsTo` 1185-1186, 1706 | Ported. |
| Per-step deadline | Application (a step is bounded; user decision 2026-10-23) | 1707, 2962-2964 | Ported, enforced by the controller's poll. |
| Environment: `R_LIBS`, one thread, `TMPDIR`, worker flag | Application: the same library; `TMPDIR` is what the controller deletes after an uncatchable stop (`R/utils_workflow_dia_state.R:917-920`, 1951-1955); the worker flag guards the head (107, 1508) | 1178-1190 | Ported, same values. |
| Admission cap as a number | Application: the handlers size their batches from it and refuse an estimate above it (`R/utils_workflow_partitions.R:33-41`, `:108-109`, `:162-163`; `request$worker_cap`, 1700) | `workflowAdmissionCap` 485-495 | Ported as a number: host free memory from `ps`. |
| Admission cap enforced (`MemoryMax`, `OOMPolicy=kill`) | Application safety, and one status (`oom-kill` → `memory_guard`, 1733) | 1184 | Not ported. |
| Slice `MemoryMax` 16 GiB, swap off | Application safety envelope; the application refuses to run without it (1006-1008), but no result depends on it | 2116-2118 | Not ported. |
| Monitor guards: `memory_guard` (15.5 GiB), `host_headroom` (2 GiB free), `nonworker_reservation` | Application safety; they surface as a step's failure reason | 1351-1352, 1369, 1437-1444 | Not ported. |
| Monitor observation: `max_anon`, `max_file`, sampling interval | Harness (`measurement_valid`) | 945-948, 1982-1983 | Not ported; `measurement_valid = FALSE`. |
| Accounting and the peak window | Harness: `peak_bytes`, `whole_service_peak_bytes`, `accounting_delta` in the receipts. The application never decides on them, but today it fails if they cannot be read (1019-1021) | 934-949, 1966-1975 | Not ported; `NA` in the receipts. |
| Controller and monitor `RuntimeMaxSec` | Harness and journey bound | 1139, 1149-1151 | Not ported; the shell that started R owns the application's lifetime. |
| Log file | Diagnosis only | 1177 | Ported (one file per unit). |
| Working directory | Nobody | not set | Set explicitly to the unit's scratch directory. |
| `--wait`, supervisor binding, browser membership, guard injection, journey record | Harness | 1061-1070, 1085-1087, 1182-1183, 1390-1400 | Refused under the plain launcher. |

## 3. Worker or in process: how the code decides, and what Windows does today

The rule is one function. `workflowSessionApplication(session)` (509-512) returns the session's
`WorkflowExecutionApplication`, or NULL. Only `workflowManagedServer` (2011-2017) sets it, and only
`workflowManagedMain` builds the server that way, inside a controller unit started by `launchWorkflowApplication`.

- **With an application:** the module calls `managed_step$apply(...)`. `step$apply` also requires the on-disk state
  owner (710-711) and refuses otherwise.
- **Without one:** the same observer runs the same registered declaration in the Shiny process, for example
  `runPeptideRollupApplyObserver` (`R/mod_prot_qc_peptide_rollup.R:232-246`). The omics tabs get no runner
  (`R/utils_workflow_omics_actions.R:14`) and call their wrappers in process. The normalisation figures are drawn in
  process (`R/mod_prot_norm_workflow_helpers.R:182-183`). The DA publication runs in process
  (`R/mod_prot_da_handlers_run_analysis.R:68`). An in-process step takes the default cap of 12 GiB
  (`R/utils_workflow_dia_qc_rollup.R:1304-1305`).

On Windows today:

- `launchWorkflowApplication` refuses at its first line (1054-1056), and `checkWorkflowApplication` refuses any
  application on a system other than Linux (992-994). There is no managed application.
- So `run_app()` runs every heavy step inside the one Shiny process. What a step allocated stays with that process;
  the allocator trim is compiled for Linux with glibc only (`src/workflow_memory.c:14`, `:23-24`).
- The same is true of a plain `run_app()` on Linux, because nothing in `R/` calls `launchWorkflowApplication`.
- By reading, the on-disk owner does not get that far on Windows (section 5.2). This must be observed on a runner
  before anything is built on it.

## 4. The launcher interface

### 4.1 Shape

The launcher is a value in the frozen configuration, not a class hierarchy. `config$launcher` is absent for systemd
(so `application.rds` and every request file are what they are today) and `"process"` for the plain launcher.
The choice is made once, in `launchWorkflowApplication`:

1. `getOption("multischolar.workflow_launcher")` when set (`"systemd"` or `"process"`). This is the option the tests
   force.
2. Otherwise `"process"` when the system is not Linux, or when Linux has no `systemctl` binary.
3. Otherwise systemd, exactly as today, including today's refusal when the user manager does not answer. The
   harness therefore never lands on the plain launcher by accident.

Four unit functions carry the contract. In each, the systemd side is the existing call, unchanged:

| Function | Contract | systemd | process |
|---|---|---|---|
| `workflowLaunchService(config, unit, role, path, hash, cap, binds, deadline)` (existing name and arguments) | Start one package entry point for a hashed request file. Returns after the start is accepted; raises `workflow_execution_unsupported` if it cannot start. | lines 1167-1194 as they are | `processx::process$new` of the same `Rscript` command line and the same environment values; `binds` and `cap` are not enforced |
| `workflowUnitState(config, unit)` (new) | A list with `available`, `SubState` (`running`, `exited`, `failed`), `Result` (`success`, `exit-code`, `signal`, `timeout`), `MainPID`. A unit never started gives `available = FALSE`. Enforces the deadline when asked. | `workflowServiceProperties(unit)` | from the handle: alive, exit status, elapsed time against the deadline |
| `workflowUnitStop(config, unit)` (new) | Stop the unit and everything it started. Returns when they are gone, or after a bounded wait. Never an error for a unit that is not running. | `workflowSystemd(c("stop", unit), required = FALSE)` | kill the marked tree, then wait up to 2 s |
| `workflowServiceEmpty(application, unit)` (existing) | TRUE when no process of the unit is left. | lines 2042-2049 as they are | the handle is dead and no process carries the unit's marker |

Because `pollWorkflowAction` and `finishWorkflowExecutionJob` already read only `available`, `SubState` and `Result`,
their logic does not change; `Result` is never `oom-kill` under the plain launcher.

The handles live in a package-level environment keyed by unit name, as `.workflowToolCache` and the operation
registry do. They are meaningful only in the process that started them, which is the controller.

The existing source tests mock `workflowSystemd`, `workflowServiceProperties`, `workflowLaunchService`,
`workflowServiceEmpty`, `workflowMemoryCounter`, `workflowCgroupAccounting`, `openWorkflowPeakWindow`,
`closeWorkflowPeakWindow`, `workflowHostAvailable` and `checkWorkflowApplication` by name. The systemd side of
every new function calls those same functions, so those tests keep testing the systemd path.

### 4.2 The plain launcher, piece by piece

- **Topology.** `launchWorkflowApplication(root, mode = "gui", port)` writes `application.rds` and then runs
  `workflowManagedMain` in the calling process. The call returns the configuration when the application stops. Under
  systemd it returns after the units are started. The shell that starts R (KIT-050) must handle both.
- **Worker.** One `processx` child per action, with `cleanup = TRUE`, `cleanup_tree = TRUE`, `stdout` to
  `<root>/worker-<unit>.log`, `stderr = "2>&1"`, `wd` = the unit's scratch directory, `windows_hide_window = TRUE`,
  the environment inherited plus the eight values of section 1.1, plus one marker variable named after the unit.
- **No monitor.** The monitor's work is cgroup reading. Under the plain launcher `job$monitor_unit` is NULL and the
  monitor waits are skipped (otherwise each action would wait 2 s for a snapshot that never comes: 925-926, 1926).
- **No measurement.** The peak window is NULL, accounting is `NA`, `measurement_valid` is FALSE. The receipts are
  still written, so the tests can read status and cleanup from them.
- **Admission cap.** `workflowAdmissionCap(policy, current, host)` as today, with `host` =
  `ps::ps_system_memory()$avail` where the system is not Linux, and `current` = the controller's resident memory
  from `ps`. On Linux `workflowHostAvailable` keeps its `/proc/meminfo` code.
- **Scope.** `gui` mode only. The harness's `protocol`, `resources` and `adapters` suites, a supervisor unit and a
  journey record are refused under the plain launcher: they assert cgroup facts.

Edits in `R/utils_workflow_execution.R`, all guarded by `identical(config$launcher, "process")`:

| Function | Edit |
|---|---|
| `launchWorkflowApplication` | Choose the launcher; under process skip the Linux gate, the systemd probe, the slice and the monitor, and run the controller here. |
| `checkWorkflowApplication` | Under process check the class, the owner pid and the launcher only. |
| `workflowManagedMain` | `app$cgroup` is NULL under process. |
| `beginWorkflowGUIAction`, `finishWorkflowGUIAction` | No peak window, no accounting, no monitor snapshot wait. |
| `runWorkflowAction` | No peak window, no accounting, no monitor unit. |
| `pollWorkflowAction` | No monitor state or readiness wait; cap from `ps`; state through `workflowUnitState`. |
| `finishWorkflowExecutionJob`, `cancelWorkflowAction` | Stop through `workflowUnitStop`; no monitor wait; `NA` measurements. |
| `workflowHostAvailable`, `workflowExecutionProcess` | A branch for systems other than Linux (`ps`). Linux code unchanged. |

No new file: `Collate` and `NAMESPACE` stay as they are (`processx::` and `ps::` are called with their namespace).

### 4.3 Dependencies

`DESCRIPTION` names neither `processx`, `callr` nor `ps`. All three are already installed with the package:
`multidplyr` (Imports) imports `callr`, which imports `processx`, which imports `ps`. The proposal is to declare
`processx (>= 3.9.0)` and `ps` in Imports. Nothing new enters the bundle. Version 3.9.0 is the first with
`linux_pdeathsig`. Base R is not enough: `system2(wait = FALSE)` gives no process id, no exit status and no tree
kill on Windows.

### 4.4 What each property becomes

| Property | systemd | Windows, plain launcher | Linux, plain launcher (forced) |
|---|---|---|---|
| Lifetime bound to the application | `BindsTo` and the slice | `processx` puts the child in a Job Object with kill-on-close, so the worker and its children die when the controller exits or crashes (installed `processx` help, "Automatic GC cleanup") | The exit finalizer kills the child's process group; `linux_pdeathsig = SIGKILL` covers a crash, for the direct child only |
| Whole-tree stop | `KillMode=control-group`, 1 s grace | `ps::ps_kill_tree(marker)`: every process that inherited the marker variable is terminated. Needed even for a plain worker: `Rscript.exe` is a front end and the R session is its child (to be confirmed by the test: the pid in `worker-started.rds` differs from the launched pid) | The same call sends SIGTERM, then SIGKILL after 200 ms |
| Emptiness | `cgroup.procs` | `ps::ps_find_tree(marker)` is empty | the same |
| Per-step deadline | `RuntimeMaxSec`, enforced by systemd | The controller's poll stops the tree at the deadline and reports `Result = "timeout"`. It holds while the controller is alive and polling; if the controller dies, the Job Object applies | the same; after a controller crash only the direct child is covered |
| Exit status | `Result`, `--remain-after-exit` | `get_exit_status()`; the handle is kept until the job is finished | the same |
| Memory limit of the worker | `MemoryMax`, `OOMPolicy=kill` | None. `processx` exposes no Job Object memory limit, and writing one would be new native code. An exhausted allocation is an R error in the worker, reported as `worker_failed` with its message | None; the kernel's OOM killer gives `worker_failed` |
| Slice limit, swap off, monitor guards | yes | None | None |
| Memory, CPU, I/O accounting | cgroup v2 | None | None |
| Environment | manager's environment plus the listed values | the controller's environment plus the same values | the same |
| Log | appended, one per role | one file per unit, so no worker overwrites another's | the same |

**Explicitly not ported:** the enforced worker limit and its OOM classification, the 16 GiB slice limit, swap off,
the three monitor guards, all accounting and the peak window, the controller's own runtime limit, and the harness's
suites. **Consequence:** no memory number may be claimed for Windows. What Windows keeps is the design itself: the
application process stays small and a step's memory is returned when its worker exits.

## 5. Other Linux-only code that would break on Windows

### 5.1 On the launcher path (this ticket)

| Item | Where | Note |
|---|---|---|
| Linux gate | 993, 1054 | Relaxed under the plain launcher only. |
| `systemctl`, `systemd-run`, `busctl` | 395-405, 1192, 2120 | Behind the seam. |
| `/proc/meminfo` | 453-461 | `ps` on other systems. |
| cgroup files | 429-448, 1016-1042, 2023-2049, 1346, 1433 | Not read under the plain launcher. |
| `/proc/self/cgroup` | 1062, 1489 | Supervisor check and a diagnostic fixture; not reached. |
| `/proc/<pid>/stat`, boot id | 2060-2064 | The branch for other systems (2057-2059) identifies a controller by pid alone. Windows reuses pids quickly, so a restart after a crash could, rarely, be refused as "this project already has a managed controller" (129-131). Adding the start time from `ps::ps_create_time()` closes it; small, optional. |
| Fork | 1484-1494 (`parallel::mcparallel`) | The `allocate` fixture only; it already refuses on a non-POSIX R. |
| `Rscript` path | 1191 | `file.path(R.home("bin"), "Rscript")` has no `.exe` and starts a child R on Windows; see 4.4. |
| Signals | none | `R/` sends no signal itself; every stop goes through `systemctl stop`. |
| File locks | 133, 836-837, 1562-1563 (`filelock`) | `filelock` supports Windows. Not run there yet. |
| Symbolic link checks | 85, 966, 1954 (`Sys.readlink`) | Returns `""` on Windows; harmless. |
| Slash-joined path check | 2085 | Diagnostic fault mode only. |

### 5.2 In the work the worker does (not the launcher, but it gates a Windows lane)

| Item | Where | Note |
|---|---|---|
| **Precursor rollup is Linux x86_64 only** | `R/utils_workflow_dia_relational_rollup.R:11-13` (gate), `:15` (`inst/duckdb/linux_amd64`), `:400-401` (only the native P14 route); `src/workflow_dia_ordered_reduce.c:13`, `:227-228`; the step is in both DIA routes (`R/mod_prot_qc_peptide.R:29-33`); the in-process call is the same (`R/utils_workflow_dia_qc_rollup.R:1326-1330`) | Also needs `fastduckdb`, which CI cannot read. Blocks the DIA and DIA_LIMPA routes on Windows in any process. KIT-052's subject. |
| **Path separators in the state owner** | `R/utils_workflow_saved_entities.R:120`, `:266-267`, `:742-745`, `:1100`, `:1509`; `R/utils_workflow_saved_exports.R:617-619`; `R/utils_workflow_dia_sessions.R:1760-1761`; 39 `normalizePath()` calls without `winslash` across `R/utils_workflow_*.R` | On Windows `normalizePath()` returns backslashes and the checks compare against `paste0(dir, "/")`. By reading, `checkSavedEntityPath` refuses every payload. `winslash = "/"` is ignored on Linux, so the repair cannot change Linux. |
| GNU tool probe | `R/utils_workflow_saved_entities.R:678-695`; pipelines at `R/utils_workflow_dia_sessions.R:108-122`, `R/utils_workflow_dia_state.R:915-931` | KIT-044 gave every tool an R fallback. `bash` is accepted by presence alone (691). A Windows machine with Rtools or Git on the PATH (a CI runner) would select the GNU pipelines, which quote for a POSIX shell. The existing switch `MULTISCHOLAR_DISABLE_TOOLS` (681) selects the fallbacks; the Windows job should set it, and the bundle has none of these tools. |
| Durability | `R/utils_workflow_saved_entities.R:707`, `:746` | Without GNU `sync` the declared contract is "not guaranteed". Existing, recorded. |
| Replacing and deleting open files | `writeWorkflowExecutionMessage` 1315-1323 and the store's renames | Windows refuses to replace or delete a file another process holds open. The plain launcher has no monitor reading beside the controller, and a worker's handles close when it exits. Not run yet. |
| Native build flags | `src/Makevars:3-4`, `:8-12` (`uname`, `-Werror`, `-fPIC`); no `configure.win` | The package has not been installed on Windows in this series (KIT-130 does not install it). The first install is a finding in itself. |

## 6. Test plan

The worker is `Rscript` loading the installed package, and the request carries the installed build's identity
(`workflowExecutionBuildIdentity`, 96-101). So these tests run against an installed package, not `load_all`, and
skip otherwise. Each starts at most one worker at a time.

### 6.1 Lifecycle tests, new file `tests/testthat/test-workflow-execution-process.R`

Layer 1 needs no project store: it calls the process side of the unit functions with a trivial R command. It runs
on Linux (option forced) and on `windows-latest` from the first day.

| # | Property | Test |
|---|---|---|
| L1.1 | Start and exit status | A command that exits 0 gives `Result = "success"`; one that stops gives `exit-code`; a command that cannot start raises `workflow_execution_unsupported`. |
| L1.2 | Whole-tree stop | The unit starts a grandchild and records both pids. After `workflowUnitStop` neither is running and `workflowServiceEmpty` is TRUE. On Windows also: the recorded R pid differs from the launched pid and is gone. |
| L1.3 | Deadline | A 30 s sleep with a 2 s deadline reports `timeout` within about 3 s, and the tree is gone. |
| L1.4 | Lifetime bound to the application | A child R process acts as controller, starts a sleeping unit, and is then killed hard. The unit's pid disappears within 5 s. On the Windows runner this also exercises nested Job Objects. |
| L1.5 | Environment, working directory, log | The unit writes its `MULTISCHOLAR_WORKFLOW_WORKER`, `TMPDIR`, thread variables, `R_LIBS` and working directory to a file; a message to standard error is in the unit's log. |
| L1.6 | A unit never started | `workflowUnitState` gives `available = FALSE`; `workflowUnitStop` is silent. |
| L1.7 | Launcher choice | The option wins; a system that is not Linux gives `process`; Linux with `systemctl` gives systemd (mocked). |

Layer 2 drives the real protocol with the diagnostic fixture (`wfx_fixture`) through a plain application whose
controller is the test process. It runs on Linux now, and on Windows once the state owner works there.

| # | Property | Test |
|---|---|---|
| L2.1 | Result read back, one commit | `identity` completes; the head advances once; the receipt has `cleanup$worker_stopped`, `measurement_valid = FALSE`, `NA` peak. |
| L2.2 | Failure | `fail` gives `worker_failed` with the worker's message; `forbidden_commit` gives `worker_failed` (the worker flag reached the child); the head is unchanged. |
| L2.3 | Cancellation | `sleep` 20 s, then `cancelWorkflowAction`: the worker is gone, the head and history are unchanged, the action lock is free, the scratch directory is removed. |
| L2.4 | Per-step deadline | With `workflowActionWorkerDeadline` mocked to 2 s, `sleep` 20 s ends as `worker_failed` within about 3 s. |
| L2.5 | Application shutdown | `closeWorkflowManagedApplication` with a running worker: the worker is gone, the controller is closed, the lease can be taken again. |
| L2.6 | Controller death | A child controller is killed while its worker sleeps. A new controller on the same project finds the head unchanged, records `controller_lost` and can run an action. |
| L2.7 | Admission | `admitted.rds` holds a finite positive cap and the request carries it. |
| L2.8 | Refusals | A non-`gui` mode, a supervisor unit or a journey record under the plain launcher is refused. |

### 6.2 A lane with the worker

`tests/testthat/helper-gui-entry.R` drives a module server under `shiny::testServer` with no managed application.
A variant gives the session a plain application; after each step the driver advances the session's timers until no
managed action is active and answers the browser acknowledgement, as `test-perf107-audit.R` already does by hand.

**Assumption A1.** The comparison is the lane with the worker against the same lane in process **on the same
machine** (head, history, state digests, table hashes). It is not a comparison with the oracles frozen on Linux
under R 4.6.1: agreement across platforms is KIT-052's tolerance comparison. On CI the engine is `limpa` and the
reader `duckdb` on both sides, because the two internal Suggests are not readable there.

Which lane runs on Windows depends on decision D3 (section 8): `prot_dia_limpa` cannot pass the rollup step there.

### 6.3 What proves Linux with systemd is unchanged

| Evidence | What it covers |
|---|---|
| `tests/testthat/test-workflow-execution.R` (source cases from line 528, systemd mocked), `test-workflow-execution-launch-record.R`, `test-perf-107.R`, `test-perf107-audit.R`, `test-workflow-omics-worker.R`, `test-prot-norm-worker-figures.R`, `test-prot-norm-core.R` | The controller logic through the mocked systemd functions. |
| `test-workflow-kit-harness.R` | The harness's reading of the receipts and the launch. |
| Installed lifecycle suites `protocol` and `resources` (`tools/profiling/run_workflow_checkpoint_proof.sh --workflow-lifecycle`) | Real units: one commit, cancellation, controller death, containment, guards, peak window. |
| The four screen GUI journeys of the kit, exact against their records (last at KIT-123, `docs/qa/workflow-execution/kit/kit-123/journeys.md`) | The ticket's "Linux journeys exact". |
| `git diff --stat -- tools/profiling` empty | The harness is untouched. |
| `tests/testthat/test-gui-entry-*.R` | The in-process path is unchanged. |

The installed suites and the journeys start real units and take long. Their purpose, bindings, ceiling and attempts
are declared at that step under the kit protocol, not here.

### 6.4 Windows runner

A separate workflow file on a `kit-131` branch (`workflow_dispatch` and pushes to that branch), jobs
`ubuntu-latest` (option forced) and `windows-latest`: install the dependencies and the package, then run the new
test file against the installed package with `MULTISCHOLAR_DISABLE_TOOLS=bash,dd,sha256sum,sync,tar,cp` on Windows.
Budget to declare at the step: at most 8 workflow runs per platform, one correction pass.

## 7. Implementation plan

| # | Step | Files | Risk | Hours |
|---|---|---|---|---|
| 0 | **Windows discovery, no package code.** The workflow file; on `windows-latest` install the package, run the state-owner unit tests and one entry lane in process, keep every failure. Sizes D3 and D4 with evidence instead of reading. | `.github/workflows/kit-131-launcher.yml` (new, branch only) | Low. The first Windows install may itself fail (`src/Makevars`). | 3 |
| 1 | **The seam, systemd side only.** `workflowUnitState`, `workflowUnitStop` and the two memory helpers, each calling the existing function; call sites switched. A pure move, checked by the mocked source tests before anything else changes. | `R/utils_workflow_execution.R` | Medium: it edits the path Linux runs. | 2.5 |
| 2 | **The process side.** `processx` and `ps` in Imports; the launcher choice; start, state, stop, empty for a process; host memory and process identity where the system is not Linux. | `DESCRIPTION`, `R/utils_workflow_execution.R` | Medium. | 3 |
| 3 | **The application under the plain launcher.** The edits of the table in 4.2. | `R/utils_workflow_execution.R` | Medium to high: about 25 small guarded edits in the most delicate file. | 4 |
| 4 | **Lifecycle tests on Linux**, option forced: layers 1 and 2. | `tests/testthat/test-workflow-execution-process.R` (new) | Medium. | 3.5 |
| 5 | **Layer 1 on Windows**, with one correction pass. | the workflow file | Medium. | 4 |
| 6 | **Only if D4 says so: the state owner's paths on Windows** (`winslash = "/"` at the 39 calls, the 9 checks), then layer 2 on Windows. One correction pass, then stop and report. | `R/utils_workflow_saved_entities.R`, `R/utils_workflow_saved_exports.R`, `R/utils_workflow_dia_sessions.R`, others under `R/utils_workflow_*.R` | High: nobody has run the owner on Windows; file locking may follow the paths. | 5 |
| 7 | **A lane with the worker**: the driver variant; Linux first, then the lane D3 names on Windows. | `tests/testthat/helper-gui-entry.R`, one new test file | High: timers and acknowledgement under `testServer`. | 4 |
| 8 | **Linux unchanged**: targeted files, dependents grep, `R CMD check`, install, the installed suites, the screen journeys. | none | Low for the code; long runs. | 3 |
| 9 | **Record**: outcome, ticket, NEWS. | `docs/qa/workflow-execution/kit/kit-131/outcome.md`, `tickets/kit/KIT-131.md`, `NEWS.md` | Low. | 1 |

Total 28 hours without step 6, 33 with it. The ticket says 24. Each step is committed on its own, so a regression
is traced to one mechanism.

## 8. Decisions for the user

| # | Decision | Recommendation |
|---|---|---|
| D1 | **New declared dependencies:** `processx (>= 3.9.0)` and `ps` in Imports. | Yes. Both are already installed through `multidplyr`; nothing new enters the bundle. |
| D2 | **Public API.** None is changed: `run_app()` and the arguments of `launchWorkflowApplication` stay as they are, and the launcher is chosen by an option. Two things follow. Under the plain launcher `launchWorkflowApplication` returns only when the application stops. And a user who calls `run_app()` directly still gets no worker, on any platform; changing that would change an exported function. | Leave both as they are. KIT-050's shell calls `launchWorkflowApplication`. |
| D3 | **The exit line "the DIA_LIMPA entry lane on Windows with the worker" cannot pass**: the precursor rollup is Linux x86_64 only. | Prove the worker on Windows with a lane that has no rollup (a lipidomics or metabolomics lane; their tabs use the same runner), and leave the DIA routes on Windows to KIT-052. The alternative is to qualify the rollup on Windows first, which is a ticket of its own. |
| D4 | **Scope: the state owner's paths on Windows** (plan step 6). Without it only layer 1 runs on Windows, and no lane. | Put step 6 in this ticket, bounded to 5 hours and one correction pass. Otherwise move it and the Windows lane to KIT-052. |

Stored results: none change on Linux with systemd. Under the plain launcher the receipts carry `NA` measurements
and `measurement_valid = FALSE`; no scientific state differs by launcher.
