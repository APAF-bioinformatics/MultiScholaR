# KIT-134: the distributable, the parts without the read token (2026-11-17)

Scope: revision 2.0 of the ticket with the user's decisions of revision 2.1 (D1 to D7), the parts that need neither the read token nor KIT-131: P2 without `fastduckdb` and `fastlimpadpc`, P3's packaging side, P4, and D5. Branch `kit-134`, cut from `performance-v2` at `f7feeae8`; the spike's pieces were brought in as new files under `packaging/` (the spike branch was not merged). Workflow `.github/workflows/distributable.yml`, triggered by pushes to `kit-134` whose commit message names a stage (`[dist:snapshot]`, `[dist:library]`, `[dist:shell]`, `[dist:installer]`, `[dist:container]`). Nothing under `R/`, `tests/`, `DESCRIPTION`, `NAMESPACE` or `NEWS.md` changed; `tools/profiling/` untouched; `.Rbuildignore` gained `^packaging$`.

**What the builds do not contain.** No run used the token, so neither internal package is in any bundle. The DIA, DIA_LIMPA and enrichment lanes cannot run on what these runs produced (F1); nothing here works around that. Every artifact of these runs is token-free, so artifacts were uploaded (revision 1.2's rule binds token builds only).

## Decision record (current)

- **Snapshot date: 2026-07-30** (Posit Package Manager; `packaging/multischolar.lock.json`, `Snapshot`). Bioconductor's Windows binaries come from the same snapshot.
- **MOFA2 at the snapshot is 1.22.0**, which pins **mofapy2 0.7.3** (KIT-130 had 1.22.1 and 0.7.5); `packaging/python/requirements.txt` follows the bundled MOFA2, and the build fails if they disagree.
- **Electron 45.3.3, electron-builder 26.18.0**, pinned with `package-lock.json` (33.2.0 with 25.1.8 stays the fallback; not needed).
- **The R entry is `packaging/electron/entry.R`**: it calls `MultiScholaR:::runPackagedApplication` when the installed package has it, and until then starts the application in process with a watchdog on the shell, saying so on the `MSR-LAUNCHER` line. Every lifecycle run therefore reports `in-process` today; the expectation becomes `process` (runners) and `systemd` (laptop) with the package side.
- **The start-up probe is in the shell** (`main.js`, `probeLauncher`); the package takes its result from `MSR_LAUNCHER_PROBE` (decision U1 below).

## Runs

| # | Run | Commit | Stage | Result |
|---|---|---|---|---|
| 1 | [43940407142](https://github.com/APAF-bioinformatics/MultiScholaR/actions/runs/43940407142) | `85cbee59` | snapshot (P2 step 1) | passed: date chosen, Windows Bioconductor probe |
| 2 | [43941357316](https://github.com/APAF-bioinformatics/MultiScholaR/actions/runs/43941357316) | `3cce1a83` | library (P2) | failed on both platforms at the MOFA2 pin check, after every package had installed (375 snapshot binaries and 4 source builds on each platform) |
| 3 | [43942816664](https://github.com/APAF-bioinformatics/MultiScholaR/actions/runs/43942816664) | `21dafb4d` | shell (P2 correction pass, P3) | Linux: everything passed (library check, lifecycle on the 22.04 build runner and on 24.04 from `install.sh`). Windows: the library check failed on R's `library/translations` folder (not a package: no namespace); 383 packages and MultiScholaR had loaded; the P3 steps did not run |
| 4 | [43945108934](https://github.com/APAF-bioinformatics/MultiScholaR/actions/runs/43945108934) | `82ae9126` | installer (P3 Windows, P4 both) | RUN4 |

RUNS-EXTRA

### Failures, kept

| Run | Where | What failed | Cause | Class | After |
|---|---|---|---|---|---|
| 2 | install-library.R, both platforms | `requirements.txt disagrees with MOFA2's pins: mofapy2, , , ...` | Two things: the snapshot's MOFA2 (1.22.0) pins mofapy2 0.7.3, not the 0.7.5 written from KIT-130 (a real difference the check exists to catch); and the check lost the dependency names in `vapply` (a defect of the check) | build input plus instrumentation | requirements to 0.7.3; names kept (P2's correction pass) |
| 3 | check-library.R, Windows | `FAILED translations: package 'translations' does not have a namespace` | R for Windows ships `library/translations` with a DESCRIPTION; `installed.packages()` lists it. The other 383 packages and MultiScholaR loaded, from the moved folder, network cut, every namespace from the bundle | instrumentation | the check, the version table and the notices skip folders without a NAMESPACE. This is a second fix to P2's instrumentation after the correction pass; it changes no build input and no threshold, and is recorded here for the user's judgement |

## P1 (laptop, read only)

`installed.packages()` over the lockfile's three library paths against the lockfile, at 2026-11-17 06:2x: **414 of 414 rows the same** (409 at the locked version, the 5 never-installed Suggests still absent). The lockfile was not regenerated. The build's copy is `packaging/multischolar.lock.json`: P1's lockfile plus a `Snapshot` block (date, rule, run, Windows Bioconductor source) and a filled `Python` block; the evidence copy beside this file is unchanged. P1 step 2 (the internal clones against the pins) was not done: it concerns the token build.

| File | sha256 (at the tip) |
|---|---|
| `packaging/multischolar.lock.json` | `df12be0ed077db697f3e32f0ce6138baf79c3a24cc70f702d39d95f2ec89546e` |
| `packaging/python/requirements.txt` | `010700e537f2bcb7a7d5b82797f01be6d5dd505fc6f45081723f2b668717cd86` |

## P2: the library

### The snapshot date (run 1; `p2/snapshot-*`)

Every calendar date from 2026-04-29 (Bioconductor 3.23's release) to 2026-11-16 was scored against the lockfile's 394 repository packages (CRAN 348, Bioconductor 46; 160 in the hard closure) from the snapshot's source indexes: 202 dates, 199 with indexes. `duckdb` was 1.5.5 from 2026-07-27 to 2026-09-29 (65 eligible dates). Best: 2026-07-27 to 2026-07-30 tie at **152 of 160 hard-closure and 365 of 394 packages at the lock**; the latest, **2026-07-30**, was chosen. At that date the hard-closure packages not at the lock are `arrow` 24.0.0 to 25.0.0, `future` 1.70.0 to 1.75.0, `gh`, `httr2`, `Matrix` 1.7-5 to 1.7-6, `plotly`, `RcppParallel` 5.1.11-2 to 6.1.1, `zip` (all newer). Windows: every CRAN package and 45 of 46 Bioconductor packages have a Windows binary at the snapshot's version (GO.db is source only, as expected); a Bioconductor zip is served, so Windows takes Bioconductor from the snapshot, not from bioconductor.org.

### The build (runs 3 and 4; `p2/linux`, `p2/windows`)

Order as revision 2.0: the snapshot's packages (Linux: 375 jammy binaries, nothing compiled; Windows: 375 binaries), then GO.db, RUVIIIC (`cran/RUVIIIC@c534151f`) and the Glimma fork (`APAF-bioinformatics/GlimmaV2@743af094`) from source, then MultiScholaR from the branch's commit. The closure is DESCRIPTION's fields less the two internal packages and the four never-installed Suggests, plus MOFA2 (D7), less `packaging/r/library-exclude.txt`, recursively: 378 packages. The system libraries for jammy came from Posit Package Manager's requirements API (31 apt packages, `p2/linux/apt-packages.txt`). Library install took about 4 minutes on each platform.

**The difference table** (`version-table.tsv`; identical on both platforms but for R's own folders): 438 rows (lockfile 414). In both: 372; **same version 346, different 26** (8 in the hard closure, 16 compiled; all snapshot binaries). Of the 26, 24 are newer than the lock and 2 older: `limpa` 1.4.2 to **1.4.0** (Bioconductor's version at the snapshot; the laptop's 1.4.2 came from the proof library) and `nanoarrow` 0.9.0 to **0.8.0-1**. In the lockfile, not bundled: 37 (the excluded tools and their own dependencies, 35; the two internal packages). Bundled, not in the lockfile: MOFA2 and the 16 packages it adds (basilisk, rhdf5, HDF5Array, uwot, ...), MultiScholaR itself, and R's own folders. `duckdb` is 1.5.5 on both. No hard-closure package is missing.

**One declared minimum is not met:** DESCRIPTION's Suggests has `nanoarrow (>= 0.9.0)`; the snapshot has 0.8.0-1. Its only use is the direct Arrow writer (`R/utils_workflow_saved_table_codec.R:206`, `fast_arrow` backend), which the token-free build cannot select; the snapshot rule did not look at declared minimums. Decision U4.

**Load check** (`check-library.R`; every package in a child R of 10, from the moved folder, network cut): Linux, in a bare `ubuntu:22.04` at `/srv/any where/msr app`, read-only, network `none`: 383 of 383 load, every loaded namespace from the bundle, MultiScholaR attaches (111 s). Windows, from `...\re located\msr moved`, reduced environment, outbound blocked for the bundle's programs: run 3, 383 loaded and only `translations` refused; run 4: RUN4WINCHECK.

Sizes (MiB): Linux bundle 2,301 (R 1,410 of which library 1,383; Python 606; system libraries 131, 117 of them; pandoc 157), 48,570 files. Windows bundle 2,172 (R 1,465, library 1,351; Python 484; pandoc 223), 50,294 files. Largest R packages: BH (headers only, 156 / 123), RUVIIIC (96 on Linux, unstripped), GO.db 68, Rhdf5lib 61, arrow, duckdb, V8, magick (67 on Windows).

Notes for P5: the bare container has no time-zone database (`Failed to locate timezone database` when MultiScholaR attaches; most desktop systems have `tzdata`); `libtbbmalloc.so` is not found by `ldd` for RcppParallel (it loads; TBB's allocator is optional).

## P3: the shell (packaging side)

`packaging/electron/main.js` (from the spike): single-instance lock; the probe; a free port on 127.0.0.1; the bundled `Rscript --vanilla entry.R <port> <run folder> <shell pid>`; a waiting page, then the application at the origin only (navigation and new windows to anything else refused; `https` links to the system browser; downloads from the origin allowed); `contextIsolation`, no Node in the page, the renderer sandboxed; the window title and, in process, the non-blocking notice of revision 2.0; Help, About with the version stamp and the probe's reason; closing stops the R tree (Linux: the process group, `SIGKILL` after 5 s; Windows: `taskkill /T /F`); logs in Electron's `logs` folder (one shell log, one R log per start); run folders under `%LOCALAPPDATA%\MultiScholaR\run` or `$XDG_STATE_HOME/multischolar/run`, five kept.

The environment, one definition per platform: Linux `packaging/linux/Rscript` (unsets `R_HOME`, `R_LIBS`, `R_ENVIRON*`, `R_PROFILE*`, `PYTHONPATH`, `PYTHONHOME`, `RETICULATE_PYTHON_ENV`, `CONDA_PREFIX`, `VIRTUAL_ENV`; sets `RHOME`, `R_LIBS_USER`/`R_LIBS_SITE` to folders that do not exist, `LD_LIBRARY_PATH` to the bundle's libraries, `RETICULATE_PYTHON`, `PYTHONNOUSERSITE=1`, `RSTUDIO_PANDOC`, `PATH=/usr/bin:/bin`, `MSR_BUNDLE`, and `C.UTF-8` when the locale is not UTF-8); Windows `main.js`, `rEnvironment` (the same library and Python settings, `RTOOLS45_HOME` to a folder that does not exist, `MULTISCHOLAR_DISABLE_TOOLS=bash,dd,sha256sum,sync,tar,cp`, `PATH` System32 and Windows only). `--no-sandbox` is passed only by the Linux launcher that `install.sh` writes, and only where `/proc/sys/kernel/apparmor_restrict_unprivileged_userns` is 1. Pandoc 3.11 is bundled.

**The probe** reads the host facts of revision 2.0, 2.0.5 step 1 in the shell (under 50 ms) and passes them to R as `MSR_LAUNCHER_PROBE`. Found: the 22.04 runner answers `process` ("kernel 6.8 is older than 6.12"); **the ubuntu-24.04 runner answers `systemd`: kernel 7.0.0-1017-azure, `systemctl --user` answers, cgroup v2, the memory controller delegated.** F9's premise (no runner can run the systemd branch) no longer holds for `ubuntu-24.04`; the laptop check may not be the only systemd check available (decision U5).

**The lifecycle test** (`lifecycle-test.js`, from `smoke.js`): no bundle process before; network closed (Linux: a network namespace with only loopback, checked by a failing connection; Windows: outbound firewall rules on every program of the application); start; the window's live Shiny session renders MultiScholaR's home tab (session connected, the omic tiles, the server-rendered log panel); the `MSR-LAUNCHER` line; HTTP 200; a second start exits 0, the first window receives it, no second R; the window closed, exit 0, no process of the bundle left, the port refuses; started again and the shell killed hard, R gone within 10 s.

| Where | Result | First tab after launch | R gone after the shell was killed |
|---|---|---|---|
| ubuntu-22.04 build runner, unpacked (run 3) | **passed**, every check | 3.5 s | 1.3 s |
| ubuntu-24.04, installed with `install.sh`, `--no-sandbox` (restriction on) (run 3) | **passed**, every check | 5.5 s | 1.3 s |
| RUN4P3 | | | |

Not done: the laptop's systemd run (F9; needs the package side), and the probe's trial unit and fallbacks (package side).

## P4: the installers

RUN4P4

## D5: Windows Server Core

RUN4D5

## Deviations from revision 2.0's file table

- `packaging/electron/entry.R` (new): the R entry; replaces the `-e` expression so the shell works before and after the package side exists.
- `packaging/verify/check-library.R` and `packaging/verify/windows-clean-run.ps1` (new): the moved-folder library check and the Windows reduced environment for it (the spike's checks).
- `packaging/electron/build/icon.png` is not committed: the workflow copies `inst/shiny/www/MultiScholaR.png` at build time.
- `packaging/r/notices.R` also writes `VERSION.json`.
- The Linux archive's launcher script is written by `install.sh` (no separate file).
- Not written (P5 and the token build): `packaging/verify/compare-results.R`, `run-lanes.sh`, `run-lanes.ps1`; the e2e helper change.

## The package side of P3, to make on top of KIT-131 (not made here)

Nothing under `R/`, `tests/`, `DESCRIPTION`, `NAMESPACE` or `NEWS.md` was changed on this branch. These are the changes, against KIT-131's branch (`kit-131` at `57c8fc04`), once it is merged:

1. **`R/run_app.R`: `runPackagedApplication(arguments)`** (internal, `@noRd`, no export; `NAMESPACE` and `Collate` unchanged). `arguments` = `c(port, run folder, shell pid)`, as `packaging/electron/entry.R` passes them (the shell already calls it through `entry.R` when it exists, so the shell needs no change).
   - Reads the shell's probe from `MSR_LAUNCHER_PROBE` (JSON: `launcher` `"systemd"` or `"process"`, `reason`, `facts`); see decision U1 for whether the probe instead moves into R as `workflowLauncherProbe()`.
   - systemd: `options(multischolar.workflow_launcher = "systemd")`; `launchWorkflowApplication(root = <run folder>, mode = "gui", port = <port>)`; on `workflow_execution_unsupported`, the next step. After the launch returns, stays alive: every 2 s checks the shell (`ps::ps_is_running(ps::ps_handle(pid))`) and the controller unit (`workflowUnitState(config, config$controller_unit)`); when either is gone, stops the slice (`workflowSystemd(c("stop", config$slice))`) and returns.
   - process: `options(multischolar.workflow_launcher = "process")`; a trial unit through KIT-131's `workflowProcessStart("multischolar-trial-<id>", file.path(R.home("bin"), "Rscript"), c("--vanilla", "-e", "quit(status = 0)"), character(), <run folder>, <log>, 30)` then `workflowProcessState` until it has exited 0 and `workflowProcessStop`; then `launchWorkflowApplication(<run folder>, mode = "gui", port = <port>)`, which returns when the application stops (KIT-131's `launchWorkflowProcessApplication`).
   - in process, only when the trial fails: `shiny::runApp(run_app(), host = "127.0.0.1", port = port, launch.browser = FALSE)`.
   - Before starting, prints exactly one line `MSR-LAUNCHER {"launcher": "systemd"|"process"|"in-process", "reason": "..."}` to standard output (and `logger`), as `entry.R` does today.
   - The shell watchdog for the process and in-process branches: `later::later` every 2 s, `quit(save = "no")` when the shell's pid is gone (as `entry.R` does today). `later` is not in `Imports` (it comes with shiny): `R CMD check` would warn on `later::` — see decision U2.
   - Never calls `loadDependencies()` (F12).
2. **`R/utils_workflow_execution.R`, `workflowLaunchService` only, systemd branch** (after KIT-131's early return for the plain launcher, before the `Rscript` arguments are appended): when `nzchar(Sys.getenv("MSR_BUNDLE"))`, append `paste0("--setenv=", name, "=", Sys.getenv(name))` for each of `RHOME`, `LD_LIBRARY_PATH`, `R_LIBS_USER`, `R_LIBS_SITE`, `RETICULATE_PYTHON`, `PYTHONNOUSERSITE`, `RSTUDIO_PANDOC`, `MSR_BUNDLE` that is set (F3). Without `MSR_BUNDLE` the argument vector is identical to today's. The bundle's launcher also sets `PATH=/usr/bin:/bin` and, on a non-UTF-8 locale, `LANG`/`LC_CTYPE=C.UTF-8`; whether units get those too is decision U3.
3. **`tests/testthat/test-run-packaged-application.R`** (new): the probe reading and the order (systemd raising `workflow_execution_unsupported` falls to process; a failing trial falls to in process; the `MSR-LAUNCHER` line each time), with `launchWorkflowApplication`, `workflowProcessStart` and `shiny::runApp` mocked; `workflowLaunchService`'s `systemd-run` arguments captured with `workflowSystemd` mocked, with and without `MSR_BUNDLE` (identical without).
4. **`NEWS.md`**: the packaged entry point, the launcher order and the unit environment.
5. Then on the packaging side: `lifecycle-test.js --expect-launcher process` in the workflow (it is `in-process` today), and the laptop's one lifecycle run of the Linux archive for the systemd branch (F9).
6. Laptop checks by the test rule for 1 to 4: the targeted test files, the dependents grep, `R CMD check`, install.

## Decisions for the user

| # | Decision | What the runs show | Recommendation |
|---|---|---|---|
| U1 | Where the start-up probe lives | The host facts are read in the shell (`main.js`, under 50 ms) and handed to R; revision 2.0 put them in R (`workflowLauncherProbe()`, testable with testthat) | Keep the facts in the shell (packaging knowledge, no package API); the package reads `MSR_LAUNCHER_PROBE` and owns only the trial and the fallback. The other way costs one more R function and its tests |
| U2 | `later` for the package's watchdog | `later::` in package code needs `later` in `Imports` (`R CMD check` warns otherwise); it is installed already as shiny's dependency. A `DESCRIPTION` change is outside revision 2.0's file table | Add `later` to `Imports` with the package side |
| U3 | What the systemd units get beyond F3's list | The bundle's R runs with `PATH=/usr/bin:/bin` and, where the locale is not UTF-8, `C.UTF-8`; units get the user manager's `PATH` and locale | Pass `PATH` too (a unit should not see `/usr/local/bin` tools the application does not); the locale only if set by the launcher |
| U4 | `nanoarrow` below DESCRIPTION's minimum; `limpa` older than the lock | The snapshot has `nanoarrow` 0.8.0-1 where Suggests asks `>= 0.9.0`, and `limpa` 1.4.0 where the laptop ran 1.4.2. The snapshot rule does not look at declared minimums. `nanoarrow` serves only the direct Arrow writer of the token build; `limpa` is DIA_LIMPA's (unavailable here) | Build both at the locked versions from source in the token build (as RUVIIIC is), and add a declared-minimum check to `install-library.R`; or move the date and accept the table's other changes |
| U5 | The systemd branch on a runner | The `ubuntu-24.04` runner meets every requirement (kernel 7.0, user manager, cgroup v2, memory delegated): F9 does not hold there | Run the systemd lifecycle on `ubuntu-24.04` in the workflow once the package side exists, beside the laptop's one check |
| U6 | The written offer of the sources | `THIRD-PARTY-NOTICES.txt` offers the complete sources for three years, by an issue on the repository or by writing to the maintainer named in DESCRIPTION (Ignatius Pang). The sources archive itself is not built yet (each component's source address is listed) | Confirm the wording and the contact; build the archive as a release asset with KIT-132 |
| U7 | Size | Linux archive 1,008 MiB, 2,577 MiB installed in 48,650 files; Windows installer 766 MiB, 50,371 files. Largest avoidable parts: BH headers (156 / 123 MiB, needed only to compile), pandoc (157 / 223 MiB), Python packages' tests | Trim BH and the Python `tests` folders after the lanes pass (P5), not before |
| U8 | Windows `PATH` | KIT-130's table gave R only System32 and the Windows folder; the application needs Windows PowerShell's folder (shinyFiles lists drives with it where WMIC is gone) | Keep PowerShell's folder on R's `PATH` (the change in run 5) |
| U9 | The budget rule | P2's instrumentation needed a second fix (`translations`) after its correction pass; it changed no input and no threshold | Accept it as recorded, or count P2 as stopped at run 3 on Windows |


## Budget

BUDGET
