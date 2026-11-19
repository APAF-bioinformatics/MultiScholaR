# KIT-134: the distributable, the parts without the read token (2026-11-17)

Scope: revision 2.0 of the ticket with the user's decisions of revision 2.1 (D1 to D7), the parts that need neither the read token nor KIT-131: P2 without `fastduckdb` and `fastlimpadpc`, P3's packaging side, P4, and D5. Branch `kit-134`, cut from `performance-v2` at `f7feeae8`; the spike's pieces were brought in as new files under `packaging/` (the spike branch was not merged). Workflow `.github/workflows/distributable.yml`, triggered by pushes to `kit-134` whose commit message names a stage (`[dist:snapshot]`, `[dist:library]`, `[dist:shell]`, `[dist:installer]`, `[dist:container]`). Nothing under `R/`, `tests/`, `DESCRIPTION`, `NAMESPACE` or `NEWS.md` changed; `tools/profiling/` untouched; `.Rbuildignore` gained `^packaging$`.

**What the builds do not contain.** (Runs 1 to 6; runs 7 and 8 used the token: see "Token build".) No run used the token, so neither internal package is in any bundle. The DIA, DIA_LIMPA and enrichment lanes cannot run on what these runs produced (F1); nothing here works around that. Every artifact of these runs is token-free, so artifacts were uploaded (revision 1.2's rule binds token builds only).

## Decision record (current)

- **Snapshot date: 2026-07-30** (Posit Package Manager; `packaging/multischolar.lock.json`, `Snapshot`). Bioconductor's Windows binaries come from the same snapshot.
- **MOFA2 at the snapshot is 1.22.0**, which pins **mofapy2 0.7.3** (KIT-130 had 1.22.1 and 0.7.5); `packaging/python/requirements.txt` follows the bundled MOFA2, and the build fails if they disagree.
- **Electron 45.3.3, electron-builder 26.18.0**, pinned with `package-lock.json` (33.2.0 with 25.1.8 stays the fallback; not needed).
- **The R entry is `packaging/electron/entry.R`**: it calls `MultiScholaR:::runPackagedApplication` when the installed package has it, and until then starts the application in process with a watchdog on the shell, saying so on the `MSR-LAUNCHER` line. Every lifecycle run therefore reports `in-process` today; the expectation becomes `process` (runners) and `systemd` (laptop) with the package side.
- **The start-up probe is in the shell** (`main.js`, `probeLauncher`); the package takes its result from `MSR_LAUNCHER_PROBE` (decision U1 below).
- **W1 closed (the user, 2026-11-18: "no extra browser we gud"): no browser is bundled.** Only a Chromium snapshot build is clearly redistributable (BSD plus third-party notices); Chrome for Testing's `chrome-headless-shell` is not clearly so. And the Summary report starts no browser: its templates attach `webshot2` but print no interactive plot. Shown by run 44146578653: from the installed application, on Linux and on Windows, the lipidomics and metabolomics Summary reports render to `.docx` with no Chrome-type browser to be found. Section "W1" below.
- **The shell's probe asks the user D-Bus** before it chooses the systemd worker; on GitHub's `ubuntu-24.04` runner it gives the plain launcher with "the user D-Bus does not answer" (run 44146578653; W1, part 4).
- **Without `fastduckdb` the lipidomics and metabolomics intensity filter is refused** ("optional fastduckdb reader is unavailable"), on the narrow GUI-entry fixtures too, on both platforms (run 44146578653; W1, part 4). A token-free bundle therefore cannot pass those lanes.

- **`fastduckdb` is a declared import (the user, 2026-11-18: "fastduckdb has to be included. all the architecture depends on it.").** DESCRIPTION moves it from Suggests to Imports (`fac2bd06`); a lane run without it is an invalid configuration, not a finding. Every build now needs the read token; under revision 1.2's rule a token build makes no shell, archive or installer, so **no installer can be built by the workflow until the two repositories are public or the user changes that rule** (section "2026-11-18").
- **Linux lanes on the bundle: the seven browser e2e lanes pass with the laptop's counts** against the bundle's installed package (run 44195739758). The metabolomics GUI-entry lanes differ from the laptop's oracles from the normalisation step on, identically on the Linux and the Windows runner: not classified yet.
- **The Linux bundle does not start on the laptop** (CachyOS): its `LD_LIBRARY_PATH` puts Ubuntu 22.04's `libreadline.so.8` in front of the host's, and the host's `/bin/bash` (which R's front-end script runs under) then fails. A defect of the bundle on every distribution whose `bash` links readline dynamically.
- **Windows: the proteomics protein route (TMT, LFQ) fails at QC** ("required Q-value figures are incomplete"); the lipidomics, metabolomics and transcriptomics e2e lanes pass with the laptop's counts.

## Runs

| # | Run | Commit | Stage | Result |
|---|---|---|---|---|
| 1 | [43940407142](https://github.com/APAF-bioinformatics/MultiScholaR/actions/runs/43940407142) | `85cbee59` | snapshot (P2 step 1) | passed: date chosen, Windows Bioconductor probe |
| 2 | [43941357316](https://github.com/APAF-bioinformatics/MultiScholaR/actions/runs/43941357316) | `3cce1a83` | library (P2) | failed on both platforms at the MOFA2 pin check, after every package had installed (375 snapshot binaries and 4 source builds on each platform) |
| 3 | [43942816664](https://github.com/APAF-bioinformatics/MultiScholaR/actions/runs/43942816664) | `21dafb4d` | shell (P2 correction pass, P3) | Linux: everything passed (library check, lifecycle on the 22.04 build runner and on 24.04 from `install.sh`). Windows: the library check failed on R's `library/translations` folder (not a package: no namespace); 383 packages and MultiScholaR had loaded; the P3 steps did not run |
| 4 | [43945108934](https://github.com/APAF-bioinformatics/MultiScholaR/actions/runs/43945108934) | `82ae9126` | installer (P3 Windows, P4 both) | Linux: everything passed, P4 included (archive installed with `install.sh` on fresh `ubuntu-22.04` and `ubuntu-24.04`, lifecycle, uninstall). Windows: library check passed (383 of 383), NSIS installer built, path lengths within the limit; the lifecycle failed: the page was served but no Shiny session connected; the fresh install was skipped |
| 5 | [43949090047](https://github.com/APAF-bioinformatics/MultiScholaR/actions/runs/43949090047) | `d156dd43` | installer, Windows only (P3's correction pass) | passed: library check, lifecycle on the build runner, NSIS installed silently per user on a fresh `windows-2025`, lifecycle from the installed application |
| 6 | [43953990999](https://github.com/APAF-bioinformatics/MultiScholaR/actions/runs/43953990999) | `a5c1282d` | container, Windows only (D5) | failed: in Server Core 382 of 383 packages load (`rgl` refuses: no OpenGL there); the application answers HTTP 200 in 2 s |

### Failures, kept

| Run | Where | What failed | Cause | Class | After |
|---|---|---|---|---|---|
| 2 | install-library.R, both platforms | `requirements.txt disagrees with MOFA2's pins: mofapy2, , , ...` | Two things: the snapshot's MOFA2 (1.22.0) pins mofapy2 0.7.3, not the 0.7.5 written from KIT-130 (a real difference the check exists to catch); and the check lost the dependency names in `vapply` (a defect of the check) | build input plus instrumentation | requirements to 0.7.3; names kept (P2's correction pass) |
| 3 | check-library.R, Windows | `FAILED translations: package 'translations' does not have a namespace` | R for Windows ships `library/translations` with a DESCRIPTION; `installed.packages()` lists it. The other 383 packages and MultiScholaR loaded, from the moved folder, network cut, every namespace from the bundle | instrumentation | the check, the version table and the notices skip folders without a NAMESPACE. This is a second fix to P2's instrumentation after the correction pass; it changes no build input and no threshold, and is recorded here for the user's judgement |

Run 4, Windows lifecycle: `first-tab` reported `connected: false, rendered: false` after 120 s with the home page served (HTTP 200) and R running. Cause: `app_server` calls `shinyFiles::getVolumes()()` when a session starts; on Windows without WMIC (Server 2025; Windows 11 24H2 drops it too) shinyFiles runs `powershell`, which is not on a `PATH` of System32 and the Windows folder only, so the server function fails and the session closes. Class: a real defect of the launch environment (it would hit users), not of the test. P3's correction pass (run 5) put Windows PowerShell's folder on R's `PATH` (`main.js`, `windows-clean-run.ps1`, D5's script); the session then connected in 5.6 s. The lifecycle test now also prints the application's R and shell logs when it fails.

## P1 (laptop, read only)

`installed.packages()` over the lockfile's three library paths against the lockfile, on 2026-11-17 before the first run: **414 of 414 rows the same** (409 at the locked version, the 5 never-installed Suggests still absent). The lockfile was not regenerated. The build's copy is `packaging/multischolar.lock.json`: P1's lockfile plus a `Snapshot` block (date, rule, run, Windows Bioconductor source) and a filled `Python` block; the evidence copy beside this file is unchanged. P1 step 2 (the internal clones against the pins) was not done: it concerns the token build.

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

**Load check** (`check-library.R`; every package in a child R of 10, from the moved folder, network cut): Linux, in a bare `ubuntu:22.04` at `/srv/any where/msr app`, read-only, network `none`: 383 of 383 load, every loaded namespace from the bundle, MultiScholaR attaches (111 s). Windows, from `...\re located\msr moved`, reduced environment, outbound blocked for the bundle's programs: run 3, 383 loaded and only `translations` refused; runs 4 and 5: 383 of 383 load, every namespace from the bundle, MultiScholaR attaches (109 to 129 s).

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
| ubuntu-22.04 build runner, unpacked (run 4) | **passed** | 5.0 s | 0.8 s |
| fresh ubuntu-22.04, installed with `install.sh` (the runner's own R at `/usr/bin/R` not used) (run 4) | **passed** | 8.1 s | 0.8 s |
| fresh ubuntu-24.04, installed with `install.sh`, `--no-sandbox` (run 4) | **passed** | 6.6 s | 1.0 s |
| windows-2025 build runner, unpacked (run 4) | **failed**: no session (above) | | 0.9 s |
| windows-2025 build runner, unpacked (run 5) | **passed** | 5.6 s | 2.1 s |
| fresh windows-2025, installed with the NSIS installer (the runner's R 4.6.2 on the PATH not used) (run 5) | **passed** | 6.0 s | 1.8 s |

Not done: the laptop's systemd run (F9; needs the package side), and the probe's trial unit and fallbacks (package side).

## P4: the installers

| | Linux | Windows |
|---|---|---|
| Installer | `MultiScholaR-0.5.0-linux-x86_64.tar.gz`, 1,057,363,126 bytes (1,008 MiB; gzip, 102 s to pack) | `MultiScholaR-Setup-0.5.0.exe`, 803,600,389 bytes (766 MiB; NSIS one-click, per user; electron-builder 694 s) |
| Installed | 2,577 MiB in 48,650 files under `~/.local/share/multischolar/0.5.0` | 2,530 MiB in 50,372 files under `%LOCALAPPDATA%\Programs\multischolar` |
| Install time on a fresh runner | 22 s (ubuntu-22.04), 26 s (ubuntu-24.04), unpack and `install.sh` | 398 s, silent (`/S`) |
| Administrator rights | none (runner's user; `install.sh` refuses root) | none (per-user NSIS) |
| Shortcut | menu entry `~/.local/share/applications/multischolar.desktop` with the icon | Start-menu `MultiScholaR.lnk` (and a desktop shortcut, electron-builder's default) |
| Longest path | | 143 characters under the application folder (limit 170); 200 installed with the user `runneradmin` |
| Removal | `uninstall.sh` removed the application, the menu entry and the icon | the NSIS uninstaller (not exercised) |
| Licences and stamp | `THIRD-PARTY-NOTICES.txt` (532 components: R, every bundled R package, 117 Ubuntu libraries, CPython and its wheels, pandoc, Electron), `sources-list.tsv`, `VERSION.json`; beside the application and in `resources/` | the same (415 components: no system libraries), in `resources/` |

`VERSION.json`: release 0.5.0, the lockfile's sha256 (`df12be0e...`), the commit, the snapshot date, the build date, the runner image, R and Electron versions, and a line saying the internal packages are not included. Help, About shows it with the probe's reason. The sources archive itself (every source tarball) is not built: the notices give each component's source address and the written offer (U6).

## D5: Windows Server Core

Run 6 (`d5/servercore-run6.log`): the Windows bundle (built in the same job) mounted read-only into `mcr.microsoft.com/windows/servercore:ltsc2025` (pulled in 296 s; Docker 29.7.2 is on the `windows-2025` image), `--network none`, the shell's R environment.

- **The application answers HTTP 200** within about 2 s, started through `entry.R` as the shell starts it.
- R 4.6.1, home and library inside the bundle, network off, every loaded namespace from the bundle, MultiScholaR attaches: passed.
- **382 of 383 packages load; `rgl` fails** in `.onLoad` ("Loading failed"). Server Core has no OpenGL (`opengl32.dll`), which `rgl`'s DLL links against; `RGL_USE_NULL` does not help a DLL that cannot load. `rgl` is in the closure only through Suggests (no package code calls it). Class: environment (E), a valid result about Server Core, not a defect of the bundle on a desktop Windows, where OpenGL is part of the system; on the build runner and the fresh runner it loaded. Not corrected: nothing in the bundle should change for Server Core, and D5's second run is left unspent.

What D5 shows beyond the runners: no Visual C++ runtime, no R, no Rtools and no Python were needed from the host for R, the library and the application (the container has none of them). The Electron window itself cannot be shown in Server Core (no desktop); the manual install on a Windows 11 PC remains with KIT-052 (D5, revision 2.1).

## Token build (2026-11-17; P2 step 4, P5 step 3, U4)

Bounded item: the token part of P2 and its verification. Stage `[dist:internal]` (`distributable.yml`): the library as before plus `fastduckdb` and `fastlimpadpc` at DESCRIPTION's Remotes commits, downloaded with the secret `DEPS_READ_TOKEN` (as `GITHUB_PAT`, in the install step's `env` only), named explicitly (package to repository, subfolder and commit; the downloaded DESCRIPTION's `Package` is checked; the pins are checked against DESCRIPTION's Remotes), installed after `duckdb` and before MultiScholaR. Linux: required. Windows: attempted (D4), a failure recorded and the build continued. U4 (both platforms, every build): `nanoarrow` 0.9.0 and `limpa` 1.4.2 from source (CRAN and Bioconductor 3.23, sha256 in the log), and `check-library.R` checks every version MultiScholaR's DESCRIPTION declares. New: `packaging/verify/check-internal.R`, `packaging/verify/run-lanes.R`.

**No artifact, no cache.** Stage internal never builds the shell, archive or installer; every library, archive, installer and evidence upload also requires `internal != true`; the only upload is a summary folder of named text files checked before upload (under 20 MiB, no NUL byte; no lane records). The internal sources are unpacked under the runner's temporary folder, never in the output folder. The plan job fails if the workflow names `actions/cache` or a cache input. Read after the runs: run 44015367782 uploaded `linux-internal-summary` (433,045 bytes) and `windows-internal-summary` (634,353); run 44017206786 `windows-internal-summary` (634,605); nothing else. The repository's caches on `kit-134` are all `Ubuntu 26.04.1 LTS-R version 4.6.2` keys written by Module CI (`setup-r-dependencies`; `module-ci.yml` uses no secret), none by these runs. Text evidence: `token/`.

| # | Run | Commit | Platforms | Result |
|---|---|---|---|---|
| 7 | [44015367782](https://github.com/APAF-bioinformatics/MultiScholaR/actions/runs/44015367782) | `55d2863a` | both | Linux: library and internal packages installed, library check passed, internal check **failed** (DIA manifest, rollup), lanes **failed** (rollup refused). Windows: `fastlimpadpc` installed and loads; `fastduckdb` failed to build (LNK1104, path of 260 characters: instrumentation); library check passed |
| 8 | [44017206786](https://github.com/APAF-bioinformatics/MultiScholaR/actions/runs/44017206786) | `f9e52d67` | Windows (correction pass) | passed: both internal packages built and load |
| 9 | [44041155687](https://github.com/APAF-bioinformatics/MultiScholaR/actions/runs/44041155687) | `e9b7f8f4` | both (after T1: `fastduckdb` 0.2.0 at `2694b972`, `fastlimpadpc` at `87a90b89`) | Windows passed (both build and load). Linux: library, internal packages, DIA manifest, rollup and writer all pass; both lanes **failed at their last-but-one step**: the summary report cannot render (`webshot2` excluded from the bundle) |

| 10 | [44044598918](https://github.com/APAF-bioinformatics/MultiScholaR/actions/runs/44044598918) | `e5b23654` (`performance-v2` 9f7e371d merged, K1; the audit repairs and `webshot2`, `chromote` shipped) | both (dispatch; this item's correction pass) | **passed**: both platforms build and load both internal packages; Linux: the DIA and DIA_LIMPA lanes pass and every step is identical to the laptop's oracle |

### The correction pass (run 10): the result

- **K1.** `performance-v2` (`9f7e371d`, DESCRIPTION at `3ec13e95`) merged into `kit-134` (`e5b23654`, a merge commit, no conflict; nothing under `packaging/` or `.github/` changed by it). The Remotes check passed: DESCRIPTION and `install-library.R` both name `fastduckdb@e724f90b` and `fastlimpadpc@87a90b89`.
- **Install.** Linux: `fastduckdb` 0.2.0 (193 s), `fastlimpadpc` 0.1.0 (22 s). Windows: 0.2.0 (406 s), 0.1.0 (43 s). Only the `INTERNAL` result lines were printed (F1).
- **Versions.** 378 packages in both lockfile and bundle, 354 same, 24 different (`webshot2` 0.1.2 and `chromote` 0.5.1 now bundled at the lock). Declared minimums: 6 met, none unmet, `testthat` not bundled. Both platforms.
- **Load, moved library, network closed.** Linux bare container and Windows reduced environment: **389 of 389** packages load, every namespace from the bundle, MultiScholaR attaches; both internal packages load at their commits; the reader round trip is identical. Linux: DIA manifest `supported`, `native_query` present, the writer recorded; the rollup and the direct writer admit the build; DIA_LIMPA resolves `fastlimpadpc`.
- **The lanes** (Linux, in process against the moved bundle, network closed (checked), the laptop's locale, `fastlimpadpc` the engine, `testthat` 3.3.2 in a separate library shadowing nothing):

| Lane | Oracle steps | Pass | Fail | Error | Steps differing from the oracle | Discrete | Digest | Seconds |
|---|---|---|---|---|---|---|---|---|
| prot_dia | 29 | 44 | 0 | 0 | 0 | 0 | 0 | 302 |
| prot_dia_limpa | 24 | 38 | 0 | 0 | 0 | 0 | 0 | 211 |

**Classification under revision 1.2's rule** (relative 1e-6, absolute 1e-9; discrete exact): every step of both lanes is **identical** to its laptop-frozen oracle (heads, histories, tab statuses, notifications, refusal checks, outputs, artifact digests and state digests), so every compared value is in class **Agree**, with nothing to classify. The state digests are of values rounded to 10 significant digits, so identity there is agreement to 10 significant digits, stricter than the rule; no digest differs, so **no laptop values are needed** for these two lanes. The rule was not changed. Expected differences from `-march=native` did not appear at this precision on these fixtures.

**The report step and Chrome.** Both lanes rendered the summary report (`DIANN_report_proteomics_prot_dia.docx`, `DIANN_limpa_report_proteomics_prot_dia_limpa.docx`) with `webshot2` from the bundle. *Corrected 2026-11-18:* this paragraph said that `webshot2` "drove the runner's preinstalled Chrome". That was an inference, not an observation: each report rendered within one second and no log line names a browser; the templates attach `webshot2` and print no interactive plot (section "W1").

**What the token run uploaded.** Two artifacts only: `linux-internal-summary` (22,516 bytes) and `windows-internal-summary` (16,201 bytes); neither holds `install-library.log`; scanned: no compiler output (`Compiling`, `.rs:`, `warning:`) and no token text (`ghp_`, `github_pat_`, `Bearer`). Caches on `kit-134`: Module CI's R 4.6.2 keys only, none from this run.

### After T1 (run 9)

The user's T1 (2026-11-17): `fast_duckdb` main fast-forwarded to `2694b972` (0.2.0), `fast_limpa_dpc` main `87a90b89` (declares Imports `limpa`, `methods`); `performance-v2` repins DESCRIPTION (`daf1e341`, then `3ec13e95` with `fastduckdb` at `e724f90b`, the same package code with Cargo online when nothing is vendored). Run 9 built `2694b972` with the crates fetched first into the package's own Cargo home, as its installer at that commit does (removed again for `e724f90b`, which needs no fetch).

- **Install and versions, both platforms.** `fastduckdb` 0.2.0 (Linux 192 s, Windows 293 s), `fastlimpadpc` 0.1.0 (21 s, 41 s). Version table: 374 in both, 350 same, **24 different** (`fastduckdb` now at the lock). Declared minimums: 6 met, none unmet, `testthat` not bundled.
- **Load, moved library, network closed.** Linux bare container: 385 of 385 load from the bundle; both internal packages load at their commits; the reader round trip is identical. Windows, reduced environment: the same (385 of 385; both load; round trip identical).
- **Linux, the DIA side** (`check-internal.R`): manifest `supported`, `native_query` present (`duckdb-c-scope-v1`), the `fastduckdb` writer recorded; **the precursor rollup and the direct writer admit the build**; DIA_LIMPA resolves `fastlimpadpc`. P2 step 4's acceptance (the manifest with `native_query`) is met.
- **Windows answer.** `fastduckdb` 0.2.0 builds under Rtools45 with Rust's GNU target and loads; its reader works, which is what the lipidomics and metabolomics steps need for tables wider than 64 columns. Its `Makevars.win.in` builds no `fastduckdb_query` library and its configure marks only Linux supported, so the DIA query runtime does not exist on Windows (D4 stands).

**The lanes** (Linux; network closed, the laptop's locale, `fastlimpadpc` the engine):

| Lane | Steps in the oracle | Pass | Error | Where it stopped | Cause | Class |
|---|---|---|---|---|---|---|
| prot_dia | 29 | 1 (the enrichment record check) | 1 | `summary_report` (after import, QC, rollup, normalisation, DA, enrichment, summary save and copy ran) | `Failed to render R Markdown report: there is no package called 'webshot2'`; no `.docx`, so the step's check `rendered[[1L]]` is out of bounds | **D, a defect of the bundle**: `packaging/r/library-exclude.txt` excluded `webshot2` (and `chromote`) as development tools, but the report templates attach `webshot2` (`library(webshot2)`; corrected 2026-11-18: not "knitr needs it to render the report's htmlwidgets", the templates print none, section "W1") (pandoc was found: the launcher sets `RSTUDIO_PANDOC`). Repaired after the run: both now ship; verified by run 10 |
| prot_dia_limpa | 24 | 0 | 1 | `summary_report` | the same | the same |

Because the GUI-entry helper writes a lane's record only when the lane ends, a lane that stops writes none: **no step of either lane was compared with its oracle**, so revision 1.2's rule (relative 1e-6, absolute 1e-9; discrete exact) has still been applied to nothing, and nothing is classified beyond the defect above. The rule is unchanged. When the steps are compared, each `digest` difference (`lane-differences.tsv`) needs the laptop's values for the same lane, step and key before it can be classed: the oracles hold sha256 digests of values rounded to 10 significant digits.

Also seen: `fastlimpadpc`'s Linux `Makevars.in` at `87a90b89` compiles with `-C target-cpu=native`: a library built on a runner may use instructions an older user CPU lacks ("illegal instruction"). For `fast_limpa_dpc`, before release. Artifacts: the summary artifacts of runs 7 and 8 and run 9's Windows summary are no longer listed (removed outside this work); run 9's `linux-internal-summary` remains (it holds `install-library.log`, F1 below). Caches on `kit-134`: Module CI's R 4.6.2 keys only.

**Install and version check.** Linux: `fastduckdb` 0.1.0 (172 s) and `fastlimpadpc` 0.1.0 (19 s), rustc 1.100.0. Windows (run 8): `fastduckdb` 0.1.0 (327 s) and `fastlimpadpc` 0.1.0 (34 s), rustc 1.99.1 with the `x86_64-pc-windows-gnu` target and Rtools45. `fastlimpadpc` declares no Imports although its NAMESPACE imports `limpa` and `methods` (a defect of the internal package, reported by the CI lane); this build installs `limpa` before it, so it was not hit here. Version table: 374 packages in both, 349 same, **25 different** (was 26: `nanoarrow` and `limpa` now at the lock; **`fastduckdb` is new: locked 0.2.0, bundled 0.1.0**). Declared minimums: 7 declared; on both platforms 6 met, 0 not met, `testthat` not bundled (excluded tool). The "compiled no" for packages installed from a folder (the internal packages, the Glimma fork) is the table reading `NeedsCompilation`, which only `R CMD build` writes; not a result.

**Load from the moved library, network closed.** Linux, bare `ubuntu:22.04`, `--network none`, `/srv/any where/msr app`: 385 of 385 packages load, all from the bundle; MultiScholaR attaches; both internal packages load from the bundle at their commits; `fastduckdb`'s reader reads back a 100-column Parquet file written by duckdb identically. Windows (run 8), reduced environment, outbound blocked: 385 of 385 load; both internal packages load; the reader round trip is identical.

**Why the DIA lanes cannot pass at the Remotes pin (P1 step 2, closed by this).** DESCRIPTION's Remotes pins `fast_duckdb/fastduckdb@d7851150`, the tip of `main` (2026-09-16), which is `fastduckdb` **0.1.0**: it exports only the Parquet reader and writer (`read_parquet_fast`, `write_parquet_fast`, ...). The laptop's verified library (the lockfile, the oracles, the simulations) has **0.2.0**, installed from a local tree on 2026-10-06; its NAMESPACE is byte-identical to `feat/scoped-native-query@2694b972` ("Document and bootstrap clean installation of fastduckdb 0.2.0", unmerged). The package code needs 0.2.0: the DIA rollup (`R/utils_workflow_dia_relational_rollup.R:6-49`) requires `open_query_fast` and the other query exports, `query-manifest.rds` and `fastduckdb_query`; the direct writer (`R/utils_workflow_saved_table_codec.R:161-189`) requires `parquet_writer_capabilities`. On the runner, `check-internal.R` read: manifest `supported` TRUE, `native_query` absent, the `fastduckdb` writer absent, the rollup refused ("direct writer requires the qualified incremental fastduckdb build"); the DIA_LIMPA route resolves `fastlimpadpc`. `fastlimpadpc`'s pin `d4f9cc1d` is the tip of its `main` and matches the laptop's 0.1.0 (built 2026-10-21).

**The lanes** (`run-lanes.R`: `test_file(desc = <lane>, load_package = "installed")` with the bundle's `Rscript`, `testthat` 3.3.2 and the 7 packages it adds in a separate library after the bundle's, shadowing nothing; network namespace with loopback only, checked closed; locale `en_GB.UTF-8` with `en_AU.UTF-8` numeric, time, monetary, paper and measurement, as the laptop; default DPC engine `fastlimpadpc`):

| Lane | Pass | Fail | Error | Record | Cause | Class |
|---|---|---|---|---|---|---|
| prot_dia | 1 | 0 | 1 | none (the lane stopped) | the precursor rollup refused (`fastduckdb` 0.1.0), then every later step had no state; the lane ended at `file(con, "r")`: cannot open the connection | build input (the Remotes pin), not a result difference |
| prot_dia_limpa | 0 | 0 | 1 | none | the same | the same |

No step reached the comparison, so revision 1.2's rule (relative 1e-6, absolute 1e-9; discrete exact) was applied to nothing: there is no numeric or discrete difference to classify. The rule was not changed. For when the lanes run: the GUI-entry oracles hold sha256 digests of values rounded to 10 significant digits, not values; `run-lanes.R` lists every differing key as `discrete` (compared exactly) or `digest` (the tolerance cannot be applied to a digest; such a difference needs the values from both machines before it can be classed).

**Windows: does `fastduckdb` build, and do the lipidomics and metabolomics steps need it?** At the pin, yes: 0.1.0 builds under Rtools45 with Rust's GNU target, loads, and its reader works (run 8; run 7's failure was this script's 260-character build path, fixed by unpacking in short folders, no input changed). Whether 0.2.0 (the query runtime the DIA routes need) builds on Windows is not known: not built here (D4). The lipidomics and metabolomics partitioned steps, the defaults since KIT-089 (`R/utils_workflow_lipid_declarations.R:37`, `R/utils_workflow_metab_declarations.R:37`), stream rows through `streamWorkflowTableRows` (`R/utils_workflow_omics_declarations.R:95`, `:770`, `:788`, `:843`), which for more than 64 columns opens a `fastduckdb` read scope unconditionally (`R/utils_workflow_saved_table_codec.R:2093-2101`). So on every platform, a lipidomics or metabolomics table wider than 64 columns (about 60 or more samples) fails at the intensity filter with "optional fastduckdb reader is unavailable" unless `fastduckdb` is installed; the reader 0.1.0 has (`read_parquet_fast`) is all those steps use. The GUI-entry and e2e fixtures are narrower, which is why the token-free lanes passed. *Corrected 2026-11-18:* no token-free lane run stood behind that sentence; in run 44146578653 the intensity filter of both narrow GUI-entry lanes was refused with this message on a bundle without `fastduckdb`, on Linux and on Windows (section "W1", part 4). Not worked around.

### Independent audit of the token build (2026-11-17)

Repaired on the branch (packaging only; takes effect from the next stage-internal run, which counts as this item's correction pass):

| # | Finding | Repair |
|---|---|---|
| F1 (high) | The token summaries uploaded `install-library.log`, and the job logs printed it: it carried compiler messages quoting `fast_limpa_dpc` source lines and `fast_duckdb`'s crate layout. The text check guarded the type, not the content. Affected: the summary artifacts and job logs of runs 44015367782, 44017206786 and 44041155687 (not deleted: the user's decision) | The internal packages' `R CMD INSTALL` output goes to `$RUNNER_TEMP/internal-build-logs/<package>.log`, never uploaded or printed; only the `INTERNAL` result line is printed. `install-library.log` is no longer in either token summary. The copies under `token/` hold no build output (checked) |
| F4 | `GITHUB_PAT` stayed in the environment of every configure script, Makevars and crate build script the install started | `install-library.R` reads it once into a variable and unsets it before anything else runs; the variable is used only in the internal downloads' headers |
| F5 | The script's pins could differ from DESCRIPTION's Remotes and only print `DIFFERENT` | `fastduckdb` pinned to `e724f90b` (DESCRIPTION on `performance-v2` since `3ec13e95`; 0.2.0's package code, Cargo offline only with vendored crates, so the crate prefetch is removed); any difference from the checkout's Remotes stops a token build. **`kit-134`'s DESCRIPTION still names the old commits, so the next token run on this branch stops at that check until the branch carries `performance-v2`'s DESCRIPTION** (this branch does not edit DESCRIPTION) |
| F8 | The Linux text check (`grep -qI`) passes a file whose NUL bytes come after its first block | Every byte is scanned (the byte count with NULs removed must equal the size), as the Windows check does |

Recorded for the package side of P3 and for P5 (nothing under `R/` is changed here):

- **F3.** systemd workers started from the bundle get neither `RHOME` nor `LD_LIBRARY_PATH` (F3 of revision 2.0, still open: the `MSR_BUNDLE` change to `workflowLaunchService` is part of the package side), and the package's launcher choice does not read `MSR_LAUNCHER_PROBE`, the probe the shell already hands it (U1).
- **F6.** Downloads in the workflow are not all pinned or checksummed: Rtools45, the Rust toolchain (`rustup ... stable`), R's installers, python-build-standalone and pandoc are fetched by version or channel without a recorded sha256 check; DESCRIPTION's Remotes give `cran/RUVIIIC` and `APAF-bioinformatics/GlimmaV2` no commit (the build pins them in `install-library.R`; D3's Remotes fix is still the recorded separate change). For P5 and the release build.
- **F9.** The in-process watchdog (`later` every 2 s) cannot run while a long step blocks R's event loop, so R can outlive a killed shell for the length of that step. For the package side (the plain launcher moves heavy steps out of the application process, which bounds this).

## W1: a browser for the Summary report (2026-11-18)

Bounded item (the user's W1, 2026-11-18: "ship headless chrome ... edge is very unsafe"; one browser for both platforms): bundle one open-source Chromium headless build in both installers and hand it to chromote in the packaged launcher; the licence first, before any build. **The item stopped before its build, on the two findings below. No workflow run was started, nothing was downloaded into a build and nothing is bundled.** The decision is the user's (end of this section).

### 1. The licence finding

Read on 2026-11-18 from the archives themselves (listed over HTTP range requests; the one licence file saved and read), the projects' own pages and the maintainers' statements. Not legal advice.

| Candidate | What it is | Licence evidence | May it go into the installers? |
|---|---|---|---|
| `chrome-headless-shell` from Chrome for Testing, 157.0.8103.78 (the stable channel on 2026-11-17) | Google's build of Chromium's old headless mode, one per Chrome release. Linux zip 127,965,683 bytes (276 MiB unpacked), Windows zip 126,047,001 bytes (282 MiB) | For: the archive's only licence file is `LICENSE.headless_shell` (2.3 MB): Chromium's BSD 3-clause licence, then 313 third-party notices (FFmpeg under the LGPL among them); no EULA and no terms file. Against: the only maintainer statement on Chrome for Testing binaries (GoogleChromeLabs/chrome-for-testing, issue 21, June 2023, a Google maintainer): "Chrome for Testing is a flavor of Google Chrome, so https://www.google.com/chrome/terms/ applies." Those terms bring in the Google Terms of Service, whose software licence is "personal" and "non-assignable" and says "You may not copy, modify, distribute, sell, or lease any part of our services or software" (an open-source licence overrides where it says so). Whether the statement covers the headless shell, which joined Chrome for Testing with Chrome 120, after it, is asked in the open issue 226 and has no answer from Google | **Not clearly.** Two readings, not settled by Google |
| A Chromium snapshot build (the Chromium project's continuous builds, bucket `chromium-browser-snapshots`) | The whole Chromium browser, built from the open-source tree at one commit. Linux also has a `headless-shell.zip`, Windows has none, so one browser for both platforms means the full `chrome` | The source licence at the commit (`LICENSE`, read from chromium.googlesource.com): BSD 3-clause, "Redistribution and use in source and binary forms, with or without modification, are permitted" on reproducing the notice; the third-party parts under their own open-source licences. chromium.org's download page attaches no terms to the builds. The archives hold **no** licence or notice file: whoever passes them on adds Chromium's `LICENSE` and the third-party notices of that build (its `chrome://credits` page) and names the source commit (FFmpeg and other LGPL parts are linked in) | **Yes**, under BSD plus the third-party notices, which we would have to produce and ship |

By the item's rule, only Chromium proper is clearly redistributable. Its pins, should it be bundled: one commit on both platforms, `109167ec53e2daebc8ad3df8cd499571c44b9c29` (`refs/heads/main@{#1717304}`, built 2026-10-13, 45 commits before Chrome 157's branch point 1717349). The sha256 was computed here by streaming each archive; the md5 of the same stream equals the bucket's.

| Platform | URL | Bytes | sha256 | Unpacked |
|---|---|---|---|---|
| Linux | `https://commondatastorage.googleapis.com/chromium-browser-snapshots/Linux_x64/1717304/chrome-linux.zip` | 251,185,343 (240 MiB) | `a15d33424d043c42f6dd27455b1d2e639d383838138e71ae1a9b5060a8aa33a5` | 677 MiB in 483 files (`chrome` alone 508 MiB) |
| Windows | `https://commondatastorage.googleapis.com/chromium-browser-snapshots/Win_x64/1717304/chrome-win.zip` | 362,399,196 (346 MiB) | `04632b0d01fc20fc2d28a8403c1415b15d73c50777b5df87deba55abfba35106` | 820 MiB in 258 files (352 MiB of it `interactive_ui_tests.exe`, a test program that would be left out) |

For the record, the Chrome for Testing archives' sha256: Linux `f90329720fb61e5a76da9b6b509b6ad12a484c9c5df0359403e935177c6d153a`, Windows `2138c1583d25d6f3abec5269943b2f430207d87a539895aa9e90cc4697ec0e1f`.

What the snapshot would cost. The "+150 to 200 MB each" given to the user was an estimate made before the archives were read. Estimated from the archives, nothing built: the Linux archive grows by about 240 MiB (1,008 to about 1,250 MiB; installed +677 MiB); the Windows installer by roughly 200 MiB once the test program is dropped (766 to about 970 MiB; installed about +470 MiB). And three things beyond size:

- A snapshot is the development trunk at one commit: "built from arbitrary revisions that don't necessarily map to user-facing Chrome releases", "do not auto-update" (chromium.org). It never receives the stable channel's security fixes, so the bundle would carry a browser that is unpatched from the day it ships.
- On Ubuntu 23.10 and later a per-user install cannot give Chromium its sandbox (AppArmor restricts unprivileged user namespaces and the setuid helper needs root: F10's restriction, for a second program). It would run with `--no-sandbox` there, from a wrapper script, because chromote adds that flag by itself only in CI, in a container or for a user without an account entry (`chromote::default_chrome_args()`).
- Its notices have to be generated from the pinned build and added to `THIRD-PARTY-NOTICES.txt` and the sources list; a systemd unit needs `CHROMOTE_CHROME` passed to it (one name in `workflowLaunchService`'s `MSR_BUNDLE` list, with its test).

### 2. The premise: the Summary report starts no browser

The item rests on "the Summary report renders its interactive plots into the Word file through webshot2 and chromote". Read against the code and the runs, it does not.

- **The templates.** All eight report templates under `inst/reports/` (lipidomics, metabolomics, the four proteomics ones, transcriptomics, the shared analysis record) were parsed and every function they call was listed, chunks and inline code: 125 distinct functions. Figures enter as PNG files through `knitr::include_graphics`, tables through `knitr::kable`, the rest is text. No call creates or prints an htmlwidget and none takes a screenshot. `flextable`, `gt`, `plotly` and `webshot2` appear only in the `library()` lines of six templates (the only called names that plotly exports are the dplyr verbs it re-exports). Nothing under `R/` or `inst/` calls `webshot2`, `chromote`, `gtsave` or an image export of a widget; `htmlwidgets::saveWidget` writes HTML and needs no browser.
- **knitr.** knitr 1.51 calls `webshot2::webshot` from `html_screenshot` only, which runs when an htmlwidget, a Shiny application or an embedded URL is printed into a document that is not HTML. Attaching `webshot2` starts nothing.
- **Run 9 (44041155687)** failed one second into the render with "there is no package called 'webshot2'": the templates' `library(webshot2)`, not a screenshot.
- **Run 10 (44044598918)** rendered `DIANN_report_proteomics_prot_dia.docx` in one second (08:29:40 to 08:29:41) and the DIA_LIMPA report within one second (08:33:15); a browser does not start and take screenshots in that time, and no log line names one. The sentences in this record that said `webshot2` "drove the runner's preinstalled Chrome" were an inference; they are corrected where they stand.
- **A probe on the laptop** (2026-11-18; R 4.6.1, knitr 1.51, rmarkdown 2.31, webshot2 0.1.2, chromote 0.5.1, pandoc 3.11; `CHROMOTE_CHROME` set to a path that does not exist). A document that does what the templates do (the four `library()` lines, `include_graphics`, `kable`, text) renders to `.docx`. The same document with one printed plotly widget fails: "Functions that produce HTML output found in document targeting docx output". So the probe sees a browser need where there is one, and the templates have none.

What the report needs is the **package** `webshot2` (and `chromote`, which it imports), because six templates attach it; both ship since run 10. The comment in `packaging/r/library-exclude.txt` said otherwise and is corrected.

A complete Summary report rendered from the installed application with every browser hidden, on both platforms, was not shown when this was written; run 44146578653 shows it (part 4).

### 3. The decision, for the user

**Decided (the user, 2026-11-18: "no extra browser we gud"): option 1.** The table is kept as it was put.

| Option | What it means |
|---|---|
| **1. Bundle no browser (recommended).** | The application as it is needs none; nothing is added to either installer, no unpatched browser ships, no `--no-sandbox` second program on Ubuntu 24.04. The install notes need no browser requirement. If a later report prints an interactive plot into Word, the question returns then, and the Electron shell already carries a Chromium (with its `LICENSES.chromium.html`) that could take the picture without a second browser. |
| 2. Bundle the pinned Chromium snapshot anyway. | For a future interactive figure in the Word report. The costs above: +240 MiB and about +200 MiB of installer, a browser without security updates, `--no-sandbox` on Ubuntu 23.10 and later, the notices to generate. One installer run plus its correction pass, as budgeted. |
| 3. Ask Google about `chrome-headless-shell` first. | Half the size (122 and 120 MiB) and a stable-channel build with its licence file inside; usable only if Google says the BSD file is the whole of its terms (issue 226). |

Separate and small, package code outside this item: the six templates could drop their unused `library(flextable)`, `library(gt)`, `library(plotly)` and `library(webshot2)` lines; `webshot2` and `chromote` could then leave the bundle. Not done.

### 4. The two follow-ups and the report check: run 44146578653

The two follow-ups, as made in `6c73b1fb`:

- **The probe asks the user D-Bus** (`packaging/electron/main.js`, `probeLauncher`): after every other requirement holds, `busctl --user call org.freedesktop.systemd1 /org/freedesktop/systemd1 org.freedesktop.DBus.Peer Ping` must answer, else the plain launcher with the reason "the user D-Bus does not answer". `busctl` being installed is no longer enough. Checked on the laptop by running the function's text under plain Node: `systemd` in this login session (58 ms for the whole probe); `process` with that reason when the user manager answers and the bus address leads nowhere, which is the hosted runner's state.
- **ubuntu-24.04 expects the plain launcher with the bus reason** (`distributable.yml`, `lifecycle-test.js --expect-reason`): both Linux install legs run with the network closed again and expect `process`; the 24.04 leg also expects the reason to name "D-Bus". The `dbus-user-session` attempt of run 44059705566 is removed. The systemd worker from the bundle is tested on the laptop, in a real login session, after the KIT-124 studies.

**The run.** [44146578653](https://github.com/APAF-bioinformatics/MultiScholaR/actions/runs/44146578653), commit `d958b479`, stage installer, both platforms, 2026-11-17 18:19 to 18:52 UTC; no token, no cache. Before it, `performance-v2` (`624153a4`: KIT-061, the Module CI changes, B-01's record) was merged into `kit-134` as `1d7d1ec9` (a merge commit; one conflict, the comment in `library-exclude.txt`, kit-134's corrected wording kept). The packaging inputs match: `multischolar.lock.json`, `requirements.txt` and `install-library.R` are identical on the two branches, and `install-library.R` names DESCRIPTION's Remotes commits (`fastduckdb@e724f90b`, `fastlimpadpc@87a90b89`). **Every job passed; the correction pass was not used.** No step that should have run was skipped (read from the run's step list; the report check runs on the fresh `ubuntu-22.04` and `windows-2025` install runners, not on 24.04, by design). The log lines quoted below are kept in `w1/run-44146578653-lines.txt`.

**Lifecycle, launcher and reason** (every check of the lifecycle test passed on each):

| Where | Launcher | The probe's reason | First tab | R gone after the shell was killed |
|---|---|---|---|---|
| ubuntu-22.04 build runner, unpacked, network closed | `process` | kernel 6.8 is older than 6.12 | 5.0 s | 0.3 s |
| fresh ubuntu-22.04, `install.sh`, network closed (the runner's `/usr/bin/R` not used) | `process` | kernel 6.8 is older than 6.12 | 5.0 s | 0.5 s |
| fresh ubuntu-24.04, `install.sh`, `--no-sandbox`, network closed | `process` | **the user D-Bus does not answer** | 5.5 s | 1.0 s |
| windows-2025 build runner, unpacked, network cut | `process` | win32: no systemd | 6.6 s | 1.2 s |
| fresh windows-2025, NSIS per user, network cut (the runner's R 4.6.2 on the PATH not used) | `process` | win32: no systemd | 5.5 s | 2.2 s |

On ubuntu-24.04 the probe's facts were: kernel 7.0.0-1017-azure, `systemctl`, `systemd-run` and `busctl` present, the user manager answers, cgroup v2, controllers `cpu memory pids`, `userBus` false; 43 ms. Every earlier requirement held, so the ping is what decided, and `--expect-reason D-Bus` passed. The other side (a user bus that answers, then `systemd`) is shown only by the laptop's run of the function under Node; the bundle's systemd lifecycle on the laptop is still to do.

**The report check** (new in this run: `run-lanes.R` with `MSR_REPORT_CHECK=1`). The lipidomics and metabolomics GUI-entry lanes (`lipid_canonical`, `metab_lc`: `mod_lipidomics_server` and `mod_metabolomics_server` under `testServer`, import to the Summary tab's Generate Report) ran in process with the installed application's own R and library, after the lifecycle, with `testthat` 3.3.2 and the 7 packages it adds in a separate library. `run-lanes.R` refuses to start if a browser can still be found, prints every `rmarkdown::render` with what it wrote, and reads each lane's `summary_report` step.

| | Linux: fresh ubuntu-22.04, installed with `install.sh` | Windows: fresh windows-2025, installed with the NSIS installer |
|---|---|---|
| Browsers on the runner | `google-chrome`, `google-chrome-stable`, `chromium`, `chromium-browser`, `microsoft-edge`, `microsoft-edge-stable` in `/usr/bin`: removed for the step | Chrome and Edge installed under Program Files; none on the runner's PATH |
| What R had | PATH `/usr/bin:/bin`, none of nine browser names on it; `CHROMOTE_CHROME=/nonexistent/no-browser/chrome` (missing), which is what `chromote::find_chrome()` returns | PATH System32, the Windows folder and Windows PowerShell's (and R's two Rtools entries, which do not exist), none on it; `CHROMOTE_CHROME=C:\no-browser\chrome.exe` (missing), which `find_chrome()` returns, so the registry is not asked |
| Network | closed (namespace with loopback only; `cran.r-project.org` not reached) | cut for the application's 39 programs (not reached) |
| `lipid_canonical` | `lipidomics_report_lipidomics_lipid_canonical.docx`, 280,956 bytes, a zip container; "Report generated successfully!"; `analysis_record.docx` 16,851 bytes; lane 51 s | the same report, 278,962 bytes, a zip container; the same message; `analysis_record.docx` 16,834 bytes; lane 105 s |
| `metab_lc` | `metabolomics_report_metabolomics_metab_lc.docx`, 232,280 bytes; the message; `analysis_record.docx` 17,188 bytes; lane 45 s | 234,993 bytes; the message; `analysis_record.docx` 17,171 bytes; lane 110 s |
| Verdict | `KIT-134 REPORT CHECK PASSED` (both lanes) | `KIT-134 REPORT CHECK PASSED` (both lanes) |

**So: a Summary report renders to `.docx` from the installed application with every Chrome-type browser hidden, on Linux and on Windows, for the lipidomics and the metabolomics lane.** What that rests on: the application's own report path, driven through the module's inputs in the installed R, not clicked in the Electron window; each `.docx` checked for its size, the zip signature and the application's success message, not opened (the lane's folder is removed when the lane ends, so no document was kept). Proteomics and transcriptomics reports were not run here (the DIA reports rendered in run 10, with the runner's browsers present).

**The lanes against their oracles (not what this check decides; P5's).** On both platforms both lanes ended with 19 expectations passed and 16 failed, no error: 16 of the 20 steps differ from the laptop's oracle (Linux: `lipid_canonical` 280 discrete and 119 digest differences, `metab_lc` 226 and 79; Windows: 290 and 119, 236 and 79). The log shows one cause at the first differing step, `qc_intensity`: "Error applying lipid intensity filter: optional fastduckdb reader is unavailable" (the metabolite filter the same). The state `lipid_intensity_filtered` is never made, and the later steps' histories lack it (the first five differences printed per lane are exactly that; the difference tables were not uploaded and the rest was not read). The lanes went on without it (duplicates, finalise, normalisation, DA, summary), which is the state the reports above were rendered from.

**Finding for P5 and P6: without `fastduckdb` the lipidomics and metabolomics intensity filter fails on the narrow six-sample fixtures too, on Linux and on Windows.** This record said that only tables wider than 64 columns need it and that "the token-free lanes passed" on narrow fixtures; no token-free lane run stands behind that sentence, and this run contradicts it (corrected where it stands). The message is `newWorkflowTableReadScope`'s (`R/utils_workflow_saved_table_codec.R:509-511`); which caller asks for the `fastduckdb` reader on a narrow table was not traced. Nothing was worked around and nothing under `R/` was changed. Consequence: a bundle without `fastduckdb` is not a working lipidomics or metabolomics application, and those two lanes can pass only on a bundle that carries it (the token, or the repositories public), on Windows as well (`fastduckdb` 0.2.0 builds and loads there, runs 9 and 10).

**Also from this run.** The GUI-entry lanes run on Windows with the installed R through `run-lanes.R` (first time; all 20 steps of each lane ran). Bundle: 387 packages load from the moved folder on both platforms, network off; version table 376 in both, 352 same, 24 different; declared minimums 6 met, none unmet, `testthat` and `fastduckdb` not bundled. Linux archive 1,060,501,852 bytes (1,011 MiB), installed 2,588 MiB in 49,412 files in 25 s; Windows installer 805,855,646 bytes (769 MiB, sha256 `7e530712fc0e7a0f7e43446ec6b50523a058a78e46ad819e1c1759c9977275b1`), installed 2,544 MiB in 51,133 files in 372 s, longest path 143 characters under the application folder. Artifacts (token-free, so allowed): `linux-archive`, `windows-installer`, the two evidence folders.

## Deviations from revision 2.0's file table

- `packaging/electron/entry.R` (new): the R entry; replaces the `-e` expression so the shell works before and after the package side exists.
- `packaging/verify/check-library.R` and `packaging/verify/windows-clean-run.ps1` (new): the moved-folder library check and the Windows reduced environment for it (the spike's checks).
- `packaging/electron/build/icon.png` is not committed: the workflow copies `inst/shiny/www/MultiScholaR.png` at build time.
- `packaging/r/notices.R` also writes `VERSION.json`.
- The Linux archive's launcher script is written by `install.sh` (no separate file).
- Not written (P5): `run-lanes.ps1`; the e2e helper change. `packaging/verify/compare-results.R` came with `performance-v2` (B-01, `kit-134/b01-ci-run.md`). The token build added `packaging/verify/run-lanes.R` (an R script, in place of `run-lanes.sh`) and `packaging/verify/check-internal.R`; the secret is `DEPS_READ_TOKEN`, not revision 2.0's `MSR_INTERNAL_READ_TOKEN`.
- Run 44146578653 added: `packaging/verify/test-tools.R` (testthat beside the bundle, either platform; the token stage keeps its own inline copy, not touched because that stage did not run); `run-lanes.R` finds a lane in any `test-gui-entry-*.R` and has the report mode (`MSR_REPORT_CHECK=1`); `windows-clean-run.ps1` takes `-NoBundleArgument`, which is how `run-lanes.R` runs on Windows. The two install jobs check out the whole tree (the lane helpers read it) and have longer time limits (60 and 90 minutes); the Windows lifecycle step no longer depends on the test-tool step before it.

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

Revision 2.0's budgets are workflow runs per platform per part, one correction pass per part. A run counts against the part of its stage; the earlier parts it rebuilds are repeats, not attempts.

| Part | Linux runs (budget) | Windows runs (budget) | Correction pass |
|---|---|---|---|
| P2 | 3: runs 1, 2, 3 (6) | 2: runs 2, 3 (4) | used (run 3); Windows also needed the `translations` instrumentation fix after it (U9) |
| P3 | 1: run 3 (6) | 2: runs 4, 5 (4); run 3 counted nothing here, its P3 steps did not run | Linux none; Windows used (run 5) |
| P4 | 1: run 4 (4) | 2: runs 4, 5 (4) | none needed beyond P3's |
| D5 | | 1: run 6 (2) | not used |
| Token build (P2 step 4, P5 step 3) | 1: run 7 (1 plus a correction pass) | 2: runs 7, 8 | Windows used (run 8, the build path); Linux not used: its failure is the Remotes pin, which needs a decision, not a correction |
| Token build after T1 (fresh budget) | 2: runs 9, 10 | 2: runs 9, 10 | used (run 10: the audit repairs, K1, `webshot2` and `chromote`); passed |
| P3 package side (2026-11-17) | 2: runs 44055833237, 44059705566 | 2: the same runs | used (44059705566: `dbus-user-session` on ubuntu-24.04, no user bus all the same) |
| W1 and the two follow-ups (2026-11-18; one installer run plus a correction pass) | 1: run 44146578653 | 1: the same run | not used: every job passed (W1, part 4) |

Six workflow runs in all, about 3 hours of runner time across the jobs. Each push to `kit-134` also starts the repository's existing `Module CI` workflow (not part of this record). Nothing heavy ran on the laptop: the P1 table (seconds of R), and three light local checks of the new scripts (two snapshot dates, the plan resolution, the version table against R's own library; no package installed; `package-lock.json` written with `npm install --package-lock-only --ignore-scripts`, no `node_modules`).

## What is left for KIT-134

- W1 and its installer run: done (run 44146578653; W1, part 4). No browser is bundled; the probe's user-bus check and the ubuntu-24.04 expectation hold; the lipidomics and metabolomics Summary reports render to `.docx` from the installed application with every browser hidden, on both platforms.
- The package side of P3: done (`f154988b`, `1f228360`; runs 44055833237 and 44059705566). Left of it: the systemd lifecycle of the bundle on the laptop, in a real login session (hosted runners have no user D-Bus).
- P5 (the lanes and the suite from the installed bundle): the e2e helper's installed surface, the e2e lanes and the suite through a runner, the comparison under the rule (`compare-results.R` is on the branch since the merge; B-01 is closed on `performance-v2`). **The lipidomics and metabolomics lanes need `fastduckdb` in the bundle** (run 44146578653: the intensity filter of the GUI-entry lanes is refused without it; the e2e lanes use the same fixtures and the same filter, by reading, not run), as DIA, DIA_LIMPA and enrichment do (F1): of the seven lanes only transcriptomics and TMT and LFQ are left as candidates for a token-free bundle, and neither has been run against one. So P5's lanes, and P6's on Windows, need a bundle built with the token and verified in the job (no artifact), or the two repositories public. Whether the intensity filter should ask for `fastduckdb` on a narrow table is package code and a decision outside this item.
- The token build: done (run 10). Both internal packages build and load on both platforms; on Linux the DIA and DIA_LIMPA GUI-entry lanes pass against the bundle, identical to the laptop's oracles. Open: W1 (section "W1": the report needs no browser as it stands; the user's decision), the `fastlimpadpc` `target-cpu=native` build flag (its repository), the enrichment e2e lane and the suite (P5).

### Decisions needed (token build)

| # | Decision | What the runs show | Options |
|---|---|---|---|
| T1 | Which `fastduckdb` the bundle names | Decided by the user (2026-11-17): 0.2.0 on `fast_duckdb` main; `performance-v2` DESCRIPTION at `3ec13e95` | done |
| K1 | `kit-134`'s DESCRIPTION | Decided (the coordinator, engineering): `performance-v2` merged into `kit-134` (`e5b23654`) | done; run 10 |
| W1 | Chrome for the report | Repaired on the packaging side (the coordinator's instruction): `webshot2` and `chromote` are no longer excluded, so they ship (`library-exclude.txt`). The user chose (b) on 2026-11-18 ("ship headless chrome ... edge is very unsafe"). *Corrected 2026-11-18:* this row said that `webshot2` "still needs a Chrome" and that the report fails without one; the templates take no screenshot, so (c) is already the case (section "W1") | **Closed** (the user, 2026-11-18: "no extra browser we gud"): nothing is bundled; shown by run 44146578653 |
- P6 on Windows, after KIT-131's step 6: the GUI-entry lanes now run there with the installed R (`run-lanes.R`, run 44146578653); the lanes themselves need `fastduckdb` in the Windows bundle (above); the e2e lanes and the suite are not started.
- The sources archive and the release (KIT-132).

## 2026-11-18: `fastduckdb` in Imports; the lanes on the bundle (P5, P6); the laptop's systemd lifecycle; Module CI; the notices

Five bounded items, each with one run and one correction pass where a run was needed. What was run: one workflow run (44195739758, stage internal, both platforms), one laptop test file, and the lifecycle test on the laptop (one invalid start, one valid). Everything else below was read.

### 1. Module CI on `kit-134` (read from the logs; nothing changed)

Every failing Module CI run of `kit-134` since the `ci-token` merge fails in one job, `module-ci impacted (transcriptomics)`, and the summary job that follows it. Run 44146578640 (`d958b479`): eight expectations, `test-module-ci-trx-da.R:7`, `:28`, `test-module-ci-trx-norm.R:13`, `:30`, `test-module-ci-trx-qc.R:64`, `test-module-ci-trx-summary.R:40`, `:52`, `:65`. These are the transcriptomics digest expectations of `b01-ci-run.md` (rows 3 to 10; row 11, the e2e journey, was not in this run's impact set): the values agree with the laptop under the 1e-6 rule and the stored digests are of the laptop's `-march=native` build. `performance-v2` fails the same job with the same first expectation whenever it is impacted (run 44111160249 at `624153a4`); its three latest runs are green only because documentation commits impact no job. They stay red until KIT-129's harness revision; no test was changed. Run 44195739804 (`fac2bd06`) additionally lost its metabolomics job to a download ("limpa 1.4.2: no source at bioconductor.org", the U4 source build; the other jobs of the same run fetched it): transient, and audit F6 (unpinned downloads without a mirror) again.

### 2. `fastduckdb` from Suggests to Imports (`fac2bd06`)

- `DESCRIPTION`: `fastduckdb (>= 0.2.0)` is under Imports; the floor and the Remotes pin (`e724f90b`) are unchanged. `NAMESPACE` is unchanged and needs no change: every use is `fastduckdb::`. `NEWS.md` has the entry. `roxygen2::document()` was not run.
- `packaging/r/install-library.R` stops an install without the read token before anything is installed ("DESCRIPTION requires fastduckdb (an internal repository)"), instead of failing at MultiScholaR's own install after 380 packages. The lockfile is not edited: its `Declared: Suggests` and `Closure: suggests-only` for `fastduckdb` describe the library as exported on 2026-11-17 and are now stale in those two fields (as its `DescriptionRemotes` already was).
- **Shown by run 44195739758:** with the import declared, MultiScholaR installs and loads from the moved bundle on Linux and on Windows (389 of 389 packages, network off, every namespace from the bundle; `fastduckdb` 0.2.0 at `e724f90b` and `fastlimpadpc` 0.1.0 at `87a90b89` on both; declared minimums met; Linux: DIA manifest, rollup and writer admit the build).
- **Jobs that change.** `distributable.yml`: stages library, shell, installer and container now stop at the library install without the token; only stage internal builds. `module-ci.yml` (through `.github/actions/snapshot-library`): unchanged with the secret; a fork or Dependabot pull request (no secret) now stops at the library install instead of running without the internal packages (audit F2's path is gone). `kit-131-launcher.yml` (branch `kit-131` only): its `fastduckdb=?ignore` can no longer work. None of the three files was edited for this.
- **On the laptop** `pkgload::load_all()` from the default library now stops ("dependency 'fastduckdb' is not available"): `fastduckdb` lives in the proof library, which has to be on `R_LIBS`. `test-e2e-browser-harness.R` with it: 123 passed, 0 failed, 1 skipped (no browser on the laptop).
- **Guards left in `R/`** (not removed; for a later cleanup): `R/utils_workflow_dia_da.R:153` (falls back to `duckdb`), `R/utils_workflow_saved_table_codec.R:165` (the direct writer's admission) and `:509` ("optional fastduckdb reader is unavailable"), and the reader argument check at `R/utils_workflow_saved_entities.R:654-657`. The tests' `skip_if_not_installed("fastduckdb")` (12 files) and DESCRIPTION's SystemRequirements sentence ("Cargo and rustc ... only for the optional fastlimpadpc engine"; `fastduckdb` needs them too) are in the same list.
- **Release consequence (the user's call; not worked around).** While `fast_duckdb` is internal, an installer contains it and so may not be uploaded as an artifact (revision 1.2, IR-B2), and stage internal builds no shell or installer at all. Either the repository goes public, or the user changes the rule, before any installer can be built again, the release candidate included.

### 3. P5 and P6: the lanes on the bundle (run 44195739758, `fac2bd06`, stage internal)

[44195739758](https://github.com/APAF-bioinformatics/MultiScholaR/actions/runs/44195739758), 2026-11-17 23:02 to 23:46 UTC. Budget: one run and one correction pass; **the correction pass is not used**, because nothing that failed is a fault of the run (below). New in the job: `packaging/verify/run-e2e.R` (the browser lanes with `MULTISCHOLAR_E2E_INSTALLED=1`, so the application the browser drives is `library(MultiScholaR)` from the bundle: F2, `tests/testthat/helper-e2e-browser.R`), `shinytest2` 0.5.1 in the separate test-tool library (shadowing nothing), the runner's Chrome. Uploaded: `linux-internal-summary` (38,141 bytes) and `windows-internal-summary` (33,166 bytes), text only, scanned (no compiler output, no token text); no cache written. Evidence: `token/linux-run44195739758/`, `token/windows-run44195739758/`.

**The seven browser e2e lanes, Linux** (moved bundle, network closed, the laptop's locale; laptop counts of KIT-118's exit, `e94b92b1`):

| Lane | Passed | Failed | Errors | Laptop | Seconds |
|---|---|---|---|---|---|
| lipidomics canonical | 88 | 0 | 0 | 88 | 71 |
| metabolomics LC and GC | 190 | 0 | 0 | 190 | 161 |
| transcriptomics screen | 16 | 0 | 0 | 16 | 50 |
| proteomics DIA | 42 | 0 | 0 | 42 | 159 |
| proteomics DIA_LIMPA (`fastlimpadpc`) | 86 | 0 | 0 | 86 | 115 |
| proteomics enrichment | 98 | 0 | 0 | 98 | 342 |
| proteomics TMT and LFQ | 153 | 0 | 0 | 153 | 161 |

What this rests on: the bundle's R, library and installed MultiScholaR in the job that built them (F8: a token build cannot go to a fresh runner), not the installed archive and not the Electron window; the laptop's counts are of an earlier tree (no laptop e2e run was started).

**The four e2e lanes without DIA, Windows** (moved bundle, the shell's reduced environment, network cut):

| Lane | Passed | Failed | Errors | Laptop |
|---|---|---|---|---|
| lipidomics canonical | 88 | 0 | 0 | 88 |
| metabolomics LC and GC | 190 | 0 | 0 | 190 |
| transcriptomics screen | 16 | 0 | 0 | 16 |
| proteomics TMT and LFQ | 136 | 0 | 1 | 153 |

The TMT and LFQ lane stops at `Download from '#proteomics_workflow-session_summary-download_report' failed`.

**The GUI-entry lanes against their oracles** (the stricter comparison: every step's states as digests of values rounded to 10 significant digits):

| Lane | Linux | Windows |
|---|---|---|
| prot_dia | 44 passed, identical (29 steps) | not run (D4) |
| prot_dia_limpa | 38 passed, identical (24 steps) | not run (D4) |
| prot_tmt, prot_lfq, prot_lfq_fragpipe | 37 passed each, identical (21 steps each) | error each: QC refused, "required Q-value figures are incomplete" |
| lipid_canonical | 35 passed, identical (20 steps) | 27 passed, 8 failed: every state and artifact identical; one output, the normalisation log text, is recorded as a digest where the oracle holds the text |
| metab_lc, metab_gc, metab_combined | 25 passed, 10 failed each: `metab_normalized` differs from the `norm` step on | the same 10 steps, **the same state digests as the Linux runner**, plus the log-text difference |
| trx_screen | 22 passed, 7 failed: one artifact, `top_de_genes.tsv`, in 7 steps | error after 3 s: "module session has been destroyed" at the design save |

**Classification under 2.0.4, and what is not classified.**

| Difference | Class | Evidence |
|---|---|---|
| Six Linux lanes identical | Agree | no digest and no discrete difference |
| trx_screen, `top_de_genes.tsv` (Linux) | not classified here; expected Agree | the table is edgeR's `topTags` written as text; `b01-ci-run.md` row 4 measured that table on the same library within 5.1e-9 relative with identical order and calls, on another fixture. The values of this lane were not compared |
| metabolomics, `metab_normalized` and what follows (both runners) | **open, not explained** | the two runners agree with each other to 10 significant digits and differ from the laptop. GC assay: mean sample correlation 1.000 on the runners, 0.972 in the laptop's oracle, which is far outside the rule. The fixtures keep 3 metabolites by 6 samples after filtering and normalise by cyclic loess; the lipidomics lane (8 and 9 lipids) is identical. A loess over 3 points is the likely ill-conditioned step, and the laptop (R built with `-march=native`) is the one machine of three that differs; that is a reading, not a measurement. To classify: the normalised matrix of one metabolomics lane from the laptop and from a runner, compared with `compare-results.R` |
| lipid_canonical and metabolomics, normalisation log text (Windows) | not classified; likely L or E | every state identical; the recorder stores a long text as its digest, so the Windows text is longer or differs (a path, a line ending). The text itself was not read |
| Windows protein route: QC refused | **D, a defect on Windows** | `R/utils_workflow_saved_exports.R:1045`: after the Q-value figure renderer ran, the twelve expected `.png` and `.pdf` files are not all in the staged folder. First at "protein accession cleanup" in all three protein lanes; the e2e TMT and LFQ lane fails at the report for want of those outputs. Package code; not changed |
| Windows trx_screen entry lane | not classified | the browser lane of the same fixture passes 16 of 16 on Windows; the entry lane's design save loses its module session under `testServer` there. Not traced |

So: **Linux meets the lane half of 2.0.9's Linux criterion in the form a token build allows; Windows does not meet P6** (TMT and LFQ fail; D). The suite half is not run: 2.0.4 compares it with a laptop check-mode table of a tree with the same `R/`, `src/` and `tests/`, the last such table is of KIT-123's tree (before KIT-131, KIT-061 and this ticket's package side), and P5 may not start a laptop suite. It needs the next checkpoint suite on the laptop, then one sharded run.

### 4. The systemd lifecycle on the laptop (P3 step 4, F9): failed, with its cause

The Linux archive of run 44146578653 (token-free, `d958b479`; sha256 checked) was downloaded into the session's scratch folder and installed there with `install.sh` under `XDG_DATA_HOME` (nothing under the home folder), with Electron's and the run folders under scratch `XDG_CONFIG_HOME` and `XDG_STATE_HOME`; `lifecycle-test.js --expect-launcher systemd`; removed afterwards with `uninstall.sh`. Kernel 7.3.3, cgroup v2, a login session.

- First start: invalid. The session's own environment carried `ELECTRON_RUN_AS_NODE=1` (the editor's), so the shell ran as plain Node and exited in 0.5 s; the test observed nothing (`p3/laptop-systemd-lifecycle-laptop-invalid-run-as-node.log`). An instrumentation failure, kept.
- The run, with those variables removed (`p3/laptop-systemd-lifecycle-laptop.log`): the shell started (Electron 45.3.3, sandbox on); **the probe chose `systemd`** ("every host requirement holds": kernel 7.3.3, the three tools, the user manager, cgroup v2, controllers `cpuset cpu io memory pids`, the user bus answers; 21 ms); R was started and **exited with status 127 after 18 ms**: `/bin/bash: symbol lookup error: /bin/bash: undefined symbol: rl_trim_arg_from_keyseq`. No unit and no slice was started; nothing was left behind.
- Cause, confirmed by one command outside the test: with `LD_LIBRARY_PATH` at the bundle's `syslib`, the host's `/bin/bash` and `/bin/sh` fail with that message. The bundle's launcher (`packaging/linux/Rscript`) exports `LD_LIBRARY_PATH=<bundle>/syslib` for R and everything R starts; `syslib` holds Ubuntu 22.04's `libreadline.so.8`; R's front-end (`R/bin/R`) is a `#!/bin/bash` script, and this host's `bash` links readline dynamically and needs a newer one. `systemctl`, `busctl`, `systemd-run` and `ps` still run under that path. On Ubuntu and Debian `bash` carries its own readline, which is why every runner passed.
- **Class D: the Linux bundle cannot start on a distribution whose `bash` links the system readline (Arch and its derivatives at least; others not checked), whatever the launcher.** D1's "glibc 2.35 or later" is not a sufficient description of the supported hosts. The same path also shadows the host's `libstdc++`, `libssl`, `libcurl` and `libz` for any host program R starts.
- Not corrected: the repair is a change to how the bundle finds its libraries (for example a run path inside the bundle's own binaries in place of the exported variable, so that host programs see the host's libraries), it needs a new archive, and no archive can be built until the repository question of item 2 is answered. The lifecycle's remaining checks on the laptop (the slice while running, gone after close) are therefore not shown.

### 5. The notices (U6), made and not run

`packaging/r/notices.R`: the three-year written offer and the maintainer's name and address are gone; the notices say that the complete corresponding source of every component is published as one archive beside the installers of the release (`<repository>/releases`) and point to the Issues page. `VERSION.json`'s `internal_packages` line now reports what the bundle holds instead of a fixed "not included". Parsed only: the script runs in the shell stage, which cannot run now (item 2). The sources archive itself is KIT-132's and is not built.

### Budget, 2026-11-18

| Item | Runs | Correction pass |
|---|---|---|
| Module CI diagnosis | none (logs read) | |
| `fastduckdb` in Imports; P5 Linux lanes; P6 Windows lanes | 1 run for both platforms: 44195739758 | not used |
| Laptop systemd lifecycle | 1 (after one invalid start that observed nothing) | not used: the repair needs a new archive |
| Notices | none | |

## In-job build and test (2026-11-18, decision M20); the Linux library path (defect 2)

M20 (the user, 2026-11-18): `fast_duckdb` becomes public only at the public build. Until then nothing built with the read token may leave the job that built it. `distributable.yml` is restructured accordingly, and the Linux bundle finds its libraries by run path. Budget declared in `b835b7f1`: one run (stage full, both platforms) and one correction pass. Both are used.

| Run | Commit | What ran | Result |
|---|---|---|---|
| [44226281203](https://github.com/APAF-bioinformatics/MultiScholaR/actions/runs/44226281203) | `b835b7f1` | stage full: ubuntu-22.04, ubuntu-24.04, windows-2025 | Windows complete (below). Linux, both runners: the library loads with the run paths (389 of 389, bare ubuntu:22.04, no `LD_LIBRARY_PATH`), then the internal check failed: "direct writer requires the qualified incremental fastduckdb build"; the later Linux steps were skipped (a fault of the new workflow's conditions) |
| [44232555102](https://github.com/APAF-bioinformatics/MultiScholaR/actions/runs/44232555102) | `cbc46e8f` | the correction pass, Linux only | every Linux step ran; open: the host check in the Arch Linux and Fedora containers, and the three metabolomics entry lanes (the fixture, defect 4) |

Neither run has an artifact (0 and 0 by the API) and neither wrote a cache (the eight caches on the branch are of 2026-11-17, R 4.6.2, another workflow). The job logs hold the install output of the public packages, the check lines and the tables below; the internal packages' build output stays in the runner's temporary folder as before.

### 1. The workflow

- Stages: `snapshot` (token-free, unchanged), `library`, `full`, `container`. Every build step that installs the library has the token; no other step does.
- Stage `full`, one job per platform: library, checks from the moved folder, shell, Linux archive or NSIS installer, install on the same runner (Linux: `install.sh` under a scratch `XDG_DATA_HOME`; Windows: the silent per-user install into the runner's `%LOCALAPPDATA%\Programs`), the build's other copies and the archive or installer deleted, then from the installed application: the lifecycle test, the GUI-entry lanes, the browser e2e lanes (network closed), and on Linux the host check and `uninstall.sh`.
- Linux has two jobs: ubuntu-22.04 (everything) and ubuntu-24.04 (to the lifecycle and the host check). The 24.04 job cannot receive an archive, so it builds its own bundle in an `ubuntu:22.04` container by the same steps (`packaging/linux/build.sh`, one phase per step; the token passed by name). Its archive is therefore a second build of the same commit, not the same file (sha256 `8fed632a...` against `85b08e03...`; 117 against 118 bundled libraries, the runner image having one more than the bare container).
- Nothing is uploaded by a build job. The plan job fails if the file names a cache, a download of an artifact, or an upload outside the snapshot job. Lane and e2e logs go to files on the runner; the log shows their tables and, only where no table was written, the last 40 lines.

**What building and testing in one job weakens** (not worked around):

| Criterion | As written | What the single job shows instead |
|---|---|---|
| P4: "installs without administrator rights on a fresh runner" | the installer goes to a runner that never held the build | installed on the build runner, from the archive or installer, after the build tree, the moved bundle and the unpacked shell were deleted. The runner still has R's build dependencies (the `-dev` system libraries, compilers, Rust, on Windows Rtools and R 4.6.2): a library missing from the bundle but present on the build runner would not be seen. Linux mitigation, partial: every package loads in a bare `ubuntu:22.04` container with the network off (before packing), and the installed R starts in Arch Linux and Fedora containers. Windows has no such check in stage full (D5's Server Core check is stage `container`, last run token-free) |
| P3 and 2.0.9: lifecycle on `ubuntu-24.04` | the 22.04-built archive installed on 24.04 | a bundle built on the 24.04 runner inside `ubuntu:22.04`: same sources, same steps, a different file from the one tested on 22.04 |
| 2.0.9 lanes "from the installed archive" on a clean machine | fresh runner | installed application, build runner (same reservation as P4) |
| Release candidate: "fresh-runner installs on both platforms" | unchanged | not attempted; needs M20 lifted |
| Transfer integrity (download, checksum on another machine) | artifact and `sha256sum -c` | checksum verified in the same job only |

### 2. The library path (defect 2)

**Repair, packaging only.** `packaging/linux/runpath.sh` writes a run path to `syslib/`, relative to each object's own place (`$ORIGIN/...`, `patchelf` 0.14.3), into R's binary, its modules, every package library and the bundled libraries themselves (326 objects). The launcher (`packaging/linux/Rscript`) no longer exports `LD_LIBRARY_PATH` and removes any it inherits. `assemble.sh` then checks, with only R's own `lib` folder on the path, that nothing but the C library resolves outside the bundle (none did).

**Why this mechanism.** R starts every child through the C library's shell (`system`, `popen`) or `processx`, all of which inherit R's environment, and R's own front-end is a bash script: any exported path reaches host programs, the shell first. A run path is read only by the object that carries it, so the bundled libraries are visible to the bundled R and to child R processes (which load the same objects) and to nothing else. The alternatives: `R_HOME/etc/ldpaths` would keep the path in R's environment, so `system()` children would still get it; restricting what is bundled would leave R without its readline and ICU on hosts that carry other versions. `patchelf` was not on the build before; it is one build-time tool (already on the 22.04 runner image), not a dependency of the bundle.

**One consequence found by the first run.** MultiScholaR's install records the sha256 of the binaries of `duckdb`, `duckplyr` and `fastduckdb` (`src/build_dia_rsum.R`), and the DIA rollup and the direct writer refuse a build whose binaries changed afterwards. Writing the run path after the install changed them. So `install-library.R`, only in the bundle build (`MSR_BUNDLE_RUNPATH=1`), writes the run path into those three packages before MultiScholaR is installed, and `assemble.sh`'s pass leaves objects that already carry it alone ("3 objects", then "323 objects (3 had it already)"). After that: manifest supported, rollup admitted, writer admitted, on both runners. No package code changed.

**Hosts (run 44232555102, the installed application's R; containers with the network off):**

| Host | bash under the old exported path | Children inherit a path into `syslib` | `bash -c true`, `sh -c true`, a pipe of host programs from R | `processx` bash; child R loading stringi, xml2, curl, MultiScholaR; the entry (plain launcher, trial unit) answering HTTP 200 |
|---|---|---|---|---|
| ubuntu-22.04 runner | runs | no | pass | pass (8 of 8) |
| ubuntu-24.04 runner | runs | no | pass | pass (8 of 8) |
| `archlinux:latest` (glibc 2.44, `/bin/sh` is bash) | **fails, exit 127: the defect reproduced** | no | **pass** | **not shown**: R's `utils` did not load in the container |
| `fedora:latest` (Fedora 45, glibc 2.44) | runs | no | pass | not shown, same cause |

So the defect itself is shown repaired on Arch's userland: the bundled R starts through the host's bash and starts `bash` and `sh`, where the old path stops bash. The rest of the check did not run there, for a cause outside the library path: both minimal images have no `/usr/bin/which`, R's `Sys.which()` runs that fixed path, and `utils`' load hook fails on it ("error in running command"; R's startup then warns that `utils` and `stats` were not found). **Finding: the bundled R needs `/usr/bin/which` on the host.** Desktop installs normally have it; a minimal one does not, and then the application cannot start. The budget is spent, so the check was not rerun with `which` added to the images.

Not shown on any host: the laptop itself and the systemd worker from the bundle (P3 step 4, still open: it needs an archive on the laptop, which M20 forbids); Debian; any host on glibc 2.35 to 2.38 other than Ubuntu 22.04.

**Recorded, not changed.** R's own `etc/ldpaths` puts `/usr/local/lib` and `/usr/lib/x86_64-linux-gnu` in front of everything, as it did before the repair. On Debian and Ubuntu hosts a host library of the same file name is therefore taken before the bundled one (on the 24.04 runner R used the host's `libreadline.so.8.2`, `libcurl`, `libgfortran`; on Arch and Fedora that folder does not exist and the bundle's are used). Children inherit those two folders and a Java folder, none inside the bundle. Whether the bundle should set `R_LD_LIBRARY_PATH` to its own `lib` only is a separate decision: it changes which BLAS and Fortran runtime a Ubuntu 24.04 user computes with.

### 3. Lifecycle and lanes from the installed application

| | ubuntu-22.04 | ubuntu-24.04 | windows-2025 |
|---|---|---|---|
| Run | 44232555102 | 44232555102 | 44226281203 |
| Archive or installer | 1,020 MiB, installed 2,645 MiB, 21 s | 1,019 MiB, 2,645 MiB, 21 s | NSIS 772 MiB, installed 2,563 MiB, 326 s; longest path under the folder 143 (limit 170) |
| Library from the moved folder, network off | 389 of 389 | 389 of 389 | passed |
| Internal packages, DIA manifest, rollup, writer | passed | passed | passed (no DIA, D4) |
| Lifecycle, installed, network closed | passed, 16 checks; launcher `process` (kernel 6.8) | passed, 17 checks; `process`, "the user D-Bus does not answer"; `--no-sandbox` (AppArmor 1) | passed, 15 checks; `process` |
| Browser e2e lanes | 7 of 7 with the laptop's counts (88, 190, 16, 42, 86, 98, 153) | not run (lifecycle only) | 3 of 4 (88, 190, 16); TMT and LFQ 136 of 153, the report download (defect 3) |
| GUI-entry lanes | 7 identical (prot_dia, prot_dia_limpa, prot_tmt, prot_lfq, prot_lfq_fragpipe, lipid_canonical, trx_screen); 3 metabolomics lanes differ from `norm` on, the same digests as run 44195739758 (defect 4, the fixture) | not run | as run 44195739758: three protein lanes error (defect 3), trx_screen errors under `testServer`, metabolomics as on Linux, lipid log text |
| Uninstall | passed | passed | not in the job |

`trx_screen` on Linux is identical to its oracle here (29 passed) where run 44195739758 had 7 failures on `top_de_genes.tsv`. Nothing in this change targets it and the cause was not looked for: read it as one observation, not as closed.

### Acceptance (2.0.9), 2026-11-18

| Criterion | State | Evidence, and what is weaker |
|---|---|---|
| Built by `distributable.yml` from the lockfile at the snapshot date | met for these builds | version table check passed in each job; installers are not kept (M20) |
| Linux: seven e2e lanes with the laptop's counts, from the installed archive | met in the weaker form | 44232555102, installed application on the build runner |
| Linux: the suite per file | not run | needs a laptop check-mode table of the same tree, then one sharded run |
| Windows: four non-DIA lanes and the suite | not met | TMT and LFQ fail (defect 3, `kit-052-win-paths`); suite not run |
| Lifecycle on ubuntu-22.04, ubuntu-24.04, windows-2025, network closed | met in the weaker form | installed in the build job; 24.04 on its own container-built bundle |
| Lifecycle, systemd branch, laptop | not met | defect 2 repaired in packaging and shown on the runners and on Arch's shell; the laptop run needs an archive there |
| No token run uploaded a library or installer; no cache | met for these two runs | 0 artifacts, no new cache |
| Release candidate | not attempted | needs M20 lifted |

Defects: #2 repaired (packaging), shown as above; #3 not mine (once `kit-052-win-paths` is merged into `kit-134`, one `[dist:full] [dist-windows]` run is enough: the job already runs prot_tmt, prot_lfq, prot_lfq_fragpipe and the TMT and LFQ e2e lane from the installed application); #4 waits for the fixtures' re-record.

## Windows figure paths, Linux hosts, M25 and M26 (2026-11-18, run 44241609555)

One run, stage full, both platforms: [44241609555](https://github.com/APAF-bioinformatics/MultiScholaR/actions/runs/44241609555) on `fca957ae` (after the merge of `kit-052-win-paths`, `390a0f7c`), 05:04 to 06:44 UTC. Budget declared in `fca957ae`: this run and one correction pass; **the correction pass is not used** (every failure below is a finding, not a fault of the run). 0 artifacts; no cache written (the eight caches on the branch are of 2026-11-17). Evidence is the job logs. Package code under `R/` and `tests/` was not changed beyond the merge; no oracle changed.

New verification files: `packaging/verify/figure-paths.R` (one case: a device probe, KIT-052's tests, the managed DA routes, the prot_tmt entry lane), `figure-trace.R` (traces every PNG/PDF device, `file.copy` and `checkWorkflowFigureFiles` call and records the path's length and whether the file is in its folder afterwards), `windows-figure-cases.ps1` (one case per project-root length, set through `TMPDIR`, under which every fixture makes its project). Instrumentation note: the `prot_worker_da` item's test filter did not narrow `test-omic-da-module.R` to its one test, so the whole file ran (1,126 expectations); that gives the lipidomics, metabolomics and transcriptomics worker routes as well, and is reported as such.

### 1. Windows (windows-2025, the installed application's R, the shell's environment, network cut)

**Lanes now.**

| Lane | Before (44226281203) | Now |
|---|---|---|
| e2e TMT and LFQ | 136 of 153, the report download | **153 of 153** (laptop count) |
| e2e lipidomics, metabolomics, transcriptomics | 88, 190, 16 | 88, 190, 16 |
| prot_tmt, prot_lfq (entry) | error at QC, "required Q-value figures are incomplete" | **QC and every later step pass**; 35 passed, 2 failed: the `session_summary` text at `summary_report` and `summary_export` holds the report's absolute path (`C:/Users/runneradmin/...`, or `D:/a/_temp/...`) where the oracle has `<ROOT>`: the recorder's root substitution does not catch the Windows form of the path (harness, not classified further) |
| prot_lfq_fragpipe (entry) | error at QC | QC passes; 31 passed, 6 failed: the same summary text, and `da_proteins_groupKO-groupWT_long_annot.tsv` differs by digest in 6 steps (not classified; values not compared) |
| lipid_canonical, metab_* | as before | as before (log text; the fixture, defect 4) |
| trx_screen (entry) | "module session has been destroyed" | unchanged |

**The diagnosis, by observation (device probe, folder 24 characters, names under 150).** `exists` is from the folder's listing.

| Path length | grDevices::png (windows) | png (cairo) | ragg::agg_png = ggsave's PNG | pdf / ggsave pdf | writeLines | file.copy |
|---|---|---|---|---|---|---|
| 239, 254, 258, 259, long paths on | yes | yes | yes | yes | yes | yes |
| 260, 269, 299, long paths on | yes | yes | **no, and no error or warning** | yes | yes | yes |
| 239 to 259, long paths off | yes | yes | yes | yes | yes | yes |
| 260, 269, 299, long paths off | no ("unable to start png() device") | no (error) | no, silent | no (error) | no (error) | no (warning) |

**Confirmed:** `ggsave` draws PNG with `ragg` (installed in the bundle), and `ragg` cannot write a path of 260 characters or more on Windows **even with long paths enabled**, and says nothing. In the prot_tmt lane from the runner's temporary folder (project root 98 characters) the Q-value figures are now drawn in the short temporary folder (107 to 121 characters, all present) and copied into the staged folder (244 to 258 characters with the short names, all present); with the earlier names (35 characters longer) they were drawn directly at 279 to 293, which is the failure of runs 44195739758 and 44226281203. `test-workflow-figure-paths.R`: 65 passed with long paths on; with them off 58 passed and 1 error (its first test makes a folder over 248 characters, which needs long paths).

**The worker's request folder, observed.** Figures drawn by the worker in `<root>/scripts/<omic>/.dia-native-<14>/.workflow-entities/.workflow-export/wfx-<32>/omics-publication/<role>/...` (transcriptomics MCI-031.2, managed route, worker entry point in process):

| Case | Project root | Figure | Path length | Device | Exists |
|---|---|---|---|---|---|
| long paths on | 63 | `publication_graphs_dir/normalization_comparison_grid.png` | 260 | ragg | **no** |
| long paths on | 63 | `Volcano_Plots/volcano_plot_Sera_vs_RPMI.png` (and `RPMI_vs_Sera`) | 270 | ragg | **no** |
| long paths on | 63 | `NumSigDeMolecules/summary_barplot_Sera_vs_RPMI.png` (and the other) | 277 | ragg | **no** |
| long paths on | 63 | the same figures as PDF; the `.xlsx` copy | 260 to 277; 274 | pdf; file.copy | yes |
| long paths on | 63 | every other figure of the run | up to 259 | ragg, pdf | yes |
| long paths on | 102 | 12 PNG figures, QC (`gene_qc_dir/trx_*_pca.png`, `_rle.png`, 268 to 269), canonical correlation (296), normalisation grid (297), volcanoes (307), bar plots (314) | 268 to 314 | ragg | **no** |
| long paths on | 102 | the PDFs and copies | up to 314; 311 | pdf; file.copy | yes |
| long paths off | 63 | normalisation grid and volcano, PNG **and** PDF | 260, 270 | ragg, pdf | **no** |
| long paths off | 63 | `canonical_correlation_optimization_summary.xlsx` | 274 | file.copy (FALSE) | **no** |
| long paths off | 100 | the run stops earlier (the test's "Could not save" and file checks fail); every figure it reached was 192 or shorter and written | | | |

With long paths on, the test passes (8 of 8) while those PNGs are missing: **the worker route loses figures silently** (nothing checks them, unlike the Q-value route after KIT-052). The same appears on the lipidomics and metabolomics worker routes of `test-omic-da-module.R` at roots 126 and 127 (volcano and bar-plot PNGs 279 to 291, long paths on; heatmap PNGs 283 to 289, long paths off). The proteomics managed DA fixture (KIT-112 W3) stores its request folder at `<root>/.workflow-entities/...` (no `scripts/proteomics/.dia-native-<14>`, about 47 characters shorter than the application's): every figure written, longest 259 at root 93. In process (the prot_tmt lane): every figure written at root 118 with long paths on (longest PNG 253, copies up to 310); with long paths off it passes at root 80 and at root 120 QC refuses: all 48 Q-value staged copies fail (266 to 280 characters).

**LongPathsEnabled** on the runner: 1. Set to 0 with `Set-ItemProperty` for the later step of the same disposable job; R processes started afterwards read 0 and every writer failed at 260 or more, so the change took effect.

**Safe project-path length on Windows, as observed:**

| | Long paths on (the runner) | Long paths off (Windows' default) |
|---|---|---|
| Limit per path (every case) | 259: `ragg` PNG fails silently from 260; everything else works beyond | 259 for every writer |
| In-process steps (the prot_tmt lane: QC, Q-value figures, DA) | works at 118 (the largest root tried) | works at 80; fails at 120 (Q-value staged copy) |
| Worker route (managed DA and omics publications; the route the plain launcher uses on Windows) | **PNG figures missing at 63**, silently | figures missing at 63, the step stops at 100 |
| Safe root, observed | **below 63** (no root between 52 and 62 was tried on the worker route) | **below 63** |

Calculated from the observed 259 and this run's longest worker figure (`summary_barplot_<contrast>.png`, 214 characters after the root with a 12-character contrast): 45 characters, less for longer contrast names. The user documentation should not state a number until the worker route draws in a short folder and copies (the KIT-052 repair, applied to the request folder) or the request folder is moved; then the limits of KIT-052's table apply.

### 2. Linux hosts; M25 (`which`) and M26 (etc/ldpaths)

**M25.** `packaging/linux/Rscript` stops with exit 69 and a message naming `/usr/bin/which` and the install command per distribution; `install.sh` refuses without it; `main.js` shows the same advice for exit 69 (not exercised: no window in a container); the notices carry a Linux requirements paragraph (checked in the install step). In the published `archlinux:latest` and `fedora:latest` images: exit 69, the message, R not started (both runners).

**Host check** (the installed application's R; containers with network none, `which` added to Arch and Fedora; `debian:stable` = Debian 13 has it): 8 of 8 on every host, before and after M26, on both runners: ubuntu-22.04, ubuntu-24.04, Arch Linux (glibc 2.44; bash under the old exported path still exits 127, the defect, shown repaired), Fedora 45 (glibc 2.44), Debian 13 (glibc 2.41): child R, `processx`, the entry answering HTTP 200 on the plain launcher, `bash -c true`, `sh -c true`.

**M26, measure, change, measure.** The change: `assemble.sh` rewrites R's default in `etc/ldpaths` from `${R_HOME}/lib:/usr/local/lib:/usr/lib/x86_64-linux-gnu` to `${R_HOME}/lib`; nothing from `syslib/` is exported. The "before" rows are the same installed copy with R's original file put back (`$OUT/ldpaths.original`), then the build's.

| Host | readline | BLAS | LAPACK | curl | gfortran | quadmath | host objects mapped in R |
|---|---|---|---|---|---|---|---|
| ubuntu-22.04, before | host 8.1 | **host** OpenBLAS 0.3.20 (`openblas-pthread/libblas.so.3`) | **host** | host 4.7.0 | host | host | 52 |
| ubuntu-22.04, after | bundle | bundle (`syslib/libblas.so.3`, `libopenblas.so.0`) | bundle | bundle | bundle | bundle | 4 (libc, libm, libresolv, the loader) |
| ubuntu-24.04, before | host 8.2 | bundle | bundle | host 4.8.0 | host | none | 50 |
| ubuntu-24.04, after | bundle | bundle | bundle | bundle | bundle | bundle | 4 (libc, libm, libresolv, the loader) |
| Arch, Fedora, before and after | bundle | bundle | bundle | bundle | bundle | bundle | 5 (C library, nss_systemd) |
| Debian 13, before | bundle | bundle | bundle | bundle | bundle | bundle | 17 (libstdc++, libssl, libcrypto, ...) |
| Debian 13, after | bundle | bundle | bundle | bundle | bundle | bundle | 4 |

The child R process shows the same copy as its parent on every row. The bundled R's BLAS is the bundled OpenBLAS (`La_library()` = `syslib/libopenblas.so.0`).

**Lanes after M26** (installed application, network closed, the laptop's locale): ubuntu-22.04 and ubuntu-24.04 each: prot_dia 44, prot_dia_limpa 38, prot_tmt, prot_lfq, prot_lfq_fragpipe 37 each, lipid_canonical 35, trx_screen 29: **identical to their oracles** (digests of values rounded to 10 significant digits); metab_lc, metab_gc, metab_combined differ as before (the fixture, defect 4) with **the same state digests as run 44232555102 before the change**, on both runners. So no value moved: nothing for `compare-results.R` to judge. ubuntu-22.04 e2e lanes: 88, 190, 16, 42, 86, 98, 153 (the laptop's counts). Lifecycle passed on both runners (plain launcher).

### Acceptance (2.0.9), rows that changed

| Criterion | State | Evidence |
|---|---|---|
| Windows: four non-DIA lanes and the suite | **lanes: met** (e2e 88, 190, 16, 153 with the laptop's counts); suite not run | 44241609555; entry lanes: protein QC passes, summary-text and fragpipe DA-table differences open (above) |
| Linux: seven e2e lanes | met in the weaker form, again after M26 | 44241609555, ubuntu-22.04 |
| Linux hosts (defect 2's "not shown") | **shown**: Arch, Fedora, Debian 13, full check, with `which` (M25) | 44241609555 |

Defects: #3 (Q-value figures on Windows) **repaired** by `kit-052-win-paths` and shown; **new, #5: the worker route's figures over 259 characters are lost silently on Windows** (package code, KIT-052's scope; not changed). #2 extended by M26 (bundle libraries on every host). #4 unchanged.

## Windows figure paths after the worker-route repair; the FragPipe table (2026-11-18, run 44275615462)

One stage-full run, both platforms (the Linux jobs serve as the reference for the FragPipe table; a second dispatch on this ref would cancel the first): [44275615462](https://github.com/APAF-bioinformatics/MultiScholaR/actions/runs/44275615462) on `7eb6fc73` (merge of `kit-052-win-paths` at `b0f4719c` as `72383d0e`, then the instrumentation). Budget: this run and one correction pass; **the correction pass is not used**. Every job ended "failure" because of a lanes step whose results are findings, not faults of the run: Windows lanes as in 44241609555, Linux metab lanes (defect 4), and on ubuntu-22.04 trx_screen's `top_de_genes.tsv` (7 steps, digest; seen before in 44195739758 but not in 44241609555; not classified here). 0 artifacts, no cache. Package code under `R/` and `tests/` unchanged beyond the merge; no oracle changed.

Instrumentation added (packaging/verify only): `figure-trace.R` also traces `file.rename` (the final name of a published file) and `rlang::abort` of class `workflow_figure_unwritten` (RAISED lines, with the message), and pairs every drawn figure with its later copies and renames of the same name (FIGPAIR lines: drawn path and length, final path and length, present and not empty). `figure-paths.R` item `omic_worker_da`: the named lipid managed test of `test-omic-da-module.R`; testthat's `desc` filter also runs the file's top-level loops ((b) in process and (b-managed) lipid, metab, trx), so the item covers the lipid and metabolomics worker routes, not the whole file. A short TMPDIR (`D:\o`) gives the lowest root these fixtures can have. Fixture roots: trx `.../project/kit097` (TMPDIR + 43); lipid/metab: the trace's root is the store prefix up to `.dia-native-<hex>`, so it reads about 25 longer than the fixture's project folder.

### 1. Figures (windows-2025, installed application's R, network cut)

"Figures" = every PNG or PDF a device drew in the item; "final ok" = every copy and rename of it present and not empty. The worker's figures are drawn in `tempdir()/wfx-<hex>/...` (the longest drawn path of any worker item: 239).

| Case (LongPathsEnabled as R reads it) | Root | figure tests | trx managed DA (MCI-031.2) | lipid/metab worker routes | trx in the omic file | Longest final path |
|---|---|---|---|---|---|---|
| on (1), root 63 | 62-63 | 88 passed, 2 skipped | 42 of 42 final ok, 8/0 | 129 of 129 | 170 of 170 | 275 |
| on, root 102 | 101-102 | 88/0 | 42 of 42, 8/0 | 129 of 129 | 170 of 170 | 315 |
| on, root 124 | 123-124 | 88/0 | 42 of 42, 8/0 | 129 of 129 | 170 of 170 | 336 |
| off (0), root 47 (the lowest possible) | 47 | 70 passed, 3 errors (tests that make 248+ folders) | **38 of 42**: both bar plots, PNG and PDF, final 261, refused | 129 of 129 (longest 237) | 168 of 170: 2 bar plot files refused (261) | 261 |
| off, root 63 | 61-62 | 70/3 errors | **32 of 40**: volcanoes 267, bar plots 274, refused; the `.xlsx` 271 | 129 of 129 (longest 250) | 164 of 168: 4 refused (268, 275) | 275 |

`omic_worker_da` has 9 failures in every case with long paths on: frozen-capture differences of the `*_da_results_<DATE>.csv` sha and the trx row digest (as in 44241609555, before the repair; not figures).

**The refusal (long paths off), as raised**, class `workflow_figure_unwritten/saved_entity_outputs_pending`: "the worker's publication files could not be copied into its request folder | missing or empty: 'summary_barplot_RPMI_vs_Sera.pdf', 'summary_barplot_RPMI_vs_Sera.png', 'summary_barplot_Sera_vs_RPMI.pdf', 'summary_barplot_Sera_vs_RPMI.png', 'volcano_plot_...' | in '<root>/scripts/transcriptomics/.dia-native-<hex>/.workflow-entities/.workflow-export/wfx-<32>/omics-publication/publication_graphs_dir/...' | the longest of these paths is 271 [261 at root 47] characters: the project folder's path is too long for Windows with long paths disabled, where a path may have at most 259 characters. Move the project to a folder with a shorter path, or enable Windows long paths." The test fails (trx 5 passed, 3 failed): a loud failure, not a silent loss and not a pass.

**prot_tmt entry lane, long paths on** (the roots reached were 120 and 125, not 118 and 124: the fixture folder's hex suffix varies by one or two):

| Root | Figures | Lane | Missing |
|---|---|---|---|
| 120 | 121 of 121 (longest drawn 255, longest copy 280) | 35 passed, 2 failed (the summary-text root, as before) | none |
| 125 | **120 of 121** | 31 passed, 6 failed | `publication_graphs_dir/NumSigDaMolecules/da_proteins_num_sig_da_molecules.png`, drawn in place by ragg at **260**: no file, **no refusal raised** (the in-process DA publication does not check its figures) |

### 2. Safe project-path length on Windows after the repair

| | Calculated (KIT-052 note) | Observed here |
|---|---|---|
| Long paths on | 124 (in-process proteomics DA number-significant PNG drawn in place) | worker routes: no limit up to root 124 (trx final 336; lipid/metab 311); in process: 120 safe, **125 loses exactly that PNG (260), silently**. So 124, as calculated; the loss at 125 is silent (a remaining defect of the in-process DA route, not of the worker route). |
| Long paths off | 41 (metab `num_sig` xlsx, 218 after the root); 43 for the trx bar plot at a 12-character contrast | the fixtures cannot go below root 47: there the trx bar plots' copies (214 after the fixture root, 261) are refused loudly; at 63 the step stops loudly. Not shown passing at 41. Consistent with the calculation (fixture's trx limit 45 = app layout 43, the `.dia-native` suffix 2 shorter in the fixture); the metab xlsx bound (41) is not exercised by these fixtures, whose metab copies stay at 116-131 characters. |

### 3. e2e lanes (Windows, installed application)

lipidomics-canonical 88, metabolomics-lc-gc 190, transcriptomics-screen 16, proteomics-tmt-lfq 153: the laptop's counts, nothing regressed. ubuntu-22.04: 88, 190, 16, 42, 86, 98, 153. GUI-entry lanes on Windows as in 44241609555 (prot_tmt and prot_lfq 35/2, prot_lfq_fragpipe 31/6, lipid_canonical 27/8, metab 25/10, trx_screen error).

### 4. prot_lfq_fragpipe `da_proteins_groupKO-groupWT_long_annot.tsv`: classified numeric, within the rule

Each copy the recorder read (6 per platform, one distinct content per platform) was kept in its job and fingerprinted (FPRINT lines): bytes, line endings, header, shape, row order, per column a text digest and digests of the values rounded to 15/12/10/8/6 significant digits, two sample cells.

- **Not textual, not structural**: header identical (41 columns), 400 rows, row order identical, LF only on every platform (no CRLF, no BOM); every text column (17, including `annotation_score`, ids, GO terms) identical on Windows, ubuntu-22.04 and ubuntu-24.04.
- **Numeric**: 3 statistics columns and 20 intensity columns differ.

| Columns | Windows vs Linux (either) | ubuntu-22.04 vs 24.04 | Sample |
|---|---|---|---|
| `fdr_qvalue`, `raw_pvalue` | equal at 10 significant digits, not at 12 | equal at 10, not 12 | 1.2396646615107782e-23 (Linux) vs 1.2396646615105722e-23: 1.7e-13 relative |
| `log2FC` | equal at 8 (vs 22.04) or 10 (vs 24.04) | equal at 8, not 10 | 3.1801726079068473 vs 3.180172607906851: 1.2e-15 |
| `log2norm.*`, `raw.*` (14 of 20) | equal at 12, not 15 | identical text | 12.397269757248953 vs 12.397269757248955: 1.6e-16 (last bit) |
| `log2norm.*`, `raw.*` (6 of 20) | equal at 15 | identical | |

Bounds under the fixed rule (relative 1e-6, absolute floor 1e-9): same rounding at s significant digits bounds the relative difference near 10^(1-s), so the largest relative difference is below about 1e-9 for the p and q values, below about 1e-7 for `log2FC` (some value near zero, where the floor applies), and about 1e-11 or less for the intensities: **every column within the rule**. The exact maximum needs the two tables side by side, which this run deliberately does not allow. The Linux hosts differ from each other by as much and both still match the oracle: the recorder's mask (`helper-gui-entry.R`, `guiEntryMask`) replaces every run of 10 or more digits with `<HEX>`, so long mantissas never reach the oracle's digest. Windows' digest difference must therefore come from a number whose printed form changes length (for example a value printed as an integer on one platform and with a long fraction on the other); inferred, not shown. Correction to the run's commit message: the oracle's digest covers the whole masked text (`guiEntryMask` collapses the lines before `digest(serialize = FALSE)`), not the first line only; the "first-line digest" in the FPRINT lines is not comparable with it, and the classification above does not use it.

### What remains

- The in-process proteomics DA publication draws `da_proteins_num_sig_da_molecules.png` in place and does not check it: from root 125 with long paths on (and earlier with them off) the PNG is lost **silently**. Package code (KIT-052's scope; not changed here).
- Long paths off: the limit (about 41 to 43) is set by the worker's request-folder names (the pending harness decision); not shown at 41 because these fixtures cannot have a root shorter than 47.
- The FragPipe table's exact maximum difference (only bounds here); the oracles' table digests do not see digits past the ninth (the mask), a weakness of the recorder, not changed.

## Windows figure paths after the in-process and folder-name repairs; trx_screen `top_de_genes.tsv` (2026-11-19)

**Declared before the run.** Purpose: confirm on windows-2025, with the installed application, `kit-052-win-paths` at `4b6770fd` (merged into kit-134 as `e7cabee5`, no conflict): 463abfb9 (in-process proteomics DA publication drawn in a short folder, copied and checked), fe6b7f54 (worker request folder `p/<role letter>`, old layout still read), 4b6770fd (normalisation, RUV covariates, protein QC in-process publications through the same helper); and classify on the Linux jobs the open difference of trx_screen's `top_de_genes.tsv` (7 steps, digest, ubuntu-22.04 in 44275615462). Bindings: kit-134 at the commit that carries this section, stage full, both platforms (the Linux jobs serve the table's classification; ubuntu-24.04 is the reference). Package code under `R/` and `tests/` unchanged beyond the merge; no oracle changed; instrumentation in `packaging/verify` and the workflow only. Ceiling: 300 minutes per job as before; each figure step 120 minutes (was 90; more cases). Budget: this run and one bounded correction pass, used only for an instrumentation fault; an unfavourable valid result is a finding, not a reason to rerun. Decision after the run: each KIT-052 Windows row confirmed or reported with its cause; the table's difference within the fixed rule (relative 1e-6, absolute floor 1e-9, discrete exact) or real.

Cases (root lengths are approximate until the run states them; trx root = TMPDIR + 43, prot_tmt root = TMPDIR + 59 to 60 in 44275615462):

| Long paths | Root | Items | Expected |
|---|---|---|---|
| on | 102 | figure tests, trx managed DA (MCI-031.2), lipid/metab/trx worker routes, old-layout DA check | every figure present |
| on | 124, 130 | trx managed DA, worker routes | every figure present |
| on | 125, 127, 130 | prot_tmt entry lane (in process: QC, normalisation, RUV covariates, DA) | every figure present, including `da_proteins_num_sig_da_molecules.png` (lost at 125 before) |
| off | 63 | figure tests, trx managed DA, worker routes, old-layout DA check | every figure present |
| off | 70, 77 | trx managed DA, worker routes | every figure present (calculated safe: trx 77 at a 15-character contrast; the fixture's is 12) |
| off | 85, 95 | trx managed DA, worker routes | the step stops with the classed message naming the files; nothing lost silently |
| off | 90 | prot_tmt entry lane | every figure present (below the Q-value staged copy's calculated 99) |
| off | 118 | prot_tmt entry lane | stops loudly (expected at QC's Q-value staged copy, before normalisation and RUV covariates can be reached); nothing lost silently |

Also: the four non-DIA e2e lanes on Windows (expected 88, 190, 16, 153). Table instrumentation (`run-lanes.R`): the recorder's copies of `top_de_genes.tsv` stay in the job; FPRINT lines give bytes, line endings, the oracle's own form (the recorder's `guiEntryMask`, parsed from the helper, then sha256) against the oracle's stored digests, digests of the body in order and sorted (raw and masked), per column text and masked-cell digests, the printed-digit range and digests of values rounded to 15, 13, 12, 11, 10, 9, 8, 7 and 6 significant digits.

**Run 1: [44368107499](https://github.com/APAF-bioinformatics/MultiScholaR/actions/runs/44368107499) on `1542823e` (18:17 to about 19:30 UTC, 2026-11-18).** Windows: **instrumentation/build fault, nothing tested.** Both internal packages failed to build (`fast_duckdb@e724f90b`, `fast_limpa_dpc@87a90b89`, "R CMD INSTALL failed"; their build output is deliberately not printed), so MultiScholaR could not install. The pins were the same as in 44275615462, where both built. What changed: runner image 20261103.298.1 to 20261116.315.1, and Rust stable (unpinned, audit F6) 1.99.1 to 1.100.0. The Linux jobs built both packages with Rust 1.100.0, so the failure is specific to Windows (the toolchain on the GNU target, or the image). The cause between the two was not isolated. Linux (both runners): lanes as before (metab fixture, defect 4; everything else identical to its oracle); ubuntu-22.04 e2e 88, 190, 16, 42, 86, 98, 153. **trx_screen identical to its oracle on both runners**: one distinct copy of `top_de_genes.tsv` (10 read by the recorder), the same bytes on 22.04 and 24.04 (sha256 `0de193f5614262e0`, 2,975 bytes, 20 lines, LF, 19 rows by 10 columns), its masked form `605111e573e3ed73` = the oracle's. That is the reference fingerprint for any differing copy. In 44275615462 the differing copy (ubuntu-22.04 only) had digest `8c3a4a5f380990c6`, from step `da` on (7 steps on `da_genes/`, 3 on `Publication_tables/`). So the difference does not recur on every run of the same commit and runner type (seen in 44195739758 and 44275615462; absent in 44232555102, 44241609555 and here).

**Correction pass (the single one), declared before launch.** Purpose: the same Windows cases (unchanged), with the Windows Rust toolchain pinned to 1.99.1 (the one that built these pins in 44275615462), and one more chance at a differing `top_de_genes.tsv`: trx_screen runs eight times per Linux runner (labelled `trx_screen_r2` .. `_r8`; about 45 s each). Both platforms, stage full, on the commit carrying this paragraph. Ceilings as declared above. If Windows fails to build again, or the table does not differ, that is the result: no further run.

**Correction pass: [44379470655](https://github.com/APAF-bioinformatics/MultiScholaR/actions/runs/44379470655) on `075a2640`.** The budget is now spent. With Rust 1.99.1, Windows built both internal packages again (image 20261116.315.1 unchanged), so the run-1 failure was Rust 1.100.0 on the Windows GNU target. Every Windows step ran. 0 artifacts, no cache. Root lengths are those the trace reached (trx = the fixture's project folder; prot_tmt = the lane's project). "Final ok" = every copy and rename of a drawn figure is present and not empty.

| Long paths (as R read it) | Root reached | Route | Drawn / final ok | Refused loudly | Lost silently |
|---|---|---|---|---|---|
| on (1) | 101 to 110 | figure tests 116/0 (2 skipped); trx managed DA 8/0; worker routes; old-layout DA check 25/0 | 8/8, 42/42, 299/299, 5/5 | none | none |
| on | 124 | trx managed DA 8/0; worker routes | 42/42, 299/299 (longest final 301) | none | none |
| on | 130 | trx managed DA 8/0; worker routes | 42/42, 299/299 (longest final 307) | none | none |
| on | 126, 128, 131 | prot_tmt entry lane (in process) 35/2 each (the summary-text root, harness) | 121/121 each; `da_proteins_num_sig_da_molecules.png` present (final 261, 263, 266) | none | **none** (was lost at 125 in 44275615462) |
| off (0) | 62 to 72 | figure tests 84 passed, 7 errors; trx managed DA 8/0; worker routes; old-layout DA check 25/0 | 42/42, 299/299, 5/5 (longest final 240) | none | none |
| off | 69 | trx managed DA 8/0; worker routes | 42/42, 299/299 (245) | none | none |
| off | 76 | trx managed DA 8/0; worker routes | 42/42, 299/299 (252) | none | none |
| off | 85 | trx managed DA 5 passed, 3 failed; omic file +1 failure | final ok 38/42 and 297/299 | **yes**: the 4 and 2 trx bar plots (261, 262), named, class `workflow_figure_unwritten/saved_entity_outputs_pending`, "the project folder's path is too long for Windows with long paths disabled" | none |
| off | 94 | trx managed DA 5/3; omic file +6 | final ok 32/40 and 293/297 | **yes**: the `.xlsx` (266, 267), volcanoes (262, 263), bar plots (269, 270), named | none |
| off | 91 | prot_tmt entry lane (QC, normalisation, RUV covariates, DA in process) | 89/89 final ok (DA PNG final 226) | none | none among the figures; **but the lane errors at `summary_report`: no Summary report `.docx` in `results_summary`** (see below) |
| off | 117 | prot_tmt entry lane | 61 drawn, 13 final ok | **yes**: QC's Q-value staged copies, 4 refusals naming 48 files (263 to 271) | none |

Every final file missing at 85, 94 and 117 is named in a refusal (FIGPAIR against RAISED). Lipid and metabolomics worker routes: everything present at every root up to 94, long paths on or off. The baseline 9 failures of `test-omic-da-module.R` (frozen-capture `*_da_results_<DATE>.csv` sha and trx rows, as in 44241609555) are present in every case.

- **Observed safe project root, long paths off, worker routes:** passes at 76, refused at 85. The fixture's longest trx file is 177 characters after its root, so its limit is 82. The calculation gives 77 for a 15-character contrast in the application's layout; the fixture's contrast has 12 characters (+3) and its `.dia-native-` suffix is 12 hex, not 14 (+2): 82. **The observation agrees with the calculation.** Lipid 85, metabolomics 78 and proteomics 80 are not bounded by these fixtures (their files stay shorter).
- **In process, long paths off, from above 115:** the lane stops loudly at QC (the Q-value staged copy, calculated 99), so normalisation and RUV covariates are not reached there. They are shown under long paths off at 91 (every figure present) and on at 126 to 131. In `test-workflow-figure-paths.R` (with doubles) they pass on Windows with long paths on.
- **New finding, long paths off, root 91:** the prot_tmt lane reached `summary_report` but no report `.docx` was written (the helper's `rendered[[1L]]` is empty; 32 fewer drawn figures than with long paths on, so the report's own figures were never drawn). At root 80 in 44241609555 (before these merges) the lane passed, so the Summary report's limit with long paths off lies between about 80 and 91. Whether the application said so to the user is not captured (the lane's notifications are in its record, not printed). Report rendering is not one of the routes KIT-052 repaired; package code; not changed.
- **Figure tests, long paths off:** 7 errors (3 in 44275615462). Only 4 of the non-passes are printed: folder creation over 248 characters ("cannot create dir"; the `.pending` stage folder under a 190-character project) and one "invalid export publication request" error after it. These are fixtures that need long paths, as before. The other 3 are not shown.
- **Old layout on Windows:** `test-workflow-dia-products.R`'s DA worker handler test (which publishes a request staged as `da-publication/<role>`) 25/0 with long paths on and off; the omics old-layout test inside the figure tests passes with long paths on.
- **Windows e2e:** lipidomics-canonical 88, metabolomics-lc-gc 190, transcriptomics-screen 16, proteomics-tmt-lfq 153: the laptop's counts. GUI-entry lanes as before (prot_tmt and prot_lfq 35/2, fragpipe 31/6, lipid 27/8, metab 25/10, trx_screen error). ubuntu-22.04 e2e: 88, 190, 16, 42, 86, 98, 153.

### trx_screen `top_de_genes.tsv`: classified numeric, within the rule

In this pass ubuntu-22.04 gave the differing table in all 8 repeats: 2,979 bytes, sha256 `6460b0a3b12cc950`, masked form `8c3a4a5f380990c6`, the digest of 44275615462. ubuntu-24.04 gave the oracle's table in 8 of 8 (`0de193f5614262e0`, masked `605111e573e3ed73`), as both runners did in run 1. So within one build the table is deterministic; between builds of the same commit on ubuntu-22.04 runners it is not (inferred cause: the runner's CPU, through the bundle's OpenBLAS kernel choice; not shown).

- **Not structural, not textual:** header identical, 19 rows by 10 columns, row order identical (first-column order digest the same), LF only, no BOM; every text column (`GeneID`, `GeneName`, `DisplaySymbol`, `regulation`, `comparison`) identical, so gene membership and the up/down calls are the same.
- **Numeric**, 22.04 against the reference (24.04 = oracle):

| Column | Equal at (significant digits) | Not at | Relative bound | Sample |
|---|---|---|---|---|
| `logCPM` | 15 | (text differs in the 16th or 17th digit only) | below 1e-14 | 14.866874080741063 (identical sample) |
| `logFC` | 12 | 13 | below 1e-11 | 1.8590001946067354 vs 1.8590001946067227: 6.8e-15 |
| `PValue` | 10 | 11 | below 1e-9 | 2.0058636769723173e-9 vs 2.0058636769655497e-9: 3.4e-12 |
| `FDR` | 10 | 11 | below 1e-9 | 3.7509650759382335e-7 vs 3.750965075925578e-7: 3.4e-12 |
| `F` | 9 | 10 | below 1e-8 | 141.48639598552006 vs 141.48639598564907: 9.1e-13 |

(Equal after rounding to s significant digits bounds the relative difference by 10^(1-s).) **Every column is within the fixed rule** (relative 1e-6, absolute floor 1e-9; the discrete columns exact). The oracle sees it only through `FDR`: of the per-column masked digests only FDR's differs, so at least one FDR cell prints with a digit run shorter than the mask's 10 characters. The printed-digit range of `logFC` (14 to 17 against 15 to 17) accounts for the 4 bytes.

### What remains

- The Summary report is not written at a project root of about 91 with long paths off (limit between about 80 and 91); whether the user is told is not shown. Package code (KIT-052's scope or the report's), next bounded item if wanted.
- The worker routes' limit with long paths off is 82 in the fixture's layout, which is 77 in the application's for a 15-character contrast, as calculated. Next limits (store `.pending` 95, Q-value staged copy 99) belong to the harness revision (v1.11 names).
- Windows Rust is now pinned at 1.99.1 in this workflow. Rust 1.100.0 does not build the internal packages on the Windows GNU target (cause inside the packages' build, not printed by rule). This needs a decision or a fix in the internal repositories before Rust is moved (F6).
- trx_screen's oracle digest depends on the ubuntu-22.04 runner build; the values agree within the rule. The recorder's mask is the weakness, as with FragPipe (not changed).

## Release-candidate run on runners (2026-11-20, run 44598642342, branch run-k134-rc)

Declaration and full results: `rc-run-declaration.md` (this folder, branch `run-k134-rc`). From `performance-v2`
`8b4ff691`, stage `rc`: 2 Linux main, 1 Windows main, 13 suite shards; 0 artifacts, no cache; correction pass unused.

- **The suite on the bundle ran for the first time** (installed package, the bundle's R and library): Linux 383/413 ok,
  every difference from the runner reference classified (known state; E installed-vs-load_all; numeric under the
  bundle's OpenBLAS); **no class D on Linux**; open: OpenBLAS digest values not compared. Windows 357 ok, 53 failed,
  3 timeouts: mostly E (checkout line endings, GNU tools, D4, path forms); **one class D (D-1)**; a few open.
- **D-1 (Windows, package code):** the R-fallback native archive (`runDIANativeTarR`) leaves its staging folder in TMP
  (unlink before the working directory is restored). Not changed.
- NIPALS "bit-identical to mixOmics" holds under reference BLAS only; under the shipped OpenBLAS it agrees to 1e-15.
- Hosts: 8/8 on ubuntu-22.04, ubuntu-24.04, Arch, Fedora, Debian, before and after M26; M25 stop shown on Arch and Fedora.
- Windows: e2e 88, 190, 16, 153; Summary `.docx` at roots 90 and 94 with long paths off; figure-paths test on/off only
  D-1 plus long-folder fixtures; workers draw figures under TMP (shown); the session archive's TMP staging is not
  reachable on a Windows lane (no native archive outside DIA), shown only by the unit test.
- MOFA2 in the installed bundle, double precision, both platforms: trains with the bundle's Python, offline; 206
  iterations, ELBO -24285.131276 on all three machines, factors agree to 5e-15 (the declared 186/-27000.32 belong to
  another corpus revision).
