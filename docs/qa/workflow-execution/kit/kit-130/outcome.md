# KIT-130: packaging spike, a bundled R, library and Python that run from any folder (2026-11-17)

Scope: ticket revision 1.1. GitHub Actions on `ubuntu-22.04` and `windows-latest` (Windows Server 2025), branch `kit-130-dist-spike`, sources under `tools/dist-spike/`, workflow `.github/workflows/dist-spike.yml`. No package code changed.

**Result: yes for every piece on both platforms.** R 4.6.1, limpa, RUVIIIC, MOFA2 and a standalone CPython with mofapy2 run from one folder that was moved, with nothing from the host but the C library (Linux) or the operating system (Windows); an Electron shell shows a Shiny page served by that R and stops R when the window closes; the AppImage and the NSIS per-user installer start on a fresh runner. Two workflow runs per platform were used, of 15 allowed.

| Run | Commit | Result | Link |
|---|---|---|---|
| 1 | `5a2c3488` | pieces 1, 3, 4, 5 passed on both platforms; piece 2 stopped on a defect of the check script; two instrumentation failures in the added pieces | https://github.com/APAF-bioinformatics/MultiScholaR/actions/runs/43922649323 |
| 2 | `9dc7bbcd` | everything passed (the one correction pass) | https://github.com/APAF-bioinformatics/MultiScholaR/actions/runs/43925336690 |

## The pieces

Evidence is run 2 unless stated; the number is the job in that run. Every check prints `KIT-130 CHECK PASSED: <name>` after its assertions.

| Piece | Linux | Windows | How it was shown |
|---|---|---|---|
| R 4.6.1 in a folder | yes (job 133503161621) | yes (job 133503162340) | Linux: R installed at `/opt/R/4.6.1`, copied into the bundle, `/opt/R/4.6.1` deleted, folder moved, then run inside a bare `ubuntu:22.04` container (no R, no Python, none of R's shared libraries), mounted read-only at another path, network off. Windows: folder moved, PATH reduced to `System32`, every R and Python variable removed. Both: `--vanilla`, R home and every library path inside the folder, no other R on the PATH, linear algebra, a graphics device and compression work. |
| limpa, RUVIIIC, MOFA2 | yes | yes | Same conditions. `limpa::dpc` on a simulated set, `RUVIIIC::RUVIII_C` (compiled path, with missing values), `MOFA2::create_mofa` and `prepare_mofa` plus an HDF5 round trip; all 78 loaded namespaces came from the bundle. limpa 1.4.2, RUVIIIC 1.0.19, MOFA2 1.22.1 (Bioconductor 3.23). |
| Python for MOFA2 | yes | yes | CPython 3.12.16 (python-build-standalone 20261110) with mofapy2 0.7.5 and MOFA2's own dependency pins. `MOFA2::run_mofa(use_basilisk = FALSE)` trained the package's example model (15 factors) with the network cut: container network `none` on Linux; outbound firewall rules on the bundle's 12 programs on Windows. A probe showed on both that the same R sees the network when it is up and not when it is cut. Python, its library, mofapy2, numpy and h5py all resolved inside the folder; basilisk was not used; on Windows no basilisk, reticulate or uv cache folder existed afterwards. |
| Shiny in Electron | yes | yes | The packaged, unpacked application on the build runner: the bundled `Rscript` started the app on a free port, the window fired `did-finish-load`, a live session rendered text from the server (R 4.6.1, R home and library inside the application's resources), the port answered 200, and after the window closed the application exited 0, the R processes were gone and the port refused connections. First render 1.5 s (Linux), 3.1 s (Windows). |
| Installer | yes (job 133506460503) | yes (job 133508702744) | Artifact downloaded on a fresh runner, same smoke test. Linux: the AppImage; the runner has no FUSE 2, so it ran with `--appimage-extract-and-run`. Windows: `MSRSpike-Setup-0.0.1.exe /S` installed per user under `%LOCALAPPDATA%\Programs\msr-dist-spike` and the installed application started. The runners' own R is 4.6.2; the page reported 4.6.1 from inside the application. |
| Dependency closure (rev 1.1) | resolved (job 133503161732) | resolved (job 133503161296) | See "Dependency closure". |

Extras, run in the same two runs:

| Extra | Linux | Windows |
|---|---|---|
| The three checks from a path with spaces | passed (`/srv/any where/msr app`) | passed (`D:\a\_temp\re located\msr moved`; short 8.3 names are disabled on that volume, so R had no short path to fall back on) |
| AppImage mounted read-only through FUSE 2 (after adding `libfuse2` to the runner) | passed, first render 2.5 s | not applicable |
| `R CMD INSTALL` of MultiScholaR itself (asked for during the spike; a separate job with the standard setup actions, hard dependencies only, not part of the bundle) | not run | installs (job 133503161620); see below |

## Sizes and times

| | Linux | Windows |
|---|---|---|
| Bundle folder | 1,361 MiB, 34,610 files | 1,056 MiB, 36,767 files |
| of which R with its library | 639 MiB (library 613 MiB, 133 packages) | 572 MiB (library 457 MiB) |
| of which Python | 606 MiB | 484 MiB |
| of which bundled shared libraries | 117 MiB, 92 libraries | none needed |
| Installer | AppImage 586,251,452 bytes (559 MiB) | NSIS 363,011,934 bytes (346 MiB) |
| Installed | not installed; mounted, or 1.3 GiB extracted to `/tmp` at every start without FUSE 2 | 1,325 MiB, 36,910 files |
| Build job, end to end | 321 s | 539 s |
| of which R packages | 152 s (103 binaries from Posit Package Manager, RUVIIIC compiled) | 135 s (binaries from CRAN and Bioconductor, RUVIIIC compiled) |
| of which packaging | 65 s (AppImage about 46 s) | 268 s (NSIS about 205 s) |
| Install on the fresh runner | none | silent install 262 s (279 s in run 1) |
| Start to first rendered session, fresh runner | 7.5 s extracted (13.0 s in run 1); 2.5 s mounted | 3.1 s (3.6 s in run 1) |

Largest parts: BH 156 MiB (headers, needed only to compile), RUVIIIC 96 MiB on Linux (the library compiled on the runner; probably unstripped, not checked), Rhdf5lib 61 MiB; scipy 103 MiB, pandas 75 MiB, scikit-learn 54 MiB. Nothing was trimmed.

## Failures, kept

| Run | Where | What failed | Cause | After |
|---|---|---|---|---|
| 1 | Piece 2, both platforms | `Error in asNamespace(ns, base.OK = FALSE): operation not allowed on base namespace` | The check script asked the base namespace for its path. The three toy calls had already passed (their output is in the run 1 log). A defect of the check, not of the bundle. | Corrected (base excluded; its location is covered by piece 1); passed in run 2. |
| 1 | The path-with-spaces extra, both platforms | the step failed | Only through the same defect; R and the MOFA2 training passed from the path with spaces in run 1 as well. | Passed in run 2. |
| 1 | Dependency closure, Windows | `trying to use CRAN without setting a mirror` | `setup-r` leaves CRAN unset on Windows; the step that installs the JSON reader had no repository. Instrumentation. | Repository given; resolved in run 2. |
| 1 | `R CMD INSTALL` job, Windows | `Can't find package called fastduckdb, fastlimpadpc` | Inside `r-lib/actions/setup-r-dependencies`: its pandoc check resolves `.` on its own, without the two ignores, when pandoc is not on the runner. The dependencies themselves had installed. Instrumentation. | `install-pandoc: false`; passed in run 2. |

No piece needed more than its one correction pass.

## Dependency closure (rev 1.1)

A pak lockfile for the whole `DESCRIPTION` (Depends, Imports, LinkingTo, Suggests) with `fastduckdb=?ignore` and `fastlimpadpc=?ignore`; nothing installed; R 4.6.1 on both runners. The branch is cut from `performance-v2` at `5534aa97`, so the `DESCRIPTION` is that one. Tables: `closure-linux-all.csv`, `closure-windows-all.csv` and the two `-no-binary.csv` files in this folder. Nothing was unresolvable on either platform.

| | Linux (jammy) | Windows |
|---|---|---|
| Dependencies resolved (including base-recommended packages already installed) | 416 | 413 |
| With a binary, as pak resolves by default | 356 | 411 |
| Without a binary | 60 | 2 |
| of which need compilation | 28 | 1 |

**Windows, no binary:** `RUVIIIC` 1.0.19 (archived on CRAN in 2023; from `cran/RUVIIIC`; needs compilation) and `GO.db` 3.23.1 (an annotation data package, source only, nothing to compile).

**Linux, no binary as pak resolves by default:** `RUVIIIC`, the `Glimma` fork (`APAF-bioinformatics/GlimmaV2`, 2.15.2, from GitHub) and 58 Bioconductor packages, which pak takes as source from bioconductor.org. Posit Package Manager does serve jammy binaries for Bioconductor 3.23: the spike's own build installed 19 of the 58 from it as binaries, and a request for each of the 58 the way R 4.6.1 asks (`ppm-bioconductor-binary-probe.txt`) returned a binary for 57. With that repository configured, the Linux list is three: `RUVIIIC`, the `Glimma` fork and `GO.db`.

One difference between the platforms needs a decision. On Windows pak resolved `Glimma` to Bioconductor's 2.22.1 binary, not to the fork named in `Remotes` (2.15.2); on Linux it resolved to the fork. A `Remotes` entry is a candidate, not a pin, and the higher-versioned binary wins on Windows. Glimma is a Suggests.

## MultiScholaR on Windows (report only)

`R CMD INSTALL` of the package at `9dc7bbcd` (the `performance-v2` tree plus the spike files) on `windows-latest`, R 4.6.1, Rtools45 (gcc 14.3.0), hard dependencies: **exit status 0; the package installs and loads.** The four C files compile, including the two built with `-std=c11 -O2 -fPIC -fno-fast-math -Wall -Wextra -Werror`; the DLL links. The log is the artifact `windows-package-install-log`.

One thing does not happen there: R prints `WARNING: this package has a configure script. It probably needs manual configuration`. On Windows R runs `configure.win`, not `configure`, and the package has only `configure`. So `src/build_dia_rsum.R` does not run, and neither the DIA exact-sum DuckDB extension nor its manifest under `inst/duckdb/linux_amd64` is produced (the script would mark the platform unsupported anyway). What the package does at run time without that manifest was not tested. Nothing in the package was changed.

## Design consequences

### What the launcher must set

| | Linux (KIT-050) | Windows (KIT-052) |
|---|---|---|
| Finding R | `RHOME=<bundle>/R` for `Rscript`, whose build prefix is compiled in; the front-end script `R/bin/R` patched once at build time to compute its own home (it carried `/opt/R/4.6.1` in nine places) | nothing: R for Windows finds its home from the executable |
| Library paths | `R_LIBS_USER` and `R_LIBS_SITE` set to folders that do not exist; `R_HOME`, `R_LIBS`, `R_ENVIRON*`, `R_PROFILE*` unset; `--vanilla`. Without the first, an existing user library joins the search path even under `--vanilla`. | the same |
| Shared libraries | `LD_LIBRARY_PATH=<bundle>/syslib` | none |
| Python | `RETICULATE_PYTHON=<bundle>/python/bin/python3`, `PYTHONNOUSERSITE=1`, `PYTHONPATH` and `PYTHONHOME` unset | `RETICULATE_PYTHON=<bundle>\python\python.exe`, the rest the same |
| PATH | `/usr/bin:/bin` | `System32` and the Windows folder |
| Stopping R | R is one process (the launcher, `Rscript` and the front-end all `exec`); started in its own process group and the group signalled | `Rscript.exe` runs `Rterm.exe` as a child: `taskkill /T /F` on the tree; ending the parent alone leaves R running |
| Serving | `127.0.0.1`, a port the shell picks | the same |

The Linux launcher is `tools/dist-spike/linux/Rscript`; the Windows settings are in `main.js` (`rEnvironment`) and repeated in `checks/windows-clean-run.ps1`. The real bundle should have one definition per platform.

### What must be built per platform

- **R.** Linux: Posit's r-builds package for the oldest distribution to be supported, installed at its prefix, packages installed there, then copied and patched (`build/linux-assemble.sh`). Windows: the CRAN installer run silently into the folder (`/VERYSILENT /DIR=...`, no registry or icon tasks); 22 to 35 s. R 4.6.1 was available on both; no substitute was needed.
- **Shared libraries, Linux only.** Everything `ldd` resolves outside the bundle for R, its modules and the package libraries, except glibc: 92 libraries (OpenBLAS, gfortran, ICU 70, readline, the curl and OpenSSL chain, cairo, pango, X11, Tcl/Tk). The bundle then needs from the host only glibc, at least the build distribution's (2.35 for jammy; not tried on any other distribution), and a shell with the core utilities.
- **Packages.** Binaries on both platforms except RUVIIIC, compiled from the CRAN archive source on both; that needs Rtools45 on the Windows build machine (already on the runner image) and a compiler on Linux, at build time only. CRAN's Windows binaries follow the current patch release: `RSpectra`, `progress` and `shiny` warn `was built under R version 4.6.2` on every load. Harmless here, but the versions float with the build day; the real bundle should install from a dated snapshot.
- **Python.** The python-build-standalone archive unpacked, then `pip install --only-binary=:all:` of `mofapy2==0.7.5` and the pins MOFA2 1.22.1 carries for its own basilisk environment (numpy 1.26.4, scipy 1.12.0, pandas 2.2.1, h5py 3.10.0, scikit-learn 1.4.0, dtw-python 1.3.1); all exist as wheels for both platforms. `run_mofa` checks the mofapy2 version against MOFA2's pin, so the build reads the pin from the installed MOFA2 (`build/mofapy2-requirements.R`). reticulate found the standalone build's `libpython` without help.
- **Shell.** electron-builder 25.1.8 with Electron 33.2.0; the bundle goes in as `extraResources`, outside the asar archive.

### What was hard, and what will bite

Nothing blocked. These are the points to carry into KIT-050 and KIT-052:

- **AppImage and FUSE.** The runner image (Ubuntu 22.04) has no FUSE 2. Without it the AppImage only starts with `--appimage-extract-and-run`, which unpacks 1.3 GiB into `/tmp` at every start (7.5 to 13 s, against 2.5 s mounted) and removes it at exit. Read-only locations are not a problem: in the container the bundle was mounted read-only for every check, and in the mounted AppImage R served the Shiny page from the read-only mount.
- **Windows install time.** The silent install of 36,910 files took 262 to 279 s. File count, not size, is the cost; trimming (BH headers, tests inside Python packages, pip, R's own `unins000.exe`, translations) is the lever.
- **Windows path length.** The longest installed path is 202 characters with the user `runneradmin`; 58 characters remain before the 260 limit for a longer user name or product folder.
- **Signing (KIT-052).** electron-builder visits every `.exe` under the bundle for signing (R's, Python's, Tcl's, pip's launchers); unsigned here.
- **R for Windows puts `c:/rtools45/...` at the front of its own PATH**, also under `--vanilla`. If the user has Rtools there, the bundled R sees its tools first. Not a failure here; the launcher cannot remove it through PATH.
- **The Linux patch is of R's own front-end script.** A new R version needs the patch re-checked; the build fails if the build prefix is still named after patching. `Rscript` and `etc/Makeconf` still carry the build prefix (the first is overridden by `RHOME`; the second matters only for compiling).

### Not shown by this spike

- **A clean Windows.** The fresh Windows runner has R 4.6.2, Python 3.12.10, Rtools45 and the Visual C++ runtimes. Pieces 1 to 3 ran with those off the PATH and every steering variable removed, and the installed application reported its own R, but no run was made on a Windows image without them (the Linux pieces were run in a bare container). A Windows Server Core container would be the equivalent test.
- **Other Linux distributions**, and Electron's own needs on Linux (a desktop's GTK and NSS libraries; the sandbox on Ubuntu 24.04 and later). Only `ubuntu-22.04` was run.
- **R when the shell dies.** R stops when the window closes. If the Electron process is killed, nothing here stops R; that needs a watchdog on the R side (KIT-131's ground).
- **MultiScholaR in the bundle.** By the ticket's declaration it is not there. The spike's library is 133 packages; the closure is 413 to 416, so the real bundle will be larger and slower to install than these figures.

## For the user to decide

1. **Oldest Linux to support.** A bundle built on jammy needs glibc 2.35 or later; an older floor means building on an older base.
2. **AppImage as the Linux format**, given the FUSE 2 requirement and the cost of the fallback, or a format that installs files.
3. **Glimma on Windows.** The resolver picks Bioconductor's 2.22.1 over the fork in `Remotes`; if the fork is the intended one, the bundle build must name it explicitly.
4. **`configure` on Windows.** The package installs, but the DuckDB extension step is skipped silently; KIT-052 needs a `configure.win` or a decision that the step is Linux only.
5. **Whether a clean-Windows run is wanted** (Server Core container) before KIT-052.
6. **A dated package snapshot** for the bundle, so the Windows binaries do not float with CRAN's current patch release.

## Notes

- The `R CMD INSTALL` job was added during the spike at the coordinator's request; it is outside the ticket's text ("MultiScholaR itself is not installed in the spike"), which holds for the bundle.
- The check scripts and the hello-world app are spike scripts, not package code: single-file app without modules or bslib, `cat()` for check output.
- `Module CI` also ran on each push of the branch; its runs are not part of this record.
- Artifacts (AppImage, NSIS installer, closure lockfiles, install log) are kept by GitHub for 14 days from the run.
