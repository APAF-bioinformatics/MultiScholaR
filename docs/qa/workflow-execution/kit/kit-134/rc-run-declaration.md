# KIT-134 release-candidate run on hosted runners (declared 2026-11-20, before launch)

## Purpose

Run on GitHub-hosted runners what KIT-134 still lacks and what runners can show, from `performance-v2` at
`8b4ff691` (KIT-134's packaging, KIT-052's Windows path repairs, KIT-131's Windows worker fix `d5c4e8d1`):

1. Both bundles built by `distributable.yml` (Linux archive, Windows NSIS installer), installed on the build runner,
   lifecycle with the network closed (as stage `full`).
2. **The package's test suite ON THE BUNDLE** (criterion 2.0.9: "the suite per file ... no class D and nothing
   unexplained"): the bundle's own R, library and installed MultiScholaR (`testthat::test_file(...,
   load_package = "installed")`), not `pkgload::load_all`; Linux and Windows.
3. Linux host checks where a runner or a container can stand in: M25 (`/usr/bin/which` absent: plain stop; present:
   full check) and M26 (R resolves the bundle's libraries, not the host's) on ubuntu-22.04, ubuntu-24.04, Arch Linux,
   Fedora, Debian (the existing `hosts` step, unchanged).
4. Windows confirmations left open: the session archive staged under `TMP` on the worker lane (`d5c4e8d1`); the
   Summary `.docx` at project roots of about 90 and 96 with long paths off; `test-workflow-figure-paths.R` on Windows
   (long paths on and off); the four e2e lanes (88, 190, 16, 153).
5. MOFA2 trained in double precision inside the **assembled, installed** bundle on both platforms (the smoke lanes used
   a bundle-shaped folder). `performance-v2` has no integration model yet, so MOFA2 is called directly with KIT-062's
   `diagnose.R` (copied as `packaging/verify/mofa-train.R`): the KIT-137 tier-1 corpus, the package's view
   construction, `use_float32` off.

## Bindings

| | |
|---|---|
| Branch | `run-k134-rc` (run-only), cut from `origin/performance-v2` at `8b4ff691` |
| Commit run | the commit after this declaration's; it changes only `packaging/verify/` and this branch's copy of `.github/workflows/distributable.yml`. Nothing under `R/`, `tests/`, `DESCRIPTION`, `NAMESPACE`, `src/`, `inst/`; no oracle; no harness |
| Workflow | `KIT-134 distributable`, `workflow_dispatch` on `run-k134-rc`, stage `rc` (new: stage `full` plus suite shards, MOFA2, the Windows cases below), platform `both` |
| Jobs | Linux main on ubuntu-22.04 and ubuntu-24.04 (24.04 builds its bundle in an `ubuntu:22.04` container, as before); 5 Linux suite shards (ubuntu-22.04); Windows main (windows-2025); 8 Windows suite shards. Every job builds its own bundle (M20: nothing built with the read token leaves its job) |
| R, Rust, Python | R 4.6.1; Rust stable on Linux, 1.99.1 on Windows (as `8b4ff691`'s workflow; F6); the bundle's CPython and mofapy2 0.7.3 |
| Read token | `DEPS_READ_TOKEN`, exactly as the existing workflow: given only to the library install step; no cache; **no artifact** (the plan job's check is unchanged and still enforces it); results are read from the job logs |
| Instrumentation (new or brought in) | `run-suite.R`, `run-suite-file.R`, `suite-seconds.tsv` (new: the suite driver); `watch-scratch.ps1` (new: the TMP watcher); `mofa-train.R` (KIT-062's `diagnose.R` plus printed MOFA lines); `lane-worker.R`, `worker-probe.R`, `figure-paths.R`, `figure-trace.R`, `windows-figure-cases.ps1` as on `run-win-k131` (`5fb86c18`) |

## The suite on the bundle: how it runs and how it is judged

- **How.** `run-suite.R` deals the 413 files of `tests/testthat` to the shards (longest first, by the seconds of the
  runner reference); each file runs in its own process of the bundle's R: `test_file(<file>, package =
  "MultiScholaR", load_package = "installed")`, helpers sourced from the checkout, `NOT_CRAN=true`. Test tools
  (testthat, sessioninfo, as the run-tests library) in a separate library after the bundle's, checked to shadow
  nothing; `shinytest2` is not installed for the suite (the reference skips the browser lanes for its absence; the
  browser lanes run separately in the main jobs). The checked bundle is moved to a path without spaces for the suite.
  Network as the reference had it (open); the runner's own locale. Per-file timeout 45 minutes; no file starts after
  200 minutes of the step. One `SUITE` line per file and one `SUITEFAIL` line per failing expectation in the log.
- **Known state (Linux).** run-tests run 44438735444 (head `791cc5f3`, the full suite with `load_all`, 20 shards,
  reference BLAS): 413 files, 384 ok, 29 failed (36,850 passed, 800 failed, 11 errors, 559 skipped); the failures are
  the stored-hash and digest expectations recorded on the laptop and a few harness tests. `8b4ff691` differs from
  `791cc5f3` only by `kit-131-win-worker` (run-tests 44468417325 and 44479976606: no new failure; figure-paths 171/0)
  and two records. A file failing on the bundle with the same failing tests as in the reference is "known", not a
  regression.
- **Classification of every other difference** (revision 2.0, 2.0.4): E (environment: installed package instead of
  `load_all`, e.g. a test reading `R/` sources or calling `pkgload`; a tool or package present on one side only), V, C,
  T, L (as 2.0.4), or D (defect). The bundle computes with its own OpenBLAS (M26), the reference with reference BLAS:
  a numeric difference behind a stored digest is expected and classified from its message, not called a regression.
  What the log alone cannot classify is listed as **open**, not passed.
- **Windows.** There is no Windows runner reference (run-tests refuses Windows). Windows files are set against the
  Linux reference; platform causes (GNU tools, `bash`, `systemd`, path forms) are class E with the message as evidence;
  anything else is listed open.

## Windows cases (main job; roots set through `TMPDIR` as before)

| # | Item | Pass |
|---|---|---|
| W-A | `lipid_canonical` on plain workers (`lane-worker.R`), `TMP`/`TEMP` = a fresh folder watched by `watch-scratch.ps1` | the lane runs to its end; identical to the in-process lane of the same job except the store's entity counts (accepted at Linux step 8) and the run-dependent `norm_log`/`session_summary` text (44472789530); **`archive-stage-*` folders seen made directly in the watched `TMP` (depth 0) during the lane** = the workers' session archive staged under `TMP`. If none appears at depth 0, not met (if they appear only at depth 1, the archive was staged in process) |
| W-B | the four browser e2e lanes | 88, 190, 16, 153 |
| W-C | prot_tmt entry lane, long paths off, roots about 90 and 96 (TMPDIR 31, 36) | the Summary `.docx` rendered (`RENDER` with a `.docx`, "Report generated successfully!"), every figure present |
| W-D | `test-workflow-figure-paths.R` (item `figure_tests`), long paths on (root about 102) and off (about 63) | on: 0 failed, 0 errors; off: 0 failed; errors only from fixtures that make folders over 248 characters (as 44457268709: 92 passed, 7 errors); every printed error of that kind |
| W-E | MOFA2 in the installed bundle | trains with the bundle's Python (no basilisk), double precision; 186 iterations, ELBO -27000.322787, planted correlations 0.9937 / 0.9936 / 0.9892 as KIT-062's double-precision runs (44466163909); factors within the 1e-6 rule of those runs (from the printed hex values) |

Linux main jobs: lifecycle, hosts, GUI-entry lanes, e2e lanes (22.04) as stage `full`, plus MOFA2 (pass as W-E).

## Runtime ceiling and attempts

- Every job 300 minutes (the workflow's); suite step 220 minutes (files stop starting at 200); MOFA2 step 30; worker
  step 60; figure steps 120 each.
- **One run plus one bounded correction pass**, the pass only for an instrumentation fault (a script error, a job that
  never reached its step, a case that never started). An unfavourable valid result is a finding, kept, and not rerun.

## Decision after the run

Per criterion of KIT-134 2.0.9: met, met in the weaker single-job form, not met, or not runnable on runners (with why);
the suite's per-file table with every difference classified or listed open; the Windows cases met or not met; the
MOFA2 result. Results below this line, in `kit-134/outcome.md` on this branch, and as dated Status lines in the
tickets. Nothing is ticked.
