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

## Results: run [44598642342](https://github.com/APAF-bioinformatics/MultiScholaR/actions/runs/44598642342) on `72fca08c` (2026-11-19 19:28 to 21:47 UTC)

Read from the job logs (no artifact; 0 uploads, no cache). Conclusion "failure" comes from three main jobs, each failing
only at its GUI-entry lanes step, whose differences are the known ones (below). **No instrumentation fault; the
correction pass is not used.** All 13 suite shards and every other step passed. Builds took about 11 minutes per job.

### Main jobs

| | ubuntu-22.04 | ubuntu-24.04 | windows-2025 |
|---|---|---|---|
| library loads from the moved folder, network off; internal packages | passed | passed | passed |
| shell, archive / NSIS installer, install (no admin), path lengths | passed | passed | passed |
| lifecycle, installed, network closed | passed, plain launcher | passed, plain launcher (D-Bus reason) | passed, plain launcher ("win32: no systemd") |
| hosts: M25 without `which` (Arch, Fedora: exit 69, message, R not started); full check before and after M26 on the runner, Arch, Fedora, Debian | 8/8 on all four hosts, both phases | 8/8 on all four hosts, both phases | - |
| GUI-entry lanes | 14 identical (prot x5, lipid, trx_screen x8); metab_lc/gc/combined 10 steps differ (defect 4, the 3-feature cyclic-loess fixture; unchanged) | the same | prot_tmt, prot_lfq 35/2 (summary text root, harness); fragpipe 31/6; metab as Linux; lipid 27/8 (norm_log text); trx_screen error under testServer: all as in 44379470655 |
| browser e2e lanes | 88, 190, 16, 42, 86, 98, 153 = laptop | not run (by design) | **88, 190, 16, 153 = laptop** |
| MOFA2, double precision, installed bundle, network closed | trained | trained | trained |

### MOFA2 in the assembled, installed bundle (W-E)

| | iterations | ELBO | planted abs r | Python, MOFA2 |
|---|---|---|---|---|
| ubuntu-22.04 (AMD EPYC 9V74) | 206 | -24285.131276 (`-0x1.7b74866d265eep+14`) | 0.9979 0.9987 0.9893 | the bundle's `python3`; MOFA2 1.22.0 from the bundle's library; mofapy2 0.7.3; `use_float32` off; BLAS the bundle's OpenBLAS |
| ubuntu-24.04 (AMD EPYC 7763) | 206 | -24285.131276 (`...265efp+14`) | same | same |
| windows-2025 (Xeon 8573C) | 206 | -24285.131276 (`...265fp+14`) | same | the bundle's `python.exe` |

Factors (20 x 5) across the three: max abs difference 5.3e-15, max rel 3.9e-14; ELBO equal to 1 to 2 ulp. **MOFA2 trains in
the assembled bundle on both platforms, in double precision, with the bundle's Python and no network, and the three
machines agree to about 1e-14.** The declared numbers (186 iterations, ELBO -27000.322787, 0.9937/0.9936/0.9892) are
**not comparable**: they are of the corpus regenerated on branch run-mofa-smoke-k062's package commit (`e1ab6482`, KIT-062
prerequisite b), while `8b4ff691` carries KIT-137's corpus `5fd98172` (20 files differ). That was an error in the
declaration's reference, not a result; the rule's own test (same inputs, same model across machines) is met.

### Windows cases

| # | Result |
|---|---|
| W-A worker lane | lipid_canonical on plain workers ran to its end: 9 managed waits, 38 worker processes, 0 left after close; 19 worker launches including `lipid_session_export`. Against the in-process lane of the same job: only the store's entity counts and the run-dependent `norm_log`/`session_summary` text (as 44472789530). **Watcher:** 20 `wfx-*` folders made directly in TMP (depth 0, 28 to 30 characters): the workers draw their figures under TMP (d1fa3ada) - shown. **No `archive-stage-*` folder appeared at depth 0 or 1.** Because such folders are left behind on Windows (D-1 below), one could not have been missed: the lipidomics session export makes no native tar archive, so the d5c4e8d1 archive path is **not reachable on a Windows lane** (the native archive belongs to the DIA routes, D4). It is shown only by `test-workflow-figure-paths.R` ("in a worker the session archive stages its members in the system's temporary folder"): the members are copied into the TMP staging folder (expectations pass), but the folder is not removed (D-1) |
| W-B | e2e 88, 190, 16, 153: met |
| W-C | long paths off: root 90: Summary `.docx` written, "Report generated successfully!", 121/121 figures, lane 35/2 (the two known mask failures). Root 94 (target 96): `.docx` written, 121/121 figures; the lane has 21 failures, all `enrichment_analysis-stringdb_plot` text digests (not classified; not a Summary or figure failure). **Met for the Summary at 90 and 94** |
| W-D | long paths on (root 102): 148 passed, 1 failed, 5 skipped; off (root 63): 116 passed, 1 failed, 7 errors (fixtures making folders over 248 characters, and the known "invalid export publication request"), 5 skipped. The 1 failure in both is D-1 |

**D-1 (class D, Windows, package code, not changed):** `runDIANativeTarR` (`R/utils_workflow_dia_sessions.R`) registers
`on.exit(unlink(staging))` before `setwd(staging)` / `on.exit(setwd(previous))`, so the unlink runs while the working
directory is still the staging folder. Windows cannot delete the current directory: every R-fallback native archive
(every DIA session export on Windows, and the session archive refresh of the Summary save) leaves a full copy of the
archive's members in TMP. Linux deletes it. Repair: restore the directory before the unlink (order of the two `on.exit`).
Present since 767ad7d6 (tempdir staging), not introduced by d5c4e8d1.

### The suite on the bundle (413 files per platform, every file run; installed package, not load_all)

**Linux** (5 shards, ubuntu-22.04, the bundle's OpenBLAS): 383 ok, 30 failed, 0 crashed, 0 timeout; 36,621 passed, 847
failed, 13 errors, 565 skipped (reference 791cc5f3 with load_all and reference BLAS: 384 ok / 29 failed; 36,850 / 800 /
11 / 559). Per file against the reference: 27 files differ in counts.

| Class | Files | Evidence |
|---|---|---|
| Known state (same failing tests as the reference) | prot-norm-core (272), workflow-saved-da (171), trx-kernels, trx-publication, workflow-omics/trx-declarations, gui-entry-metabolomics, module-ci-trx-*, omic-summary-module (+1: trx_technical RUV-III matrix digest), omic-enrich-module, omics-summary-core, trx-design, kit-harness, checkpoint-harness, gui-diagnostic, perf101, sleep-guard, tech-reps-limma, saved-exports, filemgmt-load-dependencies, ruv-pearson | messages as in the reference's failures.txt |
| E, installed package instead of load_all: srcref-guarded expectations do not run | omics-enrich-core 241 to 206, omics-summary-core, omics-qc-core, omic-qc-module, omic-norm-module, omic-gui-config, omics-da-core, trx-norm-core, omics-norm-core, prot-design/import/qc-core, perf107-audit (fewer passes, no failure) | the tests' own guard "source references are kept under load_all only" |
| E, installed namespace is locked | general-heatmap-support-direct-shared, prot-07d-da-model-stats-shared | "cannot add bindings to a locked environment" |
| E, helper not visible from namespace code | prot-import-quantity-text | `workflowCheckpointSequences` (tests/testthat/helper-workflow-checkpoint-fixtures.R) not found |
| Numeric under OpenBLAS: digests of DA tables (values not compared here) | omic-da-module +20: lipid/metab no_friendly and managed, `run rows` digest and `de_*_long_annot.tsv` sha | the B-01 family; reference BLAS on the runners reproduces the laptop's digests, the bundle's OpenBLAS does not. Expected; to be value-compared under 2.0.4 before acceptance |
| Numeric, within the rule: last-bit differences | workflow-nipals-identity 40 failures | MultiScholaR's NIPALS against mixOmics on the same machine: values equal to about 1e-15 relative. **The "bit-identical to mixOmics" contract (KIT-087) holds under reference BLAS only**, not under the shipped OpenBLAS; a note for the methods ledger, not a defect |
| Improved (fail in the reference, pass here) | gui-entry-proteomics, gui-entry-transcriptomics, optional-packages, ruv-pearson (1 error fewer) | |
| Count changes from the tree (8b4ff691 vs 791cc5f3) | workflow-figure-paths 142 to 171, workflow-execution-process 185 to 191, dia-products 191 to 192 | kit-131-win-worker |

**No class D on Linux.** Open: the OpenBLAS digest differences have not had their values compared.

**Windows** (8 shards): 357 ok, 53 failed, 3 timeout; 33,241 passed, 1,048 failed, 49 errors, 567 skipped. No Windows
reference exists. 25 files fail as in the Linux reference (same families, some counts lower). 28 files fail on Windows
only, and the 3 timeouts:

| Class | Files |
|---|---|
| **D-1** | workflow-figure-paths (1) |
| E, checkout line endings (git autocrlf on the Windows image: fixture and corpus file hashes; inferred, not shown) | module-ci-fixture-integrity (161), e2e-manifest (57), workflow-omics-import (18), integration-corpus, trx-design (fixture sha), gui-entry-lipidomics (norm_log text) |
| E, GNU tools absent in the shell's environment (`cp`, `bash`, `git`; MULTISCHOLAR_DISABLE_TOOLS) | workflow-lipid-declarations, workflow-metab-declarations, workflow-analysis-record, perf107-audit, ci-impact-source-map, gui-passive-diagnostic |
| E, D4: no DIA query runtime or direct writer on Windows | perf-106, perf107-consumers, perf106-repair-v2, workflow-partitions, workflow-core-harness (worker phase hook); timeouts checkpoint-harness, perf101-diagnostic, checkpoint-dia-state (inferred: DIA checkpoint harness) |
| E, Windows path forms in expected values | workflow-omics-sessions, run-packaged-application (TMPDIR form), metab-00d-import-module-characterization, omic-norm-module (path in a message, digests) |
| E, installed package (as Linux) | general-heatmap-support-direct-shared, prot-07d-da-model-stats-shared, prot-import-quantity-text |
| Open, not classified | workflow-native-session ("native archive creation needs an output path" from the R fallback; one connection error), workflow-execution ("missing platform counters" expected class not raised; one connection error), module-ci-prot-enrich (zip cannot be opened), omic-da-module (+3 over Linux), gui-entry-transcriptomics (testServer session destroyed, known) |

### Correction pass

Not used. No job failed from an instrumentation fault. The missing evidence (W-A's archive staging on a lane) cannot be
produced on Windows by any instrumentation, and MOFA2's comparison failed only through the declaration's reference corpus.
