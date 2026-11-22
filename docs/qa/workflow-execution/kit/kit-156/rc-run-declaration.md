# KIT-156 bundle runs on run-k156-rc (declared 2026-11-22, before any run)

## Purpose

The test build for data testing: `performance-v2` 12d68ca7 plus KIT-156 (`kit-156-data-testing` 0c4c2abf: the D-1
repair, the test setting `MULTISCHOLAR_TEST_FEATURES=1` in the packaged entry, the lifecycle run twice). On both
platforms: build, install, lifecycle (without and with the test setting), the installer suite on the installed
package, the end-to-end journeys, and with the test setting a headless smoke of the RUV design panel and the
integration module end to end on the synthetic corpus (import, match, train, results, record). Compared with KIT-134's
classified results of run 44598642342 (`run-k134-rc`, `rc-run-declaration.md` there).

## Bindings

| | |
|---|---|
| Branch | `run-k156-rc` (run-only), cut from `kit-156-data-testing` 0c4c2abf (= `origin/performance-v2` 12d68ca7 + 3 KIT-156 commits) |
| Instrumentation | 141e7b45 (= d7e51f94 cherry-picked: stage suite, `MSR_SUITE_FILES`, `MSR_SUITE_TMP`) and 71669d61 (KIT-134 rc's verify files from 72fca08c; stage rc without its MOFA2 step; part smoke; stage smoke; the vector-instruction count). Only `packaging/verify/` and this branch's workflow copy |
| Workflow | `KIT-134 distributable`, `workflow_dispatch` on `run-k156-rc`, platform `both` |
| Read token | `DEPS_READ_TOKEN` exactly as the existing workflow (library install step only); no cache; no artifact (the plan job's check unchanged); results read from the job logs |
| R, Rust, Python | R 4.6.1; Rust stable on Linux, 1.99.1 on Windows; the bundle's CPython and mofapy2 0.7.3 (as `performance-v2`'s workflow) |
| Smoke files | `packaging/verify/k156/smoke-files.txt`: `test-test-features.R` (the entry's own test-setting block evaluated, then the home page and each omic's normalisation tab rendered), `test-integration-record-mofa.R` (import, match with an uploaded table, MOFA2 training, results, record, report, replay), `test-integration-results-mofa.R`, `test-integration-gui.R`, `test-omic-norm-ruv-design.R`, `test-kit129-run-path.R`, `test-trx-ruv-designs.R`; each in its own process of the bundle's R against the installed package, `MULTISCHOLAR_TEST_FEATURES=1` |

## Runs, ceiling and attempts

1. **Lifecycle run** (instrumentation proof): stage `smoke`, input `files` = the smoke files plus
   `packaging/verify/lifecycle/test-lane-lifecycle.R` (one pass, one failure, one error, one skip on purpose, and a
   leftover folder). Two jobs (one per platform), each builds its own bundle. Pass: on each platform the probe's four
   kinds and its leftover are reported (`SUITE`, `SUITEFAIL`, `SUITETMP`), every smoke file reports a `SUITE` line,
   `test-test-features.R` sees the variable. Ceiling: job 300 min, suite step 220 min.
2. **Counted run**: stage `rc`: Linux main on ubuntu-22.04 and ubuntu-24.04 (build, install, lifecycle off and on,
   hosts, GUI-entry lanes, e2e lanes on 22.04, vector-instruction count), 5 Linux suite shards, Windows main
   (build, installer, install, lifecycle off and on, GUI-entry lanes, e2e lanes, KIT-131 worker lane, figure cases
   with long paths on and off), 8 Windows suite shards, one smoke job per platform. Ceiling: job 300 min; suite step
   220 min (no file starts after 200); 45 min per file.
3. One correction pass of the lane (instrumentation faults only). An unfavourable valid result is a finding, kept and
   not rerun.

## Judgement (fixed now)

- The suite and the lanes per file against run 44598642342's classified table (Linux: against its reference, as
  there; Windows: against its Windows table). Every new difference gets a class of KIT-134 2.0.4 (Agree, E, V, C, T,
  L, D); a class D small and inside KIT-156's files is repaired, else listed.
- D-1 on Windows: `test-workflow-figure-paths.R` with no failure in the suite shard, and in the figure cases (W-D),
  where run 44598642342 had 1 (D-1); the new KIT-156 test passes.
- Test setting: the lifecycle passes off and on on every main job (on: R reports both options TRUE and the served home
  page's Integration tile is live; off: no report and the placeholder tile).
- Smoke: every smoke file runs with its expectations passing on both platforms except failures that also fail on the
  same file in the plain suite shards (then classified there).
- Portability: the counts of AVX-512 (zmm), AVX/AVX2 (ymm) and FMA instructions in `fastlimpadpc`, `fastduckdb`
  and MultiScholaR's compiled libraries, with the runner's CPU; a non-zero AVX-512 count in `fastlimpadpc` means the
  Linux bundle stops with "illegal instruction" on CPUs without AVX-512.

Results below this line after the runs.

## Lifecycle run 44952988256 (stage smoke, a10f24fe): met on both platforms

The probe reported each kind (1 pass, 1 failure, 1 error, 1 skip) and its leftover (`SUITETMP`), every smoke file gave a
`SUITE` line, and `test-test-features.R` saw the variable (41/0 on both). Smoke results (Linux / Windows):
test-test-features 41/0 / 41/0; test-integration-gui 395/0 / 395/0; test-integration-record-mofa 27/0 / 27/0;
test-integration-results-mofa 46/0 / 46/0; test-omic-norm-ruv-design 1153/0 / 1153/0; test-trx-ruv-designs 344/6 /
344/6; test-kit129-run-path 687/6/0 err / 662/6/4 err. The 12 failures on both platforms are one expectation, the
transcriptomics corpus manifest's sha256 against the stored value (to be set against the run-tests reference of the full
suite); the 4 Windows errors are `system2("cp", ...)` in the test's replay helper (E: GNU tools absent on Windows). Every
file's temporary root was empty afterwards. Instrumentation proven; the counted run follows (stage rc).

## Results: counted run [44956911834](https://github.com/APAF-bioinformatics/MultiScholaR/actions/runs/44956911834) on 6fd6cfa (2026-11-21/22 UTC)

Read from the job logs (0 artifacts, no cache). Conclusion "failure" comes only from the three main jobs' GUI-entry lane
steps (below), as in run 44598642342. Every suite and smoke job passed its step. No instrumentation fault; the
correction pass is not used.

### Main jobs

| | ubuntu-22.04 | ubuntu-24.04 | windows-2025 |
|---|---|---|---|
| library, internal packages, shell, archive / installer, install | passed | passed | passed |
| lifecycle, test setting off (no report; placeholder tile) | passed, process ("kernel 6.8") | passed, process ("D-Bus") | passed, process ("win32") |
| lifecycle, test setting on (both options TRUE; tile live) | passed | passed | passed |
| hosts (M25, M26; Arch, Fedora, Debian) | passed | passed | - |
| browser e2e lanes | 88, 190, 16, 42, 86, 98, 153 = laptop | not run (by design) | 88, 190, 16, 153 = laptop |
| GUI-entry lanes | every lane differs from its oracle in the steps after the omic writes its integration input (`integration/<omic>/*`, KIT-147..153) and, trx_screen, from the raw state on (KIT-128 technical replicates); metab as before (fixture, defect 4) | the same | the same families plus the earlier Windows-only ones (summary text root, fragpipe DA table, lipid norm_log text, trx_screen testServer error) |
| KIT-131 worker lane (W-A) | - | - | as run 44598642342 (identical counts) |
| figure paths W-D | - | - | on 159/0/0 (was 148/1/0: D-1); off 127/0/7 (was 116/1/7: D-1; the 7 errors the known fixtures over 248 characters) |
| Summary at roots 90, 96, long paths off (W-C) | - | - | `.docx` written at both; lane 33/4 (the 4 = integration-input artifacts + summary text root) |

The GUI-entry differences are performance-v2's own: run-tests 44953187101 (load_all, reference BLAS) fails the same files
with the same counts on base 12d68ca7 and head (gui-entry lipid 4, metab 30, proteomics 22, transcriptomics 13). They are
oracles not re-recorded since KIT-147..153 and KIT-128 (rows for the single re-record, REC-049); not bundle differences.

### Native builds (item 4)

| runner CPU | fastlimpadpc.so zmm (AVX-512) / ymm / FMA | fastduckdb.so | MultiScholaR.so |
|---|---|---|---|
| AMD EPYC 9V74 (avx512f yes), ubuntu-22.04 | **1522** / 629 / 0 | 26 / 208 / 4 | 0 / 0 / 0 |
| AMD EPYC 7763 (avx512f no), ubuntu-24.04 | 0 / **2730** / 0 | 26 / 208 / 4 (identical) | 0 / 0 / 0 |

`fastlimpadpc` takes the instruction set of whichever runner built it (`-C target-cpu=native`, src/Makevars.in:39 at
87a90b89): the 22.04 bundle needs AVX-512 and stops with "illegal instruction" on most desktop CPUs (Intel consumer
chips, AMD before Zen 4); the 24.04 bundle needs AVX2. `fastduckdb`'s counts are the same on both builds (no CPU-specific
flag; its AVX-512 paths are DuckDB's own run-time dispatch). Windows (Makevars.win.in) has no such flag.

### Suite on the bundle (installed package), against run 44598642342

Linux 447 files: 405 ok, 42 failed (was 383/30 of 413). Windows 447: 384 ok, 60 failed, 3 timeout (was 357/53/3 of 413).

| Class | Files (L = Linux, W = Windows) |
|---|---|
| Fixed | workflow-figure-paths (W: D-1 gone, 159/0), integration-corpus (W: 111/0, was 39/2/1) |
| Known state of performance-v2 (fail identically in run-tests 44953187101 on base and head) | gui-entry-lipidomics, -proteomics, -transcriptomics; omic-summary-module; omics-summary-core; omic-da-module; workflow-omics-owner; workflow-operation-registry; module-ci-trx-design and -trx-qc (error); dia-limpa-ruv-rules, kit127-omics-ruv-factors, prot-ruv-factors-in-model (new files); both platforms |
| Numeric, B-01 family (bundle OpenBLAS; stored sha256 of the generated transcriptomics corpus manifest, whose truth.json doubles move in the last bits; passes under reference BLAS) | integration-input-transcriptomics (3), kit129-run-path (6), trx-ruv-designs (6); L and W; as run-win-k153 |
| E, GNU tools absent (Windows) | kit129-run-path 4 errors (`system2("cp")` in the test's replay helper) |
| Open, not classified (Windows only; not KIT-156 code) | omic-enrich-module +3 over run 44598642342 (kernel_gprofiler and refusals text digests; a second run's display text differs from the first) |

No class D. No new difference attributable to KIT-156's changes.

### Smoke with the test setting (both platforms, identical to the lifecycle run)

test-test-features 41/0; integration-gui 395/0; integration-record-mofa 27/0 (import, uploaded match, MOFA2 training,
results, record, report, replay); integration-results-mofa 46/0; omic-norm-ruv-design 1153/0; trx-ruv-designs 344/6 and
kit129-run-path 687/6 (Windows 662/6/4 errors): the B-01 manifest checks and the `cp` errors above. Temporary roots empty.
