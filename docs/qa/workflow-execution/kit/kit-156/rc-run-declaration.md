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
