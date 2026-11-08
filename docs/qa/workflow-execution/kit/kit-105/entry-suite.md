# KIT-105 entry-point suite

## Simulation Setup: two lanes alone on the merged tree (before run 2)

Tree `96f80fd8` (W2 to W5 merged plus W6's wiring; `98da5939` adds only a test block to `test-omics-norm-core.R`). R 4.6.1, `load_all`, `MULTISCHOLAR_STRICT_SKIPS=1`, `R_LIBS` = kit_104 (xocrNB). Each lane was run alone through `guiEntryExpectLane()` under testthat, with counting `trace()`s on the new code (script: [lifecycle/lifecycle.R](lifecycle/lifecycle.R)).

| Lane | Expectations | Failed / errors / skips | Elapsed | Traced calls |
|---|---|---|---|---|
| `metab_combined` (two assays through the frame's in-process path) | 35 | 0 / 0 / 0 | 116 s | `mod_omic_norm_server` 1, `runWorkflowOmicsInProcessNormalisation` 1, `runWorkflowOmicsInProcessCorrelation` 1, `collectWorkflowOmicsNormSession` 1, `saveWorkflowOmicsNormSessionArtifacts` 1 |
| `prot_dia_limpa` (F-20's removal; the moved pass-through) | 37 | 0 / 0 / 0 | 257 s | `runProtNormTabEntryWorkflow` 23 (the F-20-free entry, once per tab change), `runDIAReferenceLimpaNormalisationPassThrough` 1 |

So the suite drives the new code through completion, not just its start. Each lane's expectation count equals its share of the file's (`test-gui-entry-metabolomics.R` 105 = 3 x 35; the prot lane's 37 as in KIT-104). The ERROR lines in the logs are the lanes' own refusal steps (import refusal, ITSD absent, session load before export), as in every frozen run.

**Harness record (instrumentation, not results):**
- The first attempt ended in 0 s for both lanes: the script did not `setwd()` into `tests/testthat`, so `test_path()` failed, and the silent reporter hid it. The script was fixed and both lanes were rerun once; the first attempt's logs were overwritten by the rerun.
- In the rerun, the script's printed `identical_to_oracle=FALSE` is a reading defect: `with_reporter()` returns the reporter object, so `isTRUE()` was always FALSE (shown with a one-expectation passing test). The summary reporter's output in each log is the result: 35 and 37 passes, no F, E or S marks, "DONE" with no failure section ([lifecycle/metab_combined.txt](lifecycle/metab_combined.txt), [lifecycle/prot_dia_limpa.txt](lifecycle/prot_dia_limpa.txt)). The script now derives the flag from a `ListReporter`. The lanes were not rerun for it.
