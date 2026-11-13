# Checkpoint suite after KIT-121, KIT-117 and KIT-114 (2026-11-13)

Full suite, both modes, one process per file, six in parallel, 1800 s per file, on `da5b9c82` (branch `gui-journey-20261113`), build `browser-smoke.ZzEGOQ`. Baseline: KIT-120's exit (`kit-120/exit/`).

| Mode | Passed | Failed | Error blocks | No result |
|---|---|---|---|---|
| `load_all`, KIT-120 exit | 36,769 | 0 | 2 | 0 |
| `load_all`, this run | 36,775 | 0 | 2 | 0 |
| check, KIT-120 exit | 31,710 | 149 | 237 | |
| check, this run | 31,716 | 149 | 237 | |

`kit-105/skipcmp.py` per file: **no difference** in either mode beyond the pass counts. The +6 are KIT-121's assertions (`test-workflow-omics-owner.R` +5, `test-workflow-gui-memory-lifetimes.R` +1). `R CMD check --no-tests`: 3 WARNINGs, 3 NOTEs (baseline). No regression.
