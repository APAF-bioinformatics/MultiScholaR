# KIT-118 and KIT-119 exit: full suite at `e94b92b1` (2026-11-15)

Same protocol as KIT-120's exit: `tests/testthat` (405 files), one process per file, six in parallel, 1800 s per file, both load modes. Build `browser-smoke.cScT7A` at `e94b92b1` (all of KIT-118 M1, M3, M4 and KIT-119), installed from `browser-smoke.encpwV` (`kit/scripts-k118/verify.sh`). The baseline is KIT-120's exit ([../../kit-120/exit/suite.md](../../kit-120/exit/suite.md), `9a16bfa3`).

`R CMD check --no-tests`: 3 WARNINGs, 3 NOTEs, the same six as the baseline.

| Mode | Files | Passed | Failed | Error blocks | No result |
|---|---|---|---|---|---|
| `load_all`, KIT-120 exit | 404 | 36,769 | 0 | 2 | 0 |
| `load_all`, this run | 405 | 35,005 | 0 | 2 | 2 |
| `load_all`, with the two files rerun alone | 405 | 36,905 | 0 | 2 | 0 |
| check, KIT-120 exit | 404 | 31,710 | 149 | 237 | |
| check, this run | 405 | 31,839 | 149 | 238 | |
| check, after the correction below | 405 | 31,839 | 149 | 237 | |

## Every difference (`kit-105/skipcmp.py`, per file)

- **Two `load_all` files gave no result in the suite run; both pass alone on the same tree and library** (`rerun/`). These were failures of the run, not of the code:
  - `test-applets-direct-shared.R`: the first file of the run compiled the package's C code while five other processes did the same (`ld: cannot find workflow_memory.o`). Alone: 49 passed, 0 failed, as in the baseline.
  - `test-prot-norm-core.R`: stopped at the 1800 s ceiling; the machine was also running the KIT-122 agent's test files. Alone: 1,851 passed, 0 failed, as in the baseline.
- **`test-prot-import-quantity-text.R`, check mode: one new error block, corrected.** KIT-118 M1's new test sources `tools/profiling/workflow_perf107_synthetic.R`, which is not in the installed package. It now skips outside the source tree, as the other corpus tests do (`96129e0c`). Rerun in check mode: 367 passed, 0 failed, 1 skipped, 0 errors; `load_all`: 372 passed.
- **`test-prot-design-save-preflight.R`:** new file (KIT-118 M4, P-10), 8 passed in both modes.
- No other file has fewer blocks run, fewer passes, more failures, more errors or a new skip. The two `load_all` error blocks are the diagnostic ones, as before.

## e2e lanes (`e2e.txt`, one attempt each, same build)

| Lane | Passed | Failed |
|---|---|---|
| lipidomics canonical | 88 | 0 |
| metabolomics LC and GC | 190 | 0 |
| transcriptomics screen | 16 | 0 |
| proteomics DIA | 42 | 0 |
| proteomics DIA_LIMPA | 86 | 0 |
| proteomics enrichment | 98 | 0 |
| proteomics TMT and LFQ | 153 | 0 |

**Decision:** no regression. One defect in a new test found and corrected.

Tables: `skips-load_all.tsv`, `skips-check.tsv`, `header-load_all.tsv`, `rerun/`; `check.log`; `e2e.txt`.
