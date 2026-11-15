# Final tree of the GUI journey (KIT-118, KIT-119, KIT-122): full suite at `ed99a0a8` (2026-11-16)

Same protocol as KIT-118's exit ([../../kit-118/exit/suite.md](../../kit-118/exit/suite.md), `e94b92b1`), which is the baseline here. Build `browser-smoke.y5ordS` at `ed99a0a8`. Nothing but KIT-122, its correction and the check-mode skip of one KIT-118 test (`96129e0c`) lies between the two.

`R CMD check --no-tests`: 3 WARNINGs, 3 NOTEs, the same six.

| Mode | Files | Passed | Failed | Error blocks | No result |
|---|---|---|---|---|---|
| `load_all`, KIT-118 exit (with its two reruns) | 405 | 36,905 | 0 | 2 | 0 |
| `load_all`, this run | 405 | 36,909 | 0 | 2 | 0 |
| check, KIT-118 exit | 405 | 31,839 | 149 | 238 | |
| check, this run | 405 | 29,973 | 149 | 237 | 1 |
| check, with the one file rerun alone | 405 | 31,842 | 149 | 237 | 0 |

## Every difference (`kit-105/skipcmp.py`, per file)

- **`test-prot-norm-core.R`, check mode: no result in the suite run; alone on the same build 1,869 passed, 0 failed, 1 skipped** (`rerun/`). Its log ends 80 s before the run finished: the 1,800 s per-file ceiling. KIT-122 added two full drives of the tab to this file. It now takes 1,612 s in `load_all` under six processes and 1,570 s alone in check mode, so it will keep hitting the ceiling under load. Not changed here.
- **`test-prot-05b-norm-module-contracts.R`, both modes:** 476 passes to 439, skips 14 to 22. Eight blocks on the deleted functions are quarantined (Q-187); assertions on them removed from live blocks (Q-188).
- **`test-prot-05d-norm-module-characterization.R`, both modes:** 100 passes to 82 (Q-188).
- **`test-omic-norm-module.R`, check mode:** one new skip. The KIT-122 F6 test reads the package sources and runs under `load_all` only, where it passes.
- **`test-prot-import-quantity-text.R`, check mode:** the error block of KIT-118's exit is now a skip (`96129e0c`).
- No new failure and no new file-level skip. The two `load_all` error blocks are the diagnostic ones, as before.

## e2e lanes (`e2e.txt`, one attempt each, same build)

All seven pass with the counts of KIT-118's exit: lipidomics 88, metabolomics 190, transcriptomics 16, DIA 42, DIA_LIMPA 86, enrichment 98, TMT and LFQ 153; no failure.

**Decision:** no regression. One file is at the suite's time ceiling.

Tables: `skips-load_all.tsv`, `skips-check.tsv`, `header-load_all.tsv`, `rerun/`; `check.log`; `e2e.txt`.
