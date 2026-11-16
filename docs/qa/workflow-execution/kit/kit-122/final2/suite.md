# KIT-122 revision 1.1 (the composite figure on the worker): full suite at `ea236689` (2026-11-16)

Same protocol as the final tree of the GUI journey ([../final/suite.md](../final/suite.md), `ed99a0a8`), which is the baseline here. Build `browser-smoke.26ynxx` at `ea236689`. Between the two: revision 1.1 (the composite built by the worker, the worker-figures drives moved to a file of their own) and evidence files.

`R CMD check --no-tests`: 3 WARNINGs, 3 NOTEs, the same six.

| Mode | Files | Passed | Failed | Error blocks | No result |
|---|---|---|---|---|---|
| `load_all`, baseline | 405 | 36,909 | 0 | 2 | 0 |
| `load_all`, this run | 406 | 34,907 | 0 | 2 | 2 |
| `load_all`, with the two files rerun alone | 406 | 36,954 | 0 | 2 | 0 |
| check, baseline (with its rerun) | 405 | 31,842 | 149 | 237 | 0 |
| check, this run | 406 | 31,887 | 149 | 237 | 0 |

## Every difference (`kit-105/skipcmp.py`, per file)

- **`test-prot-norm-worker-figures.R`, new file, both modes:** 63 passed. It holds the two worker-figures drives that were in `test-prot-norm-core.R`, with the composite's assertions.
- **`test-prot-norm-core.R`:** 28 blocks to 26; 1,872 passes to 1,854 (alone). Passes in both modes: +63 in the new file, -18 here, +45 in total, which is the whole change in each mode's count.
- **`test-prot-norm-core.R` and `test-gui-entry-proteomics.R`, `load_all`: no result in the suite run; alone on the same tree 1,854 and 193 passed, 0 failed** (`rerun/reruns.tsv`). Two simulator pilots (4 cores) ran beside the suite's six processes; the same files gave results in check mode in the same run. A load effect, not an unfavourable result; both attempts are kept.
- No new failure, no new skip. The two `load_all` error blocks are the diagnostic ones, as before.
- The check-mode aggregation stopped at its header (`git rev-parse` in the check directory, as in every earlier run); the per-file table was written before it.

## e2e lanes (`e2e.txt`, one attempt each, same build)

All seven pass with the baseline's counts: lipidomics 88, metabolomics 190, transcriptomics 16, DIA 42, DIA_LIMPA 86, enrichment 98, TMT and LFQ 153; no failure.

**Decision:** no regression.

Tables: `skips-load_all.tsv`, `skips-check.tsv`, `header-load_all.tsv`, `rerun/`; `check.log`; `e2e.txt`.
