# KIT-105 run 8: the skip comparison

The per-file runner (`tools/testing/test-baseline.R suite`) over the whole of `tests/testthat` on `98da5939`, in both load modes. `load_all` uses the pinned source worktree. Check mode uses the installed `kit_105` check build. The comparison is against KIT-104's `skips-load_all.tsv` and `skips-check.tsv`, by [skipcmp.py](skipcmp.py), with full output in [skip-compare.tsv](skip-compare.tsv). Rows: [skips-load_all.tsv](skips-load_all.tsv), [skips-check.tsv](skips-check.tsv).

| Mode | Files | Passed | Failed | Error blocks |
|---|---|---|---|---|
| `load_all`, KIT-104 | 390 | 31,530 | 0 | 2 |
| `load_all`, KIT-105 | 393 | 30,645 | 0 | 2 (the same two known diagnostics) |
| check, KIT-104 | 390 | 22,726 | 167 | 394 |
| check, KIT-105 | 393 | 22,560 | 149 | 393 |

**Every change against KIT-104 is accounted for:**
- **Recorded quarantines:**
  - `test-lipid-12` is a whole-file quarantine (Q-95). It skips at the file level under `load_all`. Under check it already had a file-level skip, its pre-existing source-file guard, and now also carries the Q-95 category.
  - `test-lipid-12a`, blocks 4 to 8 (Q-96).
  - `test-metab-03b`, whole file (Q-105). The explicit line-1 skip replaces its silent symbol loader.
  - `test-metab-03c`, whole file (Q-106). In check mode this also removes its 4 error blocks.
  - `test-metab-03d`, blocks 9 to 11 (Q-107).
  - The 14 comparison blocks of the new `test-omics-norm-core.R` (Q-86 to Q-89, Q-108 to Q-111).
- **New files:**
  - `test-omic-norm-module.R`: 525/0/0/0 under `load_all`, 520/0/0/0 under check.
  - `test-omics-norm-core.R`: 62/0/0/14 under `load_all`, 56/0/0/14 under check.
  - `test-prot-norm-core.R`: 57/0/0/0 under `load_all`, and 35 passed, 3 errors, 1 skip under check. The check skip is its source-reference block, which runs under `load_all` only, as `ca35c462` set for KIT-104's new files. The 3 check errors are the known source-tree read: the shared helper `helper-workflow-checkpoint-fixtures.R:12` reads `tests/testdata/workflow-checkpoints/scenarios.json` from the repo root ("cannot open the connection", reproduced). Its source file `test-workflow-checkpoint-dia-state.R` has the same errors in check mode at KIT-104 and KIT-105 (14 each). The new files' lower check-mode pass counts come from expectations guarded to the source layout. No check-mode failure or skip is new in them.
- **A retargeted block:** `test-prot-05b` goes from 504 to 502 passes. This is W3's F-20 block. Its two cases for `pending`/`in_progress` on `peptide_count_filtered` cannot occur since KIT-104's P-16, so they were dropped. The ticket specifies the remaining cases, which are kept. Nothing is skipped or hidden.
- **Improved:** in check mode, `test-omic-qc-module.R` (17 failures → 0) and `test-omics-qc-core.R` (1 → 0). These are KIT-104's `ca35c462` guards, which postdate KIT-104's recorded check rows.

**Decision:** no file has a new file-level skip, more hidden blocks, fewer tests run or a new skip category without a Q-entry, in either load mode. The one count drop without a Q-entry is the retargeted F-20 block, which the ticket specifies. The check-mode aggregate stopped at its header step (`git rev-parse` inside the check directory, as in KIT-102). Every per-file row was written, and the TSV above is assembled from them.
