# KIT-216 R2e: predictions, committed before the run (VAL-372 (1); AUDIT-KIT-216-impl-20261126)

The run compares two commits; there is no R CMD check, because there is no `R/` change:
- **Head:** 392efa32, the two edits.
- **Base:** `kit-216-r2c-base` 923cc370 (tests-only plus staging 80f26158).

Two files, one shard per side. Ceiling 60 job-minutes; expected about 29.

| File | Base | Head | Verdict |
|---|---|---|---|
| `test-kit213-study-parameters-labels.R` | 56/0/0/0 | 56/0/0/0 | `same` (E4's printed lines identical) |
| `test-omics-enrich-core.R` | 221/20/0/1 | 222/20/0/1 | `message-changed` (see below) |

Detail for `test-omics-enrich-core.R`:
- The only changed message is the W2 `kernel_clusterprofiler` reader's run `states`: the base digest on base,
  `"<unchanged>"` on head.
- F-16 test 4 goes from 14 to 15 passes (declared).
- Per case, passes are as on base: gprofiler 21, clusterprofiler 19, kernel_gprofiler 20, kernel_clusterprofiler 14.

A miss is a result: record it and hand back.
