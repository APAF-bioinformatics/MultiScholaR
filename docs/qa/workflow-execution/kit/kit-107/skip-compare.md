# KIT-107 run 9: the skip comparison

The set is run 1's 20 files, in both load modes. The baseline is KIT-105's full-suite rows overlaid with KIT-110's and KIT-106's, compared with `kit-105/skipcmp.py`. Files outside the set are filtered out. The full comparison is in [skip-compare.tsv](skip-compare.tsv).

- **`load_all`** ([skips-load_all.tsv](skips-load_all.tsv)) is run 1, under strict skips: 7,476 passed, 0 failed, 0 error blocks.
- **Check mode** ([skips-check.tsv](skips-check.tsv)) ran on the `kit_107` check build (`66sPmE`, `8b92d7ac`).

**Every change against the baseline is one of these:**
- **A recorded quarantine:**
  - Q-163 and Q-164: whole files, each with a line-1 skip;
  - Q-165: MCI-010.6, so MCI-010 runs 92 → 79 expectations;
  - Q-166: the comparison side in the new core test.
- **A new file:** `test-omic-enrich-module.R` (975/0) and `test-omics-enrich-core.R` (223/0/1 under `load_all`; 175/0, 1 error, 1 skip under check).
- **Check mode only, the known source-tree read:**
  - `test-omics-enrich-core.R` D6 test 4: "cannot open the connection" on a repo-root fixture;
  - `test-workflow-analysis-record.R` 16 → 17 error blocks: the new F-16 replay block sources `tools/profiling/workflow_kit.R`, which the check layout does not ship. The file's existing 16 errors are the same reads.

  Both run under `load_all`, where they pass.

**Also:** check-mode failures in this set fall from 149 to 0.

**Decision:** no file has a new file-level skip, more hidden blocks, fewer tests run or a new skip category without a Q-entry.
