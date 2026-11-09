# KIT-111 run 9: the skip comparison

**The set:** run 1's 49 files, in both load modes.
**The baseline:** KIT-105's full-suite rows overlaid with KIT-110's, KIT-106's, KIT-107's and KIT-108's (KIT-108's `skips-*.tsv` laid over KIT-108's own overlaid baseline), compared with `kit-105/skipcmp.py`. Files outside the set are filtered out. The full output is in [skip-compare.tsv](skip-compare.tsv).

- **`load_all`** ([skips-load_all.tsv](skips-load_all.tsv)) is run 1 under strict skips: 15,753 passed, 0 failed, 0 error blocks.
- **Check mode** ([skips-check.tsv](skips-check.tsv)) ran on the `kit_111` check build (`MlQroA`, `c568dfa5`): 10,373 passed, 0 failed, 213 error blocks. Against the baseline, failures go from 149 to 0 and error blocks from 392 to 213.

**`skipcmp.py` reports two changes. Each is accounted for:**
1. **`test-workflow-execution.R`, `load_all`:** 106 → 98 passed, 10 → 11 skips.
   - **Cause: an instrumentation difference.** Run 1 left `NOT_CRAN` unset (chain A's environment), so "KIT-082: the monitor's read survives real units starting and stopping" took its `skip_on_cran()`. The baseline ran with `NOT_CRAN=true`. The file is unchanged since `6c04767f`.
   - **Rerun once,** with `NOT_CRAN=true` and otherwise run 1's tree and environment, after 111n (the block starts real units): 106 passed, 0 failed, 10 skips. That equals the baseline ([run9b/wfx-not-cran.txt](run9b/wfx-not-cran.txt)).
2. **`test-prot-norm-core.R`, check mode:** 3 → 14 error blocks.
   - Every one of the 14 is "cannot open the connection" ([run9b/norm-core-check.log](run9b/norm-core-check.log), each block's message). Three of them reach it through `rollUpPrecursorToPeptide`'s argument.
   - This is the file's known check-mode source-tree read: the shared DIA checkpoint fixture reads `tests/testdata/workflow-checkpoints/scenarios.json` from the repository root, which the check layout does not ship (KIT-105 run 8).
   - The 3 baseline blocks and KIT-111's 11 new blocks (W1 (a), W2, W3) all build on that fixture. The ticket recorded them as that category in advance (experiment budget, run 9).
   - All 14 pass under `load_all` (run 1: 727 passed, 0 failed).

**Decision:** no file gains a file-level skip, hidden blocks, fewer blocks run, or a new skip category without a Q-entry. No new check-mode error has a cause other than the known source-tree read.
