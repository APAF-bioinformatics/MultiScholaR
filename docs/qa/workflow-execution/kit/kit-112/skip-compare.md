# KIT-112 run 9: the skip comparison

**The set:** run 1's 54 files, in both load modes.

**The baseline:** KIT-111's overlaid baseline (`kit-111-lead/baseline-*.tsv`) with KIT-111's `skips-*.tsv` laid over it, compared with `kit-105/skipcmp.py`. Files outside the set are filtered out. The full output is in [skip-compare.tsv](skip-compare.tsv).

- **`load_all`** ([skips-load_all.tsv](skips-load_all.tsv)): this is run 1 under strict skips with `NOT_CRAN=true`. 16,994 passed, 0 failed, 0 error blocks. **No change.**
- **Check mode** ([skips-check.tsv](skips-check.tsv)): run on the `kit_112` check build (`OftnRb`, `5ecc9fd1`). 10,845 passed, 0 failed, 221 error blocks. Against the baseline, failures went from 149 to 0 and error blocks from 403 to 221.

**`skipcmp.py` reports one change, and it is accounted for:**
- **The change:** `test-workflow-saved-da.R` in check mode, error blocks 0 → 8.
- **Which blocks:** the 8 are KIT-112's new blocks:
  - W1 (a): three blocks;
  - W2 (W1 (b)): three blocks;
  - F112-03: two blocks.
- **Same message in every block:** "cannot open the connection". In two of the blocks it arrives through `rollUpPrecursorToPeptide`'s argument. The per-block messages are in [run9b/saved-da-check.log](run9b/saved-da-check.log), from one rerun of the file in check mode that only lists them.
- **Cause: the known check-mode source-tree read.**
  - Their owners build on `readWorkflowCheckpointRecipe("small")`, the shared DIA checkpoint fixture.
  - That fixture reads `tests/testdata/workflow-checkpoints/scenarios.json` from the repository root, which the check layout does not ship.
  - This is the cause KIT-105 found in run 8, and the one KIT-111 recorded for `test-prot-norm-core.R`.
- **The same blocks pass under `load_all`:** the file has 19 blocks, 861 passed, 0 failed and 0 error blocks in run 1.

**Decision:** the skip rule holds.
- No file gains a file-level skip or hidden blocks.
- No file runs fewer blocks.
- No new skip category appears.
- The one new check-mode error set has no cause other than the known source-tree read.
