# KIT-106 run 9: the skip comparison

The set is run 1's 49 targeted and dependent files, in both load modes. The baseline is KIT-105's full-suite rows overlaid with KIT-110's rows for its 45 files; the comparison script is `kit-105/skipcmp.py`. Files outside the set show as "gone" and are filtered out. Full output: [skip-compare.tsv](skip-compare.tsv).

- **`load_all`** ([skips-load_all.tsv](skips-load_all.tsv)) is run 1 itself, under strict skips: 9,441 passed, 0 failed, 0 error blocks, on `9152e613`.
- **Check mode** ([skips-check.tsv](skips-check.tsv), [skips-check.sh](skips-check.sh)) used the corrected build's check directory (`J57eHO`, `e1a5b960`), one process per file. In the first pass, the new `test-omics-da-core.R` failed 16 times. Each failure was its source-file check, which the installed package cannot answer. `e30ca0b2` applies KIT-104 `ca35c462`'s rule, so the check runs under `load_all` only. The file's row was rerun once in check mode (113/0/0/2) and replaced.

**Every change against the baseline, in both modes, is one of:**
- **A recorded quarantine:** Q-130, Q-131, Q-133 to Q-135, Q-143 to Q-148, Q-151 to Q-154. The four whole-file quarantines skip at the file level with the Q-id on line 1.
- **A new file:** `test-omic-da-module.R` (1,250/0 at run 1, 1,260/0 after the correction's block) and `test-omics-da-core.R`.
- **MCI-016's retarget:** 142 to 138. The 4 checks on the deleted resolver became checks on the load's status and the configured file name ([findings.md](findings.md), W4 note).
- **Check mode only:** one more error block in `test-workflow-kit-harness.R` (53 to 54). The new D10 block reads `tools/profiling/workflow_kit_driver.R` from the source tree, which the check layout does not ship. This is the same source-tree read behind the file's other 53 check-mode errors.
- **Improvements:** fewer error blocks (`test-lipid-02b` 1 to 0, `test-metab-02aa-…characterization` 4 to 0), and 149 to 5 failures overall in check mode for this set.

The trx MCI check-mode skips stay as declared (MCI-030 3, MCI-031 3, MCI-032 4).

**Decision:** no file has a new file-level skip, more hidden blocks, fewer tests run or a new skip category without a Q-entry or a recorded retarget.
