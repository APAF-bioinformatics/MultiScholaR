# KIT-110 run 9: the skip comparison

The targeted and dependent set (run 1's 45 files) in both load modes. The comparison is against KIT-105's `skips-load_all.tsv` and `skips-check.tsv`, using `kit-105/skipcmp.py`. Files outside the set are reported as "gone" by the script and filtered out ([skip-compare.tsv](skip-compare.tsv)).

- **`load_all`** ([skips-load_all.tsv](skips-load_all.tsv)): this is run 1, with strict skips. 5,881 passed, 0 failed, 0 error blocks. The only changes are the two new files, `test-trx-norm-core.R` (58/0/0/0) and `test-trx-norm-module.R` (89/0/0/0).
- **check** ([skips-check.tsv](skips-check.tsv)): the installed `kit_110` check build, run one process per file ([skips-check.sh](skips-check.sh)).
  - **Changes:** again only the two new files, `test-trx-norm-core.R` (49/0/0/0; its source-reference expectations run under `load_all` only) and `test-trx-norm-module.R` (89/0/0/0).
  - **Declared skips:** the blocks that need the KIT-040 references skip as at KIT-105: MCI-030 3, MCI-031 3, MCI-032 4.
  - **Known failures:** the 9 failures and 10 error blocks are in the same files, at the same counts, as at KIT-105 (the known source-tree reads).
  - **First attempt:** an instrumentation failure. The script did not export `P` into its `sh -c` workers, so every file exited 2 before any test ran. It was fixed and rerun once.

**Decision:** no file has a new file-level skip, more hidden blocks, fewer tests run, more failures or errors, or a new skip category, in either mode. No Q-entry was needed (Q-120 to Q-126 unused).
