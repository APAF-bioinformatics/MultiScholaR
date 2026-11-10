# KIT-113 run 9: the skip comparison

**The set:** run 1's 59 files, in both load modes.

**The baseline:** KIT-105's full-suite rows overlaid in ticket order with KIT-110's, 106's, 107's, 108's, 111's and 112's `skips-*.tsv`, compared with `kit-105/skipcmp.py`. Files outside the set are filtered out.
- **The baseline was rebuilt after the scratch overlay was lost** in a session restart. Before use, it was checked by reproducing KIT-112's recorded comparison (`kit-112/skip-compare.tsv`) exactly, in both modes: the totals and the one change.
- Output: [skip-compare.tsv](skip-compare.tsv).

**The runs:**
- **`load_all`** ([skips-load_all.tsv](skips-load_all.tsv)): run 1 under strict skips with `NOT_CRAN=true`. 19,924 passed, 0 failed, 0 error blocks. **No change.**
- **Check mode** ([skips-check.tsv](skips-check.tsv)): on the `kit_113` check build (`cmBMvS`, `95600989`). 13,675 passed, 0 failed, 226 error blocks. Against the baseline, failures went from 149 to 0 and error blocks from 411 to 226. **No change** in any targeted file.
- **Interrupted, then completed:** the check-mode run was cut off once by a session restart before it started any file. It was rerun once, alone, with the e2e lanes (run 8) kept from the interrupted chain, which had completed.

**Decision:** the skip rule holds in both modes.
