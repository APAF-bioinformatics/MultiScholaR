# KIT-106: entry-lane changes

Every re-frozen entry oracle in KIT-106, the decision it is attributed to, and every capture that changed.

## D3 (the user's unification, 2026-11-06): the interactive volcano is off by default and built on demand

The `volcano_interactive` box is unticked by default for every omic, and `volcano_glimma` is `<req>` until it is ticked (W3 `43517efa`, lipid and trx; W4 for metab). An entry lane never ticks the box. Every capture where the lane's volcano output used to show the default-ticked Combined banner now reads `<req>`.

### lipid_canonical (re-frozen once at `39cfe7af`, the W3 merge)

- **Verify before the freeze** ([refreeze/lipid-verify-before.txt](refreeze/lipid-verify-before.txt)): 29 passed and 6 failed against the oracle. The record's JSON diff against the oracle shows exactly 6 values changed, all `lipidomics_workflow-de-volcano_glimma`, at the steps `da`, `summary_refusal`, `summary_save`, `summary_copy`, `summary_report` and `summary_export`. Each changed from `<ui> d1968b0f12a6b4ca08dcf08bc4b7e6bb225a387b0558e2285e7bbf2d0ec601dd` (the Combined banner) to `<req>`. Nothing else differs.
- **The freeze** (`MULTISCHOLAR_GUI_ENTRY_FREEZE=1`, once; [refreeze/lipid-freeze.txt](refreeze/lipid-freeze.txt)): the oracle's diff ([refreeze/lipid_canonical.diff.txt](refreeze/lipid_canonical.diff.txt)) is those 6 values plus the header fields `frozen` and `source_sha`.
- **Verify after** ([refreeze/lipid-verify-after.txt](refreeze/lipid-verify-after.txt)): 35 passed, 0 failed.
- `trx_screen` does not change: trx's volcano was never default-ticked in the lane (W3: its record's diff against the oracle is empty).

## W1 unit captures (not entry oracles)

Five `test-omic-da-module.R` values, in two cases, move with D3 and D4. The lead decided to re-freeze them once each, with attribution: [w3/w1-refreeze.md](w3/w1-refreeze.md).
