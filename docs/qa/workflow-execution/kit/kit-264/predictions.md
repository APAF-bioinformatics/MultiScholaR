# KIT-264 R2 predictions (written before R2 is sent)

R2: head 317c3630 (the fix) against base d68ec438 (tests-only), 67 files (66 by the ticket's grep at 988512e6, plus
the new file), 5 shards per side, check both, MULTISCHOLAR_TRX_MANIFEST_NOTE=1. Ceiling 980 job-minutes.

- New file `test-kit264-trx-export-current.R`: verdict `only-on-head` is not expected (the file is on both sides); base
  70 passed / 15 failed (R1-2, run 45687693688), head 85 passed / 0 failed / 0 errors / 0 skipped; verdict `fixed`.
  Every E row (E1, E2, E3, E5, E6, E7, E8) passes on head; L1, L2, G1, G2, G5 to G8, E8 control pass on both; X1 prints
  one line on both and `succeed()`s.
- Every other file: `same` (equal counts, identical failure messages). Known failing on both sides, predicted `same` at
  the latest trains' counts (passed/failed): test-prot-norm-core.R 1548/306, test-workflow-saved-da.R 755/196,
  test-omic-enrich-module.R 1011/59, test-omic-summary-module.R 1675/28, test-omics-summary-core.R 676/25,
  test-omics-enrich-core.R 221/20, test-workflow-trx-declarations.R 313/31, test-workflow-operation-registry.R 1095/2,
  test-module-ci-trx-da.R 52/1, test-module-ci-trx-norm.R 63/1, test-module-ci-trx-summary.R 87/1,
  test-module-ci-trx-design.R 66/0 with 1 error, test-workflow-kit-harness.R 734/1, test-workflow-omics-binding-steps.R 58/1.
  `test-omics-norm-state-guards.R` 651/0 and `test-kit254-omics-correlation-skip.R` 379/0 (train 13 values).
- `test-kit257-norm-own-files.R` 133/0 `same` (no KIT-257 expectation reads the trx list).
- `R CMD check`: identical on both sides (no findings only on head).
- No stored reference moves. Any other verdict is stop rule 3.
