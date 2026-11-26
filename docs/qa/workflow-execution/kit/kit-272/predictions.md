# KIT-272 predictions, committed before each run is sent

Core-hours = job-minutes x 4 / 60. Ticket budget: expected 12.06, ceiling 60.33 core-hours.

## R1: tests-only commit, failing first (ticket section 7, stop rule 1)

- Purpose: show the new file's rows failing BY VALUE on the old templates and the guards passing. Binds the commit
  that adds `tests/testthat/test-kit272-prot-box-panel-text.R` (nothing under R/ or inst/). Base = origin/v0.9-staging
  7f17547f9 (file absent there: verdict `only-on-head`). One shard, check false. Ceiling 45 job-minutes (3.00 core-h),
  expected 20 (1.33).
- Decides: stop rule 1 (a base value other than the stated one means the reading is wrong: stop and report).

Test blocks (four templates LFQ, TMT, DIANN, DIANN_limpa unless stated):

| Row | Prediction on the old templates | Failing expectations |
|---|---|---|
| T1 to T3, T10 | FAIL by value (new bullet and PCA bullet absent, both `expect_false` fail) | 12 per template, 15 for limpa (its extra stages-not-known case): 51 |
| T5 | FAIL by value (every fact absent; the two `expect_false` pass) | 48 per template: 192 |
| T6 | FAIL by value (new phrase count 0 against 1; "density boxplots" present) | 2 per template, 4 for limpa (two routes): 10 |
| T7 | FAIL by value ("density plot after log2 transformation") | 4 |
| T4 | FAIL by value (3 matching lines against 0, twice per template) | 8 |
| T8, G1', T9 | PASS | 0 |
| T11 | FAIL by value (new bullet count 0 against 2; old count 2 against 0), once for lipid and once for metab | 4 |

File totals predicted: 21 failing test blocks, 12 passing; about 269 failing expectations, about 128 passing, 0 errors,
0 skipped (counts of expectations are an estimate; the pass/fail pattern per block is the prediction). An ERROR in any
block (not a failed expectation) is an instrument defect, recorded as such (VAL-54), not a result.

The declared edits of test-kit227-report-ruv-sentences.R (G1), test-kit236-report-sentences.R (G1) and
test-kit235-lipid-metab-report-parser.R (T6) are NOT in this commit: edited against the old templates they cannot pass
(they substitute the new text on the base side). They go in the fix commit (see outcome.md, deviation 1).
