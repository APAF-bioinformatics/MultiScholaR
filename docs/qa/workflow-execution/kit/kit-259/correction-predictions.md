# KIT-259 correction pass: predictions (written before the run)

The coordinator's ruling of 2026-11-26 on R2's stop (run 45629589384): `test-workflow-saved-exports.R:730-731`
required ggplot2's removed-rows warning from the cut-off bar of the A1 defect. It is test code, not a stored reference,
and becomes a declared edit (ticket r1.3, section 6): the expectation now asserts that saving the correlation product
raises no removed-rows warning (`expect_no_warning(..., message = "Removed [0-9]+ rows?")`); the rest of the block
is unchanged. This is the ticket's one correction pass.

Bindings (both with origin/v0.9-staging 42846c0e merged):

- base: branch `run-k259-c-base` 545c4557 = tests-only e0cb5a99 + the edit (3a0f3c04, cherry-picked from b2703094)
  + staging; R/ as staging (the old drawing code).
- head: branch `kit-259` (this commit) = the fix 1d6737c1 + the edit b2703094 + staging (merge 1fc916b2) + this file;
  base and head differ only by R/ (the fix) and records.

One run: three files, base against head, no check (R2's check stands; the merge train repeats it); ceiling 2 jobs x
45 minutes.

| File | Base (predicted) | Head (predicted) | Verdict |
|---|---|---|---|
| `test-workflow-saved-exports.R` | as R2's base (271 passed, 1 error: lobstr absent) except the edited expectation fails ten times (the warning is still raised): 261 passed, 10 failed, 1 error | 271 passed, 0 failed, the same 1 error | `fixed` |
| `test-kit259-pearson-histogram-figure-failure.R` | 39 passed, 32 failed (R1) | 85 passed, 0 failed (R2) | `fixed` |
| `test-metab-01k-s4-plot-pearson-characterization.R` | 16 passed, 1 failed | 17 passed, 0 failed | `fixed` |

Any other verdict, or other counts, is a result: recorded, not re-run.
