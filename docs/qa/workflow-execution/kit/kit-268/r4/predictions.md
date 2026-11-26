# KIT-268 R4: the one test-only correction (VAL-351; AUDIT-KIT-268-impl-20261126 I01), predictions committed before the run

Change: `tests/testthat/test-omic-da-module.R`, the A112-03 test only: the expected note is built from the shared text for the
fixture's one contrast, `fallback <- if (is.function(text)) paste0("[warning NULL -] ", text("B_vs_A", NA_character_)) else
"<absent>"`, and `shown` counts notes equal to it. Exactly one notice naming B_vs_A per route; none when nothing fell back. No
product code changes (nothing under R/ or inst/ in this commit).

Binds: base = origin/v0.9-staging 91713498 (merged into kit-268 without conflict); head = the commit named in the request.
One file, test-omic-da-module.R, one shard per side, no R CMD check, reference BLAS, extra packages sva, limpa, mixOmics.
Ceiling 50 job-minutes (3.33 core-hours); expected 45.

Prediction: head = 1288 passed / 81 failed / 0 errors / 0 skipped with failure names identical to the base side; the A112-03
test 5 passed / 0 failed; verdict `same`. Anything else is a result: recorded, nothing further changed, handed back to the owner.
