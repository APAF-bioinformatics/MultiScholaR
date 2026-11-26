# KIT-268 R1 (tests-only), predictions committed before the run

Purpose: show each new or edited test failing at the stated value on the unchanged code. Binds: head = the tests-only commit
(this commit; nothing under R/ or inst/); base = 03a2a66b (kit-268 after the merge of origin/v0.9-staging 15d86bd4; no R/ or
inst/ change of its own). Files: test-kit268-bh-fallback-notice.R (new), test-kit236-report-sentences.R (T4 appended),
test-omic-da-module.R (T2b, the declared edit of the A112-03 test). One shard per side, reference BLAS, extra packages sva,
limpa, mixOmics. Budget: expected 55 job-minutes (3.67 core-hours), ceiling 70 (4.67); spent so far 2.24 core-hours (G0).
Decision: a different failing value stops the ticket (stop rule 4); a failure is by assertion, errors 0 everywhere.

| file | head prediction (passed / failed / errors / skipped) | what fails |
|---|---|---|
| test-kit268-bh-fallback-notice.R (only-on-head) | 12 / 11 / 0 / 0 | T1 1 (exists guard), T2 1 (guard), T3 6 (per omic: the named notice, its type "warning"; plus the tied_pvalues kernel case, lipid and metab), T5 3 ("1/pi0", "unknown size", "not on the same scale"); T6 passes: 7 (two vectors x 3, plus the trigger) |
| test-kit236-report-sentences.R | 289 + 18 = 307 / 6 / 0 / 0 (base 289 / 0 / 0 / 0, train k run) | T4: 6 failures, one per template (lipid, metab, LFQ, TMT, DIANN, DIANN_limpa): the added sentence after the named comparisons; its 18 other expectations pass |
| test-omic-da-module.R | base 1288 / 81 (newest record on the base: train 20261126i, run 45688253564; G0 on 4d588ba3 read 1289 / 80) so head = base - 2 passed, + 2 failed: 1286 / 83 / 0 / 0 | T2b: the two `shown(...) == 1` expectations of the A112-03 test (unmanaged and managed) fail, the test's other 3 pass (it passes 5 / 0 at G0's base today, so it is not skipped); every other expectation as on the base side |

The delta (head = base side with 2 passed turned into 2 failed, the same failure names otherwise) is the binding claim for
test-omic-da-module.R; the absolute figures come from train 20261126i.
