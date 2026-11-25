# KIT-267 predictions (committed before each run)

Work item: VAL-310 -- the eight Enrichment defaults lose the "not yet confirmed" mark; the three organism defaults say
"human; change for other species"; no value changes. Branch kit-267 from v0.9-staging 4d588ba3ed5c6fee76a5a1450c028b409ba419ba.
Expected cost: 11.27 core-hours (169 job-minutes x 4 / 60); ceiling 31.33 core-hours (470 job-minutes).

## R1 (tests-only commit against base 4d588ba3; no check; expected 14 job-minutes, ceiling 25)

Purpose: show the new and moved tests failing on the old code at the stated values. Decides: a different failure value
or an unexpected failure elsewhere stops the ticket (ticket section 14, rule 1).

| File | Prediction on head (= old code, tests changed) |
|---|---|
| test-kit240-integration-default.R | T6: 9 failures (1 setequal + 8 per-setting); T10: 5 failures (old-mark count 8 vs 0; organism-mark count 0 vs 3; three label matches); every other test passes (T9 skipped without mofapy2) |
| test-kit267-enrichment-defaults.R | T3: 2 failures (8 provisional of 13 vs 0; `mark_organism` absent) |

Base column (same files on 4d588ba3): kit240 file passes as on KIT-240's runs; the new file is absent (only-on-head).

R1 observed (run 45704530515): as predicted (9, 5, 2). R1 cost 45.2 job-minutes (3.01 core-hours), above the ticket's 25 job-minute step figure (four test jobs of about 9 to 14 minutes each, plus plan and results); recorded as a deviation; the total ceiling (470) is unchanged.

## R2 (fix commit against the tests-only commit; base = the commit that holds these predictions, code identical to 8f334001; expected 155 job-minutes, ceiling 210)

Scope recounted on the base: 22 files (the ticket's 21 plus the new file); the list is in the request. R CMD check on both sides.
Predictions: every file `same` except test-kit240-integration-default.R (`fixed`: T6 (b) and T10 pass, 0 failures on head; base shows its 14 failures) and test-kit267-enrichment-defaults.R (base 1 pass 2 fail; head 3 pass 0 fail: the file exists on the base, so the verdict is `fixed`, not `only-on-head`). For the 20 other files base and head counts are identical. R CMD check: identical findings on both sides (check-compare.txt lists none only on head or only on base). KIT-240 T6 (c) passes on both sides.
Decision: any other file differing, or a T6 (c) failure, is stop rule 3 / 2 of the ticket.
