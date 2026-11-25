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
