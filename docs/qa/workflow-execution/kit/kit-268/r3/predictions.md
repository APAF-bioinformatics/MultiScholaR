# KIT-268 correction pass: run R3, predictions committed before the run

Authority: the coordinator's ruling of 2026-11-26 (one correction pass, defect 1 only). R2 (run 45721176971) stands as recorded.
Fix: `daQvalueFallbackText` iterates by position (`R/func_omics_da_shared.R`; nothing else under R/). No new test: T1's "A_vs_B"
expectation of `test-kit268-bh-fallback-notice.R` failed on exactly this defect in R2 and is the test that catches it.

Binds: base = 7f17547f99f5440f1330d7e760722be994460611 (staging, merged into kit-268; it has none of this ticket's tests or code);
head = the commit named in the request (staging + the ticket's tests, fix and records). Files (6): test-kit268-bh-fallback-notice.R,
test-omic-da-module.R, test-omic-enrich-module.R, test-kit236-report-sentences.R, test-gui-entry-lipidomics.R,
test-gui-entry-metabolomics.R; R CMD check on both sides; 3 shards per side, reference BLAS, extra packages sva, limpa, mixOmics,
file timeout 120 for the two GUI-entry files. Budget: expected 145 job-minutes (9.7 core-hours), ceiling 200 (13.3); spent so far
24.0 of the ticket's 56.33 core-hours. A miss on any declared prediction ends the pass.

| file | head prediction | verdict |
|---|---|---|
| test-kit268-bh-fallback-notice.R | 30 / 0 / 0 / 0 (R2 head: 29 / 1; the one failure was T1's names expectation) | only-on-head |
| test-omic-da-module.R | 1288 / 81 / 0 / 0 (the A112-03 test 5 / 0) | head = base side with T2b passing; the base side is staging's own count, observed |
| test-kit236-report-sentences.R | base side + 24 passed, 0 failed: 313 / 0 / 0 / 0 on head if staging's file is still 289 / 0 (R2: base side 307 / 6 = 289 + 18 + 6 failing; 24 = 18 + 6) | new-failure-free; counts differ by the new T4 test only |
| test-gui-entry-lipidomics.R | the same counts as its base side (19 / 16 on train 20261126i; the base side is observed); the only difference is in the `da` step's message (the added warning) | message-changed, `da` step only |
| test-gui-entry-metabolomics.R | the same counts as its base side (57 / 48 on train 20261126i) | same (message-changed possible, as on train i) |
| test-omic-enrich-module.R | NOT PREDICTED. Both outcomes are observations: 1011 / 59 on head means R2's six stringdb_plot failures were not reproduced (R2's figure is then an unreproduced difference for the auditor, with the findings of outcome-notes.md); 1005 / 65 again is a reproducible change of a static text this ticket should not reach: STOP after recording, change nothing further, hand back with the element's failure text from both sides | recorded as observed |
| R CMD check | identical findings on both sides: ERRORs 0, WARNINGs 3, NOTEs 3 each (R2) and no finding only on head or only on base now that the records' paths are at most 86 bytes (R2's one head-only and one base-only finding were the non-portable-names list naming those records) | identical |

The other 33 files of R2 keep their R2 verdicts: they do not load the changed function (`daQvalueFallbackText` is called only by
`notifyDaQvalueFallback`, itself called only by `warnDaQvalueFallback` and `completeOmicDaRun`), so a by-position loop over the same
groups cannot change what they observe; the merge train runs the full scope.
