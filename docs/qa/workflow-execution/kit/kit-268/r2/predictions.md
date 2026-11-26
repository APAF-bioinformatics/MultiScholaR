# KIT-268 R2 (fix against the tests-only commit), predictions committed before the run

Ruling of the coordinator, 2026-11-26: NO re-freeze in this ticket. The re-freeze run 45715855490 showed the lipid oracle's
`states` digests already differ on this base (KIT-217, train 15, ledger VAL-324; the four lipid/metab GUI-entry oracles wait for the
end-of-project re-record), so the owner's condition for O-3 (a) (diff = header + one warning) cannot be met. The `da` step's added
notification joins the stored values that wait for that re-record.

Binds: base = f0541e82a5ca671c4877b6fda7a4b1502a5e4c4e (tests-only, old code, old oracles); head = the commit named in the request
(the fix 28b88665 plus records; no change under tests/testdata). Scope: `r2/files.txt` = 39 files: the ticket's grep at the fix
tree (36, including the new file; `warnDaQvalueFallback` added to the pattern) plus test-gui-entry-lipidomics.R,
test-gui-entry-metabolomics.R, test-omic-summary-module.R. 6 shards per side, check on both sides, reference BLAS, extra packages
sva, limpa, mixOmics, file timeout 120 for the two GUI-entry files. Budget: expected 270 job-minutes (18.0 core-hours), ceiling 320
(21.3); spent so far 5.2 of the ticket's 56.33. A missed prediction stops the ticket (stop rules 4 to 6).

| file | base (newest record on the base: train 20261126i, run 45688253564, unless stated) | head | verdict |
|---|---|---|---|
| test-kit268-bh-fallback-notice.R | 12 / 11 / 0 (R1 run 45712438665) | 23 / 0 / 0 / 0 | fixed |
| test-kit236-report-sentences.R | 307 / 6 / 0 (R1) | 307 / 0 / 0 / 0 | fixed |
| test-omic-da-module.R | 1286 / 83 / 0 (R1; train i had 1288 / 81) | 1288 / 81 / 0 / 0 (the A112-03 test back to 5 / 0) | fixed |
| test-gui-entry-lipidomics.R | 19 / 16 / 0 / 0 (train i; its 16 failures are the 16 steps qc_intensity to summary_export whose `states` digests drifted, the same 16 steps as in refreeze/diff-vs-stored.txt) | 19 / 16 / 0 / 0: the same 16 expectations fail; the one at step `da` (`expect_identical(steps[["da"]], oracle$steps[["da"]])`, helper-gui-entry.R:414) fails with a message that now also shows the added warning in `notifications`; the 19 passing expectations (as on base) still pass | message-changed (only the `da` step's message) |
| test-gui-entry-metabolomics.R | 57 / 48 / 0 / 0 (train i) | 57 / 48 / 0 / 0 (G0: no metabolomics lane falls back, so no new notification) | same (a message-changed row is possible: train i showed one on metab_gc between its own base and head; the counts are the claim) |
| test-omic-summary-module.R | 1675 / 28 / 0 / 0 | 1675 / 28 / 0 / 0 (G0: its captures never reach a fallback) | same |
| the other 33 files of files.txt | as on the base side of this run | the same counts | same |
| R CMD check | | identical findings on both sides | identical |

Not asserted: any number, table or file of a result (stop rule 5): a difference in a numeric output in any file stops the ticket.
