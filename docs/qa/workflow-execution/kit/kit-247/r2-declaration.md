# KIT-247 R2 declaration (written before the request is sent)

- Purpose: regressions, base against head on the same runner image, and R CMD check.
- Base: `526ec2aa` (the tests-only commit on v0.9-staging e0747fc1). Head: the commit carrying this file on `kit-247`
  (fix `6bd2c194`, fix 2 `0da49542`, the departure-3 test `de744acf`, E13's mechanics `1d412824`, this declaration).
  The coordinator named head `0da49542`; the head moved by the two test commits and this records commit, none under `R/`.
- Files: `r2-files.txt` (56): the ticket's grep over `tests/testthat` (helpers excluded) at e0747fc1, widened by the
  summary's functions (`save_workflow_args`, `removeOmicIntegrationInput`), plus the kit-212, kit-222, kit-225 tests,
  `test-prot-norm-input-current.R`, the four GUI-entry lanes and the two new files.
- Shards: 4 per side; `check: true`; ceiling 8 x 70 + 2 x 90 = 740 job-minutes (ticket 11).
- Decides: A3 (every file `same` except the new files, which go from failing to passing; check: nothing only on head).
  Stop rule 3: any G row failing on head, or any other file moving: stop and report.
