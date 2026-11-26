# KIT-270 outcome: the figure-paths guards match "259 characters"

Ticket KIT-270 r1.1. Branch `kit-270` from v0.9-staging 4d588ba3. Test only.

## Commits
- `fc3657f9` declaration, predictions and `guards-check.R` (no test, no `R/` file); the base of R1.
- `ad6314dc` the edit: `"259"` to `"259 characters"` at `:99`, `:100`, `:471` of `tests/testthat/test-workflow-figure-paths.R` (3 insertions, 3 deletions, nothing else). The package commit of R1.
- `98971967` `guards-check.R` instrument fix (see below); no `R/` file, no test file.

## Runs (`runs.tsv`)
- R1, run 45707283924: file on head and base (W1); the three script entries all exited 1 after 6 s. Log: `Error in identical(self$running[[self$current_file]]$context, context): attempt to use zero-length variable name` (`r1/scripts-1-console.log`). The script had driven testthat's `ListReporter` without a file name. No result was produced or judged: an instrument defect (brief, conditions i to iv: no expected value, threshold or stop rule changed; recorded here; ticket ceiling unchanged; instrument run 1 of 3). Fix: `reporter$start_file("guards-check")` and `end_file()` around the evaluation. The `results` job of R1 is red because of the failed script entries, as the README says.
- Instrument run 1, run 45709348632: scripts job only (`tests` empty, base `none`), package commit `98971967`, with the script reading the test file at `fc3657f9` and `ad6314dc` through `git show`. Valid because `98971967` changes no `R/` file (A1 precondition). Scripts exit 0 in 1359 s, 1364 s and 8 s. Tables in `i1/`.

## Predictions against observation (OBSERVED)

| Id | Prediction | Observed |
|---|---|---|
| S1 base | all three guards fail (message quotes `file259x`, no hint) | all three fail; folder holds 259 and no hint in the message, as predicted (`i1/seeded-fc3657f9.tsv`) |
| S1 head | all three pass | all three pass (`i1/seeded-ad6314dc.tsv`) |
| M1 base and head | all three fail | all three fail on both, hint in the message (`i1/mutant-*.tsv`) |
| P1 | other expectations pass | 0 other failures or errors in `seeded` and `mutant`, both revisions |
| F1 base (observation, no gate) | by pid class | pid 19638 = hex `4cb6`: ends in neither 2 nor 25, no 259: class 1, expected about 17 (99 % range 7 to 29). Observed **20** (guards 1, 2, 3: 10, 10, 10 failing iterations): inside the range |
| F1 head | 0 | 0 of 5,000 (pid 19824, hex `4d70`) |
| F2 | 0 | 0 on base and head |
| W1 | 181 passed, 0 failed both sides | base 181/0/0, head 181/0/0; `compare.tsv` verdict `same` (`r1/`) |

Stop rules: 1 not met (S1 passes on head: no; M1 fails on head: yes; S1 fails on base, per the coordinator's added clause: yes). 2 not met (F1 base 20, not 0). 3 not met. 4 not reached.

## Acceptance criteria
- A1 met: `git diff ad6314dc~1 ad6314dc` is the three tokens; no `R/` file.
- A2 met: S1 above.
- A3 met: M1 above.
- A4 met: F1 recorded against pid class 1 (20, within 7 to 29); head 0; F2 0.
- A5 met: W1 same counts.
- A6 met: this file, `runs.tsv`, the ledger row `VAL-355` after VAL-313 in `docs/methods-decision-ledger.md`. NEWS not touched, as the ticket says.

## Cost
Job-minutes from the jobs' start and end times: R1 about 33.5 (scripts 11.1, head 12.0, base 9.5, plan and results about 0.9); instrument run 1 about 57.7 (scripts 56.9, plan and results 0.8). Total about 91 job-minutes = about 6.1 core-hours (expected 3.20, total ceiling 31.33). The repeat loops ran at about 270 ms per iteration pair (estimate was 80 ms); each entry took 1,360 s, within its 30-minute timeout.

## Notes for the coordinator (re-check R02; not applied to the ticket text)
- Section 3 (`:61`): "to the second order; digits are not exactly uniform ..." should read "an average over fresh processes; within one process the rate depends on the pid's last hexadecimal digits, see the table below".
- Section 15 (`:188`): "section 3's p moves by a few percent with them" should read "section 3's per-process rates depend on the pid's last hexadecimal digits (the table in section 3); `RAND_MAX` is assumed 2^31 - 1".
- E-3 (`:183`): "before the next train" should read "before the next checkpoint".
- R01 (stop rule 1 with S1 on base) was applied when judging, per the coordinator's ruling.
- The ticket's per-iteration time estimate (80 ms) was low by about 3 times; the ticket's expected cost (3.20 core-hours) is below the actual (about 6.1), still far inside the ceiling.

## For the change report
Not applicable: a test-only change, no product behaviour.

## Not checked
- No R on this machine; everything above is from the runner's tables.
- Windows (runners are Linux only).
- The F1 count is one process, one pid class; the other two pid classes were not observed.
- `RAND_MAX` and the runner's C library were assumed (glibc, 2^31 - 1).
- Patterns built from variables over path-bearing messages (ticket section 2).
- The change was not run in a merge train or other test file; no file sources it.
