# KIT-270 R1: declaration and predictions (committed before the run)

Ticket: KIT-270 r1.1 (sections 7, 9, 10). Test-only: three patterns in `tests/testthat/test-workflow-figure-paths.R`
(`:99`, `:100`, `:471`) change from `"259"` to `"259 characters"`. No `R/` file, no stored reference, no NEWS.

**Purpose.** Show that the edited guards keep their meaning (S1, M1) and that the old ones match a temporary folder
holding 259; record the old flake's rate in one process (F1, an observation).

**Binding.** One request R1: `package_commit` = the edit commit, `base_commit` = its parent (the commit that adds this
declaration and the script; it changes no test file and no `R/` file). Three jobs: the file on head, the file on base,
the scripts job (`blas` reference; `shards` 1; `check` false). Scripts, 30 minutes each: `repeat` on base (5,000),
`repeat` on head (5,000), then `seeded,mutant` on base and head.

**Ceiling.** 130 job-minutes = 8.67 core-hours (3 jobs: 2 x 30 + 70). Expected 48 job-minutes = 3.20 core-hours.
Ticket ceiling in total 470 job-minutes = 31.33 core-hours (one correction pass; at most three instrument re-sends of
the scripts job at its 70-minute ceiling). Core-hours = job-minutes x 4 / 60.

**Predictions** (ticket section 7; N = 5,000).

| Id | Case | Base file | Head file |
|---|---|---|---|
| S1 | `seeded`: guards 1 to 3 (`:99`, `:100`, `:471`) | all three fail (message quotes `file259x`, no hint) | all three pass |
| M1 | `mutant`: the three guards on messages that carry the hint | all three fail | all three fail |
| P1 | all other expectations in the two blocks, in `seeded` and `mutant` | pass | pass |
| F1 (observation, no gate) | `repeat`: iterations failing on a 259 folder, by the process's pid in hex: neither ending in 2 or 25 nor containing 259: about 17 (99 % range 7 to 29); ending in 2 or 25: about 100 (75 to 126); containing 259: 5,000 | by pid class | 0 |
| F2 | `repeat`: iterations failing in any other way | 0 | 0 |
| W1 | the whole file by the normal mechanism | 181 passed, 0 failed | 181 passed, 0 failed |

**Decision.** S1 or M1 failing on head stops the work (stop rule 1). F1 base = 0 in 5,000 AND S1 on base not failing
stops it (stop rule 2); F1 alone stops nothing and a count outside its pid class's range is reported, not re-run. Any
unpredicted W1 verdict: one correction pass, then stop. Instrument defects (no result produced): the scripts job only,
conditions (i) to (iv) of the brief, at most three.
