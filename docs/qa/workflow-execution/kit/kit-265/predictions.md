# KIT-265: runs declared before they are sent (2026-11-26)

Ticket r1.3 (O-1 answered: row 51 (a) refuse, ledger VAL-304). Core-hours = job-minutes x 4 / 60. Expected total
9.37 core-hours (140.5 job-minutes), ceiling 64.33 (965 job-minutes): the ticket's section 11.

## Step 0 (ticket text) and the gate

The ticket's step-0 minors (R01 to R08 of the plan re-check) were already applied in r1.2; r1.3 records the owner's
answer (O-1) only. No gate precedes R1 in the ticket (A1: R1 is sent once O-1 is answered).

## R1 = G0: the tests-only commit against nothing (base "none")

- **Purpose:** show every E row and the three declared edits failing at their base values by assertion on the code
  before the fix, and every G row passing (stop rules 1 and 2). It is the lifecycle observation of the test harness.
- **Binding:** the tests-only commit on branch `kit-265` (R/ is byte-identical to origin/v0.9-staging 82a5df02).
- **Files:** `test-kit265-organism-routes.R`, `test-kit250-organism-dialog.R`, `test-prot-01c-import-module-contracts.R`.
  One shard, `check: false`, no scripts.
- **Ceiling:** 45 job-minutes (3.00 core-hours). **Expected:** 12 job-minutes (0.80 core-hours).
- **Decides:** A1; whether the fix may be written (stop rule 2: a G row failing or an E row passing means the reading
  is wrong and the work stops).

Predictions (counts derived by hand from the assertions; the pass/fail pattern per test is the prediction, a count off
by a few with the same pattern is not a wrong reading):

| File | Predicted on the code before the fix |
|---|---|
| `test-kit265-organism-routes.R` | 15 tests with failures, 0 errors; 52 failed, 67 passed. Failing by assertion: E1, E1b, E1c (message, `duration`, `id`: 3 each, the old text with a 5 s default and no id); E2d (the same three); E2 (refusal, message, one warning shown: 3); E2c (3); E3b (3: one warning shown, no refusal); E3, E3a (2 each: no refusal, silent); E4 (2); E7 (2: the sentence is absent); E1d (2: no notice under the id); E5, E6, G6 (7 each: no error, the import completes). Passing: G3, G4 |
| `test-kit250-organism-dialog.R` | 1 failure (G3, the declared edit 1: the old text), 0 errors; every other test passes, including E8 |
| `test-prot-01c-import-module-contracts.R` | 2 failures, 0 errors: the "warns when no valid organism choices can be built" test (the old text) and the "only runs when analysis is eligible" test (the ticked case with no FASTA table returns FALSE instead of an error); every other test passes |

Row E3b is a prediction (ticket P11): either the empty-distribution route (one warning) or an error in counting (one
warning); both give the same failures above. If the base shows an error there, it is recorded, not a stop.

Instrument risks named now: the observer tests (E5, E6, E1d, G3, G6) build a session as KIT-250's E8 does; G6's
non-owner session (a plain `WorkflowState`, Spectronaut peptide) is new mechanics. An error in the session's set-up
that occurs before any judged assertion, on the code before the fix, is an instrument defect (VAL-54: fix the mechanics
only, no expected value changes, at most three reruns, ceiling unchanged); an unexpected pass or failure of a judged
assertion is a result.

## R2: the fix against the tests-only commit (to be declared before it is sent)

Declared in this file when R1 has returned and the fix is written.
