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

## R1 result (run 45698567707, 2026-11-26): as predicted

All 52 predicted failures and 67 passes of the new file, 1 of the kit250 file (G3), 2 of the 01c file, 0 errors, by
assertion at the stated values (the old text with a 5 s default and no id; no refusal; the import completing). G3 and
G4 pass. Stop rules 1 and 2 did not fire. About 13.0 job-minutes (0.87 core-hours) against 12 expected and 45 ceiling.
The per-test failing set is identical to the prediction (checked from `failures.txt`).

## R2: the fix against the tests-only commit

- **Purpose:** show the new rows, the two declared-edit files passing on the fix, and no regression elsewhere.
- **Binding:** base 2269d03ad00746e0db3be900819c8240a2d1a52a (tests only), head = the records commit on `kit-265` that
  follows the fix c40435efe (it adds NEWS, ledger, README-free documents only: no further change under R/).
- **Files (the ticket's section 10.3, recounted at the branch point 82a5df02: 9 files by the stated `git grep`, plus
  the new one):** `test-kit265-organism-routes.R`, `test-e2e-test-mode-integration.R`, `test-kit250-organism-dialog.R`,
  `test-module-ci-prot-import.R`, `test-prot-01c-import-module-contracts.R`, `test-prot-01d-import-module-characterization.R`,
  `test-prot-import-core.R`, `test-workflow-gui-memory-lifetimes.R`, `test-workflow-server-shared.R`,
  `test-workflow-ui-shared.R`. One shard per side, `check: true`.
- **Ceiling:** 340 job-minutes (22.67 core-hours). **Expected:** 113.5 job-minutes (7.57 core-hours).
- **Decides:** A2 and A3.

| File | Base (tests only) | Head (fix) |
|---|---|---|
| `test-kit265-organism-routes.R` | 67 passed, 52 failed (R1) | 119 passed, 0 failed (`fixed`) |
| `test-kit250-organism-dialog.R` | 118 / 1 (R1) | 119 / 0 (`fixed`) |
| `test-prot-01c-import-module-contracts.R` | 436 / 2 (R1) | 438 / 0 (`fixed`) |
| the seven other files | any state | `same` |
| R CMD check | | findings identical on both sides |

Stop rules: an unpredicted verdict stops the work; a head failure or a differing check is the one correction pass.

## R2 result (run 45700788798): as predicted

Verdicts: `fixed` for the three files above with the predicted counts, `same` for the seven others (the 11 skips of
`test-workflow-gui-memory-lifetimes.R` are on both sides); `R CMD check` findings identical (0 errors, 3 warnings,
3 notes on both sides; none only on head, none only on base). About 56.5 job-minutes (3.77 core-hours) against 113.5
expected (the files ran in 10 to 16 s each, not the 1.75 minutes the ticket assumed, and the check jobs took 13 to 16
minutes) and the 340 ceiling.

## H1 (browser): declared before it is sent

- **Purpose:** the fix in a headless browser on a runner (A4). Observation only.
- **Binding:** the commit under test is 4a8915ac (R/ identical to the fix c40435efe); the request's package commit is the
  records commit that adds `h1.R`, `g0.R` and `k233-journey.R` (KIT-250's driver, byte for byte, with the cases and
  predictions replaced) under this folder only; `h1.R` checks that the difference is limited to this folder.
- **Settings:** one scripts job, `shinytest2`, test mode, network namespace (loopback only); script timeout 50.
- **Ceiling:** 60 job-minutes. **Expected:** 15 (four sessions may take 25).
- **Predictions:** written in `h1.R` (`K250_PREDICTIONS`) before the run: sessions A (headers without organisms:
  persistent notice naming Mus musculus (taxon 10090), import completes, Design tab opens), B (ticked, empty FASTA
  file: the (iv) error, window closed, nothing imported; then unticked in the same session: completes), C (unticked,
  the same empty file: completes), D (valid FASTA naming organisms, none of the data's proteins: the (c) error).
- **Instrument risk named now:** an empty (0-byte) FASTA upload, and the error-toast count, are new mechanics. If the
  upload or the toast reading fails before a judged value is read it is an instrument defect (VAL-54); an empty file
  that does not fail the FASTA step (the import then continues on route (i)) is a result, not a defect.

## H1 run 1 (run 45703724770, 2026-11-26): 4 of 4 sessions attached; 37 as predicted, 2 differ, 2 not observed

About 13.5 job-minutes (0.90 core-hours). Every judged value of the fix was as predicted: session A's notice names
"Mus musculus (taxon 10090)" and the import completes (setup_import complete, design_matrix pending, one notice under
the fixed id, still there after 12 s); B's refusal (the (iv) text, window closed, `setup_import` incomplete, Design tab
disabled, the 10 s toast on the page, no success notice); C unticked completes as before; D's (c) refusal text.

Two instrument defects, no result produced for them (VAL-54, rerun 1 of 3; no expected value, threshold or stop rule
changes; the ticket's ceiling is unchanged):
1. Session A's Design-tab readings (`design_tab_active`, `data_available`, "not observed"): the driver stopped with
   "could not find function k250DataAvailable" (a helper that lives in KIT-250's `h1.R`, not in `g0.R`). Fix: the
   helper is defined in this `h1.R`.
2. Session B after the advised untick (`import_completed` FALSE and `design_data_available` "(not completed)",
   recorded as DIFFERS): the driver's wait returned at once because the refusal's 10 s error toast was still on the
   page, so it did not wait for the second import's completion line. The same session's other readings show that the
   second import completed (`setup_import` complete, `design_matrix` pending, "Data import successful!", no error).
   Fix: a refusal is an error toast that is new since the press. The values are re-observed in the rerun, against the
   same predictions.

## H1 rerun 1 (run 45706090379): 41 of 41 as predicted

4 of 4 sessions attached; 0 differ, 0 not observed. About 13.0 job-minutes (0.87 core-hours). H1 in all: 26.5
job-minutes (1.77 core-hours) against 15 expected and the 60 ceiling.

## Correction pass (AUDIT-KIT-265-impl I01; coordinator's ruling 2026-11-26): declared before it is sent

The ticket's one correction pass. Staging 15d86bd4 was merged first (a14f956c4; NEWS and ledger kept both sides,
staging's first; no conflict under R/, tests/ or inst/; the ticket's scope grep still gives 9 files, 10 with the new
file). The notice follows the form's Taxonomy ID, which alone sets the UniProt query and the enrichment organism:
name and taxon present: as before; taxon present and the name empty: "This analysis will use the Taxonomy ID entered
in the form: <id> (the form gives no organism name). "; taxon NULL or NA: "The form has no Taxonomy ID, so UniProt
annotation will have no organism to use. " (the register's suggestion; the enrichment clause is dropped because
enrichment's behaviour is not checked). Session A's notice (both present) does not change, so no H1.

Tests first, commit 6eb7dd239 (tests only, nothing under R/): row E1e; the new blank-form text in E1b, E1c and in the
two declared edits (KIT-250 G3, the 01c test). The fix is the next commit.

### Comparison run (R2's form): base = 6eb7dd239 (tests only), head = the fix on it

- **Purpose:** E1e and the moved texts fail on c40435ef's code (the base side) for the stated reason and pass on the
  fix; no regression. Stop rule: a result that misses its prediction stops the work (no second correction).
- **Files (the 10, recounted on the merged tree):** `test-kit265-organism-routes.R`, `test-e2e-test-mode-integration.R`,
  `test-kit250-organism-dialog.R`, `test-module-ci-prot-import.R`, `test-prot-01c-import-module-contracts.R`,
  `test-prot-01d-import-module-characterization.R`, `test-prot-import-core.R`, `test-workflow-gui-memory-lifetimes.R`,
  `test-workflow-server-shared.R`, `test-workflow-ui-shared.R`. One shard per side, `check: true`, reference BLAS.
- **Expected:** about 60 job-minutes (R2 took 56.6) = 4.0 core-hours. **Ceiling:** 120 job-minutes (8.0 core-hours),
  inside the ticket's correction line of 385.

| File | Base (c40435ef's code, new tests) | Head (fix) | Verdict |
|---|---|---|---|
| `test-kit265-organism-routes.R` | 126 passed, 3 failed (E1b, E1c, E1e: one assertion each, the message: the old blank-form sentence, and for E1e the same sentence instead of the Taxonomy ID's) | 129 / 0 | `fixed` |
| `test-kit250-organism-dialog.R` | 118 / 1 (G3, the message) | 119 / 0 | `fixed` |
| `test-prot-01c-import-module-contracts.R` | 437 / 1 ("warns when no valid organism choices can be built", the message) | 438 / 0 | `fixed` |
| the seven other files | their R2 counts | the same | `same` |
| R CMD check | | | findings identical (0 errors, 3 warnings, 3 notes both sides) |

The failures on the base side are by assertion at the stated values and only the message differs (the notice's
`duration`, `id`, type and the rest already hold on that code). Counts are derived by hand from the assertions.
