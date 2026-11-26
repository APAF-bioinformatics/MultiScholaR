# KIT-228 predictions (committed before the run they concern)

## G0 (inventory, ticket section 5): done by reading, no run requested

`omicSummaryW1Frozen` (`tests/testthat/test-omic-summary-module.R`): every case has one Save step, or (`trx_main`) the
Save before the RUVIII export (the gate, on a fresh folder) followed by the export and a successful Save; the `*_no_label`
steps never reach the Save; `omicSummaryW1Start` leaves `integration/` empty. No case holds `integration/<omic>/` at the
start of a Save that ends without a new input. GUI-entry lanes (`tests/testdata/gui-entry/*.json`): `summary_refusal`
precedes `summary_save` in all ten and lists no `integration/` artifact (ticket section 6, re-read at this head: the
lane files were not edited since). Prediction: no listed case; stop rule 2 does not fire. R2 (base against head on the
same files) is the check that nothing moved; the G0 run (2.0 core-hours expected) is not requested: a by-reading inventory
the comparison run then tests is cheaper and stronger than an instrumented base run. Recorded as a departure.

## R1 (tests-only commit 96f32e91, R/ equal to staging 0bdf7213; one file, one job; ceiling 60 min = 4.0 core-hours)

Purpose: show each new failing-first row failing on the old code for its stated reason (stop rule 1). Decides: whether
the fix may be written. The E rows fail at the folder's presence and the missing notice; G rows and the class-lost checks
pass; counts of passing expectations are not predicted. Predicted failing expectations (22) by test:

| Test | Failing on the old code | Passing |
|---|---|---|
| R01 | 1: the message names "kept" with no previous input | the "could not be published at" match, the second half |
| E4, G1, G3, G4 (proteomics) | E4: folder present (1), notice absent (1); G3: no warning naming the folder (1), no T4 (1) = 4 | first Saves, G1, fallback file, G4 |
| E1 (trx) | folder present (1), T1 absent (1) = 2 | preconditions, `failed` class, G1 |
| E2, G4 (trx gate) | folder present (1), T2 absent (1) = 2 | gate warning, "saved successfully", G4 |
| E5, E6 (trx managed) | E5: present (1), notice absent (1); E6: present (1) = 3 | class lost, `done`/`failed` counts |
| E3, E5, E6, G1, G2, G4 (lipid) | E3: 2; E5: 2; E6: 1 = 5 | class-lost, G1, G2 text, G4 |
| the same (metab) | 5 | as lipid |

G3 is skipped only when the runner user can write into a read-only folder (not predicted on a runner). If an instrument
defect shows (a precondition fails, a fixture helper errors before a result), it is repaired under the instrument-defect
rule and sent again; an E row passing on the old code stops the work (stop rule 1).

## Comparison run R2 (after the fix; base = staging 0bdf7213, head = the fix commit)

Files: the new file; the 27-file grep set of ticket section 6 (R2) at this head plus the 14 module-level Summary files by
name; KIT-212's, KIT-216's, KIT-219's (when it exists on staging; its branch is not mine), KIT-184's new test files;
`test-kit240-integration-default.R` (the merge check after the retired identity script); `check: true`. Predicted:
every file `same` except `only-on-head` for the new file (all rows pass on head, G3 included); no stored reference moves;
no `new-failure`; `check` identical. Observation to record: base-side counts that already differ from the last record
(checkpoint findings of train 35) are "unchanged" when identical on both sides.

## R2b (instrument repair of R2; declared before sending)

R2's only red item: `test-kit247-design-resave.R` hit the file timeout I set (3600 s) on base and on head, so it gave no result on either side (an instrument defect of my request: no result seen, no expectation touched). It matters to KIT-228: its E23 asserts a refused Save removes the folder AND the top-level file, and E21 the design cut. Request R2b: that one file, head only (a7e3a582), `file_timeouts` 150 min, ceiling 160 job-min = 10.7 core-h; running total before it 37.0 of 92.0. Prediction: all expectations pass (the start removal keeps the folder-and-top-level-file removal); any failure there is a result, to be judged against base by a second request only if needed.
R2b attempt 1 (45794764741): infrastructure failure at library setup (apt 404 for libtiff, 16 min, no test started); rerun of the failed job only, instrument rerun 1 of 3, no change to request, ceiling unchanged (running total 37.0 + 0.3 core-h before it)

## Correction pass (coordinator's ruling of 2026-11-26; I04 to I07), declared before the run

Budget (core-hours = job-min x 4 / 60). Used so far, recounted: R1 0.55, R2 36.37, R2b attempt 1 1.07, R2b rerun 4.53 = 42.52 of the ticket's 92.0. This pass: tests-only run C1 declared at ceiling 4.0 (the coordinator's cap is 6.0; I declare 4.0 so that the plan fits: 42.52 + 4.0 + 45.0 = 91.52 under 92.0), then one comparison C2 at ceiling 45.0.

Design (tests 7e9fadcc; code not yet written):
- I04: a T1 variant without "the reason is in the message above" when no reason notice was shown (an in-process error outside the named classes); E4's expected text is that variant (DRAFT, VAL-41).
- I05: after a "publication pending" failure, the re-press must not delete the input the first request's worker wrote. Repair chosen: the start removal is skipped while the omic's runner is retrying (the resume finishes the first request, whose input is on disk). Why this one and not "judge a TRUE write by its manifest" (the register's second option): the manifest judgement would also change the stored W1 captures whose managed double writes no input and shows "Saved Integration S4 Object" (lipid_managed, metab_managed, trx_managed), which may not move; and it would still have deleted the worker's published input. The decided rule (a Save that ends without a new input leaves no earlier input) stays true on every path: the input kept on a retry is the new input of the accepted first request, not an earlier one.
- I06: T4 versus T1/T2/T3 is chosen from the state at the end of the Save: T4 only when the removal failed AND the folder is still there; otherwise the end state decides (folder or manifest gone: T1/T2/T3).
- I07: a row asserting T1 in the stale-head branch with an earlier folder.

C1 (tests 7e9fadcc on R/ = a7e3a582; one file, one job; base none). Predicted failing expectations on this code: 9. E4 (proteomics) 1 (the variant text); I06 proteomics 2 (T3 absent, T4 shown); lipid 3 and metab 3 (I05: the manifest is gone after the second press, and a T3 note appears in the first press's notes; I06 (ii): T3 absent). Pass on both: I07 (T1 already emitted in that branch on a7e3a582: a coverage row), every earlier row (E1, E2, E3, E5, E6, G1 to G4, R01). Expected after the code: 0 failures, all rows.

C2 (comparison): base = the tests-only tree 7e9fadcc, head = the code commit; the R2 file list (53 files) plus `test-kit264-trx-export-current.R`, `check: true`; `test-kit247-design-resave.R` with a 7,200 s file timeout; ceiling 45.0 core-hours = 675 job-minutes. Predicted: the new file `count-changed`/`fixed` (base 9 failures, head 0), every other file `same`, check identical. A miss stops the pass.

## Further correction (VAL-393, option (a); audit re-check C01 to C05), declared before the run (2026-11-27)

Ceiling 135.0 core-hours (VAL-393); used 84.2 (the auditor's recount: R1 0.55 + R2 36.37 + R2b 1.07 + 4.53 + C1 0.73 + C2 41.0). Plan: small run D1 (tests-only tree, ceiling 4.0), then one comparison D2 (ceiling 45.0): 84.2 + 4.0 + 45.0 = 133.2 under 135.0. Staging ac9a173a (KIT-219 and others) was merged into the branch once, before the tests (e5498e69; conflict only in the ledger, resolved by keeping staging's rows and mine after them), so the tests-only tree is the comparison's base.

Code to be written (nothing else): drop the `retrying()` skip, so every Save removes the earlier input at its start; one notice on the retry path (the runner is retrying and an earlier folder was removed: DRAFT "This omic's earlier integration input (integration/<omic>/) was removed by this save. If this save only finishes an earlier save's publication, that does not write it again: save once more before integrating."), shown just before the runner is handed the export (so a Save refused earlier shows its own T1/T2, not this).

D1 (tests ae9ee18d, R/ = 67e51e68's plus staging): predicted failing expectations on the current code (the skip): 16. Lipid 7 and metab 7: I05 2 (folder not absent after the second press; retry notice absent), C02 1 (folder present after the dropped resume), C01 failed 2 (folder present; T1 absent), C01 refused 2 (the same). Transcriptomics gate with a busy runner (C01 gate) 2 (folder present; T2 absent). Everything else passes on both, including the earlier rows (E1 to E6, G1 to G4, I04, I06, I07, R01), I05's first-press and "Saved Integration" expectations and the first-press manifest expectations of C02. Expected after the code: 0 failures.

D2 (comparison): base = the tests-only tree ae9ee18d, head = the code; the C2 file list (54 files), `check: true` on both sides, `test-kit247-design-resave.R` with a 7,200 s file timeout; ceiling 45.0 core-hours. Predicted: the new file `fixed`, every other file `same` (the C2 table with staging's partner changes in the base on both sides), check identical. A miss stops the pass.
