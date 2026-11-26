# KIT-228 outcome (implementer's hand-back; revised after AUDIT-KIT-228-impl-20261126: correction pass)

Branch `kit-228` (from staging 0bdf7213). Correction pass (coordinator's rulings, I01 accepted; I04 to I07 coded; record fixes I02, I03, I08 to I13) in commits 7e9fadcc (tests first), ca0c19a7 (predictions), 67e51e68 (code and records), and the records commit after C2. Commits: 96f32e91 tests only; e1105640 step 0 and predictions; a7e3a582 the fix, NEWS, ledger row `VAL-<merge>`, R1 records; then records of R2 and R2b. Records: `step0.md`, `predictions.md`, `runs.tsv`, `r1-*`, `r2-*`, `r2b-tests.tsv`.

## What changed (observed by diff)

- `R/utils_workflow_summary.R` (Save path only): `saveWorkflowSummaryParameters` removes the omic's earlier integration input once, as the first statement inside its `tryCatch`, before the head guard (`removeOmicIntegrationInput(earlierDir, omicType)`, errors named to the user, the Save continues). The two per-class removals in `summaryFailed` are deleted. Notices T1 (error or refusal, in both branches of `summaryFailed`), T2 (gate; the gate's `done` carries attribute `gate`), T3 (no table; Re-check 1 R02 text), T4 (the removal failed; R03 text), added only when an earlier folder existed and the Save ends with no manifest. Two text helpers `workflowSummaryInputRemovedText`, `workflowSummaryInputNotRemovedText` (DRAFTS for the owner's wording pass, VAL-41).
- `R/utils_workflow_integration_input.R`: the failed-publish text says "the previous input is kept" only when a previous input existed (Re-check 1 R01).
- `NEWS.md` (one bullet under the unreleased heading), `docs/methods-decision-ledger.md` (row `VAL-<merge>`). No version, NAMESPACE or man change.
- Copy and Report paths untouched (KIT-219).

## Departures and rule slips (stated plainly)

1. **Start removal uses `removeOmicIntegrationInput` (folder plus the top-level `<omic>_*` file), not `removeIntegrationInput` (folder only) as the ticket text says.** Reason: since the ticket was written, KIT-247 (VAL-239 (iv)) made the stale-head removals take the top-level file too, and KIT-247's E23 asserts it; deleting those calls with a folder-only start removal would have regressed it. Consequence: on every Save an earlier top-level RDS/TSV with another label is removed (the Save rewrites its own). Needs the coordinator's ruling (step0.md). R2b (E23 passes, 188/0) is the evidence the choice keeps KIT-247.
2. **G0 was done by reading, no G0 run** (predictions.md). R2's base-against-head on the same files is the check that no stored reference moved.
3. **G2 covers only KIT-207's class (lipid, metab) with no earlier folder**; KIT-206's and KIT-212's classes are covered by their own test files in R2 (all `same`), not by a row of mine.
4. T5 (workbook) dropped: no workbook ships (KIT-234). The `test-integration-workbook*.R` files were not added to R2 for that reason; the 38-file grep set was.
5. **R2 ceiling exceeded by about 5 job-minutes** (545.5 incl. plan and results against 540), and its `test-kit247-design-resave.R` hit MY file timeout (3600 s) on both sides, no result: an instrument defect of the request (it needs about 3560 s). Repaired by R2b (head only, 150 min limit); R2b's first attempt died at library setup (apt 404, 16 min, no test started) and was re-sent once (instrument rerun 1 of 3 for that item). No expectation changed.
6. R1's evidence line "manifest sha before/after" for E1 differs on base only because my G1 Save sits between the first Save and the refusal; the E1 failure (folder present) is the row's value. Not a defect of a result.

## Cost (core-hours = job-minutes x 4 / 60)

R1 8.2 min test job (0.6); R2 about 545 min (36.4); R2b 16 + 68 min (5.6). Total about 42.6 of ceiling 92.0 (expected 19.0 in the ticket: over the expected figure because R2 ran 8 jobs near its ceiling and R2b was needed). Running total never above the ceiling.

## Acceptance criteria and evidence

- **Tests first, failing for the stated reason (stop rule 1): MET.** R1 (45774845920) on 96f32e91 with R/ equal to staging: 22 failed, 95 passed, exactly the 22 predicted by test and row (R01 1; proteomics 4; trx E1 2, E2 2, managed 3; lipid 5; metab 5). Folder present on base, no notice, class lost on the worker route as stated.
- **Fix: the same rows pass: MET.** R2 head: the new file 117 passed, 0 failed (G3 ran, not skipped).
- **No regression, no stored reference moved (stop rules 2 and 3): MET.** R2 (45776244327), base 0bdf7213 against a7e3a582, 53 files: 52 `same`, 1 `only-on-head`, 0 new-failure, 0 message-changed, 0 count-changed (232 failures on both sides, identical; one file timed out on both, repaired by R2b); check identical (6 findings both sides, 0 new). Stored references: none moved.
- **KIT-247 E21/E23 (the top-level file goes with a refused Save): MET** by R2b: 188 passed, 0 failed on head.
- **G1 (successful Save after successful Save: no notice, folder valid): MET** in all four omics (lipid, metab, trx, proteomics rows). The cross-commit file-digest comparison from the `KIT-228 G1` log lines (base against head) was NOT made; the R2 logs of the head hold them, the base runs the old test-less tree, so only R1's old-code digests exist (`r1-log.txt`). Not checked below.
- **G4 (no notice and no folder when there was no earlier folder): MET** (proteomics, trx gate, lipid, metab).
- **Merge check after the retired identity script:** `test-kit240-integration-default.R` is in R2 and `same`; no identity script run (step0.md).
- **R01 text: MET** (R01 row passes on head, failed on base).

## For the change report

What was wrong: when an omic's Summary Save was refused or failed for any reason other than a stale table, or the transcriptomics Save had no RUVIII export, the input from an earlier Save stayed in `integration/<omic>/`, beside newer study parameters and an analysis record, and the integration imported it as current; on the background-worker route (lipidomics, metabolomics, transcriptomics) the refusal also lost its class. Who: anyone who integrated after such a Save without a later successful Save. Code (base 0bdf7213): `R/utils_workflow_summary.R:260-290` (only two classes removed the input), gate `:353-360`; fix a7e3a582: the start removal in `saveWorkflowSummaryParameters`. Results: no computed value changes; such an input is removed and the omic is absent from the integration until saved again. Evidence: R1 45774845920 (base values), R2 45776244327, R2b 45794764741. Limit: changes after the last successful Save reach the integration only after the next Save. Failed-publish message corrected (R01).

## Overlap with KIT-277 (VAL-353)

KIT-277 removes the integration input at an organism Save by its own path (`removeOmicIntegrationInput`, KIT-247's helper), not through the Summary Save; the two do not share code I changed. After this ticket the next Summary Save also removes the input at its start, so KIT-277's "the next Summary Save writes it again" holds unchanged. Not implemented here.

## Unrelated discoveries (left alone)

- The base-side failure count of 232 in the R2 file set is unchanged by this ticket (identical lists); not examined.
- `test-kit247-design-resave.R` needs about 3560 s, close to a 3600 s limit: any request including it needs a larger `file_timeouts` entry.

## Not checked

No R locally; the retry path (a Save pressed again after a "publication pending" failure, R02) and whether the Save button can be pressed while a worker runs; whether a worker killed at its deadline can still publish after `failed`; unlink semantics and locks on Windows (G3 ran on Linux with a read-only parent: the folder's contents were removed and the folder stayed); KIT-206 and KIT-212 classes with no earlier folder as rows of this file; the G1 base-against-head digest comparison; the wording of T1 to T4 (drafts); the full suite (folded into the next checkpoint, KIT-228 closes only when that run id is recorded).


## Correction pass (after the implementation audit)

Rulings applied: I01 accepted (the start removal uses `removeOmicIntegrationInput`; nothing in R/ or inst/ reads the files beside the folder; NEWS now says they are removed too, including an earlier label's file on a later Save; the owner is told by the coordinator). The one correction pass of the ticket was unused (R2b was an instrument re-send) and is used here for I04 to I07.

- **I07 (test).** Row added in the lipid/metab test: a Save, a Revert to QC complete, a Save into the same folder: T1 once, folder absent. Passes on both (T1 was already emitted in that branch): a coverage row.
- **I05 (code).** Row: a runner double reproduces the auditor's sequence (the worker writes the input; the step fails as publication pending; Save pressed again, which resumes the first request). On the code of a7e3a582 the re-press deleted the input (C1: manifest absent) and the first press's notes held both "Saved Integration S4 Object" and a T3 notice. Repair chosen in the first correction pass (SUPERSEDED by the further correction below, VAL-393): the start removal was skipped while the omic's runner is `retrying()`. Why not the register's other option (judge a TRUE write by its manifest): it would change the stored W1 captures whose managed doubles write no input and show "Saved Integration S4 Object" (stored references may not move), and it would still have deleted the worker's published input. (That reasoning did not hold: audit re-check C01 and C02.)
- **I06 (code).** Rows for the two cases the register names (the start removal fails and a later removal succeeds: proteomics no-source writer, lipid/metab no-source writer reporting TRUE). T4 ("could not be removed") is now shown only when the removal failed and the folder is still there; otherwise the end state decides (T1/T2/T3).
- **I04 (code, wording).** `kind = "failed_unnamed"`: T1 without "(the reason is in the message above)" when the error's reason was not shown (an in-process error outside the named classes; a failure that arrives through the runner is the managed step's own notice, so T1 keeps the parenthesis). E4's expected text is that variant (a DRAFT for the wording pass, VAL-41). Not a weakening: the expectation names a different, specific text.

C1 (45812278778, tests 7e9fadcc on the old code): 145 passed, 9 failed, exactly the 9 predicted (E4 1, I06 proteomics 2, lipid 3, metab 3); I07 passed. 10.9 job-min = 0.73 core-h (declared ceiling 4.0).

C2 (45814381840): base = the tests-only tree 7e9fadcc, head = 67e51e68; 54 files (R2's 53 plus `test-kit264-trx-export-current.R`) and `check`: 53 same, the new file `fixed` (145/9 to 154/0), 0 new-failure, 0 message-changed, 0 count-changed; check identical (6 findings both sides, none new). `test-kit247-design-resave.R` 188/0 on both sides at a 7,200 s limit (the base count also stands from train 35, run 45762830948, 188/0 on 761b4301, whose R/, tests/ and inst/ equal 0bdf7213's); `test-kit264-trx-export-current.R` 85/0 both sides. 615 job-min = 41.0 core-h (ceiling 675 / 45.0). Stored references: none moved.

## Accounting (audit I10), corrected

Core-hours = job-min x 4 / 60. R1 0.55; R2 36.37 (545.5 job-min against its 540 / 36.0 row: +5.5 job-min, +0.37 core-h; the coordinator records the overrun); R2b attempt 1 (apt 404, no test) 1.07 (not 0.3); R2b rerun 4.53 (R2b fits no declared row: charged to the unused "Correction" row, which it replaced); C1 0.73; C2 41.0. Total 84.2 core-hours at the first correction pass (the auditor's recount; not 83.7) against the then ceiling 92.0 (expected 19.0). The R2 kit247 timeouts were avoidable (a 60-minute limit against a file that takes about 3,560 s; train 35's record said up to 67.2 minutes).

## Record fixes

- I08: G1 as the ticket defines it (base-against-head file digests, stop rule 3) was not made as a comparison; G1 on the limpa route was not made (proteomics owner is LFQ without RUV; the transcriptomics RUV route shows `removed_` and `corrected_` parquet files in the folder instead). By reading: a successful Save's bytes cannot depend on the start removal (same writer, same inputs; `writeIntegrationInput` only finds no previous folder). R2 and C2 show every writer test file equal base against head.
- I12: step0.md's date corrected to 2026-11-26. G2 covers KIT-207's class (lipid, metab) only; after KIT-216 SD-4 KIT-206's class has no producer on the Save path, and KIT-212's class takes the same code as KIT-207's with no earlier folder. E1's evidence line compares with the first Save's manifest though G1 re-saved in between (the masked digests are equal); the failing-first value is the folder's presence.
- I13: 5 of the 14 named module-level Summary files in R2 are quarantined at file level (`test-lipid-14a-...`, `test-lipid-15-...`, `test-metab-05-...`, `test-prot-12-...`, `test-prot-12a-...`: 0 expectations on both sides); P05's evidence rests on the other nine.
- I11: G0 by reading is accepted by the audit; R2 and C2 answer stop rule 2.
- Ledger row `VAL-<merge>` and NEWS corrected (I02, I03): the file beside the folder, the earlier-label effect, "or one with no final table"; VAL-148's corrected reading is offered in the row, qualified by VAL-201 (2); the proteomics clause fixed; dated 2026-11-26.
- Not checked, added: whether `retrying()` is true on exactly the sequence the double reproduces in the application (by reading `R/utils_workflow_omics_actions.R:73-82`, `:93`); T1 to T4 as shown in the running application.


## Further correction (VAL-393, option (a); audit re-check C01 to C05), 2026-11-27

What the re-check found: the first pass's skip of the start removal keyed on the runner's `retrying()`, which is true while ANY managed step of the omic has a publication pending (C01: a refused, gated or failed Save then left the earlier input with no notice), and a resumed request dropped as superseded left the older head's input (C02). Decision (owner, VAL-393; ceiling raised to 135.0): option (a), fail closed.

What the code now does (the only change): the `retrying()` skip is gone, so every Save removes the omic's earlier input at its start (as at a7e3a582); when such a Save removed an earlier folder and is then handed to a runner that is retrying, it says once (DRAFT, VAL-41): "This omic's earlier integration input (integration/<omic>/) was removed by this save. If this save only finishes an earlier save's publication, that does not write it again: save once more before integrating." (`workflowSummaryInputRetryText`, shown in `startWorkflowSummaryIntegrationExport` just before the runner is handed the export, so a Save refused earlier shows its own T1/T2 only). The resumed request's own "Saved Integration S4 Object" notice still appears after a normal retry although the re-press removed the input: a stated limit for the owner's wording pass. Statements elsewhere (NEWS, the ledger row) were corrected to this (C03). Staging ac9a173a was merged once before the tests (e5498e69; conflict only in the ledger).

Tests first (ae9ee18d), with predictions committed before (d215f6f2): row I05 now expects the folder absent after the second press and the retry notice exactly once; new rows C01 (another step's publication pending, the runner double reports `retrying()` true: a Save that fails, a Save refused by the stale-head guard, and the transcriptomics gate: the earlier input is removed and the Save's own T1/T2 shown) for lipid and metab (failed, refused) and transcriptomics (gate), and C02 (the resume dropped as superseded leaves no input); C05 fixed (the evidence line logs the second press).

- D1, 45864156449 (tests-only tree on the previous code): 163 passed, 16 failed, exactly the 16 predicted (lipid 7, metab 7, transcriptomics gate with a busy runner 2). 12.1 job-min = 0.81 core-h (ceiling 4.0).
- D2, 45866659197 (base = ae9ee18d, head = 27336c99; the C2 list of 54 files, check on both sides, the kit247 file at 7,200 s): 53 same, the new file `fixed` (163/16 to 179/0), 0 new-failure, 0 message-changed, 0 count-changed; check identical (6 findings both sides). 596.7 job-min = 39.8 core-h (ceiling 45.0). Stored references: none moved.
- Running total: 84.2 + 0.81 + 39.8 = 124.8 core-hours of the ceiling 135.0 (C04: 84.2, not 83.7, for the earlier passes).

Not checked, as before plus: the notice texts in the running application; whether `workflowOmicsManagedActions(omicType)` hands `start...` the same runner object as the Save's `runner$retrying()` reading (the code reads it from the runner given to the export, by reading `R/utils_workflow_omics_actions.R:93, :173-177`); the full suite (folded into the next checkpoint).
