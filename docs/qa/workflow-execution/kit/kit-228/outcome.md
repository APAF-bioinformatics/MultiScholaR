# KIT-228 outcome (implementer's hand-back; for its implementation audit)

Branch `kit-228` (from staging 0bdf7213). Commits: 96f32e91 tests only; e1105640 step 0 and predictions; a7e3a582 the fix, NEWS, ledger row `VAL-<merge>`, R1 records; then records of R2 and R2b. Records: `step0.md`, `predictions.md`, `runs.tsv`, `r1-*`, `r2-*`, `r2b-tests.tsv`.

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
