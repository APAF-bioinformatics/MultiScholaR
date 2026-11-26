# KIT-277 outcome (implementer): STOPPED AT STOP RULE 3; not ready for the implementation audit

Code commit 53e4c745 on branch kit-277 (tests-only trees 6aeb21e6, 55aab0cc, 5e4f86dc, f16edb13 before it).
The one correction pass is UNUSED. Core-hours used: about 71 of 200.00 (R1 7.73, re-send 2.87, R1c about 3.0, SG1 about 1.9, R2
about 55; expected total was 68.95 without the extra runs).

## What was observed (run ids in runs.tsv; tables in tables/)
- R1 45777669520 (old code): new file 47 passed, 93 failed, 10 instrument errors; KIT-247 file 16 failures, all in the edited
  W288a x2 and W288b. Re-send 1 45789773919: 0 errors; F4, P1 to P6, R1b and P8 all fail by value on the old code.
- R1c 45797066259 (old code, new lines): the Copy-press lines and the P4 head rule (prot_lfq) fail by value; G4 passes on both sides.
  **The P4 head-rule expectation passes on the old code for prot_tmt, prot_dia and prot_dia_limpa (their lanes stop at the cleanup
  step, so the head is already not a DA state): it is NOT EVIDENCE for those three lanes and is discriminating only on
  prot_lfq. Not adjusted.**
- SG1 45797086832 (code): new file 113 passed 0 failed 17 skipped (S2/S3); globals 14 passed 0 failed. Held on the first run.
- R2 45804909404 (head 53e4c745, base f16edb13; 134 files, 4 shards per side, check both): 132 files `same`; test-kit247-design-resave.R
  `fixed` (base 170 passed 16 failed, head 186 passed 0 failed: W288a x2 and W288b); R CMD check identical (0 ERRORs, 3 WARNINGs, 3
  NOTEs; none only on one side); no stored value moved. Failures in other files are identical on both sides (unchanged).
- **MISS against the prediction "the new file passes every row on the head": 241 passed, 3 failed (verdict `message-changed`).**
  1. O4 (prot_dia, after an organism Save and the DA and enrichment runs again): `copy_after` is "refused" (expected "allowed") and the
     Copy press after the re-run raises "Copy error" (expected none). The log's message for it is "The differential analysis was run on an
     earlier table, not the current one" (the DA status is "stale" after the re-run). On the old code both passed. Cause NOT diagnosed (no
     R locally, the run keeps no state). Candidates read: the order of the test (an enrichment press that is refused, then the DA run, then
     enrichment) or a real defect of the head after the organism Save (a re-commit of the head by the DA run or the Summary Save changing
     the DA's source table components). It must be treated as a possible product defect until a diagnostic shows otherwise.
  2. P4 prot_lfq: `record_da` is TRUE after the FASTA cut (expected FALSE). By reading `buildWorkflowAnalysisRecord`
     (R/utils_workflow_analysis_record.R:556-566) lists the DA whenever `snapshot$da_results` is set, with no currency check, and the FASTA
     cut (`selectWorkflowStateReference`) leaves `da_results` in place (the DA is "stale", which the same row asserts and which passed). Re-check
     2 states both "the DA becomes stale" and "the record lists no DA"; they cannot both hold without either clearing `da_results` at the cut
     (then the status is "none", not "stale") or a currency check in the record. **Needs a ruling.**

## Expectation edits (old text -> new text)
- W288a (prot_lfq, lipid_canonical): "a formula-only re-save asks nothing, resets nothing and says nothing changed" -> "asks, sets the
  differential analysis back and the next run fits the new formula" (question asked and confirmed; notice "Formula saved."; DA status
  "pending"; no DA reference; the DA run's recorded formula `~ 0 + group + replicates`).
- W288b: "asks nothing, applies nothing and says nothing changed" -> "asks, and Save applies the new organism and removes the earlier
  annotation" (taxon 10090, "Mus musculus", NULL annotation, "Organism saved.", history and head unchanged).
- P8 (new file, instrument repair of R1): old `expect_error(runOmicsDesignResaveGate(NULL, "unchanged", TRUE, ..., limited = limited),
  class = "omics_design_resave_unconfirmed")`; new: the same call captured by `tryCatch(..., error = function(e) e)` and
  `expect_s3_class(condition, "omics_design_resave_unconfirmed")`. The class asserted is the same and any other error still fails the
  row; on the old code it errored ("unused argument") and now fails by value. Nothing is weaker.
No test was removed.

## Departures and rulings recorded
- Four shards per side (the job limit is 300 minutes), not six. Instrument repair 1 of 3 used (mechanics only). The organism rule found
  defective by row G4 (asked on a name-only difference) was fixed before any code run, not a correction pass.
- Coordinator rulings applied: the Copy press itself, the head rule (`assertWorkflowSummaryHeadFinal`, the function the Summary Save
  (:297), Report (:524) and Export Session State (:622) call) asserted by class and text, G4, P8's form.
- Overlap: kit-228 removes the integration input on a Save that writes none; this ticket removes it at the organism Save (carryProtDesignOrganism).

## For the change report
Before: a design saved again with only the model formula, organism, taxonomy id or FASTA file changed was answered "nothing changed" and applied
nothing, so the differential analysis kept fitting the earlier formula, enrichment and the analysis tables kept the earlier organism's annotation,
and an ordinary FASTA file was ignored whenever the project already held a FASTA table; the analysis record could name a FASTA that was not used
on any Import Existing Design into a project holding a table (R/mod_prot_design_import_helpers.R hydrate, staging :80-135, record :437-441).
Now (code 53e4c745): the Save asks, carries the change through (formula committed on the head's settings or placed on the DA tab's box; organism
applied with the annotation removed everywhere and the DA and enrichment set back; FASTA applied and recorded as the file applied; after Accession
Cleanup the history goes back to the step before it), and every full import with a new `.fasta` into a project holding a table now applies that
file (cleanup results can differ on label-free, TMT and DIA with limpa). Evidence: the runs above; limits: the three misses above are open.

## Not checked
- F2b: no real project reopen (the row asserts `omicDaDesignFormula` and the head marker). P5: the worker path is not drivable under testServer
  (the row asserts the project copy the worker reads equals the applied table and the cleanup state's `fasta_sha256`;
  R/utils_workflow_execution.R:2842-2859 read). An existing reopen or worker fixture might carry both: for the auditor.
- The Report, Summary Save and Export Session State presses themselves (the shared head rule is asserted by value; lane prot_lfq only).
- Lipidomics and metabolomics session export with a NULL DA. A session opened from a file cannot see an import's formula change.
- The enrichment tab's `enrichment_da_stale` text names a changed table as the reason; left for the wording pass.
- No stored reference was edited; R2 found none moved.

## Coordinator rulings after the stop (2026-11-27), recorded
- R2's failure output is committed under tables/r2/ (failures-head.txt, compare-failures.txt, compare.tsv, new-file-head.log). The
  log shows which refusal text O4 raised: BOTH of O4's Copy refusals on prot_dia (after the organism Save, and after the DA and
  enrichment ran again) are the "files unmade" text (`summary_da_files_unmade`: "The results folder ... holds differential analysis
  files that no differential analysis of this session made (da_proteins/da_proteins_groupKO-groupWT_long_annot.xlsx, ...Volcano_Plots/...
  and 6 more)"); the "earlier table" (stale) text quoted in the first outcome belongs to the P4 prot_lfq Copy press (the expected
  refusal after the FASTA cut), not to O4. The files-unmade text is raised whenever the DA status is not "current", so after the
  re-run the status is not "current" (none or stale) or the re-run's files do not match: the diagnostic below settles it.
- Miss 1 (O4): a diagnostic outside the correction pass (declared in predictions.md; purpose: read the store at six points; nothing is
  judged by it; ceiling 5.0 core-hours; one send, base and head). The cause is to be read from its log; no repair yet.
- Miss 2 (P4 prot_lfq `record_da`): ruling (K): the product is right as built (the DA is stale after the FASTA cut, like KIT-247's design
  cut; no writer can reach a record that lists it: the head rule refuses the Summary Save, Report and Export there). The expectation
  `expect_false(out$record_da)` was mis-specified (the ticket's 4.3 sentence and Re-check 2's R15). It is replaced IN THE CORRECTION PASS
  (not yet) by the user-level guarantee: the Save press after the cut is refused and writes no analysis record (digest of
  analysis_record.json unchanged), failing by value on the base and passing on the head. The owner is told.
- Departures stated: R2 ran 4 shards per side under the 300-minute job limit (the register's 6 was for 80 minutes; inside the budget).
  The coordinator's ruling R12 = (b) was applied: W288a and W288b are EDITED, not deleted, and the KIT-247 file ran with a
  120-minute file timeout (R1 and R2).
- Outside this ticket, recorded: a session export at a DA-state head with a stale DA still packs the DA folder
  (R/utils_workflow_dia_sessions.R:1782, `collectDIANativeGuiFiles`, condition `!is.null(snapshot$da_results)`).

## Reading of the O4 diagnostic (run 45825080961; tables/diag/head.log, base.log; no repair made)
Head (code 53e4c745) at the six points A to F (after the organism Save; after the first Save Workflow Arguments and Copy; after the refused
enrichment press; after the DA run; after the enrichment run; after the second Save and Copy):
- `da_results` entity "<none>" and `workflowDaResultStatus` "none" at EVERY point, A to F. The DA status is never "stale" and the three
  tests of the status (owner, membership, component difference) are never evaluated ("<no DA reference>"). The "files unmade" text is
  raised because the status is "none" (R/utils_workflow_dia_da.R `workflowDaFilesUnmade`, called by `assertWorkflowSummaryDaCurrent`
  with files = TRUE), as read.
- At A the set-back is as designed (da_results NULL, enrichment marker FALSE). The break is at D: the DA press on the head produced
  no accepted manifest (da_results still "<none>", publication marker "complete", unmade files "<no DA>").
- `commitDIADaArguments` committed nothing: the head entity (d1912228...) and the head's `args$daAnalysisUI` (timestamp 13:16:54, the
  lane's own DA) are identical at A to F. The annotation fetch DID run at D (O2b and O2 pass in R2: the bus table is the new organism's).
- The DA press added two notifications (the list grows 3 -> 5 between C and D); their texts are NOT in the log: the diagnostic
  printed them with `str`, which truncates, and the DA handler's own messages are `cat()` output captured by the lane harness.
Base (the same file on staging's R) at the same points: `da_results` present and "current" at A to C (owner, membership true, component
difference FALSE, ids equal), a NEW da_results entity after D (c3350c25...), the head entity changes at C, D and E (the earlier enrichment
and the DA re-run commit), unmade files empty, publication marker "complete", copy gate "allowed" throughout.
Pointer: the divergence is between the set-back (correct at A) and the DA press at D; by reading the press passes the null check, runs
the fetch (R/mod_prot_da_handlers_run_analysis.R:71) and then does not reach `acceptDIADaManifest` (R/utils_workflow_dia_da.R, accept half),
i.e. an error inside the handler's `tryCatch` (`report_error`: "Error in DA analysis: ...") is the candidate. The log does not show which.
Needed to settle it: one more diagnostic send printing the FULL text of every notification at D (`paste(notes, collapse = " || ")`, a
one-line change to `kit277Diag`) and the error class. NOT done: awaiting the coordinator's word.
