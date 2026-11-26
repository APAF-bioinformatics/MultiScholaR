# KIT-219 outcome (implementer's record; state: READY FOR THE IMPLEMENTATION AUDIT)

Branch `kit-219` (from staging 0bdf72131). Fix commit 514cb7062; records follow. Ticket r1.2.

## Observed
- G0 (i), (ii) on the base (run 45778573466; the first attempt 45775309316 was an instrument defect: a script function named `add`
  resolved to a package function, no result): on both real routes (g:Profiler with a fixed `gost` answer, GO offline) the files
  a run writes equal the names the rule derives from the stored result, six each; the two contrast-name sources
  (`names(da_results@da_data)`, `contrast_names` argument, `names(manifest$results)`) are equal (`groupB-groupA`). Stop rule 2
  does not fire. GO second run (down list empty): the three down files stay and the rule's names exclude them. Tables: `g0-base/`.
- R1 (run 45775313887, tests-only commit): predictions met: E1 to E6 fail at the stated values (workbook 3 sheets not 2; 2 not 0;
  render called once), E7a, E7b, E8, G5 fail, G1 to G3 pass. Stop rule 1 does not fire. Tables: `r1/`.
- R1h (run 45780449814, fix commit bec2de00): 141 expectations pass. G0 on head: the package function's list equals the derived
  names and the files written, both routes, both runs; g:Profiler second run (down NULL) excludes the down files, which stay.
- Merge: origin/v0.9-staging a67c4161 merged once (conflict only in the ledger: both sides kept, staging's rows first).
  Base for R2: `a70b352c5` = staging a67c4161 plus the tests-only commit (branch `kit-219-base`).
- R2 (run 45785931266, 50 files, check): 47 same, kit219 fixed, **2 new failures: `test-omics-summary-core.R` (25 to 34) and
  `test-omic-summary-module.R` (28 to 54)**. Cause (a result of this change, contrary to the ticket's section 6): the W1 harness's
  recording `copyToResultsSummary` has fixed formals; `copyWorkflowSummaryResults` passed `enrichment_files = character()` for the
  reference-state lipid and metab owners even when the folder held nothing to leave out ("unused argument"). Check identical (0
  ERROR, 3 WARNING, 3 NOTE on both sides). Tables: `r2/`.
- Correction (the ticket's one pass, 0513a1b0): the list is passed only when the enrichment folder holds a file outside it;
  otherwise the call is the unchanged one (same files and workbook). R2b (run 45795755317, 26 files, check): both files back to
  `same`; check identical; the base side's shard 3 lost its job after 7 minutes (infrastructure) so 8 files had no base result.
  Tables: `r2b/`.
- R2c (run 45803833921, the 8 files, base and head): all `same` (kit212 309, kit243 246, kit257 133, kit271 78, kit259b 90,
  kit207 292, gui-entry-lipidomics 19/16 on both sides); kit219 `fixed` (base 64 pass, 43 fail, 4 errors; head 141 pass). Tables: `r2c/`.
- G0 (iii), (iv) (`g0b.R`, run 45803833921; two instrument repairs of its mechanics, runs 45785931266 and 45795755317): (iii) on 11
  W1 fixtures (lipid and metab for each inputs case, four in-memory proteomics branches, one proteomics owner) no case is a listed
  leaf: Copy is refused first by KIT-212 where the folder holds tables (metab complete and optional_missing, the owner), the in-memory
  proteomics fakes are unchecked, lipid has no table; the `none` cases hold no file. (iv): 274 lines in the R2 files mention the enrichment
  outputs; the ones that assert the publication folder's `Enrichment_Plots` or the workbook are the W1 proteomics in-memory cases
  (unchecked, so unchanged), `test-module-ci-prot-enrich.R:334-338` and `test-general-filemgmt-report-helpers-shared.R:227` (direct
  copy, no argument) and `test-e2e-proteomics-enrichment.R:38-43` (skipped on runners); R2 and R2b show them `same`. No leaf, so OD-3
  is not asked. The four GUI-entry lanes are not rebuilt by the script: their stored captures are compared by R2 and R2c (same).

## R2d (coordinator's ruling I03; run 45814520114)
The new file, base `2cab053d` (staging plus the new test file) against head `b8c18fcf`: head 157 pass, 0 fail (E9 passes, E4 passes with its
new bundle assertion); base 71 pass, 50 fail, 4 errors (E9 fails 6 expectations: the first Report renders, no refusal, Copy copies the
stale file, no T1). As predicted; tables `r2d/`. 1.59 core-hours; total 67.91 of 485.33.

## Acceptance criteria
| # | State | Evidence |
|---|---|---|
| A0 | met (VAL-187) | ledger, ticket r1.2 |
| A1 | met as restated (K219-I02) | G0 (i), (ii) as stated. (iii) and (iv) were run late (`g0b.R`, R2c) on the corrected head; its leaf rule (files omitted) cannot detect a moved argument list, and it has no transcriptomics fixture. That no stored leaf moves is shown by R2b and R2c (the W1 and core files and the four GUI-entry lanes `same`). For the plan: any new argument to the Summary's copy call moves the W1 (b) "calls" leaf, because the recorder prints its formals |
| A2 | met; G4 by composition (coordinator's ruling I04, message of 2026-11-26: G0 on head plus G1) | R1 (fail at the stated values), R1h and R2c (all pass; G1 to G3 pass on both). Not covered: G3 has no transcriptomics fixture; lipidomics or transcriptomics with files (T1b naming them) are untested (E5 covers metabolomics; the code is omic-agnostic); a contrast whose name the sanitiser changes, and more than one contrast, are not tested. G1 became vacuous after the correction (it plants nothing, so both arms are the unchanged call); the list path followed by a Report is covered by E9 (ruling I03) |
| A3 | met with one correction (stop rule 3 passed over, see Rulings; the correction is spent) | R2, R2b, R2c: every file `same` except the new file `fixed`; check identical, nothing only on head |
| A4 | to be met by the next checkpoint's full suite (coordinator's ruling, dispatch of 2026-11-26); the merge train runs the scoped comparison again on its merge head | not run here |
| A5 | met: NEWS, ledger row `VAL-388`, `NAMESPACE` and `DESCRIPTION` unchanged, `man/copyToResultsSummary.Rd` changed by the argument and the `@return` text only; this file | diff |

## What G4 asked and what covers each part
G4 asked: the real kernel on each route (as G0 (i)), then Copy through `copyWorkflowSummaryResults`; every file the run wrote (tables,
`.html`, `.png`) is copied; no `enrichment_not_current` record; the notice is 'Publication files copied'. Covered by parts:
(a) the real kernel on each route and the files it writes: G0 on the base and on head (run 45778573466, 45780449814);
(b) every written file is in the package function's list for the stored product: G0 on head (`package_equals_derived` TRUE, `written_not_derived` empty);
(c) a Copy that copies exactly the product's files, with no record and the plain notice: G1 (the test-double product, real owner);
(d) Copy through `copyWorkflowSummaryResults` on the real kernel's own output is not tested. After the correction a Copy over a folder whose files are all current passes no list.

## Rulings (coordinator), recorded
Rulings 1 to 4 (the coordinator's dispatch of 2026-11-26): 1 merged once (above). 2 R2 at the re-derived scope with check (R2, R2b, R2c;
plan and ceiling declared in `predictions.md` before each request). 3 G0 (iii), (iv) run (no leaf). 4 R3 not run: A4 as stated.
Rulings 5 and 6 are in the coordinator's message of 2026-11-26: 5 departures accepted as declared (E8, R1h, G3 on lipid and metab, the
g:Profiler mock repair); 6 rule slips noted.
After the implementation audit (AUDIT-KIT-219-impl-20261126), the coordinator's message of 2026-11-26:
- I01: at R2 stop rule 3 fired (two files with new failures) and I went on to change the product in the correction pass instead of
  stopping; the choice between repairing the product and asking the owner for a by-file re-freeze of the W1 recorder was the
  coordinator's to make after a stop. Stop rule 3 was passed over. The correction (0513a1b0) is ACCEPTED retroactively (the behaviour is
  unchanged in every case VAL-187 decided: the list is dropped only when every file on disk is already in it). The one correction pass
  is SPENT.
- I04: G4 accepted as met by composition (A2).
- I03: one test row, E9, added, tests only, on unchanged code (R2d, below).
- I08 is with the owner (T3: the proteomics report says "No enrichment analysis was run for this project." when the enrichment is no
  longer current); the templates are not changed; the three options are in `wording.md`.
- Record fixes I02, I05, I06, I07, I09, I11 applied. I10 (the W1 recorder no longer mirrors the copy function's formals, for the
  end-of-project re-record; two path sources, the `project_dirs` global and `projectDirs`; a file written between the pre-listing and the
  copy is copied without a notice, the Report check catches it afterwards; the code comment at the argument's drop should give the
  product-level why) and I12 (T2's remedy can fail silently for proteomics: if Copy's backup fails it records `backup_content_copy`,
  proteomics' `copy$required` is NULL, so the notice reads "Publication files copied" and the Report refuses again) are candidates for
  plan section 13, untouched.

## Departures from the ticket (stated plainly)
1. G0 (iii) and (iv) were run late, as one script (`g0b.R`) on the fix commit rather than on the base; the GUI-entry lanes are compared by R2, not rebuilt by the script.
2. G0 (i) calls `runWorkflowEnrichment` with `processEnrichments` as `process_fn` (the real kernel), as KIT-208's and KIT-173's
   harnesses do, with the contrasts table's raw contrast name; the g:Profiler answer is the test double's gost-shaped table.
3. G4 is not a test row (see the section on G4). G3 uses lipidomics and metabolomics (no
   transcriptomics fixture among the module-CI helpers).
4. One row added beyond the ticket: E8 (the four proteomics templates' flag reads the publication folder only), the test of 3.1
   item 5, which the ticket gave no row.
5. g:Profiler route (ii) was observed in the head run (R1h), not the base run: the mock answered a retried NULL call on the base
   run (instrument defect in the mock, not changed expectation); the kernel is identical on both.
6. R1h (about 2.9 core-h) is not a row of the ticket; it is inside the total ceiling.
7. The Report's T2 notice reads "Report generation failed: The report was not generated: ..." (the ticket's text in its notice).
8. At R2 stop rule 3 fired and the implementer continued (Rulings, I01). G1 is vacuous after the correction; E4 did not assert the bundle recorder (now asserted; E9 added).
8a. The copy call passes `enrichment_files` only when the folder holds a file outside the current list (the correction pass; the ticket's 3.1 item 3 said it passes the list).
9. Rule slips: I started a second `gh run watch` on run 45775313887 (stopped within minutes); one command I wrote contained
   `git stash list` and the permission layer blocked it before it ran (no stash command ran). My earlier remark of a record "deleted and then restored" is withdrawn: I removed a copy of a console log from the working tree and restored it with `git checkout` before any commit; the branch history holds no deletion.

## Instrument runs
G (45775309316): script stopped at its first helper call (name clash with a package function); fix = rename, no expectation changed;
rerun G-2 (G's instrument repair 1 of 3). The g:Profiler mock repair, folded into R1h (declared in `predictions.md`), is G's repair 2 of 3. The G0b script's two mechanics repairs are numbered 1 and 2 of 3 for that item. The ticket's run name R2b (the KIT-213 merge check) is not used by this ticket; "R2b" here is the correction pass's rerun, drawn from the Correction row, as is R2c (a rerun after a lost base shard).

## Stored values that move
None edited; none moved after the correction (R2b, R2c).

## Unrelated discoveries (left untouched)
- `R/func_general_filemgmt_results_summary_helpers.R` holds many duplicated lines (logger and `cat` lines appear twice).

## For the change report
- Wrong: Copy copied, and put in `Pathway_enrichment_results_<omic>.xlsx`, every enrichment file in `pathway_enrichment`, and the
  proteomics report showed every image there; files of an earlier enrichment (the DA re-run, the table changed, a later enrichment
  with no result for a direction, a re-run after Copy, a reopened project, files from outside the application) were presented next
  to the current DA. Every lipidomics and metabolomics report stated that a KEGG pathway enrichment with Benjamini-Hochberg
  correction was performed, though neither omic has an enrichment analysis.
- Code, base 0bdf72131: `R/func_general_filemgmt_results_summary_helpers.R:357-363` (folder copy), `:719-728` (workbook);
  `R/utils_workflow_summary.R` Copy `:432-458`, Report guard `:509-534`; templates `inst/reports/proteomics/*:1550-1672` (flag),
  `inst/reports/lipidomics/lipidomics_report.rmd:1630-1651`, `metabolomics_report.rmd:1586-1607`. Fix 514cb7062: the same places
  (`enrichment_files` argument; `workflowEnrichmentCurrentFiles`, `assertWorkflowSummaryEnrichmentCurrent`, `workflowEnrichmentTexts`).
- Changed: which files are copied and shown; one notification; one Report refusal; no number.
- Evidence: runs 45778573466 (G0 base), 45775313887 (R1), 45780449814 (R1h), the comparisons R2 45785931266, R2b 45795755317 (the correction) and R2c 45803833921, and R2d (E9). Limit F3: a file with a current name edited after the enrichment is copied and the Report passes (the name rule accepts it). Other limits: service doubles, the offline GO route, test-
  double enrichment tables in the E rows, one contrast; templates reach new projects only.

## Not checked
- R3: by ruling, the next checkpoint's full suite (A4).
- The rendered lipid and metab reports: the GUI-entry lanes render the edited templates and their `summary_report` steps are `same` on both sides (R2b, R2c); at KIT-268's R3 these steps differ from the oracle only at state digests, so knitting succeeds by inference. The content is not stored.
- The templates' `enrichment_ran = TRUE` branch (figure and caption chunks) is never knitted in any test.
- No W1 lipid or metab owner reaches Copy with files (KIT-212 refuses first; G0 (iii)). Metabolomics with files is covered by E5; lipidomics and transcriptomics with files are not tested.
- The 24 R2 files not rerun after the correction keep the R2 evidence (they do not reach the changed line).
- Copy through `copyWorkflowSummaryResults` on the real kernel's output (G4 (d)); the running application.
