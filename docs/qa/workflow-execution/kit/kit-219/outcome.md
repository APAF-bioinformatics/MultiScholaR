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

## Acceptance criteria
| # | State | Evidence |
|---|---|---|
| A0 | met (VAL-187) | ledger, ticket r1.2 |
| A1 | met | G0 (i), (ii) base and head; (iii) and (iv) by `g0b.R` (no leaf); GUI-entry lanes by R2/R2c |
| A2 | met | R1 (fail at the stated values), R1h and R2c (all pass; G1 to G3 pass on both) |
| A3 | met with one correction | R2, R2b, R2c: every file `same` except the new file `fixed`; check identical, nothing only on head |
| A4 | to be met by the next checkpoint's full suite (coordinator's ruling of 2026-11-27); the merge train runs the scoped comparison again on its merge head | not run here |
| A5 | met: NEWS, ledger row `VAL-<merge>`, `NAMESPACE` and `DESCRIPTION` unchanged, `man/copyToResultsSummary.Rd` changed by the argument and the `@return` text only; this file | diff |

## What G4 asked and what covers each part
G4 asked: the real kernel on each route (as G0 (i)), then Copy through `copyWorkflowSummaryResults`; every file the run wrote (tables,
`.html`, `.png`) is copied; no `enrichment_not_current` record; the notice is 'Publication files copied'. Covered by parts:
(a) the real kernel on each route and the files it writes: G0 on the base and on head (run 45778573466, 45780449814);
(b) every written file is in the package function's list for the stored product: G0 on head (`package_equals_derived` TRUE, `written_not_derived` empty);
(c) a Copy that copies exactly the product's files, with no record and the plain notice: G1 (the test-double product, real owner);
(d) Copy through `copyWorkflowSummaryResults` on the real kernel's own output is not tested. After the correction a Copy over a folder whose files are all current passes no list.

## Rulings of 2026-11-27 (coordinator), recorded
1 merged once (above). 2 R2 at the re-derived scope with check (R2, R2b, R2c; plan and ceiling declared in `predictions.md` before each
request). 3 G0 (iii), (iv) run (no leaf). 4 R3 not run: A4 as stated. 5 departures accepted as declared: E8, R1h, G3 on lipid and
metab, the g:Profiler mock repair (instrument repair 1 of 3 for G). 6 rule slips noted.

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
8. The copy call passes `enrichment_files` only when the folder holds a file outside the current list (the correction pass; the ticket's 3.1 item 3 said it passes the list).
9. Rule slips: I started a second `gh run watch` on run 45775313887 (stopped within minutes); one command I wrote contained
   `git stash list` and the permission layer blocked it before it ran (no stash command ran).

## Instrument runs
G (45775309316): script stopped at its first helper call (name clash with a package function); fix = rename, no expectation changed;
rerun G-2 (instrument rerun 1 of 3). The g:Profiler mock repair was folded into R1h (declared in `predictions.md`).

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
  `R/utils_workflow_summary.R` Copy `:432-458`, Report guard `:497-540`; templates `inst/reports/proteomics/*:1550-1672` (flag),
  `inst/reports/lipidomics/lipidomics_report.rmd:1630-1651`, `metabolomics_report.rmd:1586-1607`. Fix 514cb7062: the same places
  (`enrichment_files` argument; `workflowEnrichmentCurrentFiles`, `assertWorkflowSummaryEnrichmentCurrent`, `workflowEnrichmentTexts`).
- Changed: which files are copied and shown; one notification; one Report refusal; no number.
- Evidence: runs 45778573466 (G0 base), 45775313887 (R1), 45780449814 (R1h). Limits: service doubles, the offline GO route, test-
  double enrichment tables in the E rows, one contrast; templates reach new projects only.

## Not checked
R2 and R3 (the regression and full-suite comparisons; stored W1, GUI-entry and hand-written leaves); G0 (iii), (iv); the rendered
lipid and metab reports; whether the lipid/metab owners of W1 now copy no enrichment file as predicted; the running application.
