# KIT-219 outcome (implementer's record; state: STOPPED before R2 for the coordinator's re-price)

Branch `kit-219` (from staging 0bdf72131). Fix commit 514cb7062; records follow. Ticket r1.2.

## Observed
- G0 on the base (run 45778573466; the first attempt 45775309316 was an instrument defect: a script function named `add`
  resolved to a package function, no result): on both real routes (g:Profiler with a fixed `gost` answer, GO offline) the
  files a run writes equal the names the rule derives from the stored result, six each; the two contrast-name sources
  (`names(da_results@da_data)`, `contrast_names` argument, `names(manifest$results)`) are equal (`groupB-groupA`). Stop rule 2
  does not fire. GO route second run (down list empty): the three down files stay on disk and the rule's names exclude them.
  Tables: `g0-base/`.
- R1 (run 45775313887, tests-only commit 6cc29e6d): predictions met: E1 to E6 fail at the stated values (workbook 3 sheets
  not 2; 2 sheets not 0; render called once), E7a, E7b, E8, G5 fail, G1 to G3 pass. Tables: `r1/`. Stop rule 1 does not fire.
- R1h (run 45780449814, fix commit bec2de00): 141 expectations pass, 0 fail. G0 on head: the package function's list equals the
  derived names and the files written, both routes, both runs; g:Profiler second run (down list NULL) excludes the down
  files, which stay on disk. Tables: `r1h/`.

## Acceptance criteria
| # | State | Evidence |
|---|---|---|
| A0 | met (VAL-187) | ledger, ticket r1.2 |
| A1 | partly: per real route files against derived names met (G0); the per-leaf survey (iii) and assertion survey (iv) were not run | G0 tables; see Not checked |
| A2 | met | R1 (fail at the stated values), R1h (all pass; G1 to G3 pass on both) |
| A3 | NOT DONE: R2 not sent | R2 scope is 49 files (list below), above the ticket's 40: the ticket requires the coordinator to re-price first |
| A4 | NOT DONE: R3 not sent | after R2 |
| A5 | partly: NEWS, ledger row (`VAL-<merge>`), `NAMESPACE` and `DESCRIPTION` unchanged, `man/copyToResultsSummary.Rd` changed by the argument and the `@return` text only; this file | diff |

R2 scope (grep of the ticket's formula at the branch, helpers dropped, plus the ticket's named additions): 39 files by the grep
(`test-e2e-proteomics-enrichment.R`, `test-enrichment-significance.R`, `test-general-enrichment-shared.R`,
`test-general-filemgmt-path-contracts.R`, `test-general-filemgmt-report-helpers-shared.R`, `test-go-propagation.R`,
`test-kit205-enrichment-background.R`, `test-kit207-summary-head-guard.R`, `test-kit208-gprofiler-identifiers.R`,
`test-kit212-da-current.R`, `test-kit213-study-parameters-labels.R`, `test-kit215-not-tested-lists.R`,
`test-kit216-prot-session-after-enrichment.R`, `test-kit219-enrichment-files-current.R`, `test-kit226-gprofiler-columns.R`,
`test-kit243-summary-current-outputs.R`, `test-kit257-norm-own-files.R`, `test-kit259b-second-correction.R`,
`test-kit271-limpa-table-name.R`, `test-kit273-limpa-ruv-delivery.R`, `test-metab-05-summary-module-characterization.R`,
`test-module-ci-metab-da.R`, `test-module-ci-metab-summary.R`, `test-module-ci-prot-enrich.R`, `test-module-ci-prot-summary.R`,
`test-multiomics-enrich-shared.R`, `test-omic-design-module.R`, `test-omic-enrich-module.R`, `test-omics-enrich-core.R`,
`test-omics-summary-core.R`, `test-omic-summary-module.R`, `test-prot-11a-enrichment-module-characterization.R`,
`test-prot-11-enrichment-module-contracts.R`, `test-prot-12-summary-module-contracts.R`, `test-workflow-analysis-record.R`,
`test-workflow-dia-products.R`, `test-workflow-gui-diagnostic.R`, `test-workflow-memory-mechanism-member.R`,
`test-workflow-omics-figures.R`), plus `test-gui-entry-{proteomics,lipidomics,metabolomics,transcriptomics}.R`,
`test-prot-summary-report-templates.R`, and the template-naming files `test-lipid-14a-summary-module-characterization.R`,
`test-lipid-15-summary-module-seams.R`, `test-module-ci-lipid-summary.R`, `test-module-ci-sentinels.R`,
`test-omic-gui-config.R` (49 distinct). Needs the coordinator's price (kit-212/216/243 files are long owner journeys).

## Departures from the ticket (stated plainly)
1. G0 (iii) (per W1 case and GUI-entry lane, whether Copy omits files) and (iv) (assertions on `Enrichment_Plots` contents) were
   not run as a script; they rest on the ticket's reading (6) and on R2's base-against-head comparison, which is not yet sent.
2. G0 (i) calls `runWorkflowEnrichment` with `processEnrichments` as `process_fn` (the real kernel), as KIT-208's and KIT-173's
   harnesses do, with the contrasts table's raw contrast name; the g:Profiler answer is the test double's gost-shaped table.
3. G4 (the real kernel then Copy through `copyWorkflowSummaryResults`) is not a test row: the real kernel's files against the
   package function's list is G0 on head (R1h), and the copy mechanics are E1/G1. G3 uses lipidomics and metabolomics (no
   transcriptomics fixture among the module-CI helpers).
4. One row added beyond the ticket: E8 (the four proteomics templates' flag reads the publication folder only), the test of 3.1
   item 5, which the ticket gave no row.
5. g:Profiler route (ii) was observed in the head run (R1h), not the base run: the mock answered a retried NULL call on the base
   run (instrument defect in the mock, not changed expectation); the kernel is identical on both.
6. R1h (about 2.9 core-h) is not a row of the ticket; it is inside the total ceiling.
7. The Report's T2 notice reads "Report generation failed: The report was not generated: ..." (the ticket's text in its notice).
8. Rule slips: I started a second `gh run watch` on run 45775313887 (stopped within minutes); one command I wrote contained
   `git stash list` and the permission layer blocked it before it ran (no stash command ran).

## Instrument runs
G (45775309316): script stopped at its first helper call (name clash with a package function); fix = rename, no expectation changed;
rerun G-2 (instrument rerun 1 of 3). The g:Profiler mock repair was folded into R1h (declared in `predictions.md`).

## Stored values that move
None edited. By reading, none predicted to move (ticket 6); R2 would show.

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
