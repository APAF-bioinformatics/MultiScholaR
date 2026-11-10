# GUI workover: remediation of the 2026-11-11 audit

Closure record for [gui-workover-audit-20261111.md](gui-workover-audit-20261111.md); the audit itself is unchanged. Base commit `ec7db20e`.

## User instruction (2026-11-11)

1. Individual approval of KIT-101's F-23, F-24, F-26, F-27, F-28: "i think so? i cant remember." Recorded as **unverified**; nothing is changed on that basis.
2. Interactive (Glimma) volcano default: "no - i want it on." The off-by-default unification of KIT-106 is reversed.
3. Full suite: "go once all fixes are in." One run, both modes, after the last package below merges.
4. "i want everything fixed according to your recommendations." The audit's decision matrix (section 3) is the scope. Items the matrix defers (measure first, at the next re-record, after KIT-109, own ticket) stay deferred.

## Work packages

One bounded item per row. Each is test-first, runs its targeted files under `load_all`, and is reviewed line by line before it merges. No harness, protocol or `tools/profiling` file changes. No oracle is re-frozen without a named attribution here.

| Package | Items | Status |
|---|---|---|
| R1 DA record | A112-01 (F112-05) clear the pending provenance; A112-02 drop the masks; A112-03 (F112-01) fallback notice on the managed path and the applied method per contrast in the analysis record; A112-04 (F112-04) refuse all-NA or one-level formula terms; A112-06 annotation identity of the reused Glimma page | merged; A112-06 stopped |
| R2 managed runner and proteomics normalisation | A111-01 stale source of a queued action; A112-05 and A113-06 guards and notices (P-39); A111-03 worker error text; P-35 repeated or mode-switched RUV covariates; A111-08 (F111-05) image paths | merged |
| R3 enrichment | P-43 revert the worker placement; P-07 define `clean_acc`; A107-02 applied flag; P-28 declare the organism filter; A107-03 effective correction; A107-04 database version and failed calls; A107-06 input identities; A107-08 swallowed commit failure | merged; A107-06 not done |
| R4 import, design, standards | P-05 refuse non-numeric quantity text at import, with a short vendor-token list; P-11 metab import preflight; P-19 manual standards by candidate id on both paths | merged |
| R5 DA display and report | volcano default on; P-26 proteomics results download; P-27 drop `GeneID`; P-12 if the trx lane is re-frozen; A108-04 (F108-01) silent save failure; P-30 route B through `reclaimWorkflowMemory()`; A108-06 | merged; P-27 stopped |
| R6 check layout | A111-06: resolve the checkpoint fixture root from the test path | merged |
| R7 records | parked rows corrected and answered; ticked boxes; outcome rewordings; P-34 raised; lipid/metab zeros parked; KIT-114 wording; KIT-109 re-authored with one queue table | **for the implementer** |
| R8 verification | install, `R CMD check`, full suite in both modes against `kit-105`'s tables; P-30 measured once | **for the implementer** |

## Deferred by the matrix (not in this remediation)

- P-01, P-02, P-03 and the lipid/metab zeros: measurement, KIT-114.
- P-09, P-14, P-18: at the next re-record of the affected oracles.
- P-10: after the series.
- P-22, P-24, P-29: one harness revision after KIT-109.
- P-25: wording first, with KIT-109.
- P-32, P-33, P-37: accepted as they are. P-40: accepted as "no benefit shown".
- The no-plot-objects enrichment redesign: its own ticket.
- F108-04, F112-06's remainder, A103-06, the quarantined-comparison freezes: KIT-109.
- Worktree cleanup: after KIT-113 closes.

## P-15, answered 2026-11-11

The GUI uses only its own design step, but `RunApplet("designMatrix")` is exported and seven shipped workbooks call it. User: "Repair it: give it its own session store and a real run table so it works again for notebook users ... IN THE FUTURE i want to support people to learn code with this." Package R9 (in progress): the applet creates its own session owner over a temporary store and builds a real run table; `RunApplet`'s signature is unchanged; the GUI's design step is unchanged.

## Closure per item

| Item | Repair | Evidence |
|---|---|---|
| A111-06 (R6) | `tests/testthat/setup-here-root.R` declares here's root from the test directory. The fixture helper is bound by the journey harness and is not edited. | `ab843f5f`; `test-prot-qc-core.R` 18/0 under `load_all`; root resolves in a copy of the check layout. Check-mode proof is R8's. |
| P-15, A103-01, A103-02 (R9) | The applet takes the caller's `data_tbl` (with `Run`) and `config_list`, runs the design tab's builder on its own owner over a temporary store, writes `design_matrix.tab`, `data_cln.tab`, `contrast_strings.tab` and the `.rds`, assigns the four variables in the caller's environment, and removes the store on exit. Proteomics only, as the original was; other omic types are refused with a classed message. `RunApplet`'s signature is unchanged. Reviewed line by line. | `60e84761` (test, fails on the old applet: "No sample files"), `a26d5d11`; `test-applets-direct-shared.R` 19/0 → 49/0, `test-omic-design-module.R` 158/0 unchanged, in the main checkout. No real browser session was run. |
| A112-01, A112-02 (R1) | `completeDIADaPublication` clears the pending step provenance. The W2 masks on `pending` and `provenance` are removed (one kept: `pending` on `dia_failure`, D6); three pending assertions flipped to NA. No frozen value changed. | `4d8716ee`, `8868faf0`, `5434d703` (spacing, review); failing first: head entry `da_publication` instead of `pearson-selection`; `test-workflow-saved-da.R` 861/0 → 917/0 |
| A112-03 (R1) | (a) the managed delivery shows the in-process BH-fallback warning; (b) each proteomics `DA-contrast-output` result in the analysis record carries `fdr_method` and, on fallback, `fdr_fallback_reason` (generic: the kernel stores only that a fallback happened). The 11 frozen `record` hashes are unchanged, computed without the new fields; ten new `fdr` entries are frozen beside them. Visible in `analysis_record.json` only. | `37f8bd93`, `c1abd709`; `test-workflow-analysis-record.R` 274/0 → 283/0, `test-module-ci-prot-da.R` 218/0/1 → 236/0/1, `test-omic-da-module.R` 1339/0 → 1344/0 |
| A112-04 (R1) | `assertProtDaDesignTerms()` refuses an all-NA or one-level formula variable before the arguments commit (class `prot_da_design_term_invalid`). | `f48f05c6`; `test-prot-07b-da-handlers-compat.R` 35/0/2, `test-omics-da-core.R` 129/0/2 unchanged |
| A112-06 (R1, follow-up) | The saved Glimma page is reused only when it was built with this run's annotations: the annotation digest (`diaAnnotationDigest`, equal for the table in memory and the worker's read-back) is recorded in the page's factors, not its identity; a page built with other annotations is saved as the next revision and the latest revision is the page (`diaInteractiveVolcanoPage`). `diaDaPublicationIdentity` adds the digest only when there is a table, so a rerun after the annotations change republishes. Without annotations every identity, key and entity is as before; no frozen value changed. Limits: A to B and back to A republishes A's files but the tab shows B's page; old revisions are never pruned. | `383dd37e`; failing first: 6 failures, 2 errors (no gene names in the republished page, in process and managed); `test-workflow-saved-da.R` 951/0, `test-omic-da-module.R` 1369/0, `test-module-ci-prot-da.R` 250/0/1, `test-workflow-dia-products.R` 189/0 |

| P-19, A105-01 (R4) | Reproduced: IS_B ticked, a second Run without Reset used IS_A. The in-process path resolves the ticks from `normData$itsd_candidate_ids`, as the managed path. Separate commit: a manual selection naming features the state no longer holds is refused on both paths (`omics_itsd_selection_stale`), so a second Run without Reset with ticked standards now asks for a Reset; droppable on its own. First-run captures unchanged. | `1d1a5ebd`, `f095ec32`; `test-omics-norm-itsd-selection.R` 22/0 (base: 10 pass, 12 fail), `test-omic-norm-module.R` 547/0 unchanged |
| P-11, A103-04, A103-05 (R4) | Metab `design_import$preflight` is TRUE; a new `preflight_contrasts` flag (lipid TRUE, metab FALSE) keeps the file's contrasts for metab so valid imports are unchanged. An unknown contrast term is now refused for metab too. One live refusal test for lipid and metab. | `f93e38f7`; `test-omics-design-import-refusal.R` 39/0 (base: 28 pass, 11 fail); `test-omic-gui-config.R` 1851/0; lipid and metab entry lanes 35/0 and 105/0 unchanged |
| P-05, A102-01, A102-02 (R4) | `checkProtImportQuantityColumns` in all five proteomics readers: `NaN`, `Filtered`, `n.q.` read as NA with a counted notice; other text refuses the import naming the column (`prot_import_non_numeric_quantity`), including vroom "number" columns (decimal commas were read 100 times too large, silently). Empty cells and "NA" are read as NA without a notice (departure from the brief: the frozen `prot_tmt_pd` fixture has 97 empty cells). The managed design-step conversion is unchanged. Frozen reader digests and lanes unchanged; two pinned expectations in `test-prot-import-core.R` (:160, :202) changed with attribution. | `7bc24709`; `test-prot-import-quantity-text.R` 313/0 (base: 47 pass, 6 fail, 7 error), `test-prot-import-core.R` 38/0, `test-gui-entry-proteomics.R` 193/0 |

| A111-01 (R2) | `runner$run` accepts its actions as a function resolved when they start; the correlation and RUV-covariates dispatchers use it, so a click that waited reads what the action before it left. | `2db1c3e5`; failing first: the queued actions ran on the previous normalisation's final. Proven at the queueing level with a mocked dispatch, not on a real worker. |
| A112-05, P-39 (R2) | DA click guards: no commit while another tab's publication waits for a retry (notice names the step); a retry says it uses the earlier run's settings; a pending publication on a superseded state is dropped (`workflow_execution_stale`). `mod_omic_enrich.R` not edited: its guard went with the P-43 revert. | `f2e6bd3f`; `test-omic-da-module.R:5068` gains the retry notice |
| A111-03 (R2) | A failed worker's own error text reaches the tab's notice (`workflowWorkerFailureReason`); the test's notice mask now compares the reason. | `81c726ee` |
| A111-08, F111-05 (R2) | `publishDIAReferenceRuvCovariates` remaps the QC image paths to the published files. | `f6e2667a`; `test-prot-norm-core.R:1518` staged-path mask removed |
| P-35, F111-02 (R2) | `diaRuvCovariatesSource` records the covariates' input once; cleared by a new correlation, by Reset and (review) by a session restore. Re-frozen with attribution in `test-prot-norm-core.R`: the two `*_ruv_repeat` captures (:977, :1107) now hold the first click's values; the first-click captures (:941, :1071) gain `ruv_covariates_input_reference` only. | `8b4bcd3d`, `f53e894c` (review) |
| P-43 (R3) | Enrichment runs in process again. The four R files are byte-identical to `132e5445` (`git diff 132e5445 2060c862 -- R` empty). KIT-113's W2/W3 managed test blocks are deleted; its W1 in-process captures are kept. | `13c9115d`; `utils_workflow_enrichment.R` 1369 → 1134 lines, `mod_omic_enrich.R` 759 → 692 |
| P-07, A107-01, A107-02, A113-07 (R3) | `clean_acc` defined as `normalizeUniprotAccession(x, remove_isoform = TRUE)` and removed from `globalVariables`; "applied" set only when a contrast table is tested. The captures that froze the error now freeze a completed filter. | `6a3cea17`, `13117f40` |
| P-28 (R3) | `organism_filter_enabled` and `organism_filter_taxon` declared only when the filter is on, so an unfiltered run's identity, product and record are unchanged; applied and rows before/after are step outcomes; replay applies and checks them. 16 frozen leaves changed, all in filter-on sessions. | `b4d5647f` |
| A107-03 (R3) | The declared correction is the one applied (BH on clusterProfiler); the selection is kept as outcome `correction_method_requested`; the user is warned when they differ. Kernel unchanged. 12 frozen leaves. | `c6d28850` |
| A107-04, A107-05 (R3) | `perform_enrichment` signals `multischolar_enrichment_service_failure` when no attempt was answered; the declaration fails the Run before any product or state. A no-terms answer still completes. The fake gost reports `meta$version`; 9 frozen leaves. Limit: an all-no-terms run still records the service as `unavailable`. | `258b018e` |
| A107-08 (R3) | A failed enrichment commit is a warning; the tab is not marked complete. | `967fd75f` |
| A107-06 (R3) | **Not done.** Hashing the annotation table and DA files on every run changes seven frozen state captures and duplicates the reader's file-name rule; needs its own decision. | for the implementer |

R2 and R3 on the merged tree, main checkout: `test-prot-norm-core.R` 1606/0, `test-omic-da-module.R` 1369/0, `test-workflow-saved-da.R` 917/0, `test-prot-05b-norm-module-contracts.R` 502/0/9, `test-workflow-omics-worker.R` 289/0, `test-module-ci-prot-dia-limpa.R` 33/0/2; `test-omic-enrich-module.R` 1060/0, `test-omics-enrich-core.R` 241/0/1, `test-workflow-operation-registry.R` 1040/0, `test-workflow-analysis-record.R` 287/0, `test-prot-01c-import-module-contracts.R` 436/0, `test-module-ci-prot-enrich.R` 79/0/2, `test-workflow-operations.R` 75/0.

R2 left alone: after the covariates, `ruv_normalized_obj` is the RUV state, so a correlation re-apply followed by covariates still compounds; a managed failure notice carries the condition's hint lines where in process shows the message only; `test-gui-entry-proteomics.R` gave 1 to 2 unattributed failures at host load 20 to 45 (base and worktrees alike) and 193/0 on rerun. R3 left alone: KIT-113's evidence scripts and a `bindings-generator.R` block still name the removed worker functions.

| Volcano default (R5) | The interactive volcano's box is ticked by default for every omic (`mod_omic_da_panels.R`); the page is still built only while ticked and shown. Before KIT-106 lipid, metab and proteomics were ticked, trx was not. Entry oracles: lipid and metab (4 lanes, 6 `volcano_glimma` values each) return to their pre-KIT-106 value; `trx_screen` gains 21 new values over 7 steps; the five proteomics oracles are unchanged. `test-omic-da-module.R` W1 captures: 18 fields, 7 back to pre-KIT-106 values, 11 new (the page is now built by default). The harness driver is not edited; ticking a ticked box is harmless there (read only). | `a61e82f4`, `8ed0944e` |
| P-12 (R5) | Trx's design builder texts say transcriptomics; done because `trx_screen` was re-frozen by the volcano change. No lane value moved; `test-omic-design-module.R:40-41` re-frozen. | `f23cc186` |
| P-26 (R5) | Proteomics' "Download All Results" writes every contrast's saved long table as `proteomics_da_results_<date>.csv`, one contrast at a time. `prot_main` download capture re-frozen (`test-omic-da-module.R:3987`). Not tried in a real browser. | `a5704be3` |
| A108-04, F108-01 (R5) | The copy's handlers assign into the function's frame, a failed save warns, the DA-read record is restored; `copyWorkflowSummaryResults` returns workbook failures as `failed`, raises an error notice and the status line says so. The copy still counts as done, so the report is not blocked (the agent's choice, for review). | `171e7205` |
| P-30 route B (R5) | `reclaimWorkflowMemory()` runs once after the copy returns or fails. **Not measured.** | `040be345` |
| A108-06 (R5) | The summary capture masks the Version Information block; 16 `study_parameters` values re-frozen in `test-omic-summary-module.R`. | `86e3648d` |
| P-27 (R5) | **Stopped.** The notebook holding the column list (`dev/distributable/plans/transcriptomics_workflow.rmd:684`) is hash-pinned in `tools/profiling/trx_notebook_reference.R:19` and the bindings. Options: GUI and `test-trx-publication.R:152` only (the notebook then differs by one column); both, with a harness revision; leave. | for the implementer, with P-22, P-24, P-29 |

All seven packages on the merged tree, main checkout (`86e3648d`), `load_all`, one process per file: `test-omic-da-module.R` 1369/0, `test-module-ci-prot-da.R` 250/0/1, `test-gui-entry-proteomics.R` 193/0, `-lipidomics.R` 35/0, `-metabolomics.R` 105/0, `-transcriptomics.R` 29/0, `test-module-ci-trx-app.R` 50/0, `test-omic-gui-config.R` 1853/0, `test-omic-design-module.R` 160/0, `test-omics-summary-core.R` 700/0/2, `test-omic-summary-module.R` 1703/0, `test-omic-enrich-module.R` 1060/0. One merge conflict (R1 and R5 each appended a test to `test-module-ci-prot-da.R`); both kept.

R5 left alone: the trx e2e block was edited and not run; the harness driver's comment about the tick is stale; the trx builder note still gives "LCMS_Pos, LCMS_Neg, GCMS" as its example; the quarantined `test-prot-da-render-direct-shared.R` still names the zip.

R4 left alone: hex text is still read as a number; DIA-NN parquet inputs and `.xlsx` error cells are not checked; the quarantined block at `test-omics-design-core.R:787` still asserts metab `preflight` FALSE.

R1 and R4 counts are from the main checkout after the merge, `load_all`, one process per file, `R_LIBS` at the kit_112 library (the default library lacks duckplyr).

**Post-publication (user, 2026-11-11: "mark it for a post-publication need"):** a standalone design applet for metabolomics and lipidomics so their workbooks work, with the limpa workbook cleanup. Not in this remediation or the KIT series; to be ticketed after the manuscript is published.

Found by R9: the applet had been broken since 2025-12-01 (`be8dd1e4`), before KIT-048. The two metabolomics workbooks call `RunApplet("designMatrix")` without an omic type and expect `metabolite_design_matrix` and `design_matrix_metabolomics.tab`, which no package code has ever produced. `DIA_workflow_limpa.rmd:39-41` removes `RunApplet` and holds a hard-coded Windows `load_all` path. `man/dot-runDesignMatrixApplet.Rd` is stale until the next `document()`.

## Handover to the implementer (2026-11-11)

State: `performance-v2` at `86e3648d` plus this record, local only, nothing pushed. R1 to R6 and R9 are merged; each diff was read line by line before merging. Branches `audit-r1` to `audit-r5` and `audit-r9` and their worktrees under `.claude/worktrees/` are left in place.

**Not verified by anyone yet:** the full suite, `R CMD check` (with tests), the check-mode run that proves R6, any journey, any e2e browser lane, any real managed worker session. Every count above is `load_all`, one process per file, with `R_LIBS` at the kit_112 library because the default library lacks duckplyr.

**The flaky proteomics entry lane, diagnosed (2026-11-11).** Three whole runs of `test-gui-entry-proteomics.R` at host load 12 to 18 gave 193/0 each; the lane times were 240 to 611 s. `guiEntryExpectLane` asserts each lane finishes under a 900 s ceiling (`tests/testthat/helper-gui-entry.R:389-391`), so at load 30 to 45 (when the earlier failures were seen, one run taking 4455 s) the slowest lanes cross it. It is a wall-time assertion under load, not a product defect. The `ERROR` lines in the log are the lanes' intended refusals. Read a full-suite failure in this file against its lane times before treating it as a regression.

**To do, in order**

1. Stopped items that need a design decision: A112-06 (annotation identity of the reused Glimma page), A107-06 (input hashes in enrichment provenance), P-27 (`GeneID`, harness-pinned notebook).
2. Review points: the stale-ITSD refusal (`f095ec32`, droppable on its own); the workbook-failure copy still counting as done (`171e7205`); empty cells and "NA" giving no import notice; the remaining RUV compounding route after a correlation re-apply.
3. R7, records: correct and answer the parked rows (P-05, P-07, P-11, P-12, P-15, P-19, P-26, P-28, P-30, P-35, P-39, P-43 done here; P-40 to P-42 moot with the revert; P-34 to raise; the lipid/metab zeros to park with KIT-114; P-01 and KIT-114 (a)/(b) reworded as measured on a group-confounded fixture; P-03 marked as blocking the real-data gate); untick or reword the boxes the audit lists (KIT-101, 102, 103, 104, 106); reword the KIT-108, KIT-111, KIT-112 and KIT-113 outcomes per the audit; close KIT-113 as reverted; re-author KIT-109 from HEAD with one queue table (F108-04, F112-06's remainder, A103-06, the quarantined-comparison freezes, the stale KIT-019 comment, the convention divergences, the diagnostic exit boundary, AX-10).
4. Bindings: `bindings-generator.R` still names KIT-113's removed worker functions; regenerate with the next build. Journeys now build the interactive volcano page at DA for every omic, so DA-step time and memory are not comparable with builds from KIT-106 to KIT-113.
5. R8, verification, once 1 to 4 are settled (user: the full suite runs once, after all fixes): install; `R CMD check` with tests; the full suite in both modes against `kit-105`'s tables; one measurement of the report step on the two large corpora for P-30 against 108v/108w.
6. Deferred by the matrix: see "Deferred by the matrix" above. Post-publication: see the note above.
