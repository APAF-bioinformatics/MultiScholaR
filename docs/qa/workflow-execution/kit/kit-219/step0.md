# KIT-219 step 0 (2026-11-27)

Base: `origin/v0.9-staging` 0bdf72131 (KIT-216 merged by train 35, VAL-377; KIT-212, KIT-243, KIT-273, KIT-274, KIT-268,
KIT-247 in). Ticket r1.2. Branch `kit-219`, worktree `/home/doktersmol/Projects/APAF/MultiScholaR-kit219`.
Bounded item: the ticket's section 3.1 (items 1 to 5) and 3.2, tests first, then the declared runs. No R was run.

## Text repairs of the plan audit's Re-check 1 (R1-01 to R1-03)

Already in the ticket at r1.2 (the coordinator applied them; closure table 17): 3.1 item 5's last sentence, NEWS's
last sentence, the F8 row, A2's "G1 to G4 pass on both; G5 fails on base", 2.2's P01 row ("the paragraph on every report
whose `Main()` succeeds ...") and 3.2's R2 sentence. Read and checked at the lines named; nothing contradicted the code.
No ticket edit was needed.

## What moved since the ticket's citations (67064baf) and what it changes

No statement of the ticket fails in a way that changes what is built or predicted. Moves, all by reading:

| Ticket citation | Now (0bdf72131) | Effect |
|---|---|---|
| `diaEnrichmentProduct` "does not exist" (KIT-216.md:137) | `R/utils_workflow_dia_state.R:965`; marker `list(state, source, da_results)` set at `:953`; the record is read with `readDIAStateObject(manager, marker$state)` (= `manager$getState`), as `restoreWorkflowEnrichment` does (`R/utils_workflow_enrichment.R:1330-1360`) | predicate exists as the ticket describes; the contrasts and directions come from the record's `enrichment_data` slot |
| Copy guard chain `R/utils_workflow_summary.R:~278` | `copyWorkflowSummaryResults` `:384-470`: one `tryCatch` of `assertWorkflowSummaryDaCurrent(files = TRUE)`, `assertWorkflowDaPublicationComplete`, `assertWorkflowSummarySaveCurrent` (KIT-243) before the fallback file; the copy call `:432`; notice `:455-458` | none of KIT-243's guards looks at enrichment files; this ticket adds nothing to Copy's guard, only to the call and the notice |
| Report guard `:365-372` | `renderWorkflowSummaryReport` `:497-540`: `refused <- tryCatch({ HeadFinal; DaCurrent(files=TRUE); DaPublicationComplete; SummarySaveCurrent })`, notice "Report generation failed: <text>" | item 4's check is appended to that chain (same catch-all) |
| Module gate `R/mod_omic_summary.R:200-205`, `:278` | `:202-221`: since KIT-243/KIT-212 K212-I02 the Report resets `files_copied` when the DA's entity id differs from the one at the last Copy (`omicSummaryDaReference`, `:250-260`) | S1 after a DA re-run is partly caught at the module; S1 by a table change, S2 and S6 are not: the core's md5 check stays needed |
| `copyToResultsSummary` `:344-357`, `:686-707`, `:719-728` | folder specs `:357` and `:359-363`; workbook `:716-737`; save `:752-759`; return `:949`. The function now also writes RUV delivery tables (KIT-273) and the limpa naming (KIT-271) | none touches `pathway_dir`. The workbook code assigns a local `enrichment_files` (`:719`), the name the ticket gives the new argument: the local is renamed in this change |
| Templates (`DIANN_report.rmd:1198-1227` etc.) | the same chunks, 440 lines further down: DIANN `:1638-1672`, limpa `:1550-1584` (STRING `:1586-1613`), LFQ `:1591-1625`, TMT `:1582-1616`; lipid `:1630-1651`, metab `:1586-1607` (content unchanged; the paragraph chunk `Pathway Enrichment Section` lipid `:930`, metab `:917`) | same edits at the new lines |
| Enrichment file names (2.1) | `generate_enrichment_plots` writes `<contrast>_<dir>_enrichment_results.tsv` only for a direction with a plot (`R/func_general_enrichment_engine_helpers.R:379-386`); plots from the sanitised `contrast_names` (`R/func_general_enrichment_orchestration.R:~540-570`); the test double writes both tables per contrast (`R/utils_workflow_dia_enrichment.R:108-122`) | as the ticket; the derived names are a superset of what a run writes (a direction with a non-NULL entry, no file when its plot is NULL) so no written file is called "not current" by the rule; G0 checks the reverse |
| `test-prot-summary-report-templates.R:43-55` | `:43-56`, unchanged | item 5 keeps the three asserted strings |
| Helpers | `goPropRunProject` (`helper-go-propagation.R:118`) now also returns `product` (`goPropProduct`); KIT-216's file builds a proteomics owner with `k216Owner`, `k216RealDa`, `k216Enrich`, `k216Report` | reused, copied with sources named, in the new test file |

Stored references: none edited. W1 and GUI-entry leaves are not touched by reading (the ticket's 6); G0 (iii) is not a
separate run here: see `predictions.md` for what the tests and the comparison decide. The ticket's addendum
(K212-I05: `test-omics-summary-core.R` also runs nine proteomics cases) is recorded: its proteomics cases use in-memory
fakes (`usesWorkflowReferenceStates` FALSE), which get `NULL` from item 2 and copy as today.

## Not changed from the ticket
Scope, rows, texts, budget (expected 111.7, ceiling 485.33 core-hours) and stop rules are as the ticket.
