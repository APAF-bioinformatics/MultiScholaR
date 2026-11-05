# KIT-097 remediation (audit_and_repair, 2026-11-06)

**Authority:** the user, 2026-11-06, "fix all recommendations", with these decisions:
1. RUVSeq and EDASeq are installed in the kit library (KIT-040).
2. KIT-095 extends the shared DA saving to record edgeR QL.
3. Existing modules are parameterised by omic type, not copied. "Unchanged" means the other omics' module CI, declaration tests and frozen dump pass unchanged.
4. RUVIII is an export step after DA, product only, with `M` recorded.
5. Transcriptomics runs under protocol v1.8 (KIT-098).

**Audit:** [KIT-097-AUDIT.md](KIT-097-AUDIT.md), left unedited.
**Edited:** `tickets/kit/KIT-097.md`, revision 1.0 to 1.1. Nothing else: no code, no other ticket, no index and no plan.
**Validator:** `next-openai-ticket.sh KIT-097` is valid, with 0 errors, 0 warnings and status consistent.

## Index-relevant changes

- **Title:** unchanged.
- **Dependencies:** unchanged (`KIT-096`).
- **Time estimate:** unchanged (8h).
- **Frontmatter description:** revised. The index's description is generic, so it needs no sync.
- **`files_to_modify`:**
  - **Added:**
    - `R/mod_transcriptomics.R`;
    - `R/mod_lipid_norm_*.R`, `R/mod_lipid_da*.R`, `R/mod_lipid_summary*.R`, `R/func_lipid_da_plotting.R` and `R/utils_workflow_omics_actions.R` (the parameterised base);
    - `R/func_general_filemgmt_results_summary_helpers.R`;
    - `inst/reports/transcriptomics/transcriptomics_report.rmd`;
    - `tests/testthat/helper-module-ci-trx-*.R`, `tests/testdata/mci/manifest.json` and `tests/testthat/test-workflow-omics-application-reads.R`;
    - `tests/testthat/helper-e2e-browser.R`, `tests/testthat/test-e2e-browser-harness.R`, `tests/testdata/e2e/manifest.json`, `tests/testdata/e2e/trx_small/` and `tests/testdata/e2e/report_templates/transcriptomics_report.rmd`.
  - **Narrowed:** `inst/reports/` to the single template file.

## Closure

"Closed in ticket" means the ticket now carries the repair and a falsifiable criterion. It does not mean the implementation is verified; that waits for the ticket's run.

| Finding | Severity | Disposition | Where |
|---|---|---|---|
| KIT097-F01 | blocking | **Closed in ticket; producer routed.** KIT-095 owns the record change (user decision 2). The DA tab shows that recorded model and offers no TREAT input; a fold-change field, if kept, is display-only. The DA-model criterion covers edgeR `glmQLFit(robust)`, `glmQLFTest`, BH, the formula with W terms and `treat_applied = FALSE`. | Outcome: DA tab; AC "DA model"; report forbidden list |
| KIT097-F02 | blocking | **Closed in ticket** (user decision 4). The state sequence ends with DA on the RUVs head. RUVIII is a product launched from the summary tab after DA, with `M` in its parameters and record outcome. The `M` choice is removed from the normalisation tab. The first criterion asserts the DA `source_state`, `M` on two runs, and DA results invariant to `M`. | Outcome: state sequence, summary; AC "Full GUI journey" |
| KIT097-F03 | should-fix | **Closed in ticket; routed.** KIT-093 names the BH `FDR` column. KIT-097 labels it "FDR (BH)" and never uses a Storey q-value as the rule. The significance criterion requires counts equal to the reference and an `identical` FDR column. | Outcome: DA tab; AC "Significance" |
| KIT097-F04 | should-fix | **Closed in ticket; routed.** KIT-095 declares the parameters and outcomes; KIT-097 exposes them. `automatic` (the default, the notebook's grid) or `manual` (`k` within the cap); no `skip`. The RUV-selection criterion records the mode, settings, `k`, percentage, number of controls and scores, with automatic runs `identical` to the reference. | Outcome: normalisation tab; AC "RUV selection" |
| KIT097-F05 | should-fix | **Closed in ticket.** The required and forbidden report statements are listed, and the UI help text says "inferred". The report-wording criterion is a grep of the rendered report, kept as evidence. | Outcome: report, norm tab; AC "Report wording" |
| KIT097-F06 | should-fix | **Closed in ticket; routed.** The output map is now in the ticket: notebook output to GUI, results summary and report. KIT-093 owns the renderers and KIT-095 the products. A missing product goes back to its owner. The output-map criterion asserts every mapped output exists. | Outcome: output map; AC "Output map" |
| KIT097-F07 | should-fix | **Closed in ticket.** `func_general_filemgmt_results_summary_helpers.R` is in scope for a transcriptomics branch fed by the map. The output-map criterion checks `results_summary_dir` and the report. | `files_to_modify`; AC "Output map" |
| KIT097-F08 | should-fix | **Closed in ticket; contract routed to KIT-095.** The summary tab launches `trx_ruviii_export` with KIT-095's contract (product, format, path). The integration-export criterion requires the matrix to be `identical` to the reference `newY`. | Outcome: summary; AC "Integration export" |
| KIT097-F09 | should-fix | **Closed in ticket** (user decision 3). The reuse map lists what is reused as is, what is parameterised (the lipidomics family as base; metabolomics not edited) and what is trx-local (the normalisation sequence only). Backward compatibility is the user's definition. The hand-off criterion lists the map as built. | Outcome: reuse; Regression Checks; Backward Compatibility |
| KIT097-F10 | should-fix | **Closed in ticket.** The allowlist is extended (see above). | Frontmatter |
| KIT097-F11 | should-fix | **Closed in ticket.** The KIT-082 probe gains a transcriptomics case with its positive control. The non-matrix values the app holds are listed. | Outcome: no whole matrix; AC "No whole matrix" |
| KIT097-F12 | minor | **Closed in ticket.** The simulation setup names the behaviour tests, each to fail for its asserted reason. | Simulation Setup |
| KIT097-F13 | minor | **Closed in ticket** (user decision 1). The e2e runs on a 12-sample corpus, one attempt, 1,800 s ceiling, against the kit library. Any skip of an acceptance test is a failure. | Experiment budget |
| KIT097-F14 | minor | **Closed in ticket.** Static and interactive volcano semantics are declared. The notebook's non-formal Glimma `lfc`/`fdr` arguments are not carried over as a rule. | Outcome: volcano semantics |
| KIT097-F15 | minor | **Closed in ticket.** Parity is defined for the notebook's first contrast; the others come from the same fit. | Outcome: contrasts |
| KIT097-F16 | minor | **Closed in ticket.** The `trx-tab-*` `data-testid`s and action ids are listed in the outcome for KIT-098. | Outcome: wiring; AC "Hand-off" |
| KIT097-F17 | minor | **Closed in ticket.** The technical-replicate option is enabled only with replicates; otherwise it shows the classed refusal. | Outcome: summary; AC "Technical-replicate option" |
| KIT097-F18 | minor | **Closed in ticket.** The reused helpers are named. The scaffold is named as the lipidomics template's structure plus the unchanged shared `analysis_record_SI.rmd`. | Outcome: report, reuse |
| KIT096-F13 (routed in) | minor | **Closed in ticket.** The e2e helper (`helper-e2e-browser.R`: `.E2E_VALID_OMICS`, `e2e_input_id`) and its harness test are in scope. The expectation at `test-e2e-browser-harness.R:120` becomes a positive transcriptomics id, and the unsupported case moves to a fictitious omic. | Outcome: e2e harness; `files_to_modify`; Simulation Setup |

## Not closable here

None. Five closures depend on producer tickets carrying the matching repair. The KIT-097 criteria will fail visibly if those owners do not deliver:

| Finding | Owner | What it must deliver |
|---|---|---|
| F01 | KIT-095 | the edgeR QL model in the manifest and record |
| F03 | KIT-093 | the BH `FDR` column name |
| F04 | KIT-095 | the k and percentage parameters and outcomes |
| F06 | KIT-093, KIT-095 | the renderers and the products |
| F08 | KIT-095 | the export contract |
