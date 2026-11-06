# KIT-097 plan (rev 1.1), written 2026-11-06 before any code

Bounded item: the transcriptomics normalisation, DA and summary tabs, built by giving the lipid norm, DA and summary
modules the `profile` argument KIT-096 introduced (`NULL` = lipidomics, unchanged), the results-summary branch and the
transcriptomics report. Exit: the acceptance criteria with test evidence, the regression checks (test first, targeted,
the lipid/metab/proteomics module CI and declaration tests unmodified, the KIT-082 probe, the dependents grep, `R CMD
check` at baseline, install) and `outcome.md`. Baseline HEAD `c07809ed`. No kit journeys (KIT-098).

## State sequence and data flow (F02)

`trx_qc_complete` -> `trx_tmm_logcpm` (product) -> `trx_ruv_grid` (product; `automatic` or `manual`) -> `trx_ruvs`
(state, the DA head) -> DA (`trx_da_contrast` per contrast on one fit, the manifest accepted on the head; tables and
publication) -> `trx_ruviii_export` (product after DA, launched from Summary & Export, `M` = `group` or
`technical_replicate`; accepted on the head under `snapshot$exports`) -> `trx_integration_export` (the RUVIII `newY` TSV).
Managed: one runner action each (`omicsManagedProduct` / `omicsManagedCommit`), as KIT-096's filter. In process: the
KIT-095 wrappers (`runTrxReferenceRuvs`, `runTrxReferenceDaAnalysis`, `runTrxReferenceRuviiiExport`) and the same
publication operations called in process (a trx-local runner that copies their staged files to the experiment paths).

## Reuse map (F09)

| What | How |
|---|---|
| `installWorkflowOmicsActions`, `omicsManagedCommit`, `omicsManagedProduct`, `workflowOmicsDisplay`, `workflowOmicsVolcanoPageFactory`, `runWorkflowOmicsManagedDa`, `reportWorkflowOmicsDaExport`, `exportWorkflowAnalysisRecordBundle`, `RenderReport` | reused as is |
| `mod_lipid_norm_ui`/`_server` | `profile`: the options panel is the profile's (TMM/logCPM, RUVs grid controls; no ITSD/log2/between-sample/RUV-III/skip); the QC tabset keeps RUV QC, PCA, RLE and Log, one assay row, the profile's stage labels; the server is the profile's (`mod_trx_norm_server`), which reuses the lipid norm state, log, image and RUV QC renderers |
| `buildLipidNormQcTabsetPanel`, `registerLipidNormStaticQcImageOutputs`, `buildLipidNormQcImageRenderer`, `buildLipidNormRuvQcUiRenderer`, `buildLipidNormLogRenderer` | `profile`/new defaulted arguments (stage prefixes, QC directory key, which tabs); lipid defaults unchanged |
| `mod_lipid_da_ui`/`_server` and the DA helpers | `profile`: the model text and realised formula in place of the formula input, "FDR (BH)" threshold label, the fold-change field labelled display-only, no session load, no heatmap; the long table, id/label/p-value/significance columns and labels; the runner omic and binding, the reference run, the contrasts (`trxContrastsTable(full_format)`), the Glimma widget (`trxVolcanoWidget` in process, the worker's page managed), the static volcano (KIT-093's `trxVolcano`), the in-process publication |
| `mod_lipid_summary_ui`/`_server` and helpers | `profile`: the RUVIII export panel (the `M` choice), the integration export (the RUVIII TSV), the final-state reader, the report template's package directory and file name, test ids |
| `func_lipid_da_plotting.R` | not edited (the trx widget and static volcano are injected) |
| `utils_workflow_omics_actions.R` | not edited unless a shared runner needs an omic parameter |
| metabolomics and proteomics modules | not edited |

## Files

| File | Change |
|---|---|
| `R/mod_trx_norm.R` (new) | the options panel, `mod_trx_norm_server` (managed and in-process sequence, RUV QC outputs, figures, reset), the in-process publication runner |
| `R/mod_trx_da.R` (new) | the DA profile entries: static volcano, widget, run hooks, contrast mapping, in-process publication |
| `R/mod_trx_summary.R` (new) | the RUVIII export panel and server, the integration writer, the final-state reader |
| `R/mod_transcriptomics.R` | the three tabs (`norm`, `de`, `summary`; test ids `trx-tab-norm`, `trx-tab-de`, `trx-tab-summary`), stepper, tab-status chain, profile entries |
| `R/mod_lipid_norm_ui.R`, `_ui_helpers.R`, `_server.R`, `_support_helpers.R` | `profile` / defaulted arguments as above |
| `R/mod_lipid_da_ui.R`, `_server.R`, `mod_lipid_da.R`, `_analysis_helpers.R`, `_results_helpers.R`, `_display_helpers.R`, `_session_helpers.R` | `profile` / defaulted arguments as above |
| `R/mod_lipid_summary.R`, `_workflow_helpers.R`, `_state_helpers.R` | `profile` / defaulted arguments as above |
| `R/func_general_filemgmt_results_summary_helpers.R` | the transcriptomics branch (the output map) and the DA workbook from `de_results_*.csv` |
| `inst/reports/transcriptomics/transcriptomics_report.rmd` (new) | the lipid template's structure (params, `results_summary_dir`, the record overlay), transcriptomics text |
| `DESCRIPTION` | Collate gains the three new files (outside the list; required to load them, as KIT-095/096) |
| `tests/testthat/helper-module-ci-trx-tabs.R`, `test-module-ci-trx-{norm,da,summary}.R` | MCI-030 norm, MCI-031 DA, MCI-032 summary |
| `tests/testthat/test-workflow-omics-application-reads.R` | the transcriptomics norm/DA/summary probe with its positive control |
| `tests/testthat/helper-e2e-browser.R`, `test-e2e-browser-harness.R`, `test-e2e-trx-journey.R` (new) | `transcriptomics` in the e2e helper (F13), the harness expectation, the one e2e journey |
| `tests/testdata/e2e/report_templates/transcriptomics_report.rmd` | not needed: the summary reads the template from the package, as lipid does |

Not changed, with reasons: `tests/testdata/mci/manifest.json` (`test-module-ci-manifest.R`, outside the file list, pins
the transcriptomics scenario set, so new scenarios would fail it: hand-off) and `tests/testdata/e2e/manifest.json` /
`trx_small/` (`test-e2e-manifest.R` pins nine lanes; the e2e reads the KIT-040 generator's 12-sample `trx_edge_12`
corpus, hash-checked as the module CI does).

## Output map as built (F06)

Each notebook output is mapped to the producer KIT-095 supplies; the summary branch copies them into
`results_summary_dir`; the report reads them there. Outputs no producer supplies (the Glimma pages as files under
`Interactive_Volcano_Plots/`) are returned to KIT-095 and named in the outcome.

## Tests first (F12)

The behaviour tests are written and run before any `R/` change and must fail on their assertions: the record's RUVIII
`M` outcome after the summary tab's export; the record's DA model fields after the DA tab; the `k`/percentage mode
fields after the norm tab; the results-summary file set; the e2e helper's transcriptomics id.

## Budget

Targeted files one process each on the kit_096 library; the e2e once, 1,800 s ceiling, `CHROMOTE_CHROME=/usr/bin/brave`
(the kit runner's), any skip a failure; one bounded correction pass.
