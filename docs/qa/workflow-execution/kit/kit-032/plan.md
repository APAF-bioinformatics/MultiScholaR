# KIT-032 plan (the planning gate): metabolomics declarations and the stage-3 closure

Authored 2026-10-31 from KIT-031's tree (e4572f01) and outcome. The metabolomics route mirrors the lipidomics one
(`kit-001/step-map-metab.json`); the metab module splits each step into finer helpers than the lipid module, so the
module sites differ while the declarations are the twins. Kernels untouched except the `ceiling()` (user decision
2026-10-28). One mechanism: the metabolomics steps as declarations on the kit; then the stage-3 closure.

## Item 1: the declarations, the module sites, record, replay, overlay

| What | Where | How |
|---|---|---|
| Shared wrapper core | `R/utils_workflow_lipid_declarations.R` | the omic-neutral helpers get neutral names with the lipid names kept as aliases (one body; the KIT-045 alias pattern): `runWorkflowOmicsReferenceOperation`, `markWorkflowOmicsStoreStep`, `commitWorkflowOmicsStoreStep`. Lipid modules unchanged |
| Twelve declarations | `R/utils_workflow_metab_declarations.R` (new) | on `MetaboliteAssayData`, `fields = metabolite_data, design_matrix`: `metab_intensity_filter` (row; per-assay quantile pre-pass; called as `metaboliteIntensityFiltering`, the module's call name, so `@args` keep the key), `metab_duplicate_resolution` (feature group; `resolveDuplicateFeaturesByIntensity` per assay as `resolveMetabDuplicateAssayData` does), `metab_qc_finalize` / `metab_norm_complete` / `metab_reset` (store), `metab_itsd_normalisation` (complete sample; mean-factor pre-pass), `metab_log2` (row), `metab_between_sample_normalisation`, `metab_ruv` (`runPerAssayRuvOptimization` + `extractBestKPerAssay`/`extractCtrlPerAssay` + `ruvIII_C_Varying`), `metab_pearson_pairs`, `metab_pearson_selection` (`filterSamplesByMetaboliteCorrelationThreshold`), `metab_da_output` (`runMetabolitesDA`, product-only; manifest `metab.da_results.v1`, long table `da_metabolites_long`). Wrappers `runMetabReference*`, `recordMetabDesignProvenance` (inputs `metab_assay_<n>` from `import_sources/`, `design_file`), `runMetabReferenceDaAnalysis`, `replayMetabolomicsStep` |
| Registry | `R/utils_workflow_operations.R` | the in-process loop over `metabolomicsWorkflowOperations()` with `omic_type = "metabolomics"` and the QC tracking labels; defaults (the shipped `[metaboliteIntensityFiltering]` keys for the filter; the shared keys elsewhere), labels; `workflowSectionOperationsByType$METABOLOMICS` (`metaboliteIntensityFiltering`, `ITSDNormalization`, `ruvIII_C_Varying`) |
| Output declaration | `R/utils_workflow_class_declarations.R` | `metab.da_results.v1` (`isMetabResultManifest`, operation `metab_da_output`, long table `da_metabolites_long`) |
| Module sites | `mod_metab_qc_intensity.R` (the apply observer), `mod_metab_qc_duplicates_server_helpers.R` (`runMetabDuplicateResolutionWorkflow` path: prepare + apply), `mod_metab_qc_s4_server_helpers.R` (`saveMetabQcS4CompletedState`), `mod_metab_norm_pipeline_helpers.R` (the ITSD, log2, between-sample steps and the RUV progress shell: one `metab_ruv` call replaces optimisation + correction on the owner), `mod_metab_norm_session_helpers.R` (correlation: `runMetabNormApplyCorrelationWorkflow` + outcome; skip `completeMetabNormSkipCorrelationState`; reset marks its provenance and keeps its `saveState` of the pre-normalisation object), `mod_metab_da_observer_helpers.R` (`runMetabDaAnalysisObserverShell`), `mod_metab_design.R` and `mod_metab_design_import_helpers.R` (design provenance and inputs) | on the owner the declared operation; a test double keeps the direct kernel path (the KIT-031 rule) |
| `ceiling()` | `R/func_metab_s4_qc_methods.R:102-105` | `quantile(...)[1]`; the effect recorded (module-CI, e2e lane, screen corpus) |
| Record / replay / overlay | `R/utils_workflow_analysis_record.R` | outcomes for `metab_ruv` / `metab_itsd_normalisation` (the lipid branches by operation regex); replay branch `omic_type == "metabolomics"`; the DA dispatch by `da$operation`; the overlay's metab sources appended to the candidate lists |
| Summary and template | `mod_metab_summary_observer_helpers.R` (save args: `analysis_record.json` beside `study_parameters.txt`; generate report: the bundle to `Study_report/`), `inst/reports/metabolomics/metabolomics_report.rmd` (the overlay after `ExtractStudyParametersNew`) | the KIT-031 pattern |
| Harness | `tools/profiling/workflow_kit.R` | `kitOperationIdentityProducts` matches `(lipid|metab)_*_science` |
| Tests | `tests/testthat/test-workflow-metab-declarations.R` + helper (the KIT-031 test's shape on a 24-sample, two-assay metab fixture), the registry lists, the touched metab module tests on their contracts | |

Item-1 gates: the metab module-CI lanes hashed before/after (identical apart from store/time-bearing files and the
ceiling's effect, recorded); the lipid declarations test, the KIT-028 replay and the LFQ/TMT/DIA_LIMPA records
unchanged; R CMD check at baseline.

## Item 2: the recorder, the exact member, the install

`kit-032/oracle-record.R` on the KIT-031 recorder's pattern: the metab import through `prepareMetabImportAssaySelectionState`
and `runMetabImportProcessing` (the driver's inputs), the design through the metab builder server driven with
`shiny::testServer` (the driver's factor, metadata, technical-replicate and contrast actions) and the module's
results flow, then the wrappers with the tab defaults (integers where the tab sends integers; the ITSD table's
preselected candidates), the DA, the report. Bound as `metab32_<step>`; ticket `032m` (kit-journey-032m, route
metabolomics, the lifecycle normalisation settings). Install `kit_032`; the recorder on the install; the member exact;
the member's record valid and its replay exact; the two KIT-068 lifecycle members.

## Item 3: the stage-3 closure (the stage's checkpoint)

On the final build: the KIT-031 member's record replayed exact under kit_032's library and the KIT-032 member's
record replayed exact (`member-record-replay.R`, generalised by omic); the KIT-028 DIA replay exact; the DIA_LIMPA,
LFQ and TMT records identical to the kit_065 oracles; R CMD check at baseline; the closure full suite on a clean
worktree of the final commit (green by the rule recorded since KIT-017); the stage-3 closure statement in
`kit-032/outcome.md` (the omics on the kit; claim limits: no memory claim for these omics, synthetic corpora,
enrichment and the report templates' shared logic excluded, the counted complete journeys with KIT-050).

## Stop-and-addendum triggers

As KIT-031's: any need for a new role, enum, partitioner or harness function. Budget: one bounded correction pass
per gate; a valid unfavourable result recorded, not retried.
