# Review of the KIT series after proteomics closed, and the stage-3 re-cut (2026-10-28)

User request: after KIT-065 closed proteomics on the kit, assess whether stage 3 (lipidomics and metabolomics; KIT-030, 031, 032) is sufficiently scoped given everything the series changed since the stubs were written, then scope and plan it. Read-only on the code; the outcome is a re-cut of stage 3 (tickets in `tickets/kit`, gitignored; this file and the plan rows are the tracked record). Same form as [the stage-2 re-cut](REVIEW-20261019-stage2-recut.md).

## What stage 3 inherits

The closure statement in [kit-065/outcome.md](kit-065/outcome.md): four proteomics routes end to end through the GUI to the report on one build, every state verified against a record or the KIT-005 oracles, the owner attached from session start (KIT-048), every step a declared operation on saved entities, the record and its replay exact (KIT-028), the capability fallbacks identical (KIT-044), the omic-agnostic names for the 25 owner, scope, commit and display primitives (KIT-045). Still carried: the kit's state-summary rule as a cross-build exactness test (KIT-050), enrichment on real data (KIT-066), the report templates' shared logic (KIT-067), the static audit's lipid/metab findings ([kit-064/AUDIT-ALIGNMENT-20261028.md](kit-064/AUDIT-ALIGNMENT-20261028.md)).

## Why the stubs could not start as written

| Assumption in the stubs (rev 1.0–1.2) | The tree on 2026-10-28 |
|---|---|
| Control `kit_036`; the journey ends at DA | Control `kit_065` (InuLC5); every journey ends at `report` (KIT-065); the driver dispatches on step kind |
| "Declarations only" on an existing kit | The kit's routes are hard-coded to `DIA`, `protein`, `DIA_LIMPA` in `kitJourneySteps` (`tools/profiling/workflow_kit.R:15`), the case binding's step-map choice (`:162-171`), the driver's tab prefixes `proteomics_workflow-*` and its import entry, the coordinator, both recorders, the three step maps and every bindings block; a new omic is two new routes through all of them |
| A corpus exists | None: the kit binds the synthetic DIA-NN corpora, the MaxQuant/PD fixtures and the KC LFQ lane; the only lipid/metab inputs anywhere are the e2e lanes (`tests/testdata/e2e`: 10 lipids or 3 metabolites × 6 samples, two groups, no batch, no technical replicates, no internal-standard rows in the files the manifest uses) and hard-coded module-CI data frames |
| The owner is general after KIT-045 | The 25 renamed primitives are; the layer around them is not: `setWorkflowType` allows `DIA, DIA_LIMPA, LFQ, TMT` only (`R/utils_workflow_dia_state.R:1383-1389`); `readWorkflowStateMetadata` reads Peptide/ProteinQuantitativeData only (`:221-224`); `createDIASessionOwner` carries the `PROTEOMICS-pending` identity and `workflowRoutesByType` four types (`R/utils_workflow_dia_sessions.R:798-879`); `buildWorkflowAnalysisRecord` and `replayWorkflowAnalysisRecord` assume `setupDirectories(..., "proteomics")`, accept proteomics input roles only and abort on any operation outside four branches (`R/utils_workflow_analysis_record.R:167, 426-453, 541`); the worker publication keys `setupDirectories(...)$proteomics` (`R/utils_workflow_execution.R:2578-2770`); the harness finds `.dia-native-*` stores only, reads `peptide_data` as the input component and gates fingerprints on the DIA checkpoint schema (`workflow_kit.R:998, 1078-1084, 1189`) |
| The inventory of 2026-10-14 | Refreshed below: 26 `saveState()` sites, globals in 26 files (nine of them kernels) |
| Not mentioned | The audit carry-overs; the intensity-filter `ceiling()` held "until their stage" (2026-10-09); the record's per-operation replay branches; KIT-066/067 boundaries; 4 h estimates |

## The two workflows as they stand (read 2026-10-28)

Twins: six tabs (`import, design, qc, norm, de, summary`), module id `<omic>_workflow` (`R/app_server.R:589`), state names one to one (`lipid_*` / `metab_*`), `WorkflowState$new()` at `R/mod_lipidomics.R:152` and `R/mod_metabolomics.R:152` handed to every submodule, no owner, no `usesDIAReferenceStates` gate. The assay is a list of wide tables, one per assay (`lipid_data` / `metabolite_data`, declared with the design in `R/utils_workflow_class_declarations.R:37-42` since KIT-003); every kernel loops over the assays.

| Pattern | Lipidomics | Metabolomics | Stage-3 treatment |
|---|---|---|---|
| `saveState(` with a full S4 | 13 sites | 13 sites | reference commits of the declared class (KIT-070) |
| `.GlobalEnv` writes in modules (`config_list`, `contrasts_tbl`, `project_dirs`) | design builder, design import, summary (4 files) | design, design import, summary (3 files) | onto the workflow bus (KIT-030) |
| `.GlobalEnv` reads / `assign` in kernels (`project_dirs`, `omic_type`, `experiment_label`, `filtering_progress_*`) | `func_lipid_qc_filtering_helpers.R:264-306`, two progress helpers | `func_metab_qc.R:298-340`, progress helpers | adapters; kernels not rewritten (KIT-030) |
| `get0(..., .GlobalEnv)` helper lookups | `func_lipid_s4_objects.R:247`, `func_lipid_s4_plotting_methods.R:125,296`, `func_lipid_s4_da_plot_methods.R:11,74`, `func_lipid_s4_da_analysis_methods.R:159` | none | resolved by namespace (KIT-030) |
| `current_s4_object` holders | DA session helpers | DA observer helpers | references (KIT-030) |
| RDS sessions (`*_filtered_session_data_*.rds`, `session_state_*.RDS`, `*_final_s4.RDS`) | norm export, DA load, summary | same | owned exports; Open adopts (KIT-071) |
| DA result | `runLipidsDA` returns a plain list; the S4 result class is never built by the GUI | same (`runMetabolitesDA`) | a declared output manifest (`lipid.da_results.v1` / `metab.da_results.v1`, the `dia.da_results.v1` pattern), not the S4 classes |
| Enrichment tab | none | none (two uncalled helpers in the multi-omics file) | KIT-066 corrected |
| Report | `lipidomics_report.rmd` by hard-coded name; `RenderReport` | `metabolomics_report.rmd` | the report step as KIT-065's; templates' shared logic KIT-067 |
| `checkParamsObjectFunctionSimplify` in kernels | 62 call sites across both | — | kernel contract, adapted to (README) |

### The step set and the unit as read (to be confirmed in the step maps KIT-068 authors)

| Step | Kernel (lipid / metab) | State | Unit |
|---|---|---|---|
| import | `importLipidSearchData` or `importLipidMSDIALData` / `importMetabMSDIALData` | none (files) | store |
| design (builder or import) | `createLipidomicsAssayData` / `createMetaboliteAssayData` | `*_raw_data_s4` | store |
| intensity filter | `lipidIntensityFiltering` (`func_lipid_s4_normalization_methods.R:868`) / `metaboliteIntensityFiltering` (`func_metab_s4_qc_methods.R:21`) | `*_intensity_filtered` | row after a `global_statistics` pre-pass: the quantile over all intensities of the assay (`ceiling()` removed, user 2026-10-28) |
| duplicates | `resolveLipidDuplicateFeaturesByIntensity` / `resolveDuplicateFeaturesByIntensity` | `*_duplicates_resolved` | feature group (the id column) |
| internal-standard QC | CV per standard, inline | none | products |
| QC finalise | validation | `*_qc_complete` | store |
| internal-standard normalisation (optional) | `normaliseUntransformedData(method = "ITSD")` | `*_itsd_norm` | complete sample after a `global_statistics` pre-pass (the mean of the per-sample factors; standards left unnormalised). This resolves KIT-031's pre-recorded first expected gap without a stage-1 addendum |
| log2 | `logTransformAssays` | `*_log2` | row (user offset; no global statistic) |
| between-sample normalisation | `normaliseBetweenSamples` (limma) | `*_normalized` | whole matrix |
| RUV-III | per-assay optimisation (`getNegCtrlMetabAnova`, `ruvCancor`, `findBestKElbow`) and `ruvIII_C_Varying` | `*_ruv_corrected` or skipped | whole matrix |
| correlation filter or skip | `pearsonCorForSamplePairs`, `filterSamplesBy*CorrelationThreshold` | `*_correlation_filtered` / `*_norm_complete` | whole matrix (DIA's determination) |
| session export | RDS writers | none | store (owned export) |
| DA load and run | `runLipidsDA` / `runMetabolitesDA` | manifest | whole matrix |
| report | `copyToResultsSummary`, `RenderReport` | none | completion-checked |

Multi-assay note: the partitioners work on one table; the per-assay loop is the kernel's. Stage 3 declares each step's unit as read and runs the partitioned units projected (`runProjected`, the `whole` variant the protein intensity filter uses today). A partitioned variant is not needed for exactness on these table sizes and is not built.

Two kernel observations recorded for the step maps, not repaired: the intensity filters and the duplicate steps take every numeric column as a sample column rather than the design's sample list (a numeric annotation column would enter the quantile); the GUI never calls the S4 `resolveDuplicateFeatures` methods (it calls the by-intensity helpers).

### Audit carry-overs (from the 2026-10-28 alignment)

BUG-10 (`mod_lipid_qc.R:134`, `mod_metab_qc.R:142`: self-invalidating observers), BUG-16 (`mod_metab_norm_render_helpers.R:86`: `1:n` composite), BUG-17 (`func_lipid_da_plotting.R:196-201`: `min(numeric(0))` is `Inf`, every row dropped when all FDR are zero; reachable), BUG-21 (the `.GlobalEnv` writes above), BUG-23 (`sapply` in subscript contexts, lipid/metab sites).

## Decisions (user, 2026-10-28)

- The lipid and metab intensity-filter `ceiling()` on the dataset-wide quantile is removed at stage 3, in the ticket that declares each filter, before any oracle is recorded (the proteomics rule of 2026-10-09).
- The kit's two routes are owned by one explicit kit-extension ticket (KIT-068) before the declaration tickets; KIT-031 remains the completeness test for the kit's functions (owner, runner, partitioners, record, store), not its tooling.
- Oracles are recorded against a new synthetic generator (`tools/profiling/kit_omics_corpus.R`, the `kit_protein_corpus.R` pattern) with batch, technical replicates and internal-standard rows; the e2e lanes remain fixtures.
- Lifecycle members only: one screen member per route through the GUI to the report, exact against the in-process record; no counted memory members (the series' memory claim is proteomics; lipid and metab tables are small). Complete-journey qualification stays with KIT-050.
- Taken here from the reading: the DA output is a declared manifest, not the S4 result classes the GUI never builds; KIT-066's lipid/metab items are re-worded (no enrichment GUI exists for either omic).

## The stage-3 cut (execution order)

| Ticket | One mechanism | Depends on | Est. |
|---|---|---|---|
| KIT-068 (new) | Stage-3 kit extension: the `lipidomics` and `metabolomics` routes in the harness, driver and coordinator; the synthetic corpus generator; the two step maps; bindings; in-memory lifecycle members that prove the instrumentation through the report | KIT-065 | 12 h |
| KIT-030 rev 2.0 | Globals onto the workflow bus and holders to references in the lipid and metab modules and kernel adapters (the KIT-013/012 patterns) | KIT-068 | 6 h |
| KIT-069 (new) | The audit carry-overs (BUG-10/16/17/23) with failing-then-passing tests, the KIT-064 pattern, before any oracle exists | KIT-030 | 6 h |
| KIT-070 (new) | The owner for lipidomics and metabolomics: the owner layer made omic-general and attached from session start in both modules; every `saveState` a reference commit; the in-memory path retired | KIT-069 | 12 h |
| KIT-071 (new) | Sessions: the RDS exports as owned exports, the DA Open adopting legacy RDS sessions, the summary's final state from the head | KIT-070 | 8 h |
| KIT-031 rev 2.0 | Lipidomics declarations (the kit completeness test), record, replay, recorder, the screen lifecycle member exact | KIT-071 | 12 h |
| KIT-032 rev 2.0 | Metabolomics declarations mirroring KIT-031; stage closure (both replays exact, check at baseline, the closure suite) | KIT-031 | 8 h |

Rule kept from the stubs: KIT-031/032 add declarations, registry entries, record branches and module request/poll code only; a kit function gap found there stops the ticket and returns as an addendum (KIT-070 is expected to have closed the known ones). KIT-040 keeps its dependency on KIT-031; KIT-050 and KIT-066 on KIT-032.

## Claim limits

No memory claim for lipidomics or metabolomics (no counted members); exactness is against records made on the same build, on synthetic corpora; the state-summary rule (KIT-050) applies to any cross-build comparison; enrichment and the templates' logic are KIT-066 and KIT-067; the two kernel observations above are recorded, not repaired.
