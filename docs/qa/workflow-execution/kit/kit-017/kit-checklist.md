# KIT-017: kit completeness checklist (DIA stage closure)

Inventory taken at `5250c88c` (the working tree of the closure; the candidate build `kit_024` is `28075aa1`, which differs from it only by KIT-034's test repairs and the KIT-024 closure record). Each row names the file and the reading; nothing here is a performance claim (those are in [outcome.md](outcome.md)).

## A. Every kit function exists, with no DIA logic in its body

| Kit function (plan table) | Where | DIA names in the body | Reading |
|---|---|---|---|
| `workflowClassDeclaration(class)`, `workflowOutputDeclaration`, `validateWorkflowClassDeclaration`, `workflowClassRequiredComponents`, `workflowDeclarationField`, `workflowAssayDimensions` | `R/utils_workflow_class_declarations.R` (146 lines) | The registry's entries name the two proteomics classes (`PeptideQuantitativeData`, `ProteinQuantitativeData`, describer `describeDIAProteinState`) and the DA manifest output (`dia.da_results.v1`, producer `runDIAReferenceDaAnalysis`, identifier `isDIAResultManifest`) | pass: those are the declarations the registry holds, not logic in the functions; the lookup, validation and dimension helpers read the declaration only |
| `saveComponentCheckpoint` / `readComponents(fields, columns)` | `checkpointDIAExecutionResult`, `readDIAExecutionComponents`, `saveDIAStateRecord`, `materializeDIAStates` in `R/utils_workflow_dia_state.R`, declaration-driven since KIT-002 (schema `dia.component_checkpoint.v2`, additive `declaration` field) | the functions keep the `DIA` prefix in their names by history; their bodies take the assay, design, table and list-of-tables slots from `workflowClassDeclaration` | pass with a recorded naming debt: the storage primitives (`readDIAExecutionComponents`, `saveDIAOwnedValue`, `saveDIATableComponent`, `rootDIAProduct`, `inheritDIAExecutionSource`) are omic-agnostic in behaviour and DIA-named. A rename is unrelated renaming under AGENTS.md and is not done here; stage 3's first consumer (KIT-030) decides whether to alias them. |
| `runSavedOperation(manager, declaration, source, parameters)`, `runProjected`, `newWorkflowOperation`, `validateWorkflowOperation` | `R/utils_workflow_operations.R` | five calls to the DIA-named storage primitives above; no DIA class, column or kernel name | pass |
| `partitionBySamples`, `partitionByFeatureGroup`, `partitionByRows` (+ `completeGroupBatchBound`, `chooseCompleteGroupBatches`) | `R/utils_workflow_partitions.R` (182 lines) | none (two roxygen lines mention "intermediates") | pass |
| `partitionByPlex` | absent (0 references in `R/`, `tools/`, `tests/`) | — | **exempt until KIT-020** (plan) |
| `registerWorkflowOperation`, `workflowOperationNames`, `workflowOperationEntry`, `assertWorkflowOperationDispatchable`, `registerBuiltinWorkflowOperations` | `R/utils_workflow_operations.R` | the built-in registrations name the proteomics operations and their declaration constructors (`dia*Operation()`); the registry functions themselves name none | pass: the registry is the one place an operation is named (KIT-007); every declaration is registered, the in-process ones with role `in_process` and dispatch refused (KIT-024) |
| request/poll from the modules | `runWorkflowAction`, `pollWorkflowAction`, `beginWorkflowGUIAction`, `installWorkflowManagedStep` in `R/utils_workflow_execution.R` | the worker-dispatched steps (KIT-008 P14 and sample filter; KIT-014/015 peptide and protein QC; KIT-019 DA publication) go through the managed step | pass |
| on-disk state manager as the default `state_manager` | `isDIANativeDefaultEnabled()` is `TRUE` (`R/utils_workflow_dia_sessions.R`); DIA-NN import attaches the owner (KIT-011) | DIA only: LFQ/TMT (KIT-020) and the list-assay omics (KIT-030) are stage 2 and 3 | pass for the DIA stage; the other omics are the stages' work, not a gap |
| `savedReactiveVal`, products, paged tables | `savedReactiveVal <- newDIASavedReactiveVal` (`R/utils_workflow_dia_state.R`); `saveDIADisplayTable`, `readDIADisplayTablePage`, `diaDisplayTableFilter` (`R/utils_workflow_dia_products.R`, KIT-018/024) | DIA-named, declaration-driven | pass; same naming debt as the storage primitives |
| harness members bound to any registered operation, per-step populations and oracles, journey driver | `tools/profiling/workflow_kit*.R`; protocol v1.4 section 4 names every member kind | — | pass |

## B. Every DIA journey step is a declaration

The step map's three `store` steps (import, design, session export) are not data kernels (`kit-001/step-map.json`, unit `store`) and are the state manager's own work. The 16 data steps:

| Step (driver id) | Declaration (registry name) | Unit | Compute role | Partitioner | Admission | Measured by |
|---|---|---|---|---|---|---|
| qvalue | `peptide_qvalue_filter` (`diaQvalueFilterOperation`) | feature_group | peptide (worker) | `partitionByFeatureGroup` through `produceDIAConfidenceTables` | streamed-page estimate, refused above the cap (KIT-022, closes KIT-014 PH-02) | KIT-014 single-1 |
| protein_peptide | `peptide_protein_evidence_filter` | feature_group | peptide | `partitionByFeatureGroup` | declared row bytes | KIT-014 single-2 |
| rollup (P14) | `precursor_rollup` (`diaPrecursorRollupOperation`, variant `p14`) | complete_sample | peptide | `partitionBySamples` (SQL aggregate, exactness-proven kernel) | declared row bytes | KIT-006, KIT-008 gui member |
| intensity | `peptide_intensity_filter` | feature_group | peptide | `partitionByFeatureGroup` | declared row bytes | KIT-014 single-3 |
| sample | `minimum_peptides_per_sample` (`diaSampleFilterOperation`) | complete_sample | peptide | `partitionBySamples` | declared row bytes | KIT-006 |
| replicate | `peptide_replicate_filter` | feature_group | peptide | `partitionByFeatureGroup` | declared row bytes | KIT-014 single-4 |
| imputation | `peptide_imputation` | feature_group | peptide | `partitionByFeatureGroup` | declared row bytes | KIT-009, KIT-014 |
| iq | `protein_iq_rollup` (`diaProteinIQRollupOperation`) | feature_group | protein | per-group working files, publication through the worker | declared row bytes | KIT-015 single-1 |
| accession | `protein_accession_cleanup` | row | protein | `partitionByRows` | declared row bytes | KIT-015 single-2 |
| protein_intensity | `protein_intensity_filter` | whole_matrix | protein | `runProjected` under the cap | projected object bound | KIT-015 single-3 |
| duplicates | `protein_duplicate_removal` | feature_group | protein | `partitionByFeatureGroup` (wide, KIT-010 sample projection) | declared row bytes | KIT-015 single-4 |
| protein_replicate | `protein_replicate_filter` (`diaProteinReplicateFilterOperation`) | row | protein | `partitionByRows` | declared row bytes | KIT-010 |
| normalisation_ruv | `normalization-ruv` (`diaNormalizationOperation`) | whole_matrix | in_process (dispatch refused) | `runProjected` | cap | KIT-005 oracle, journey members |
| correlation | `pearson-pairs`, `pearson-selection` | whole_matrix | in_process | `runProjected` | cap | KIT-005 oracle, journey members |
| da | `DA-contrast-output` (`diaContrastOperation`) + `da_publication` (publication through the worker, KIT-019) | whole_matrix | in_process + protein (worker) | `runProjected`; publication in the worker | cap | KIT-019 probe; this closure (first production dispatch of `da_publication` on the 300-sample corpus) |
| enrichment | `enrichment` (`diaEnrichmentOperation`, non-cacheable, service identity on the product) | whole_matrix | in_process | `runProjected` | cap | KIT-016 in-process check |

Every row has a registry entry (`registerBuiltinWorkflowOperations`), a declaration constructor and a unit from the step map. No kernel was rewritten; the one replacement kernel (P14) carries its CORE exactness proof.

## C. Kernel contracts adapted to, not removed (README "Kernel contracts the kit adapts to")

1. Legacy globals at publication: none (KIT-023; boundary tests assert the global environment unchanged at the figure boundary).
2. Call-stack parameter resolution (`checkParamsObjectFunctionSimplify`): the kit calls kernels with explicit arguments under the kernel's own name.
3. The active owner (`getDIAActiveManager()`: innermost `withDIAStateScope` first, the Shiny session second, KIT-022).

These are adapters to kernels that are never rewritten; they are not kit design and not gaps.

## D. Test quarantine at closure (`../quarantined-tests.md`)

- Closed by KIT-034: Q-01, Q-05, Q-09, Q-15, Q-16, Q-20; Q-04 and Q-14 are skip conditions proven on the installed `kit_024` library.
- Open at closure, pre-series fixture and library-drift debt for the stage-2 pre-step: Q-02, Q-03, Q-06, Q-07, Q-08, Q-10, Q-12 (reclassified as fixture), Q-13, Q-17, Q-18, Q-19.
- Q-11 (the DA replicate block): decided 2026-10-15, owned by KIT-035 after this closure (the oracle is re-recorded only if `journey_da` moves).

## E. Record required by the quarantine (user decision 2026-10-14)

KIT-008, KIT-011, KIT-013 and KIT-015 completed with regression checks that ran only the files they touched. The full-suite pass of KIT-024 found tests those tickets had overtaken: Q-01 and Q-16 (KIT-015: the module fixtures and the coordinator tests), Q-09 (KIT-008: the default variant), Q-15 (KIT-013: the globals moved onto the bus), Q-20 (KIT-011: the capacity-session functions removed), and the saved-QC module-edge test repaired in KIT-024. No scientific result changed; the tests had stopped describing the code. The standing rule from KIT-024 on: every ticket's regression check runs the whole suite once before completion and classifies any new failure at the previous commit.

## F. Measured at closure

Filled from stage `kit-journey-017` in [outcome.md](outcome.md) after the counted pair runs.

## G. Hand-off

- Stage 2 (LFQ/TMT, KIT-020/021): `ProteinQuantitativeData` from import (no peptide steps; `protein_s4_initial` state); import formats (MaxQuant, FragPipe, PD) adopted into the store; `partitionByPlex` for TMT; the QC tab skipping; report templates; the protein declarations inherited from KIT-015, KIT-004 (normalisation, correlation, DA) and KIT-016 (enrichment); Q-11/KIT-035's technical-replicate declaration applies to LFQ/TMT designs too.
- Stage 3 (lipidomics, metabolomics, KIT-030/031/032): the list-of-tables assay classes (`LipidomicsAssayData`, `MetaboliteAssayData`, KIT-003) on the state manager default; the naming debt in section A (DIA-prefixed storage primitives) is the first thing their consumers meet; the inventory `../omics-inventory-20261014.md` lists the kit idioms those modules already use.
