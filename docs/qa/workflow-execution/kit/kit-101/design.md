# KIT-101 design: proteomics as the pattern, the shared layer's boundary, the omic configuration contract and the family work map

KIT-101 work item W1, 2026-11-07. Written at HEAD `bb87a463` on `performance-v2`. The scoping notes in `../gui-workover-scoping/` (at `c4839d25`/`c783e97a`) were the starting point, and every file:line below was checked again at HEAD by grep, reading and R parse data. One R call was made: `pkgload::load_all` to list the bindings and the trx profile. Nothing was run against data, and no production file changed.

## 1. Verdict

**Proteomics is the pattern only in part.** Three of its pieces are the light observer layer and are the template:
- the QC sub-tab spec tables;
- the QC step functions, each of which makes one call to a core runner;
- the DA run handler: one click, one core call, and the publication as a managed step.

Proteomics' import, design, normalisation and enrichment are not light, and its GUI as a whole is not lighter than lipid's or metab's. It is the largest GUI (24,468 lines), it runs on its own core, it keeps part of its backend in the app process, and the core calls back into 24 of its module functions.

**The boundary.**
- The shared layer is built on the omics binding and its runner (`installWorkflowOmicsActions`), which lipid, metab and trx already share.
- Proteomics joins through a DIA binding adaptor (user decision). The adaptor maps the shared components' calls onto the existing `runDIAReference*`/`computeDIAReference*` runners and managed steps. Both runner kinds already sit on the same primitive, `installWorkflowManagedStep`, so the adaptor is a mapping, not a second engine.
- Proteomics' route-specific QC sub-tabs plug into the shared QC shell as configured sub-tabs.
- Enrichment is its own omic-agnostic module, with proteomics as its only consumer.
- Proteomics' in-process backend (normalisation, correlation, enrichment, the DA limma kernel) keeps its placement through the adaptor. It moves to the worker in KIT-111 to 113, not in this series.

### Evidence: the parts that are the template

| Claim | Evidence at HEAD |
|---|---|
| Each QC sub-tab's step function makes one call to a core runner | All 17 `R/mod_prot_qc_*` sub-tab files call one `runDIAReference*` or `computeDIAReference*` runner. Example: `mod_prot_qc_peptide_qvalue.R:90` calls `runDIAReferenceQvalue` (`utils_workflow_dia_state.R:2078`). Six protein sub-tabs call `computeDIAReference*` instead, for example `mod_prot_qc_protein_replicate.R:99` calls `computeDIAReferenceProteinReplicateFilter` (`utils_workflow_dia_qc_rollup.R:1537`). |
| Every sub-tab is a managed step | `installWorkflowManagedStep` (`utils_workflow_execution.R:633`) is installed at 18 module sites: 17 QC sub-tab files and `mod_prot_da_handlers_run_analysis.R:11`. With no managed application, the tab runs the same core runner in process. |
| The sub-tabs come from spec tables | `getProtQcPeptideModuleSpecs(workflowType)` (`mod_prot_qc_peptide.R:28`) and `getProtQcProteinModuleSpecs(workflowType)` (`mod_prot_qc_protein.R:70`) return `{id, title, uiFn, serverFn, needsExperimentPaths}` per route. |
| DA run: one click, one core call | `mod_prot_da_handlers_run_analysis.R:47` calls `runDIAReferenceDaAnalysis` (`utils_workflow_dia_da.R:341`). The publication is the managed step `da_publication` (`:11`). `mod_prot_da_server.R` is 59 lines and only wires six handlers. |

### Evidence: the parts that are not light

| Part | Evidence at HEAD |
|---|---|
| Normalisation | `runProtNormNormalizationWorkflow` (`mod_prot_norm_workflow_helpers.R:561-616`) calls `runDIAReferenceNormalization` (`:599`, core at `utils_workflow_dia_products.R:525`), or the DIA_LIMPA pass-through (`:594`). The pass-through, `runDIAReferenceLimpaNormalisationPassThrough`, is a core-style runner defined in the module file at `:365`. The function keeps 13 injected-default parameters it never uses (`:572-582`, `:584-585`), plus an unused `checkMemoryUsageFn` (`:568`). No `mod_prot_norm_*` file installs a managed step. |
| Core calls into the module | 24 proteomics module functions are called from the core (section 5). Example: `utils_workflow_dia_products.R:357` calls `resolveProtRuvParameters` (`mod_prot_norm_ruv_helpers.R:100`). |
| Backend work in the module | the UniProt fetch at design save: `getUniprotAnnotationsFull` (`mod_prot_design_state_helpers.R:169`). FASTA processing: `processFastaFile` (`mod_prot_design_import_helpers.R:127`). The enrichment kernel runs from the module when there is no DA manifest (`mod_prot_enrich_analysis_helpers.R:601-604`, through `:128`). Organism mapping and filtering: `resolveProtEnrichOrganismMapping` (`mod_prot_enrich_input_helpers.R:1-181`) and `applyProtEnrichOrganismFilter` (`:182`). |
| Test doubles in production code | the `ProtEnrichDeterministic*` class and the synthetic results (`mod_prot_enrich_execution_helpers.R:1-60`); the double kernel is at `:126`. |
| Import | all reading, FASTA processing, organism analysis and configuration loading run inside one observer: `registerProtImportObservers` (`mod_prot_import_server_helpers.R:181`, the `process_data` observer). The reader takes the Shiny `input` (`importProtImportDataByFormat`, `mod_prot_import_processing_helpers.R:38`). |
| Its own core | proteomics runs on `utils_workflow_dia_*.R` and `utils_workflow_prot_declarations.R`. `omicsDeclarationBindings()` (`utils_workflow_omics_declarations.R:588-591`) lists only lipid, metab and trx. |

**Where each step runs today:**

| Step | Proteomics | Lipid, metab, trx |
|---|---|---|
| QC | worker (managed step per sub-tab) | worker, through the omics runner |
| normalisation and RUV | **app process** (`runDIAReferenceNormalization`, `utils_workflow_dia_products.R:525`, which calls `runSavedOperation` at `:540`) | worker (`runWorkflowOmicsManagedNormalisation`, `utils_workflow_omics_actions.R:222`) |
| correlation filter | **app process** (`runDIAReferenceCorrelation`, `:895`) | worker (`runWorkflowOmicsManagedCorrelation`, `:315`) |
| DA kernel | **app process** ("the limma kernel runs in process before the request", `mod_prot_da_handlers_run_analysis.R:4-5`) | worker (`runWorkflowOmicsManagedDa`, `:338`) |
| DA publication | worker (`da_publication`) | worker |
| enrichment | **app process** (`runProjected` in `runDIAReferenceEnrichment`, `utils_workflow_dia_da.R:572`, `:581`) | not offered |

**Size (lines, `wc -l` at HEAD; no `mod_` file changed since the scoping):**

| Family | Proteomics | Lipid | Metab | Trx |
|---|---:|---:|---:|---:|
| import | 2,494 (10 files) | 2,251 (6) | 1,940 (7) | in `mod_transcriptomics.R` |
| design | 3,537 (10) | 3,309 (11) | 2,689 (5) | in `mod_transcriptomics.R` |
| QC | 5,186 (20) | 2,419 (7) | 2,709 (12) | 234 (`mod_trx_qc_filter.R`) |
| normalisation | 4,505 (10) | 2,665 (8) | 3,852 (7) | 235 (`mod_trx_norm.R`) |
| DA | 1,845 (9) | 2,569 (8) | 2,453 (6) | 171 (`mod_trx_da.R`) |
| enrichment | 4,882 (10) | — | — | — |
| summary | 1,544 (6) | 1,542 (4) | 1,360 (3) | 143 (`mod_trx_summary.R`) |
| shell | 475 (`mod_proteomics.R`) | 414 (`mod_lipidomics.R`) | 414 (`mod_metabolomics.R`) | 286 (`mod_transcriptomics.R`) |
| **total** | **24,468** | **15,169** | **15,417** | **1,069** |

## 2. The shared layer's boundary

### 2.1 The base

The shared components are new omic-neutral files (user decision). They call:
- the omic's declaration binding (`lipidomicsDeclarationBinding`, `utils_workflow_lipid_declarations.R:17`; `metabolomicsDeclarationBinding`, `utils_workflow_metab_declarations.R:17`; `transcriptomicsDeclarationBinding`, `utils_workflow_trx_declarations.R:702`);
- the omics runner, `installWorkflowOmicsActions` (`utils_workflow_omics_actions.R:13`), which the shell installs once per session (`mod_lipidomics.R:196`, `mod_metabolomics.R:196`, `mod_transcriptomics.R:255`);
- the in-process wrappers when there is no runner (`runWorkflowOmicsReferenceStep`, `utils_workflow_omics_declarations.R:972`, and its siblings at `:983-1094`).

**Today's omics call sequence**, which the shared components keep:
1. The tab asks for the runner: `workflowOmicsManagedActions(omic_type)` (`utils_workflow_omics_actions.R:160`).
2. With a runner, it runs `runner$run(list(omicsManagedCommit(step, config, deliver)), done, failed)` (`:61`, `:172`). The runner is a queue over one `installWorkflowManagedStep` (`:43`).
3. With no runner (no managed application, a test, or the entry-point suite's `testServer`), it calls the in-process wrapper, for example `runWorkflowOmicsReferenceStep(manager, step, config)`.
4. `deliver` receives `readManagedOmicsResult` (or the wrapper's value), and the tab renders and sets `tab_status`.

### 2.2 Proteomics through a DIA binding adaptor

The adaptor is a new binding-shaped list for proteomics. It holds only what the shared components call, and it is filled in family by family: each ticket adds the entries its component needs, and no more. KIT-102 creates it with the shell and import entries, and its validator arrives with the contract's validator. Its call sequence:

1. **The component reads `config$binding$runner`.** It is `"omics_actions"` for lipid, metab and trx, and `"dia_managed_step"` for proteomics.
2. **The component asks the adaptor for the step** by id, for example `qc/qvalue_filter`, `norm/run`, `da/run` or `enrichment/run`. The adaptor answers with today's core runner and, where one exists, today's managed operation name.
3. **The adaptor dispatches as today:**
   - **QC (a plugin).** The sub-tab server is unchanged. It keeps its own `installWorkflowManagedStep(session, input, workflow_data, "<operation>", ...)`. A managed application sends the click to the worker. With no application, the server calls the same core runner in process, for example `runDIAReferenceQvalue` (`utils_workflow_dia_state.R:2078`).
   - **Normalisation.** `runDIAReferenceNormalization` (`utils_workflow_dia_products.R:525`). On DIA_LIMPA it calls the pass-through instead, which moves from `mod_prot_norm_workflow_helpers.R:365` into the core. Then `runDIAReferenceRuvCovariates` (`:824`). The correlation filter is `runDIAReferenceCorrelation` (`:895`). All of these run in the app process, as today.
   - **DA.** `runDIAReferenceDaAnalysis(..., publish = !managed)` (`utils_workflow_dia_da.R:341`) runs in process. Under a managed application, the managed `da_publication` step follows (`mod_prot_da_handlers_run_analysis.R:11`, `:47`).
   - **Enrichment.** `runDIAReferenceEnrichment` (`utils_workflow_dia_da.R:572`), in process. The declaration is registered in process.
4. **The adaptor hands the result back** through the DIA readers: `readDIAStateMetadata` (`utils_workflow_dia_state.R:237`) and `hasDIAState` (`:1572`). It also sets the status keys proteomics uses today (`differential_expression`, `differential_abundance`, `enrichment_analysis`; `utils_workflow_dia_da.R:446-455`).

The adaptor changes no placement and no science. Moving proteomics' backend to the worker happens in KIT-111 to 113, each measured.

### 2.3 Proteomics' QC sub-tabs as plugins

The shared QC shell renders groups, and each group renders a list of sub-tab specs in the shape of `getProtQc*ModuleSpecs`. Proteomics declares two groups, `peptide_qc` and `protein_qc` (`mod_prot_qc.R:123-141`, `:165-191`). The routes give these sub-tabs:

| Route | `peptide_qc` | `protein_qc` |
|---|---|---|
| DIA | qvalue_filter, protein_peptide_filter, rollup, intensity_filter, sample_filter, replicate_filter, imputation (`mod_prot_qc_peptide.R:40-83`) | rollup, cleanup, intensity_filter, duplicate_removal, replicate_filter (`mod_prot_qc_protein.R:27-69`, `:80-84`) |
| DIA_LIMPA | qvalue_filter, rollup, log2_matrix, accession_cleanup (`mod_prot_qc_peptide.R:32-38`) | dpc_quant, normalisation, peptide_count (`mod_prot_qc_protein.R:73-79`) |
| TMT, LFQ | none: the QC tab is bypassed (`mod_proteomics.R:406-416`) | — |

The 17 sub-tab servers are not changed. Their spec names (`uiFn`/`serverFn` strings, resolved today with `exists()`/`get()`, `mod_prot_qc_protein.R:87-112`) become plugin ids in the shared layer's component registry. The contract then holds ids, not functions.

### 2.4 What stays separate, and why

| Stays separate | Why |
|---|---|
| Proteomics' 17 QC sub-tab servers | DIA-specific steps that vary by route, each with its own managed step. No other omic shares them. They plug in as configured sub-tabs. |
| Enrichment as its own module | User decision: omic-agnostic, proteomics its only consumer for now. Lipid and metab enrichment needs new kernels, outside this series. |
| Proteomics' backend placement | User decision: it moves to the worker in KIT-111 to 113. The adaptor keeps today's in-process placement. |
| The import source panel per omic | Vendor formats stay per omic (parity decision). Proteomics has its own source panel (search results, FASTA, DIA route, organism, mappings). Lipid, metab and trx share a feature-table panel. |
| Proteomics' design data model and commit | Proteomics' design acts on a long DIA-owned table (`transformDIADesignTable`, `utils_workflow_dia_state.R:996`) and commits a checkpoint per workflow type (`buildProtDesignStateCheckpoint`, `mod_prot_design_state_helpers.R:1`). It joins the shared builder through a `dia_long` data adapter and keeps its commit as a hook. |
| Per-omic ids and status keys | User decision. The e2e lanes and the frozen kit harness (protocol v1.9, `tools/profiling/workflow_kit.R:20-53`) address them. |
| Report rendering and the proteomics integration write | They stay in process inside core functions (KIT-108, D3). KIT-108's report-memory diagnosis (revision 1.1) may reopen this with the user. |

## 3. What changed since the scoping, and which scoping claims no longer hold

**What changed in `R/` since `c4839d25`:** only commit `53e66b4a` (KIT-100). It touched three declaration files:
- `utils_workflow_lipid_declarations.R` and `utils_workflow_metab_declarations.R`: `default_variants` gained `*_log2 = "partitioned"` (`:37`), and lines after `:35` moved down by 1;
- `utils_workflow_trx_declarations.R`: `default_variants` became the filter and TMM/logCPM partitioned (`:713`), and `trxRunProduct` runs the default variant in process (`:775-782`). Lines after `:711` moved down by 1, and lines after `:775` by 5.

No `mod_*`, `func_*`, `app_*` or other `utils_workflow_*` file changed. The harness moved to protocol v1.9 (KIT-100: three-digit ticket keys).

**Scoping claims that no longer hold, or were never exact:**

| Claim | Fact at HEAD |
|---|---|
| "the core calls 34 module functions (25 proteomics, 5 lipid, 4 metab)" | **32 are called: 24 proteomics, 4 lipid, 4 metab.** A comment-inclusive grep still gives 34. The two extra names occur only in roxygen comments: `handleLipidNormRunNormalization` (`utils_workflow_omics_actions.R:208`) and `runProteinDuplicateRemovalStep` (`utils_workflow_prot_declarations.R:271`). The difference is the counting method, not a code change (section 5). |
| "15 `runDIAReference*` call sites in the modules" | 19 call sites. The module files are unchanged, so 15 was an undercount, which missed names with digits such as `runDIAReferencePeptideLog2Matrix`. |
| "Each of the 17 sub-tab step functions makes one call to a core runner … `runDIAReference*`" | It holds, but 6 protein sub-tabs call `computeDIAReference*` runners, not `runDIAReference*`. |
| "lipid: msdial, lipidsearch, custom (`mod_lipid_import_ui_helpers.R:78-83`)" | Lipid offers six formats: msdial, progenesis, xcms, compound_discoverer, lipidsearch, custom (`:77-84`). The KIT-102 note had this right. |
| "That runner [`runDIAReferenceNormalization`] is defined in the module file itself, at `:365`" | `:365` defines `runDIAReferenceLimpaNormalisationPassThrough`. `runDIAReferenceNormalization` is core (`utils_workflow_dia_products.R:525`). |
| injected kernel defaults: "proteomics 84, lipid 36, metab 30" | Not reproducible, and the scoping recorded no method. By R parse data at HEAD, `*Fn` formals with a function-symbol default number: proteomics 1,147, lipid 646, metab 677, trx 1. Of these, the ones whose default names a `func_*` kernel number 56, 28 and 28. Use these figures, with their method, from now on. |
| `omicsResidencyJourney` at `helper-workflow-omics-residency.R:113-261`; `trxResidencyJourney` at `:268+` | `:112-251` and `:257-330`. The file is unchanged since KIT-096, so the scoping's lines were off. |
| "376 test files, 3,508 `test_that` blocks" | 376 files and 3,516 blocks (grep at HEAD). KIT-098 to 100 added blocks. |
| `runTrxReferenceRuviiiExport` at `utils_workflow_trx_declarations.R:810` (KIT-108 note) | `:815` (shifted by `53e66b4a`). |
| harness protocol v1.8 | v1.9 (KIT-100). The step-map lines are `tools/profiling/workflow_kit.R:20-30` (trx) and `:44-53` (lipid and metab). |

Every other file:line reference used in this design was checked at HEAD and holds.

## 4. The omic configuration contract

### 4.1 Rules

- **One plain list per omic**, built by one function per omic (proposed: `omicGuiConfig(<omic>)`; KIT-102 fixes the name). The list holds:
  - data: strings, numbers, logicals and id vectors;
  - **references to core functions** (`utils_workflow_*`, `func_*`), never server or UI functions;
  - **component ids**, chosen from the shared component library or from an omic's plugin registry. Proteomics' QC sub-tab plugins, the normalisation panels, the DA views and the summary export steps are all named this way. The library maps an id to its UI and server, and that map lives in the GUI layer, not in the configuration.
- **Nested by group:** `shell`, `binding`, `import`, `design`, `qc`, `norm`, `da`, `enrichment`, `summary`. This settles KIT-105's D2 (flat or nested): nested, one group per family ticket.
- **The binding is referenced, not repeated.** `class`, `prefix`, `feature_noun`, `long_table`, `step_labels`, `default_variants`, `qc_dir_role`, `filtering_metrics` and `da_fixed` are read from the binding. They appear below only where the GUI reads them.
- **Status keys, tab values, module ids, tabset ids and test ids stay per omic** (user decision). They are data in `shell$tabs`.
- **Parity decisions are not configuration.** Text-numeric conversion and the QC-finalise alignment check apply to every omic, so they are not fields.
- **A field whose value today is a module function** holds a core reference once its owning ticket has moved that function into the core. The tables mark these "today in `mod_…`; moves in KIT-1nn".
- **The contract is a document in KIT-101.** Its validator and first consumer arrive together in KIT-102, so KIT-101 adds no R code for it.

Type notation: `chr(1)` one string; `chr` a character vector; `lgl(1)`; `num(1)`; `ref` a reference to a core function; `id`/`ids` component or plugin ids; `list` a list with the named fields. "Req." is yes, no (optional), or derived (computed from other fields, never set).

### 4.2 `shell`: identity and the app shell (owner KIT-102)

| Field | Type | Req. | Omics | Where today's value lives |
|---|---|---|---|---|
| `omic_type` | chr(1) | yes | all | lipid `utils_workflow_lipid_declarations.R:18`; metab `utils_workflow_metab_declarations.R:18`; trx `utils_workflow_trx_declarations.R:703`, `mod_transcriptomics.R:33`; prot passed by the dispatcher (`app_server.R:576`) |
| `prefix` | chr(1) | yes | all | lipid, metab binding `:18`; trx binding `:703`; prot: none today (the test-id prefix `prot`) |
| `testid` | chr(1) | yes | all | `mod_proteomics.R:127,184`; `mod_lipidomics.R:65,118`; `mod_metabolomics.R:65,118`; `mod_transcriptomics.R:214,226` |
| `label`, `label_lower` | chr(1) | yes (prot: `label` NULL) | all | `mod_lipidomics.R:45`; `mod_metabolomics.R:45`; `mod_transcriptomics.R:33,217`; lipid's lower-case default at `mod_lipid_design_builder_helpers.R:84`; prot has no header |
| `module_id` | chr(1) | derived | all | `paste0(omic, "_workflow")` (`app_server.R:460,591`) |
| `tabset_id` | chr(1) | yes | all | `workflow_tabs` (`mod_proteomics.R:118`); `lipidomics_tabs` (`mod_lipidomics.R:55`); `metabolomics_tabs` (`mod_metabolomics.R:55`); `transcriptomics_tabs` (`mod_transcriptomics.R:219`) |
| `routes` | chr | yes | all | prot `DIA`, `DIA_LIMPA`, `TMT`, `LFQ` (`workflowTypeForImportFormat`, `utils_workflow_dia_sessions.R:804-812`); `LIPIDOMICS` (`mod_lipid_import_processing_helpers.R:258`); `METABOLOMICS` (`mod_metab_import_processing_helpers.R:284`); `TRANSCRIPTOMICS` (`mod_transcriptomics.R:33`) |
| `tabs` | list of `{status_key, value, module_id, title, icon, testid, step_name, step_icon}`, ordered | yes | all | prot tabs `mod_proteomics.R:122-182`, stepper `:380-386`; lipid `mod_lipidomics.R:59-116`, `:271-276`; metab same lines in `mod_metabolomics.R`; trx `mod_transcriptomics.R:220-225`, `:265-269` |
| `status_keys` (initial `tab_status`) | derived from `tabs` | derived | all | `mod_proteomics.R:230-238`; `mod_lipidomics.R:162-169`; `mod_metabolomics.R:162-169`; `mod_transcriptomics.R:234-235` |
| `chain` (each tab opens when the one before completes) | derived from the `tabs` order | derived | all | `mod_lipidomics.R:291-336`; `mod_metabolomics.R:291-336`; `mod_transcriptomics.R:271-283`; `mod_proteomics.R:399-470` |
| `chain_rules` | list | no | prot | the TMT/LFQ QC bypass (`mod_proteomics.R:406-416`) |
| `bus` (the initial `workflow_data` fields) | list | yes | all | `mod_proteomics.R:205-245`; `mod_lipidomics.R:139-176` (`data_type = "lipid"` `:143`, `lipid_assay_data` `:155`); `mod_metabolomics.R:139-176` (`metabolite`, `metabolite_assay_data`); `mod_transcriptomics.R:231-237` |
| `bus_registration` | chr(1): `omics` / `userData` | yes | all | `registerWorkflowOmicsBus` (`func_general_parameter_helpers.R:140`), called at `mod_lipidomics.R:180`, `mod_metabolomics.R:180`, `mod_transcriptomics.R:246`; prot `session$userData$workflow_data` (`mod_proteomics.R:248`) |
| `extras` | ids | no | prot | the peptide-checkpoint panel (`mod_proteomics.R:101-113`), its save and open (`:275-304`), the native-restore QC trigger (`:252`) |

### 4.3 `binding`: core side; the GUI never computes (owner KIT-102, extended by each family)

| Field | Type | Req. | Omics | Where today's value lives |
|---|---|---|---|---|
| `binding` | ref (declaration binding, or the DIA adaptor) | yes | all | lipid `utils_workflow_lipid_declarations.R:17`; metab `utils_workflow_metab_declarations.R:17`; trx `utils_workflow_trx_declarations.R:702`; the registry `omicsDeclarationBindings` (`utils_workflow_omics_declarations.R:588-591`); prot: the DIA adaptor (new, KIT-102) |
| `runner` | chr(1): `omics_actions` / `dia_managed_step` | yes | all | `installWorkflowOmicsActions` (`utils_workflow_omics_actions.R:13`; installed at `mod_lipidomics.R:196`, `mod_metabolomics.R:196`, `mod_transcriptomics.R:255`); prot: per-step `installWorkflowManagedStep` (`utils_workflow_execution.R:633`, 18 sites) |
| `owner` | ref | yes | all | `createWorkflowSessionOwner` (`utils_workflow_dia_sessions.R:875`; `mod_lipidomics.R:194`, `mod_transcriptomics.R:254`); prot the alias `createDIASessionOwner` (`:902`, `mod_proteomics.R:266`) |
| `read_state` | ref | yes | all | `readWorkflowOmicsState` (`utils_workflow_omics_declarations.R:673`); trx `readTrxWorkflowState` (`utils_workflow_trx_declarations.R:749`, profile `mod_transcriptomics.R:49`); prot `readDIAStateMetadata` (`utils_workflow_dia_state.R:237`) |
| `read_metadata` | ref | yes | all | `readWorkflowOmicsMetadata` (`utils_workflow_omics_declarations.R:687`, profile `mod_transcriptomics.R:48`); prot `readDIAStateMetadata` |
| `class`, `prefix`, `feature_noun`, `long_table`, `step_labels`, `default_variants`, `qc_dir_role`, `filtering_metrics`, `da_fixed` | read from the binding | derived | lipid, metab, trx | lipid binding `:18-37`, `:101-105`; metab `:18-37`, `:92-96`; trx `:703-722`, `:732-733` |

### 4.4 `import` (owner KIT-102)

| Field | Type | Req. | Omics | Where today's value lives |
|---|---|---|---|---|
| `source_panel` | id: `feature_table` / `proteomics` | yes | all | lipid panel `mod_lipid_import_ui_helpers.R`; metab `mod_metab_import_ui.R:26-329`; trx through lipid's (`mod_lipid_import_ui.R:28-35`); prot `mod_prot_import_ui.R` |
| `title` | chr(1) | yes | all | `mod_lipid_import_ui_helpers.R:359`; `mod_metab_import_ui.R:26`; `mod_transcriptomics.R:35`; `mod_prot_import_ui.R:22` |
| `formats`, `default_format` | named chr, chr(1) | yes | all | prot `mod_prot_import_ui.R:55-61` (`auto` the default); lipid `mod_lipid_import_ui_helpers.R:77-84` (`msdial`, `:93`); metab `mod_metab_import_ui.R:36-43` (`msdial`); trx `mod_transcriptomics.R:35` (`trx_counts`) |
| `route_choice` | named chr | no | prot | DIA / DIA_LIMPA (`mod_prot_import_ui.R:64-68`) |
| `detect` | ref | yes | all | `detectProteomicsFormat` (`func_prot_import_detection.R:14`); `detectLipidomicsFormat` (`func_lipid_import_detection.R:17`); `detectMetabolomicsFormat` (`func_metab_import_detection.R:17`); trx: a constant (`trx_counts`) |
| `help` | chr(1) | yes | all | `mod_lipid_import_ui_helpers.R:100`; `mod_metab_import_ui.R:50`; `mod_transcriptomics.R:37-39`; prot step headings (`mod_prot_import_ui.R:27,71,106`) |
| `assay_name`, `id_label`, `annotation_label` | chr(1) | yes (feature table) | lipid, metab, trx | lipid `mod_lipid_import_ui_helpers.R:105,168,210`; metab `mod_metab_import_ui.R:62,176,215`; trx `mod_transcriptomics.R:36` |
| `second_assay` | lgl(1) | yes (feature table) | lipid, metab, trx | lipid and metab TRUE; trx FALSE (one assay under a profile, `mod_lipid_import_ui_helpers.R:76`) |
| `id_mapping_key` | chr(1) | yes (feature table) | lipid, metab, trx | `lipid_id_col` (`buildLipidImportColumnMapping`, `mod_lipid_import_processing_helpers.R:470`; trx keeps it); `metabolite_id_col` (`mod_metab_import_processing_helpers.R:69`) |
| `aux_inputs` | ids | no | prot, lipid, metab | prot FASTA (`mod_prot_import_ui.R:71`) and organism (`:106`); lipid and metab internal-standard pattern and sanitising; lipid only `exclude_norm` |
| `validate_mapping` | ref | no | lipid, metab | `validateLipidColumnMapping` (`func_lipid_import_detection.R:317`); `validateMetabColumnMapping` (`func_metab_import_detection.R:282`) |
| `check` | ref | no | trx | the count policy: `checkTrxImportCounts` (`mod_transcriptomics.R:99`, over `validateTrxCounts` in `func_trx_counts.R`). After KIT-102 it runs after the text-numeric conversion. |
| `preview`, `process` | ref | yes | all | today in `mod_`; they move in KIT-102. Lipid `loadLipidImportAssayPreview` (`mod_lipid_import_selection_helpers.R:1`) and `runLipidImportProcessing` (`mod_lipid_import_processing_helpers.R:618`); metab `prepareMetabImportAssaySelectionState` (`mod_metab_import_processing_helpers.R:321`) and `runMetabImportProcessing` (`:519`); trx `trxImportPreview` (`mod_transcriptomics.R:82`) and `runTrxImportProcessing` (`:130`); prot `importProtImportDataByFormat` (`mod_prot_import_processing_helpers.R:38`) inside the observer (`mod_prot_import_server_helpers.R:181`) |
| `summary_file`, `count_key` | chr(1) | yes (feature table) | lipid, metab, trx | `lipidomics_import_summary.tsv` and `n_lipids` (`mod_lipid_import_processing_helpers.R:561`, `:296`); `metabolomics_import_summary.tsv` and `n_metabolites` (`mod_metab_import_processing_helpers.R:270`, `:79`); `transcriptomics_import_summary.tsv` (`mod_transcriptomics.R:144`) with `n_lipids` |
| `success_text` | chr(1) | yes | all | "Data imported successfully!" (`mod_lipid_import_processing_helpers.R:721`, `mod_metab_import_processing_helpers.R:435`); prot "Data import successful!" (`mod_prot_import_processing_helpers.R:401`) |

**Not fields (parity and unification, KIT-102):**
- text-numeric sample columns are converted for every omic. Today only metab converts (`coerceMetabImportWorkflowPayloadSamples`, `mod_metab_import_processing_helpers.R:197`). Lipid refuses them (`mod_lipid_import_processing_helpers.R:528-538`), and trx refuses them as "counts must be numeric" (`func_trx_counts.R:28`). For proteomics, KIT-102 first establishes whether any reader leaves text-numeric quantity columns; the prot lanes record today's behaviour.
- import opens the Design step for every omic. Today metab leaves `design_matrix` `disabled` (`mod_metab_import_processing_helpers.R:300-302`), where lipid sets it `pending` (`mod_lipid_import_processing_helpers.R:307-312`), and so do prot (`completeProtImportSuccessState`, `mod_prot_import_processing_helpers.R:375`) and trx (through lipid).

### 4.5 `design` (owner KIT-103)

| Field | Type | Req. | Omics | Where today's value lives |
|---|---|---|---|---|
| `module_id` | chr(1) | yes | all | `design` (`mod_lipidomics.R:75`, `mod_transcriptomics.R:221`); prot `design_matrix` (`mod_proteomics.R:52`) |
| `help` | chr(1) | yes (prot: none) | lipid, metab, trx | `mod_lipid_design_ui.R:39`; `mod_metab_design.R:59`; `mod_transcriptomics.R:41-43` |
| `formula_default` | chr(1) | yes | all | `"~ 0 + group"`: `lipidDesignResolveFormulaString` (`mod_lipid_design_validation.R:5`); metab inline (`mod_metab_design.R:224-227`) |
| `run_order` | ref, **required, no shared default** | yes | all | lipid `gtools::mixedsort` (`mod_lipid_design_builder_state_helpers.R:298,395`; default at `mod_lipid_design_builder_server.R:43`); metab `gtools::mixedsort` (`mod_metab_design_builder_server.R:237`); trx `identity` (`mod_transcriptomics.R:44`); prot `gtools::mixedsort` (`mod_prot_design_builder_state_helpers.R:348`) |
| `run_order_rule` | chr(1) | yes | all | the recorded name: `mixedsort` (lipid, metab, prot) or `input` (trx) |
| `group_levels` | ref or NULL | yes | all | trx `trxDesignGroupLevels` (`mod_transcriptomics.R:155`; read at `mod_lipid_design_builder_server.R:181`); NULL for the others |
| `validate` | ref or NULL | yes | all | today in `mod_`; moves in KIT-103. Lipid `validateLipidDesignDaPreflight` (`mod_lipid_design_validation.R:449`); metab `validateMetabDesignDaPreflight` (`mod_metab_design_validation.R:429`); trx `trxDesignCheck` (`mod_transcriptomics.R:175`); prot NULL (no save-time preflight today) |
| `build_object` | ref | yes | all | lipid, inside `runLipidDesignBuilderObserverShell` (`mod_lipid_design_builder_helpers.R:35`); metab `createMetaboliteAssayData(...)` (`mod_metab_design.R:314`); trx `trxDesignObject` (`mod_transcriptomics.R:193`); prot `buildProtDesignStateCheckpoint` (`mod_prot_design_state_helpers.R:1`) |
| `record_provenance` | ref | yes | all | `recordLipidDesignProvenance` (`utils_workflow_lipid_declarations.R:211`), `recordMetabDesignProvenance` (`utils_workflow_metab_declarations.R:193`), `recordTrxDesignProvenance` (`utils_workflow_trx_declarations.R:753`), all over `recordWorkflowOmicsDesignProvenance` (`utils_workflow_omics_declarations.R:1067`); prot `setDIAStepProvenance` (`utils_workflow_dia_state.R:742`) |
| `state_name`, `state_description` | chr(1) | yes | all | `lipid_raw_data_s4` (`mod_lipid_design_builder_helpers.R:179`); `metab_raw_data_s4` (`mod_metab_design.R:342`); `trx_raw_data_s4` (`mod_transcriptomics.R:45`); prot `raw_data_s4` / `protein_s4_initial` (`mod_prot_design_state_helpers.R:19,90`) |
| `data_adapter` | id: `wide_assay_list` / `dia_long` | yes | all | prot's design edits act through `transformDIADesignTable` (`utils_workflow_dia_state.R:996`) |
| `annotation_step` | ref or NULL | no | prot | today in `mod_`; moves in KIT-103. The UniProt fetch (`mod_prot_design_state_helpers.R:169`) and FASTA processing (`mod_prot_design_import_helpers.R:127`) |
| `design_import` | list `{reader, id_candidates, qc_baseline}` or NULL | no | lipid, metab, prot | lipid `runLipidDesignImportConfirmationShell` (`mod_lipid_design_import_helpers.R:87`); metab `createMetabDesignImportedS4Object` (`mod_metab_design_import_helpers.R:352`) and its siblings; prot `mod_prot_design_import_helpers.R`; trx NULL (no button) |
| `commit` | id: `shared` / `prot_hook` | yes | all | prot `runProtDesignBuilderSaveFlow` (`mod_prot_design_builder_helpers.R:53`) and `persistProtDesignBuilderArtifacts` (`:1`) |

**Not a field (unification, KIT-103):** a design-save failure is shown as a notification for every omic. Today metab `stop()`s inside the observer (`mod_metab_design.R:237-241`), outside the `tryCatch` that would notify. Lipid notifies and returns (`mod_lipid_design_builder_helpers.R:82-88`).

### 4.6 `qc` (owner KIT-104)

| Field | Type | Req. | Omics | Where today's value lives |
|---|---|---|---|---|
| `title` | chr(1) | yes | all | `mod_lipid_qc.R:45`; `mod_metab_qc.R:44`; `mod_transcriptomics.R:47`; `mod_prot_qc.R:43` |
| `tabset_id` | chr(1) | yes | all | `lipid_qc_tabs` (`mod_lipid_qc.R:120`); `metab_qc_tabs` (`mod_metab_qc.R:116`); `trx_qc_tabs` (`mod_trx_qc_filter.R:171`); prot nested groups (`mod_prot_qc.R:123-141`) |
| `groups` | list of `{id, title, tabs = list of {id, title, component, routes, needs_paths}}` | yes | all | lipid `intensity`, `duplicates`, `itsd`, `s4_finalize` (`mod_lipid_qc.R:59-80`, `:121-124`); metab the same (`mod_metab_qc.R:54-75`, `:117-120`); trx `filter`, `s4_finalize` (`mod_trx_qc_filter.R:171-179`); prot the spec tables (section 2.3) |
| `auto_init` | lgl(1) | yes | all | TRUE for lipid (`initializeLipidQcSubmodules`, `mod_lipid_qc.R:57`), metab (`mod_metab_qc.R:52`) and trx (`initializeTrxQcSubmodules`, `mod_trx_qc_filter.R:177`); FALSE for prot (its `qc_trigger`, `mod_prot_qc.R:150-201`) |
| `complete_state` | chr(1) or NULL | no | prot | `protein_replicate_filtered` (`mod_prot_qc.R:204-216`); NULL for the others (their finalise tab completes QC) |
| `filtering_figure` | ref | yes (lipid, metab, trx) | lipid, metab, trx | `updateLipidFiltering`; `updateMetaboliteFiltering`; `trxQcFilteringFigure` (`mod_trx_qc_filter.R:55`) |
| `itsd` | list `{name_columns, metrics, quant}` or NULL | no | lipid, metab | lipid inline in `mod_lipid_qc_itsd.R`; metab `analyzeMetabQcItsdData` (`mod_metab_qc_itsd.R:57`) |
| `finalize_alignment` | ref | yes (lipid, metab, trx) | lipid, metab, trx | lipid `validateLipidDesignAlignment` (called at `mod_lipid_qc_s4.R:383`, defined at `mod_lipid_design_validation.R:87`); metab `validateMetabDesignAlignment` (`mod_metab_design_validation.R:67`), not called at finalise today; trx a new `trxDesignAlignment` over `validateTrxDesign` |
| `finalize_zero_rows` | lgl(1) | yes | lipid, metab, trx | lipid refuses any zero-row assay (`mod_lipid_qc_s4.R:365`); metab refuses only when there are no assays (`mod_metab_qc_s4_server_helpers.R:337-338`). Kept per omic (user, "other differences stay per omic"). |
| `finalize_opens_norm` | lgl(1) | yes | lipid, metab, trx | lipid's `completeLipidQcS4TabStatus` also opens normalisation (`mod_lipid_qc_s4.R:420`); metab's does not. Kept per omic. |
| `finalize_per_assay` | lgl(1) | yes | lipid, metab, trx | TRUE for lipid and metab; FALSE for trx (one table) |

**Not a field (parity, KIT-104):** the QC-finalise design-alignment check runs for every omic. Today only lipid runs it. Metab's finalise checks only for assays and design rows (`validateMetabQcS4FinalizeState`, `mod_metab_qc_s4_server_helpers.R:320-347`), and so does trx's (`mod_trx_qc_filter.R:216`).

### 4.7 `norm` (owner KIT-105 for the frame, the library, lipid and metab; KIT-110 for trx's panel set)

| Field | Type | Req. | Omics | Where today's value lives |
|---|---|---|---|---|
| `module_id` | chr(1) | yes | all | `norm` (`mod_lipidomics.R:95`); prot `normalization` (`mod_proteomics.R:68`) |
| `panels` | ordered ids from the shared panel library, each with its options | yes | all | lipid `buildLipidNormOptionsControlPanel` (`mod_lipid_norm_ui_helpers.R:1`): `apply_itsd` `:31`, `itsd_aggregation` `:42`, `log_offset` `:59`, `norm_method` `:73` (cyclicloess, quantile, scale, none; default cyclicloess), `ruv_mode` `:93` (automatic, manual, skip), `ruv_percentage` `:186` (5), `ruv_k` `:197` (2), run `:228`, reset `:241`, `export_session` `:253`; metab the same ids by hand (`mod_metab_norm_ui.R:53-275`); trx `trxNormOptionsPanel` (`mod_trx_norm.R:25-55`); prot `mod_prot_norm_ui.R:26,69,169,180,216,229,241,250` |
| `plot_tabs` | ids | yes | all | lipid `buildLipidNormQcTabsetPanel` (`mod_lipid_norm_ui_helpers.R:285-308`): itsd, ruv, pca, density, rle, correlation, correlation_filter, log; metab the same, by hand (`mod_metab_norm_ui.R:290-658`); trx `ruv, pca, rle, log` (`mod_transcriptomics.R:57`); prot RUV QC, PCA, Density, RLE, Correlation, Post-Normalisation QC, Correlation Filtering (`mod_prot_norm_ui.R:258-470`) |
| `stage_labels` | chr(3) | yes | all | default "Post-Filtering", "Post-Normalization", "RUV-Corrected" (`mod_lipid_norm_ui_helpers.R:287`); trx (`mod_transcriptomics.R:58`) |
| `assay_slots` | int(1) | yes | lipid, metab, trx | 2 (`mod_lipid_norm_ui_helpers.R:286`); trx 1 |
| `ruv_heading`, `ruv_help` | chr(1) | yes | all | defaults at `mod_lipid_norm_ui_helpers.R:296-298`; trx (`mod_transcriptomics.R:59-61`) |
| `run` | ref | yes | all | managed `runWorkflowOmicsManagedNormalisation` (`utils_workflow_omics_actions.R:222`); in process lipid `handleLipidNormRunNormalization` (`mod_lipid_norm_workflow_helpers.R:102`) and metab `runMetabNormNormalizationPipelineShell` (`mod_metab_norm_pipeline_helpers.R:52`), which move in KIT-105; trx `runTrxNormalisation` and `runTrxNormManaged` (`mod_trx_norm.R:156`, `:130`), which move in KIT-110; prot `runDIAReferenceNormalization` (`utils_workflow_dia_products.R:525`) and the DIA_LIMPA pass-through (`mod_prot_norm_workflow_helpers.R:365`) |
| `correlation` | ref or NULL | yes | all | `runWorkflowOmicsManagedCorrelation` (`utils_workflow_omics_actions.R:315`); prot `runDIAReferenceCorrelation` (`utils_workflow_dia_products.R:895`); trx NULL |
| `pre_qc_on_entry` | lgl(1) | yes | all | TRUE for lipid, metab and prot; FALSE for trx |
| `export` | list `{prefix, counts_field, samples}` or NULL | no | lipid, metab, prot | lipid `buildLipidNormExportSessionData` (`mod_lipid_norm_workflow_helpers.R:329`) and `saveLipidNormExportArtifacts` (`:402`); metab `buildMetabNormExportSessionData` (`mod_metab_norm_session_helpers.R:499`); managed `startWorkflowManagedSessionExport` (`utils_workflow_saved_exports.R:78`); prot `export_filtered_session` and `export_s4_compatibility` (`mod_prot_norm_ui.R:241,250`); trx NULL |

Composite figure and final PCA come from the binding (`composite_figure`, `final_pca`: lipid `:103,108`, metab `:94,99`), so they are not fields. KIT-105 moves the module functions those bindings call (section 5).

### 4.8 `da` (owner KIT-106)

| Field | Type | Req. | Omics | Where today's value lives |
|---|---|---|---|---|
| `module_id`, `status_key` | chr(1) | yes | all | `de` and `differential_analysis` (`mod_lipidomics.R:105,167`; `mod_transcriptomics.R:224,235`); prot `differential_expression` (`mod_proteomics.R:76,235`). Prot's publication also sets `differential_abundance` and opens `enrichment_analysis` (`utils_workflow_dia_da.R:451-452`). |
| `model` | list `{editable, text}` | yes | all | editable formula for lipid (`mod_lipid_da_ui.R:24-30`), metab (`mod_metab_da.R:53`) and prot (`mod_prot_da_ui.R:47`); trx displays it only (`trxDaModelUi`, `mod_trx_da.R:24`; `trxDaModelText` `:33`) |
| `model_record` | from the binding | derived | trx | `da_fixed` (`utils_workflow_trx_declarations.R:722`: edgeR_QL, robust, BH, `treat = FALSE`); `trxDaModel` (`:375`) |
| `threshold` | list `{label, default}` | yes | all | "Q-value Threshold:" 0.05 (`mod_lipid_da_ui.R:36-42`, `mod_metab_da.R:64-65`, `mod_prot_da_ui.R:58-60`); trx "FDR (BH) threshold:" (`mod_transcriptomics.R:62`) |
| `lfc` | list `{label, default, tested}` | yes | all | "Log Fold-Change Cutoff:" 0, a TREAT test (`mod_lipid_da_ui.R:44-50`; metab `:72-73`; prot `:66-67`); trx "Fold-change line (display only; no test):" (`mod_transcriptomics.R:63`), `tested = FALSE` |
| `columns` | list `{id, label, logfc, pvalue, significance, significance_label}` | yes | all | lipid `lipid_id`, `lipid_name`, `raw_pvalue`, `fdr_qvalue` (`mod_lipid_da_results_helpers.R:106-107`); metab `metabolite_id`, `metabolite_name`, `assay`, `logFC`, `raw_pvalue`, `fdr_qvalue`, `significant` (`mod_metab_da_display_helpers.R:136-139`); trx `GeneID`, `GeneName`, `PValue`, `FDR`, "FDR (BH)" (`mod_transcriptomics.R:64`); prot the object's `protein_id_column`, `log2FC`, `raw_pvalue`, `fdr_qvalue` (`mod_prot_da_render_table.R:28,62`) |
| `long_table` | from the binding; chr(1) for prot | derived | all | `da_lipids_long` (lipid binding `:20`); `da_metabolites_long` (metab `:20`); `da_genes_long` (trx `:705`); prot `da_proteins_long` (`utils_workflow_dia_da.R:143`) |
| `assays` | chr(1): `multi` / `single` | yes | all | lipid and metab multi (assay selectors); trx and prot single |
| `views` | ids `{static_volcano, interactive_volcano, heatmap, table}` with options | yes | all | interactive default TRUE for lipid (`mod_lipid_da_ui.R:147`), metab (`mod_metab_da.R:174`) and prot (`mod_prot_da_ui.R:158`), FALSE for trx (`mod_transcriptomics.R:65`); heatmap tree cut for lipid (`mod_lipid_da_ui.R:310-331`) and prot (`mod_prot_da_ui.R:288-330`); metab's UI lacks the tree-cut inputs its server reads (`mod_metab_da_registration_helpers.R:90-93,165-168`) |
| `heatmap` | ref or NULL | yes | all | `generateLipidDAHeatmap` (`func_lipid_da_plotting.R:369`); `generateMetabDAHeatmap` (`func_metab_da_heatmap.R:29`); prot `generateDIAReferenceHeatmap` route; trx NULL |
| `run` | from the binding or the adaptor | derived | all | managed `runWorkflowOmicsManagedDa` (`utils_workflow_omics_actions.R:338`); in process `runWorkflowOmicsReferenceDaAnalysis` (`utils_workflow_omics_declarations.R:1094`); prot `runDIAReferenceDaAnalysis` (`utils_workflow_dia_da.R:341`) |
| `preflight` | ref | yes | all | lipid `validateLipidDaAnalysisPreflight` (`mod_lipid_da_analysis_helpers.R:84`; moves in KIT-106); metab the binding's (`utils_workflow_metab_declarations.R:76`, a no-op) with the design preflight; trx `trxDaPreflight` (binding `:724`); prot `validateProtDaModelAndContrasts` (`mod_prot_da_validation.R:108`; moves in KIT-106) |
| `contrast_entries` | ref or NULL | no | trx, prot | trx `trxDaContrastEntries` (`mod_trx_da.R:58`); prot `normaliseProtDaContrastsTable` (`mod_prot_da_validation.R:5`; moves in KIT-106) |
| `session_load` | id or NULL | no | lipid, metab, prot | lipid `mod_lipid_da_session_helpers.R`; metab `mod_metab_da_observer_helpers.R:446-639,760-930`; prot `mod_prot_da_handlers_load_session.R`; trx NULL |

**Not fields (KIT-106):**
- unification: the interactive volcano is off by default and is built on demand for every omic, as the worker's `<prefix>_volcano_page` (`computeWorkflowOmicsVolcanoPageAction`, `utils_workflow_omics_declarations.R:1309`; `workflowOmicsVolcanoPageFactory`, `utils_workflow_omics_actions.R:143`);
- parity: the metab tree-cut defect is fixed (findings F-03).

### 4.9 `enrichment` (owner KIT-107; NULL for every omic but proteomics)

| Field | Type | Req. | Omics | Where today's value lives |
|---|---|---|---|---|
| `tab` | list `{status_key, value, module_id, testid}` | yes when set | prot | `enrichment_analysis`, `enrichment`, `prot-tab-enrichment` (`mod_proteomics.R:84,169,172`) |
| `operation` | ref | yes | prot | `diaEnrichmentOperation` (`utils_workflow_dia_da.R:539-559`: `ProteinQuantitativeData` in, `EnrichmentResults` out, product `dia-enrichment-results`, not cacheable); run `runDIAReferenceEnrichment` (`:572`) |
| `parameters` | chr | yes | prot | `taxon_id`, `up_cutoff`, `down_cutoff`, `q_cutoff`, `exclude_iea`, `correction_method`, `protein_id_column`, `contrast_names`, `results` (`utils_workflow_dia_da.R:543-544`) |
| `background` | chr(1) | yes | prot | `per_contrast_da_table` (`utils_workflow_dia_da.R:545-546`) |
| `id_column`, `gprofiler_id_column` | chr(1) | yes | prot | `uniprot_acc` by default; `gene_name` when the DA table has it (`mod_prot_enrich_analysis_helpers.R:1-21`) |
| `cutoffs` | list `{up, down, q}` | yes | prot | 0, 0, 0.05 (`mod_prot_enrich.R:54-56,64-66,74-76`) |
| `correction` | list `{choices, default}` | yes | prot | gSCS, fdr; default gSCS (`mod_prot_enrich.R:145-151`) |
| `exclude_iea` | lgl(1) | yes | prot | FALSE (`mod_prot_enrich_analysis_helpers.R:41`; not an input) |
| `taxon` | list `{default, prefill_from}` | yes | prot | 9606 (`mod_prot_enrich.R:84-86`), prefilled from `workflow_data$taxon_id` |
| `services` | ids | yes | prot | gprofiler2, clusterprofiler, stringdb (a placeholder) tabs (`mod_prot_enrich.R:162-315`); supported organisms `buildProtEnrichSupportedOrganisms` (`mod_prot_enrich_observer_helpers.R:164`) |
| `annotation` | ref | yes | prot | today in `mod_`; moves in KIT-107. `resolveProtEnrichUniprotAnnotations` (`mod_prot_enrich_input_helpers.R:599`) |
| `organism_filter` | ref or NULL | no | prot | today in `mod_`; moves in KIT-107. `resolveProtEnrichOrganismMapping` (`:1-181`), `applyProtEnrichOrganismFilter` (`:182`); toggle `enable_organism_filter` (`mod_prot_enrich.R:95`) |
| `test_double` | ref | no | prot | `runProtEnrichDeterministicProcessEnrichments` (`mod_prot_enrich_execution_helpers.R:126`) |
| `state_name` | chr(1) | yes | prot | `enrichment_completed` (`saveProtEnrichCompletedState`, `mod_prot_enrich_state_helpers.R:210`) |

### 4.10 `summary`: summary, report and export (owner KIT-108)

| Field | Type | Req. | Omics | Where today's value lives |
|---|---|---|---|---|
| `module_id` | chr(1) | yes | all | `summary` (`mod_lipidomics.R:115`); prot `session_summary` (`mod_proteomics.R:92`) |
| `placeholder`, `description` | chr(1) | yes | all | `mod_lipid_summary.R:53,55`; `mod_metab_summary.R:49,51`; `mod_prot_summary_ui.R:25,27`; `mod_transcriptomics.R:68-69` |
| `report` | list `{dir, template, status}`; `template` is a string or a per-route ref | yes | all | `inst/reports/<dir>/`; lipid `mod_lipid_summary_workflow_helpers.R:359,387`; metab `mod_metab_summary_observer_helpers.R:750`; trx `mod_transcriptomics.R:50`; prot `resolveProtSummaryExpectedTemplate` (`mod_prot_summary_report_helpers.R:63-85`) |
| `final_state` | ref | yes | all | lipid and metab, head first: `workflowSummaryFinalStateName` (`utils_workflow_dia_sessions.R:1442`) over `mod_lipid_summary_state_helpers.R:322-329` and `mod_metab_summary_observer_helpers.R:315-319`; prot, the first of five states in the history (`resolveProtSummaryFinalS4State`, `mod_prot_summary_report_helpers.R:1-25`); trx `trxSummaryContext` (`mod_trx_summary.R:121`) |
| `record_when` | ref | yes | all | prot only for a `DIAStateAdapter` (`mod_prot_summary_session_helpers.R:59`); the others whenever there is a `workflow_data` |
| `copy` | list `{gate, required}` | yes | all | prot gate `assertDIAExportReady` (`mod_prot_summary_publication_helpers.R:104`); lipid and metab require "Contrasts Table", "Design Matrix", "Study Parameters" (`mod_lipid_summary_state_helpers.R:213`); prot does not classify copy failures |
| `integration` | list `{filename, operation, write, requires, missing_message}` | yes | all | `%s_%s_final_s4.RDS` (`mod_lipid_summary_workflow_helpers.R:89`, `mod_metab_summary_observer_helpers.R:421`); prot native `.msdia.tar` (`mod_prot_summary_session_helpers.R:84-85`); operation `<prefix>_integration_export` (`mod_lipid_summary_workflow_helpers.R:127`, `mod_metab_summary_observer_helpers.R:462`); trx `transcriptomics_normalized_data.tsv`, `trxSummaryIntegrationWrite`/`Ready` (`mod_trx_summary.R:136,143`), the missing message (`mod_transcriptomics.R:72`) |
| `export_steps` | list of `{id, panel, operation, choices, accept, run_reference, publications, status, figure}` | no | trx | the RUVIII export (`mod_trx_summary.R:25-118`): `trxRuviiiMChoices` `:25`, `trxHeadHasTechnicalReplicates` `:33`, panel `:43`, `runTrxRuviiiExport` `:65`; core `acceptTrxRuviiiExport` (`utils_workflow_trx_declarations.R:411`), `runTrxReferenceRuviiiExport` (`:815`) |

`copyToResultsSummary`'s per-omic branches (`func_general_filemgmt_results_summary_helpers.R:325-651`) become an internal map keyed by `omic_type`, not a configuration field, because notebooks call the exported function.

### 4.11 Example values per omic, from today's code

The blocks below show the values today's code holds, in the contract's shape. `# moves: KIT-1nn` marks a reference whose function lives in a module file today. Library and plugin ids are new names for today's components; KIT-102 to 108 fix them.

**Proteomics** (routes DIA, DIA_LIMPA, TMT, LFQ):

```r
list(
  shell = list(omic_type = "proteomics", prefix = "prot", testid = "prot", label = NULL, label_lower = "proteomics",
    tabset_id = "workflow_tabs", routes = c("DIA", "DIA_LIMPA", "TMT", "LFQ"),        # utils_workflow_dia_sessions.R:804-812
    tabs = list(                                                                         # mod_proteomics.R:122-182, 380-386
      list(status_key = "setup_import", value = "setup", module_id = "setup_import", title = "Setup & Import", icon = "upload", testid = "prot-tab-setup", step_name = "Import"),
      list(status_key = "design_matrix", value = "design", module_id = "design_matrix", title = "Design Matrix", icon = "table", testid = "prot-tab-design", step_name = "Design"),
      list(status_key = "quality_control", value = "qc", module_id = "quality_control", title = "Quality Control", icon = "chart-line", testid = "prot-tab-qc", step_name = "QC"),
      list(status_key = "normalization", value = "normalization", module_id = "normalization", title = "Normalization", icon = "balance-scale", testid = "prot-tab-normalization", step_name = "Normalize"),
      list(status_key = "differential_expression", value = "da", module_id = "differential_expression", title = "Differential Expression", icon = "chart-bar", testid = "prot-tab-da", step_name = "DE"),
      list(status_key = "enrichment_analysis", value = "enrichment", module_id = "enrichment_analysis", title = "Enrichment Analysis", icon = "network-wired", testid = "prot-tab-enrichment", step_name = "Enrich"),
      list(status_key = "session_summary", value = "session_summary", module_id = "session_summary", title = "Session Summary & Report", icon = "file-export", testid = "prot-tab-session-summary", step_name = "Summary")),
    chain_rules = list(qc_bypass = c("TMT", "LFQ")),                                     # mod_proteomics.R:406-416
    bus = c("data_tbl", "fasta_file_path", "config_list", "taxon_id", "organism_name", "design_matrix", "data_cln",
      "contrasts_tbl", "state_manager", "peptide_data", "protein_log2_quant", "protein_data",
      "ruv_normalised_for_da_analysis_obj", "da_analysis_results_list", "uniprot_dat_cln", "enrichment_results",
      "tab_status", "state_update_trigger", "processing_log"),                          # mod_proteomics.R:205-245
    bus_registration = "userData", extras = c("peptide_checkpoint")),                    # :248; :101-113, 275-304
  binding = list(binding = "<DIA adaptor, KIT-102>", runner = "dia_managed_step", owner = createDIASessionOwner,
    read_state = readDIAStateMetadata, read_metadata = readDIAStateMetadata),
  import = list(source_panel = "proteomics", title = "Setup & Import Data",
    formats = c("Auto-detect" = "auto", "DIA-NN" = "diann", "Spectronaut DIA" = "spectronaut", "FragPipe LFQ" = "fragpipe",
      "MaxQuant LFQ" = "maxquant", "Proteome Discoverer TMT" = "pd_tmt"), default_format = "auto",  # mod_prot_import_ui.R:55-61
    route_choice = c(DIA = "DIA", DIA_LIMPA = "DIA_LIMPA"),                              # :64-68 (DIA-NN only)
    detect = detectProteomicsFormat, aux_inputs = c("fasta", "organism", "mapping_files"),
    process = importProtImportDataByFormat,                                              # moves: KIT-102
    success_text = "Data import successful!"),
  design = list(module_id = "design_matrix", formula_default = "~ 0 + group", run_order = gtools::mixedsort,
    run_order_rule = "mixedsort", group_levels = NULL, validate = NULL, data_adapter = "dia_long",
    build_object = buildProtDesignStateCheckpoint,                                       # moves: KIT-103
    state_name = c(DIA = "raw_data_s4", DIA_LIMPA = "raw_data_s4", TMT = "protein_s4_initial", LFQ = "protein_s4_initial"),
    annotation_step = "uniprot_and_fasta",                                               # moves: KIT-103
    record_provenance = setDIAStepProvenance, commit = "prot_hook"),
  qc = list(title = "Quality Control & Filtering", auto_init = FALSE, complete_state = "protein_replicate_filtered",
    groups = list(
      peptide_qc = list(DIA = c("qvalue_filter", "protein_peptide_filter", "rollup", "intensity_filter", "sample_filter",
                                "replicate_filter", "imputation"),
                        DIA_LIMPA = c("qvalue_filter", "rollup", "log2_matrix", "accession_cleanup")),
      protein_qc = list(DIA = c("rollup", "cleanup", "intensity_filter", "duplicate_removal", "replicate_filter"),
                        DIA_LIMPA = c("dpc_quant", "normalisation", "peptide_count"),
                        TMT = c("cleanup", "intensity_filter", "duplicate_removal", "replicate_filter"),
                        LFQ = c("cleanup", "intensity_filter", "duplicate_removal", "replicate_filter")))),
  norm = list(module_id = "normalization",
    panels = list(between_sample = list(methods = c("cyclicloess", "quantile", "scale", "none"), default = "cyclicloess"),
      aesthetics = list(), ruv_iii = list(modes = c("automatic", "manual", "skip"), percentage = 5, k = NULL),
      actions = list(), export = list(buttons = c("export_filtered_session", "export_s4_compatibility"))),
    plot_tabs = c("ruv", "pca", "density", "rle", "correlation", "post_norm_qc", "correlation_filter"),
    run = c(DIA = "runDIAReferenceNormalization", TMT = "runDIAReferenceNormalization", LFQ = "runDIAReferenceNormalization",
      DIA_LIMPA = "runDIAReferenceLimpaNormalisationPassThrough"),                      # pass-through moves: KIT-105
    correlation = runDIAReferenceCorrelation, pre_qc_on_entry = TRUE),
  da = list(module_id = "differential_expression", status_key = "differential_expression",
    model = list(editable = TRUE), threshold = list(label = "Q-value Threshold:", default = 0.05),
    lfc = list(label = "Log Fold-Change Cutoff:", default = 0, tested = TRUE),
    columns = list(id = "<object@protein_id_column>", logfc = "log2FC", pvalue = "raw_pvalue", significance = "fdr_qvalue"),
    long_table = "da_proteins_long", assays = "single",
    views = list(static_volcano = "volcano_plot_static", interactive_volcano = list(default = TRUE),  # TRUE today
      heatmap = list(tree_cut = TRUE), table = "paged_from_store"),
    preflight = validateProtDaModelAndContrasts, contrast_entries = normaliseProtDaContrastsTable,  # moves: KIT-106
    session_load = "prot_native_archive"),
  enrichment = list(tab = list(status_key = "enrichment_analysis", value = "enrichment", module_id = "enrichment_analysis",
      testid = "prot-tab-enrichment"),
    operation = diaEnrichmentOperation, background = "per_contrast_da_table",
    id_column = "uniprot_acc", gprofiler_id_column = "gene_name",
    cutoffs = list(up = 0, down = 0, q = 0.05), correction = list(choices = c("gSCS", "fdr"), default = "gSCS"),
    exclude_iea = FALSE, taxon = list(default = 9606, prefill_from = "taxon_id"),
    services = c("gprofiler2", "clusterprofiler", "stringdb"),
    annotation = resolveProtEnrichUniprotAnnotations, organism_filter = resolveProtEnrichOrganismMapping,  # moves: KIT-107
    test_double = runProtEnrichDeterministicProcessEnrichments, state_name = "enrichment_completed"),
  summary = list(module_id = "session_summary", placeholder = "e.g., my_proteomics_analysis",
    description = "Full protein analysis workflow with config parameters",
    report = list(dir = "proteomics", template = c(DIA = "DIANN_report.rmd", DIA_LIMPA = "DIANN_limpa_report.rmd",
      TMT = "TMT_report.rmd", LFQ = "LFQ_report.rmd")),                                 # mod_prot_summary_report_helpers.R:63-85
    final_state = c("ruv_covariates", "correlation_filtered", "ruv_corrected", "protein_replicate_filtered", "imputed"),  # first in history
    record_when = "DIAStateAdapter", copy = list(gate = assertDIAExportReady, required = character()),
    integration = list(filename = c(native = "%s_%s_final.msdia.tar", other = "%s_%s_final_s4.RDS"), operation = NULL),
    export_steps = list())
)
```

**Lipidomics:**

```r
list(
  shell = list(omic_type = "lipidomics", prefix = "lipid", testid = "lipid", label = "Lipidomics", label_lower = "lipidomics",
    tabset_id = "lipidomics_tabs", routes = "LIPIDOMICS",
    tabs = list(                                                                         # mod_lipidomics.R:59-116, 271-276
      list(status_key = "setup_import", value = "import", module_id = "import", title = "Setup & Import", icon = "file-import", testid = "lipid-tab-import", step_name = "Import"),
      list(status_key = "design_matrix", value = "design", module_id = "design", title = "Design Matrix", icon = "th", testid = "lipid-tab-design", step_name = "Design"),
      list(status_key = "quality_control", value = "qc", module_id = "qc", title = "Quality Control", icon = "check-double", testid = "lipid-tab-qc", step_name = "QC"),
      list(status_key = "normalization", value = "norm", module_id = "norm", title = "Normalization", icon = "balance-scale", testid = "lipid-tab-norm", step_name = "Normalize"),
      list(status_key = "differential_analysis", value = "de", module_id = "de", title = "Differential Analysis", icon = "chart-line", testid = "lipid-tab-de", step_name = "DE"),
      list(status_key = "session_summary", value = "summary", module_id = "summary", title = "Summary & Export", icon = "download", testid = "lipid-tab-summary", step_name = "Summary")),
    bus = list(data_type = "lipid", assay_slot = "lipid_assay_data"), bus_registration = "omics"),   # :139-176
  binding = list(binding = lipidomicsDeclarationBinding, runner = "omics_actions", owner = createWorkflowSessionOwner,
    read_state = readWorkflowOmicsState, read_metadata = readWorkflowOmicsMetadata),
  import = list(source_panel = "feature_table", title = "Lipidomics Data Import",
    formats = c("MS-DIAL" = "msdial", "Progenesis QI" = "progenesis", "XCMS" = "xcms",
      "Compound Discoverer" = "compound_discoverer", "LipidSearch" = "lipidsearch", "Other/Custom" = "custom"),
    default_format = "msdial", detect = detectLipidomicsFormat,
    help = "Import one or more assay files (e.g., LCMS_Pos, LCMS_Neg, GCMS)", assay_name = "LCMS_Pos",
    id_label = "Lipid ID Column", annotation_label = "Annotation Column", second_assay = TRUE,
    id_mapping_key = "lipid_id_col", aux_inputs = c("is_pattern", "sanitize_names", "exclude_norm"),
    validate_mapping = validateLipidColumnMapping,
    preview = loadLipidImportAssayPreview, process = runLipidImportProcessing,           # moves: KIT-102
    summary_file = "lipidomics_import_summary.tsv", count_key = "n_lipids", success_text = "Data imported successfully!"),
  design = list(module_id = "design", help = "Use the tools below to define your experimental groups and contrasts for differential analysis. All assays (Pos/Neg modes) share the same design matrix.",
    formula_default = "~ 0 + group", run_order = gtools::mixedsort, run_order_rule = "mixedsort", group_levels = NULL,
    validate = validateLipidDesignDaPreflight,                                           # moves: KIT-103
    record_provenance = recordLipidDesignProvenance, state_name = "lipid_raw_data_s4", data_adapter = "wide_assay_list",
    design_import = "lipid", commit = "shared"),
  qc = list(title = "Lipid Quality Control & Filtering", tabset_id = "lipid_qc_tabs", auto_init = TRUE,
    groups = list(qc = c("intensity", "duplicates", "itsd", "s4_finalize")), filtering_figure = updateLipidFiltering,
    itsd = list(name_columns = c("Lipid name", "Name", "lipid")),
    finalize_alignment = validateLipidDesignAlignment,                                   # moves: KIT-103/104
    finalize_zero_rows = TRUE, finalize_opens_norm = TRUE, finalize_per_assay = TRUE),
  norm = list(module_id = "norm",
    panels = list(aesthetics = list(), itsd = list(aggregation = TRUE), log2 = list(offset = TRUE),
      between_sample = list(methods = c("cyclicloess", "quantile", "scale", "none"), default = "cyclicloess"),
      ruv_iii = list(modes = c("automatic", "manual", "skip"), percentage = 5, k = 2, noun = "Lipids"),
      actions = list(), export = list(button = "export_session")),
    plot_tabs = c("itsd", "ruv", "pca", "density", "rle", "correlation", "correlation_filter", "log"),
    stage_labels = c("Post-Filtering", "Post-Normalization", "RUV-Corrected"), assay_slots = 2L,
    run = runWorkflowOmicsManagedNormalisation, correlation = runWorkflowOmicsManagedCorrelation, pre_qc_on_entry = TRUE,
    export = list(prefix = "lipid_", counts_field = "lipid_counts", samples = "design")),
  da = list(module_id = "de", status_key = "differential_analysis", model = list(editable = TRUE),
    threshold = list(label = "Q-value Threshold:", default = 0.05),
    lfc = list(label = "Log Fold-Change Cutoff:", default = 0, tested = TRUE),
    columns = list(id = "lipid_id", label = "lipid_name", logfc = "logFC", pvalue = "raw_pvalue", significance = "fdr_qvalue"),
    assays = "multi", views = list(static_volcano = "volcano_plot", interactive_volcano = list(default = TRUE),  # TRUE today
      heatmap = list(tree_cut = TRUE), table = "rows"),
    heatmap = generateLipidDAHeatmap, preflight = validateLipidDaAnalysisPreflight,     # moves: KIT-106
    session_load = "omics_rds"),
  enrichment = NULL,
  summary = list(module_id = "summary", placeholder = "e.g., my_lipidomics_analysis",
    description = "Full lipidomics analysis workflow with normalization and DA",
    report = list(dir = "lipidomics", template = "lipidomics_report.rmd"),
    final_state = c("lipid_correlation_filtered", "lipid_norm_complete", "lipid_ruv_corrected", "lipid_normalized",
      "loaded_for_de", "lipid_qc_complete"),                                             # head first
    record_when = "workflow_data", copy = list(gate = NULL, required = c("Contrasts Table", "Design Matrix", "Study Parameters")),
    integration = list(filename = "%s_%s_final_s4.RDS", operation = "lipid_integration_export"), export_steps = list())
)
```

**Metabolomics.** The metab shell is identical to lipid's apart from `data_type` and the assay slot (`mod_metabolomics.R:143,155`). The differences from the lipid block are:

```r
list(
  shell = list(omic_type = "metabolomics", prefix = "metab", testid = "metab", label = "Metabolomics",
    label_lower = "metabolomics", tabset_id = "metabolomics_tabs", routes = "METABOLOMICS",
    tabs = "<as lipid, test ids metab-tab-*>", bus = list(data_type = "metabolite", assay_slot = "metabolite_assay_data")),
  binding = list(binding = metabolomicsDeclarationBinding, runner = "omics_actions"),
  import = list(title = "Metabolomics Data Import",
    formats = c("MS-DIAL" = "msdial", "Progenesis QI" = "progenesis", "XCMS" = "xcms",
      "Compound Discoverer" = "compound_discoverer", "Other/Custom" = "custom"),       # mod_metab_import_ui.R:36-42
    detect = detectMetabolomicsFormat, id_label = "Metabolite ID Column", id_mapping_key = "metabolite_id_col",
    second_assay = TRUE,                                                                 # second reader always MS-DIAL
    aux_inputs = c("is_pattern", "sanitize_names"), validate_mapping = validateMetabColumnMapping,
    preview = prepareMetabImportAssaySelectionState, process = runMetabImportProcessing, # moves: KIT-102
    summary_file = "metabolomics_import_summary.tsv", count_key = "n_metabolites"),
  design = list(validate = validateMetabDesignDaPreflight, record_provenance = recordMetabDesignProvenance,
    state_name = "metab_raw_data_s4", design_import = "metab"),
  qc = list(title = "Metabolite Quality Control & Filtering", tabset_id = "metab_qc_tabs",
    filtering_figure = updateMetaboliteFiltering, itsd = list(name_columns = c("Metabolite name", "Name", "metabolite")),
    finalize_alignment = validateMetabDesignAlignment,                                   # not called at finalise today
    finalize_zero_rows = FALSE, finalize_opens_norm = FALSE),
  norm = list(panels = list(ruv_iii = list(noun = "Metabolites")),
    export = list(prefix = "metab_", counts_field = "metabolite_counts", samples = "columns")),
  da = list(columns = list(id = "metabolite_id", label = "metabolite_name", logfc = "logFC", pvalue = "raw_pvalue",
      significance = "fdr_qvalue"),
    views = list(heatmap = list(tree_cut = FALSE)),                                      # the defect KIT-106 fixes
    heatmap = generateMetabDAHeatmap, preflight = "binding no-op plus design preflight"),
  summary = list(placeholder = "e.g., my_metabolomics_analysis",
    description = "Full metabolomics analysis workflow with normalization and DA",
    report = list(dir = "metabolomics", template = "metabolomics_report.rmd"),
    final_state = c("metab_correlation_filtered", "metab_norm_complete", "metab_ruv_corrected", "metab_normalized",
      "loaded_for_de", "metab_qc_complete"),
    integration = list(filename = "%s_%s_final_s4.RDS", operation = "metab_integration_export"))
)
```

**Transcriptomics** (today through `trxOmicProfile()` and `trxNormDaSummaryProfile()`, `mod_transcriptomics.R:32-73`):

```r
list(
  shell = list(omic_type = "transcriptomics", prefix = "trx", testid = "trx", label = "Transcriptomics",
    label_lower = "transcriptomics", tabset_id = "transcriptomics_tabs", routes = "TRANSCRIPTOMICS",
    tabs = "<as lipid; DA title 'Differential Expression'; test ids trx-tab-*>",       # mod_transcriptomics.R:220-225
    bus = list(data_type = "transcriptomics"), bus_registration = "omics"),            # :231-237, 246
  binding = list(binding = transcriptomicsDeclarationBinding, runner = "omics_actions", read_state = readTrxWorkflowState),
  import = list(source_panel = "feature_table", title = "Transcriptomics Data Import",
    formats = c("Count table (gene ID, gene name, sample counts)" = "trx_counts"), assay_name = "counts",
    id_label = "Gene ID Column", annotation_label = "Gene Name Column", second_assay = FALSE,
    id_mapping_key = "lipid_id_col",                                                     # kept: trxDesignObject reads it
    check = checkTrxImportCounts,                                                        # after the conversion (KIT-102)
    preview = trxImportPreview, process = runTrxImportProcessing,                        # moves: KIT-102
    summary_file = "transcriptomics_import_summary.tsv", count_key = "n_lipids"),
  design = list(run_order = identity, run_order_rule = "input", group_levels = trxDesignGroupLevels,
    validate = trxDesignCheck, build_object = trxDesignObject,                           # moves: KIT-103
    record_provenance = recordTrxDesignProvenance, state_name = "trx_raw_data_s4", design_import = NULL),
  qc = list(title = "Transcriptomics Quality Control & Filtering", tabset_id = "trx_qc_tabs", auto_init = TRUE,
    groups = list(qc = c("filter", "s4_finalize")), filtering_figure = trxQcFilteringFigure,  # moves: KIT-104
    finalize_alignment = "trxDesignAlignment (new, KIT-104)", finalize_per_assay = FALSE),
  norm = list(panels = list(tmm_logcpm = list(),
      ruv_grid = list(modes = c("automatic", "manual"), auto_percentage = c(min = 5, max = 30, step = 5),
        k_penalty_weight = 0.05, percentage = 10, k = 2), actions = list(reset = TRUE)),   # mod_trx_norm.R:25-55
    plot_tabs = c("ruv", "pca", "rle", "log"),
    stage_labels = c("Raw logCPM", "TMM logCPM", "RUVs-corrected logCPM (display only)"), assay_slots = 1L,
    ruv_heading = "RUVs Control Selection (canonical correlation)",
    run = "runTrxNormManaged / runTrxNormalisation",                                     # moves: KIT-110
    correlation = NULL, pre_qc_on_entry = FALSE, export = NULL),
  da = list(module_id = "de", status_key = "differential_analysis",
    model = list(editable = FALSE, text = "trxDaModelText"),                             # mod_trx_da.R:24-33
    threshold = list(label = "FDR (BH) threshold:", default = 0.05),
    lfc = list(label = "Fold-change line (display only; no test):", default = 0, tested = FALSE),
    columns = list(id = "GeneID", label = "GeneName", logfc = "logFC", pvalue = "PValue", significance = "FDR",
      significance_label = "FDR (BH)"),
    assays = "single", views = list(static_volcano = "trxDaStaticVolcano", interactive_volcano = list(default = FALSE)),
    heatmap = NULL, preflight = trxDaPreflight, contrast_entries = trxDaContrastEntries, session_load = NULL),
  enrichment = NULL,
  summary = list(placeholder = "e.g., my_transcriptomics_analysis",
    description = "Full transcriptomics analysis workflow with normalization, DE and the RUVIII integration export",
    report = list(dir = "transcriptomics", template = "transcriptomics_report.rmd"), final_state = trxSummaryContext,
    integration = list(filename = "transcriptomics_normalized_data.tsv", operation = "trx_integration_export",
      requires = "trx_ruviii_export",
      missing_message = "The integration export needs the RUVIII export: run it in the panel above, then save again."),
    export_steps = list(list(id = "ruviii_export", operation = "trx_ruviii_export",
      choices = "trxRuviiiMChoices(trxHeadHasTechnicalReplicates(manager))",          # mod_trx_summary.R:25-33
      accept = acceptTrxRuviiiExport, run_reference = runTrxReferenceRuviiiExport,
      publications = c("trx_qc_figures", "trx_composite_figure"), status = "trxRuviiiStatusText",
      figure = "ruviii_before_and_after.png")))
)
```

### 4.12 Where the 65 trx profile fields go

The trx profile (`trxOmicProfile()` has 34 fields and `trxNormDaSummaryProfile()` has 31; both counted with `load_all`) is the first draft of the contract. Each field maps to a contract field or leaves. KIT-109 checks that none remain.

| Profile fields | Goes to | Owner |
|---|---|---|
| `omic_type`, `prefix`, `testid`, `label`, `label_lower`, `workflow_type` | `shell` (`workflow_type` becomes `routes`) | KIT-102 |
| `class`, `binding`, `read_metadata`, `read_state` | `binding` (`class` is read from the binding) | KIT-102 |
| `data_format`, `import_title`, `formats`, `assay_name`, `id_label`, `annotation_label`, `import_help` | `import` | KIT-102 |
| `load_preview`, `process_import` | `import$preview`, `import$process` (core refs) | KIT-102 |
| `design_help`, `run_order`, `group_levels`, `validate_design`, `build_object`, `record_provenance`, `state_name`, `state_description` | `design` | KIT-103 |
| `qc_title` | `qc$title` | KIT-104 |
| `qc_tabs`, `qc_init`, `register_qc_summary`, `finalize` | leave: `qc$groups` ids, `qc$auto_init`, `qc$finalize_alignment` and `qc$finalize_per_assay` replace these server substitutions | KIT-104 |
| `norm_options`, `norm_tabs`, `norm_stage_labels`, `norm_ruv_heading`, `norm_ruv_help` | `norm$panels`, `norm$plot_tabs`, `norm$stage_labels`, `norm$ruv_heading`, `norm$ruv_help` | KIT-110 (trx values on KIT-105's frame) |
| `norm_server` | leaves: the frame plus trx's panel set | KIT-110 |
| `da_threshold_label`, `da_lfc_label`, `long_table`, `da_columns`, `da_static_volcano`, `da_volcano_interactive`, `da_contrasts`, `da_formula` | `da$threshold`, `da$lfc`, `da$long_table` (from the binding), `da$columns`, `da$views`, `da$contrast_entries`, `da$model` | KIT-106 |
| `da_server`, `da_model_ui`, `da_widget`, `da_with_source`, `da_export`, `da_run_reference` | leave: the shared module; the interactive volcano becomes the worker page; the run comes from the binding | KIT-106 |
| `summary_placeholder`, `summary_description`, `report_template`, `report_dir`, `summary_context`, `summary_state_object` | `summary$placeholder`, `$description`, `$report`, `$final_state` | KIT-108 |
| `integration_filename`, `integration_write`, `integration_ready`, `integration_missing_message`, `ruviii_m_choices` | `summary$integration`, `summary$export_steps` | KIT-108 |
| `summary_export_ui`, `summary_server` | leave: `export_steps` data | KIT-108 |

## 5. Core-to-module back-references at HEAD

**The count is 32, not 34:** 24 proteomics, 4 lipid, 4 metab, 0 trx. They were found with R parse data: every top-level function defined in `R/mod_*.R` that appears as a call or a symbol, outside comments and strings, in `R/utils_workflow_*.R` or `R/func_*.R`. No module function is referenced there by string (`get`, `exists`, `match.fun`). A comment-inclusive grep gives 34, which is the scoping's figure. The two extra names appear only in roxygen comments:
- `handleLipidNormRunNormalization` (defined at `mod_lipid_norm_workflow_helpers.R:102`), mentioned at `utils_workflow_omics_actions.R:208`. KIT-105 updates the comment when the function moves.
- `runProteinDuplicateRemovalStep` (`mod_prot_qc_protein_dedup.R:74`), mentioned at `utils_workflow_prot_declarations.R:271`. KIT-104 updates the comment.

**Rule for the owner:** before it deletes any module file, the owner moves the callee into the core, unchanged, and keeps its name (or re-points the caller to the moved function). A test shows the caller still works through the moved function.

| # | Callee | Defined at | Called from (caller function) | Owner |
|---:|---|---|---|---|
| 1 | `recordDIAImportMetadata` | `mod_prot_import_processing_helpers.R:114` | `utils_workflow_dia_sessions.R:1108` (`prepareDIAGUIReferences`) | KIT-102 |
| 2 | `getDIAImportMetadata` | `mod_prot_import_processing_helpers.R:145` | `utils_workflow_dia_state.R:1030` (`transformDIADesignTable`) | KIT-102 |
| 3 | `buildDIAReferenceDesignState` | `mod_prot_design_builder_server_helpers.R:205` | `utils_workflow_dia_sessions.R:1113` (`prepareDIAGUIReferences`) | KIT-103 |
| 4 | `mod_prot_design_builder_ui` | `mod_prot_design_builder_ui.R:22` | `func_general_applets.R:268` (`.runDesignMatrixApplet`) | KIT-103 |
| 5 | `mod_prot_design_builder_server` | `mod_prot_design_builder_server.R:16` | `func_general_applets.R:272` (`.runDesignMatrixApplet`) | KIT-103 |
| 6 | `.formatPeptideRollupSummary` | `mod_prot_qc_peptide_rollup.R:73` | `utils_workflow_dia_qc_rollup.R:917` (`buildDIAReferencePeptideHandoff`) | KIT-104 |
| 7 | `.peptideImputationAuditSummary` | `mod_prot_qc_peptide_impute.R:82` | `utils_workflow_dia_qc_rollup.R:869` (`buildDIAReferencePeptideHandoff`) | KIT-104 |
| 8 | `.proteinCleanupActiveKey` | `mod_prot_qc_protein_cleanup.R:99` | `utils_workflow_prot_declarations.R:160,172` (`computeDIAReferenceAccessionCleanup`) | KIT-104 |
| 9 | `resolveMetabDuplicateAssayData` | `mod_metab_qc_duplicates_server_helpers.R:1` | `utils_workflow_metab_declarations.R:54` (`metabolomicsDeclarationBinding$resolve_duplicates`) | KIT-104 |
| 10 | `buildProtNormSkippedRuvResult` | `mod_prot_norm_workflow_helpers.R:90` | `utils_workflow_dia_products.R:351` (`diaNormalizationOperation`) | KIT-105 |
| 11 | `resolveProtRuvParameters` | `mod_prot_norm_ruv_helpers.R:100` | `utils_workflow_dia_products.R:357` (`diaNormalizationOperation`) | KIT-105 |
| 12 | `validateProtNormRuvK` | `mod_prot_norm_ruv_helpers.R:26` | `func_prot_s4_norm_methods.R:205` (`ruvCancorFromInput`); `func_prot_norm_optimization_helpers.R:1114` (`computeRuvCovariates`) | KIT-105 |
| 13 | `stripProtNormQcPlotOwners` | `mod_prot_norm_qc_generation_helpers.R:32` | `utils_workflow_dia_products.R:361` (`diaNormalizationOperation`), `:770` (`diaRuvCovariatesOperation`), `:938` (`runDIAReferenceCorrelation`); `utils_workflow_dia_da.R:321` (`diaContrastOperation`); `utils_workflow_dia_sessions.R:1300` (`ingestDIANativeLegacySession`); **`utils_workflow_omics_declarations.R:295` (`omicsRuvOperation`, used by lipid and metab)** | KIT-105 |
| 14 | `buildProtNormQcPlotProduct` | `mod_prot_norm_qc_generation_helpers.R:52` | `utils_workflow_dia_products.R:936` (`runDIAReferenceCorrelation`, passed as a function) | KIT-105 |
| 15 | `buildProtNormDensityPlot` | `mod_prot_norm_qc_support_helpers.R:446` | `utils_workflow_dia_products.R:937` (`runDIAReferenceCorrelation`, passed as a function) | KIT-105 |
| 16 | `buildProtNormCorrelationPlot` | `mod_prot_norm_qc_support_helpers.R:468` | `utils_workflow_dia_products.R:937` (`runDIAReferenceCorrelation`, passed as a function) | KIT-105 |
| 17 | `resolveProtNormCorrelationResultFilenames` | `mod_prot_norm_correlation_helpers.R:1` | `utils_workflow_dia_products.R:940` (`runDIAReferenceCorrelation`) | KIT-105 |
| 18 | `saveProtNormCorrelationResults` | `mod_prot_norm_correlation_helpers.R:28` | `utils_workflow_dia_sessions.R:1466` (`restoreDIANativeLegacyArtifacts`) | KIT-105 |
| 19 | `saveProtNormExportMetadataFiles` | `mod_prot_norm_session_helpers.R:180` | `utils_workflow_dia_sessions.R:1469` (`restoreDIANativeLegacyArtifacts`), `:1605` (`refreshDIANativeSummaryArchive`) | KIT-105 |
| 20 | `buildProtNormExportSummaryContent` | `mod_prot_norm_session_helpers.R:125` | `utils_workflow_dia_sessions.R:1606` (`refreshDIANativeSummaryArchive`) | KIT-105 |
| 21 | `generateProtNormCompositeQcFigure` | `mod_prot_norm_workflow_helpers.R:283` | `utils_workflow_dia_products.R:595` (`runDIAReferenceNormalization`) | KIT-105 |
| 22 | `finalizeProtNormWorkflowState` | `mod_prot_norm_workflow_helpers.R:331` | `utils_workflow_dia_products.R:606` (`runDIAReferenceNormalization`) | KIT-105 |
| 23 | `writeLipidNormCompositeFigure` | `mod_lipid_norm_workflow_helpers.R:5` | `utils_workflow_lipid_declarations.R:103` (`lipidomicsDeclarationBinding$composite_figure`) | KIT-105 |
| 24 | `buildLipidNormCompositeFromFilesGenerator` | `mod_lipid_norm_support_helpers.R:425` | `utils_workflow_lipid_declarations.R:104` (`$composite_figure`) | KIT-105 |
| 25 | `buildLipidNormFinalQcPcaPlot` | `mod_lipid_norm_support_helpers.R:369` | `utils_workflow_lipid_declarations.R:108` (`$final_pca`) | KIT-105 |
| 26 | `runMetabNormCompositeQcFigureStep` | `mod_metab_norm_pipeline_helpers.R:632` | `utils_workflow_metab_declarations.R:94` (`metabolomicsDeclarationBinding$composite_figure`) | KIT-105 |
| 27 | `generateMetabNormCompositeFromFiles` | `mod_metab_norm_render_helpers.R:49` | `utils_workflow_metab_declarations.R:95` (`$composite_figure`, passed as a function) | KIT-105 |
| 28 | `buildMetabNormFinalQcPcaPlot` | `mod_metab_norm_session_helpers.R:836` | `utils_workflow_metab_declarations.R:99` (`$final_pca`) | KIT-105 |
| 29 | `normaliseProtDaContrastsTable` | `mod_prot_da_validation.R:5` | `utils_workflow_dia_da.R:350` (`runDIAReferenceDaAnalysis`) | KIT-106 |
| 30 | `validateProtDaModelAndContrasts` | `mod_prot_da_validation.R:108` | `utils_workflow_dia_da.R:355` (`runDIAReferenceDaAnalysis`) | KIT-106 |
| 31 | `validateProtDaFilteredSession` | `mod_prot_da_validation.R:72` | `utils_workflow_dia_sessions.R:1197` (`adoptDIANativeSession`), `:1262` (`ingestDIANativeLegacySession`) | KIT-106 |
| 32 | `validateLipidDaAnalysisPreflight` | `mod_lipid_da_analysis_helpers.R:84` | `utils_workflow_lipid_declarations.R:85` (`$da_preflight`); `func_lipid_da_model.R:413` (`runLipidsDA`) | KIT-106 |

**By owner:**

| Owner | Count | Proteomics | Lipid | Metab |
|---|---:|---:|---:|---:|
| KIT-102 | 2 | 2 | 0 | 0 |
| KIT-103 | 3 | 3 | 0 | 0 |
| KIT-104 | 4 | 3 | 0 | 1 |
| KIT-105 | 19 | 13 | 3 | 3 |
| KIT-106 | 4 | 3 | 1 | 0 |
| KIT-107, 108, 110 | 0 | 0 | 0 | 0 |
| **total** | **32** | **24** | **4** | **4** |

**Notes:**
- **Rows 4 and 5.** The callees are the exported `mod_prot_design_builder_ui`/`_server`, and the caller is the standalone design applet (`.runDesignMatrixApplet`). Exported module functions stay as thin wrappers (user decision), so these are not moved into the core. KIT-103 keeps both wrappers working for the applet and shows it with one test.
- **Row 13.** `stripProtNormQcPlotOwners` is also on the lipid and metab RUV path (`omicsRuvOperation`). Its move must leave the lipid and metab journeys exact.
- **Row 32.** `validateLipidDaAnalysisPreflight` builds on lipid's design validators, which KIT-103 may move first. KIT-106 owns this callee and consumes whatever KIT-103 leaves.
- **Outside the count.** One app-layer file calls a module function: `app_ui.R:62` calls `workflowPeptideSampleAckScript` (`mod_prot_qc_peptide_sample.R:28`). It is not core, so it is not one of the 32. KIT-104 keeps it callable.

## 6. Work map for KIT-102 to 110

**Common to every family ticket:**
- **Entry-point lanes** (KIT-101 W3): `prot_dia` (with enrichment, through the doubles), `prot_dia_limpa`, `prot_tmt`, `prot_lfq`, `prot_lfq_fragpipe`, `lipid_canonical`, `metab_lc`, `metab_gc`, `metab_combined`, `trx_screen`. Each lane runs from import to summary and report, so every lane passes through every family. "Primary" below marks the lanes whose captures the ticket changes or must prove. "Guard" marks the lanes that must stay identical.
- **A ticket's slice of the entry-point suite** is the captures at its tab's steps (committed state and history, `tab_status`, the state-content digest, outputs, artifacts) plus that tab's refusal cases, on its lanes. A ticket passes when its slice is identical to the frozen oracles, apart from the parity and unification items it owns. Each of those has a test that fails before the change.
- **Skip rule:** section 7.
- **Deletion:** each ticket deletes the old copies once nothing calls them, and quarantines their seam tests in the same ticket (user decision). Its back-references (section 5) move first.

### KIT-102: the shared app shell and import

| | |
|---|---|
| Files | `R/mod_lipidomics.R`, `R/mod_metabolomics.R`, `R/mod_proteomics.R`, `R/mod_transcriptomics.R` (shell `:200-286`, import `:76-146`), `R/mod_lipid_import_*` (6 files), `R/mod_metab_import_*` (7), `R/mod_prot_import_*` (10). New shared shell and import files, plus a core import file (KIT-102 names them; the scoping proposed `mod_omic_workflow.R`, `mod_omic_import.R` and `utils_workflow_omics_import.R`). `app_ui.R` and `app_server.R` are unchanged: dispatch is by name (`app_server.R:437,576`). |
| Functions moving to the core | the feature-table import: header read, preview, sample and column resolution, assay names, assembly, sanitising, validation, coercion, column mapping, artifacts and orchestration (scoping I1-I13; for example `runLipidImportProcessing` `mod_lipid_import_processing_helpers.R:618`, `runMetabImportProcessing` `mod_metab_import_processing_helpers.R:519`, `runTrxImportProcessing` `mod_transcriptomics.R:130`). Proteomics' readers, metadata, FASTA, organism and configuration move into the existing `func_prot_*` and `utils_workflow_dia_sessions.R` files (scoping P1-P6; for example `importProtImportDataByFormat` `mod_prot_import_processing_helpers.R:38`, `processProtImportFastaData` `:282`, `analyzeProtImportMixedSpeciesData` `mod_prot_import_organism_helpers.R:72`). |
| Config groups | `shell`, `binding` (the DIA adaptor's first entries and the contract validator), `import` |
| Back-references | #1 `recordDIAImportMetadata`, #2 `getDIAImportMetadata` |
| Parity owned | text-numeric sample columns converted for every omic, trx before its count policy; vendor formats stay per omic |
| Unification owned | metab import opens the Design step |
| Lanes | primary: all ten (each lane's import step and the shell's status chain) |
| Tests to retarget | MCI-004 (prot), MCI-012 (metab), MCI-018 (lipid), MCI-026 (trx) import matrices; `helper-workflow-omics-residency.R`; `test-workflow-server-shared.R` (stubs by namespace name) |
| Seam-test quarantine candidates | `test-lipid-00b` (file guard), `test-lipid-00c`, `test-metab-00b` ×2 (file guards), `test-metab-00c`, `test-prot-01c` (namespace guard), `test-prot-01d` |

### KIT-103: the shared design and design builder

| | |
|---|---|
| Files | `R/mod_lipid_design_*` (11 files), `R/mod_metab_design*` (5), `R/mod_prot_design*` (10), `R/mod_transcriptomics.R:149-198` |
| Functions moving to the core | validation, formula resolution, initial state and run order, design edits, save payload, object construction, commit and artifacts, design import (scoping C1-C10). Examples: `validateLipidDesignDaPreflight` (`mod_lipid_design_validation.R:449`), `validateMetabDesignDaPreflight` (`mod_metab_design_validation.R:429`), `validateLipidDesignAlignment` (`:87`), `buildLipidDesignInitialState` (`mod_lipid_design_builder_state_helpers.R:290`), `buildLipidDesignSaveResultsPayload` (`:559`), `runLipidDesignBuilderObserverShell` (`mod_lipid_design_builder_helpers.R:35`), `registerMetabDesignBuilderResultsObserver` (`mod_metab_design.R:195`), and the trx design helpers (`mod_transcriptomics.R:155-193`). Proteomics' UniProt fetch (`mod_prot_design_state_helpers.R:169`) and FASTA step (`mod_prot_design_import_helpers.R:127`) become the core `annotation_step`. |
| Config groups | `design` |
| Back-references | #3 `buildDIAReferenceDesignState`, #4 `mod_prot_design_builder_ui`, #5 `mod_prot_design_builder_server` (kept as wrappers for the applet) |
| Unification owned | metab's design-save failure is shown as a notification (`mod_metab_design.R:237-241`) |
| Lanes | primary: all ten (the design step; the fresh-design `Run` order asserted per omic) |
| Tests to retarget | MCI-027 (trx design, `helper-module-ci-trx-tabs.R`); `test-module-ci-{lipid,metab}-design.R` validator cases, with identical expectations; `test-prot-04b-design-server-shared.R` kept |
| Seam-test quarantine candidates | `test-lipid-13`, `13a`, `13b`, `14`; `test-metab-04`; the internal parts of `test-prot-04` (file guard, 182 blocks) and `test-prot-04a` |

### KIT-104: the shared QC shell with configurable sub-tabs

| | |
|---|---|
| Files | `R/mod_lipid_qc*` (7 files), `R/mod_metab_qc*` (12), `R/mod_trx_qc_filter.R`, and the proteomics orchestrators only: `R/mod_prot_qc.R`, `mod_prot_qc_peptide.R`, `mod_prot_qc_protein.R`. The 17 proteomics sub-tab servers are not changed. |
| Functions moving to the core | the counts table, the in-process QC action, duplicate detection, the resolution statistics, ITSD analysis and its figures, finalise validation (KIT-104 scoping §2.2). Examples: `buildMetabIntensityFilterStats` (`mod_metab_qc_intensity.R:173`), `handleLipidDuplicateDetection` (`mod_lipid_qc_duplicates_server_helpers.R:187`), `detectMetabDuplicateFeatures` (`mod_metab_qc_duplicates_server_helpers.R:96`), `analyzeMetabQcItsdData` (`mod_metab_qc_itsd.R:57`), `readMetabQcItsdState` (`:48`, mocked by `omicsStandInReaders`), `validateLipidQcS4FinalizeState` (`mod_lipid_qc_s4.R:347`), `validateMetabQcS4FinalizeState` (`mod_metab_qc_s4_server_helpers.R:320`), and `trxQcRawStats`, `trxQcMetricsTable`, `trxQcFilteringFigure` (`mod_trx_qc_filter.R:29,46,55`) |
| Config groups | `qc` |
| Back-references | #6 `.formatPeptideRollupSummary`, #7 `.peptideImputationAuditSummary`, #8 `.proteinCleanupActiveKey`, #9 `resolveMetabDuplicateAssayData`; and the comment at `utils_workflow_prot_declarations.R:271` |
| Parity owned | the QC-finalise design-alignment check runs for every omic (metab and trx refuse a misaligned design; lipid already does) |
| Lanes | primary: `lipid_canonical`, `metab_lc`, `metab_gc`, `metab_combined`, `trx_screen`, `prot_dia`, `prot_dia_limpa` (peptide and protein groups); guard: `prot_tmt`, `prot_lfq`, `prot_lfq_fragpipe` (protein group, the QC bypass) |
| Tests to retarget | MCI-014, MCI-020, MCI-028; `test-workflow-omics-application-reads.R` (the KIT-082 proof, which drives `mod_lipid_qc_s4_server(profile = trxOmicProfile())`); `test-omics-audit-carry-overs.R` (`countQcAutoInitialisation`) |
| Seam-test quarantine candidates | `test-lipid-01b` to `01k`; `test-metab-01o` to `01w` (`01p`, `01s`: file guards); `test-metab-qc-itsd-direct-shared`; `test-prot-02b`, `02c` (namespace guards) only if the orchestrator's internals go |

### KIT-105: normalisation core moves and the shared lipid/metab frame

| | |
|---|---|
| Files | `R/mod_lipid_norm_*` (8 files), `R/mod_metab_norm_*` (7), the shared panel library and frame (new). Proteomics: only the 13 back-referenced functions move, out of `R/mod_prot_norm_*`, into the core, plus the shared option sections (`between_sample`, `aesthetics`, `ruv_iii`) if the entry suite covers its routes. Its figure frame and run stay its plugin; placement is KIT-111 to 113. |
| Functions moving to the core | the composite figure (one copy instead of `writeLipidNormCompositeFigure`/`runMetabNormCompositeQcFigureStep`), the composite from files, the final PCA, the in-process normalisation sequence (`handleLipidNormRunNormalization` `mod_lipid_norm_workflow_helpers.R:102`; `runMetabNormNormalizationPipelineShell` `mod_metab_norm_pipeline_helpers.R:52`), and the session-export collect and save (`buildLipidNormExportSessionData` `:329`, `saveLipidNormExportArtifacts` `:402`, `buildMetabNormExportSessionData` `mod_metab_norm_session_helpers.R:499`). `getPlotAesthetics` (`mod_metab_norm_render_helpers.R:230`) is used across omics and moves to a shared file. |
| Config groups | `norm` (the frame, the shared panel library, the lipid and metab values) |
| Back-references | #10 to #28 (13 proteomics, 3 lipid, 3 metab); and the comment at `utils_workflow_omics_actions.R:208` |
| Parity and unification owned | none: every committed state, file and figure stays identical |
| Lanes | primary: `lipid_canonical`, `metab_lc`, `metab_gc`, `metab_combined`; guard: the five prot lanes (the back-reference moves, including `stripProtNormQcPlotOwners` on the omics RUV path) and `trx_screen` (trx mounts the lipid norm module until KIT-110) |
| Tests to retarget, never quarantine | `helper-workflow-omics-worker.R`, `test-workflow-omics-application-reads.R`, `-figures.R`, `-sessions.R`, `test-omics-audit-carry-overs.R`: they carry the KIT-078/082 managed-path evidence. MCI-015 and MCI-021 handlers re-pointed. |
| Seam-test quarantine candidates | `test-lipid-12` (file guard, 69 blocks), the internal parts of `test-lipid-12a`, `test-metab-03b` (silent symbol loader), `test-metab-03c`, the internal parts of `test-metab-03d` |

### KIT-106: the shared DA module

| | |
|---|---|
| Files | `R/mod_lipid_da*` (8 files), `R/mod_metab_da*` (6), `R/mod_trx_da.R`; proteomics shares the UI and the render registrations (`R/mod_prot_da_*`, 9 files), and its run, session load and init handlers stay as adaptors |
| Functions moving to the core | `mod_prot_da_validation.R` (contrast normalising, model and contrast validation, filtered-session validation); lipid's DA preflight (`mod_lipid_da_analysis_helpers.R:84`). Result reshaping for display stays in the shared module as pure view functions. The in-process interactive volcano readers (`lipidDaResultsWithSource` `mod_lipid_da_session_helpers.R:294`, `metabDaResultsWithSource` `mod_metab_da_observer_helpers.R:287`, `trxDaWithSource` `mod_trx_da.R:93`) are replaced by the worker page. |
| Config groups | `da` |
| Back-references | #29 `normaliseProtDaContrastsTable`, #30 `validateProtDaModelAndContrasts`, #31 `validateProtDaFilteredSession`, #32 `validateLipidDaAnalysisPreflight` |
| Parity owned | the metab DA tree-cut defect is fixed (UI inputs `heatmap_tree_cut_method` `none`, `heatmap_n_clusters` 4, `heatmap_cut_height` 10, `heatmap_min_cluster_size` 5, as lipid's) |
| Unification owned | the interactive volcano is off by default and built on demand for every omic |
| Lanes | primary: all ten (the DA step; the volcano default and the on-demand page; the metab heatmap) |
| Tests to retarget | MCI-009 (prot; MCI-009.4 pins `volcano_glimma` and `volcano_plot_static`), MCI-016, MCI-022, MCI-031; `test-workflow-omics-da-contrasts.R`, `test-workflow-saved-da.R` |
| Seam-test quarantine candidates | `test-lipid-02b`; `test-metab-02aa` ×2 (`-shared`: file guard); `test-module-ci-{lipid,metab}-da` internals; `test-prot-07b`, `test-prot-da-render-direct-shared`; `test-prot-07`, `09` (file guards) if their files move |

### KIT-107: the enrichment module, omic-agnostic, proteomics its only consumer

| | |
|---|---|
| Files | `R/mod_prot_enrich*` (10 files). Unchanged: the kernels `R/func_general_enrichment_*.R` and the multi-omics helpers `R/func_multiomics_enrich_*.R` (KIT-063). |
| Functions moving to the core | the analysis set (`mod_prot_enrich_analysis_helpers.R`, all 7, for example `runProtEnrichAnalysisBody` `:489`, `executeProtEnrichProcessEnrichments` `:119`); the execution set and the deterministic double (`mod_prot_enrich_execution_helpers.R`, the double at `:126`); the input set (`resolveProtEnrichOrganismMapping` `mod_prot_enrich_input_helpers.R:1`, `applyProtEnrichOrganismFilter` `:182`, `resolveProtEnrichUniprotAnnotations` `:599`); the state set (`saveProtEnrichCompletedState` `mod_prot_enrich_state_helpers.R:210`, `restoreProtEnrichResultsFromState` `:588`); `buildProtEnrichSupportedOrganisms` (`mod_prot_enrich_observer_helpers.R:164`) |
| Config groups | `enrichment` |
| Back-references | none |
| Parity and unification owned | none |
| Lanes | primary: `prot_dia` (both backends on the doubles; the organism-filter case where the W3 suite has it); guard: the other four prot lanes (the enrichment status transitions into the summary) |
| Tests to retarget | E2E-007 (`test-e2e-proteomics-enrichment.R`, both backends); the KIT-016 equality test in `test-workflow-dia-products.R`; MCI-010 |
| Seam-test quarantine candidates | `test-prot-11` (namespace guard, about 167 internal blocks); `test-prot-11a` (silent binding guard, source by path); MCI-009.6 internals |

### KIT-108: the shared summary and report; trx's RUVIII export as an omic step

| | |
|---|---|
| Files | `R/mod_prot_summary_*` (6 files), `R/mod_lipid_summary*` (4), `R/mod_metab_summary*` (3), `R/mod_trx_summary.R`; a new core summary file (scoping: `R/utils_workflow_summary.R`). `copyToResultsSummary` (`func_general_filemgmt_results_summary_helpers.R:16`) and `RenderReport` (`func_general_filemgmt_report_render_helpers.R:63`) stay exported with the same signatures. |
| Functions moving to the core | save parameters, the integration-export dispatch, the publication copy, the report, the session-state RDS, a generic export-step runner (KIT-108 scoping §2). Examples: `handleLipidSummarySaveWorkflowArgs` (`mod_lipid_summary_workflow_helpers.R:1`), `runMetabSummarySaveWorkflowArgsObserverShell` (`mod_metab_summary_observer_helpers.R:231`), `completeProtSummaryWorkflowArgsSave` (`mod_prot_summary_session_helpers.R:1`), `handleLipidSummaryGenerateReport` (`mod_lipid_summary_workflow_helpers.R:352`), `runMetabSummaryGenerateReportObserverShell` (`mod_metab_summary_observer_helpers.R:699`), `resolveProtSummaryReportTemplate` (`mod_prot_summary_report_helpers.R:131`), `runTrxRuviiiExport` (`mod_trx_summary.R:65`). The four `runLipidSummary*ObserverShell` (`mod_lipid_summary_observer_helpers.R:35,64,94,251`) have no R caller (findings F-04). |
| Config groups | `summary` |
| Back-references | none |
| Parity and unification owned | none: every behaviour difference stays per omic as configuration. The report-step memory item (revision 1.1) is its own diagnosis-then-fix. |
| Lanes | primary: all ten (the summary, report and export steps; trx's RUVIII export and the integration gate) |
| Tests to retarget | MCI-032 (trx full tabs); the about 20 `testServer` blocks in `test-prot-12` |
| Seam-test quarantine candidates | `test-lipid-14a`, `15`; `test-metab-05`; the seam parts of `test-prot-12` (namespace guard) and `12a`; MCI-011, MCI-017, MCI-023 (MCI-023 tests the dead lipid shells) |

### KIT-109: remove the profile hooks and the duplicate copies, quarantine, full suite, closure

| | |
|---|---|
| Files | what the family tickets left: the `profile` hooks in `R/mod_lipid_*` (`grep -nE 'is\.null\(profile\)\|profile = NULL\|trxOmicProfile\|trxNormDaSummaryProfile' R/` must be empty), `R/mod_transcriptomics.R`'s two profile constructors, any copy still unreferenced; `tests/testthat/test-omics-module-globals.R` (file pattern and allowlist); `tools/ci/impact-map.json` and `tools/ci/detect-impact.R` (a shared-layer rule and the trx lanes); `tools/ci/run-module-ci.R` (trx branches) |
| Config groups | none new; it checks that the 65 trx profile fields (§4.12) and the 60 lipid `profile$` reads all map or leave |
| Back-references | checks that none of the 32 remains: a parse-data search of `R/utils_workflow_*.R` and `R/func_*.R` for names defined in `R/mod_*.R` finds only the exported wrappers kept on purpose (#4, #5) |
| Owns | the Q-23 onward rows consolidated in `docs/qa/workflow-execution/kit/quarantined-tests.md`, from each family ticket's draft rows |
| Lanes | all ten, plus the full suite against KIT-101's per-file baseline under both load modes |

### KIT-110: transcriptomics' TMM/RUV-grid panel set on the shared frame

| | |
|---|---|
| Files | `R/mod_trx_norm.R` (235 lines); the trx entries `norm_*` in `R/mod_transcriptomics.R:56-61`; the lipid norm `profile` branches KIT-105 leaves |
| Functions moving to the core | `runTrxPublicationInProcess` (`mod_trx_norm.R:85`) and `runTrxNormPublicationsInProcess` (`:106`), generic over `binding$publications`; the run (`runTrxNormalisation` `:156`, `runTrxNormManaged` `:130`) becomes the trx plugin's `run`; `trxNormGridParameters` (`:57`) becomes its settings |
| Config groups | `norm` (trx values: `tmm_logcpm`, `ruv_grid`) |
| Back-references | none (trx has none) |
| Parity and unification owned | none |
| Lanes | primary: `trx_screen` (exact against KIT-098's oracle and the trx journeys under protocol v1.9); guard: `lipid_canonical`, `metab_lc`, `metab_gc`, `metab_combined` (the shared frame) |
| Tests to retarget | MCI-030 (`test-module-ci-trx-norm.R`; MCI-030.4 mounts through `mod_lipid_norm_server(profile =)`) |

### 6.1 Who owns each parity decision and unification

| Item | Kind | Owner |
|---|---|---|
| Text-numeric sample columns converted for every omic | parity | KIT-102 |
| Vendor formats stay per omic | parity | KIT-102 |
| The QC-finalise design-alignment check runs for every omic | parity | KIT-104 |
| The metab DA tree-cut defect is fixed | parity | KIT-106 |
| Metab import opens the Design step | unification | KIT-102 |
| Metab's design-save failure is shown as a notification | unification | KIT-103 |
| The interactive volcano is off by default and built on demand, for every omic | unification | KIT-106 |

Each needs a test that fails before the change and passes after. Every other difference stays per omic as configuration (user decision). That includes the summary's state choice, record writing and copy-failure handling, lipid's zero-row refusal, the finalise hand-off to normalisation, and the assay-name rules.

### 6.2 Lanes by ticket

P = primary, G = guard.

| Lane | 102 | 103 | 104 | 105 | 106 | 107 | 108 | 109 | 110 |
|---|---|---|---|---|---|---|---|---|---|
| `prot_dia` | P | P | P | G | P | P | P | P | — |
| `prot_dia_limpa` | P | P | P | G | P | G | P | P | — |
| `prot_tmt` | P | P | G | G | P | G | P | P | — |
| `prot_lfq` | P | P | G | G | P | G | P | P | — |
| `prot_lfq_fragpipe` | P | P | G | G | P | G | P | P | — |
| `lipid_canonical` | P | P | P | P | P | — | P | P | G |
| `metab_lc` | P | P | P | P | P | — | P | P | G |
| `metab_gc` | P | P | P | P | P | — | P | P | G |
| `metab_combined` | P | P | P | P | P | — | P | P | G |
| `trx_screen` | P | P | P | G | P | — | P | P | P |

A "—" means the lane does not reach that ticket's code. The full lane still runs in the ticket's checkpoint wherever the plan's test rule asks for it.

## 7. The skip rule for KIT-102 to 109 (from W4, verbatim)

> **Rule for KIT-102 to 109,** written into `design.md` and each stub: any new file-level skip, any new hidden blocks, or a drop in tests run fails, unless a quarantine entry in `quarantined-tests.md` (Q-23 onwards) names the internal that replaced it.

KIT-110 is the second half of KIT-105's split, so it follows the same rule. The comparison is per file and per load mode, against KIT-101's per-file baseline (`kit-101/test-baseline.tsv`, W4). The strict-skip sentinel (`MULTISCHOLAR_STRICT_SKIPS=1`, W2) turns the 13 file-level `R/`-path guards, the 7 namespace guards and the three partial loads into failures, so a family ticket runs its targeted files with it set.

## 8. Decisions left for the family tickets' authoring

None of these blocks KIT-101. Each is a default the authoring ticket confirms or raises with the user:
- **KIT-102:** the names of the shared files, the DIA adaptor and the config builder. Whether any proteomics reader leaves text-numeric quantity columns, which decides whether the parity change touches proteomics at all.
- **KIT-103:**
  - run order is recorded in the configuration and the evidence only, not in the `design_builder` provenance. A provenance change would alter the record the KIT-098 and v1.7 oracles cover, and needs the user.
  - Metab writes the bus before it checks `source_dir` (`mod_metab_design.R:243-252`), where lipid checks first. A single shared commit changes metab's order on that error path. This is not one of the three unifications, so the author confirms with the user, or keeps metab's order.
  - Whether to drop the `DEBUG66` messages in moved code.
- **KIT-104:** none beyond the per-omic flags in §4.6.
- **KIT-105:** proteomics takes only the shared option sections, and only with entry-suite coverage of its routes. Its figure frame and run stay its plugin.
- **KIT-106:** writing `Interactive_Volcano_Plots/` files for lipid, metab and trx (no output map needs them; default no). Heatmap clustering on the worker (new machinery; out of scope).
- **KIT-107:** whether the legacy non-manifest branch, reachable only from tests, is kept through the move or retired under the quarantine rule.
- **KIT-108:** none beyond revision 1.1's report-memory sequence. A fix that needs worker dispatch reopens D3 with the user.
