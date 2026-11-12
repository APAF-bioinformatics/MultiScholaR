# KIT-120 W0: sites and map (2026-11-12)

Read at `8ca68e93`. Reading only: no product code changed, no runs. Every line reference is to that commit.

## Summary

- **The join is feasible without touching the harness.** Every input id the journey harness sets or clicks exists in both frames or is configurable. The module id stays `normalization`.
- **It retires about 1,150 lines, not 3,900.** The ticket's figure counted every `mod_prot_norm_*.R` file. About 2,700 of those lines are action handlers (run, RUV, correlation, export, reset, QC generation) that the core and worker depend on for their field names and file products. They stay, behind declared hooks, as trx's run handler does.
- **Three findings change the plan** (section 6): a scientific default the frame does not have, the Log tab, and the image decision.

## 1. The two frames today

| | Shared frame | Proteomics |
|---|---|---|
| Files | `mod_omic_norm.R` (537), `mod_omic_norm_panels.R` (550) | ten `mod_prot_norm_*.R` files (3,927 lines) |
| UI | sections and tabs chosen by `norm$panels` and `norm$plot_tabs` | one fixed layout, `mod_prot_norm_ui.R:15-486` |
| Server | `runOmicNormServer`, `mod_omic_norm.R:59-112` | `runProtNormModuleServerShell` and `registerProtNormServerObservers`, `mod_prot_norm_observer_helpers.R:1-369` |
| Config | `omicGuiConfig(<omic>)$norm` | `norm = NULL`, `utils_omic_gui_config.R:1162` |
| Mounted by | registry `feature()` | registry `prot_norm = component("prot_norm", "tab_selected_positional")`, `mod_omic_workflow.R:225` |

The two layouts descend from the same code. Their option columns hold the same inputs with the same defaults, and their QC tabsets hold the same tabs in the same order, apart from the rows below.

## 2. Element map

Fate codes: **F** = the frame already has it; **C** = the frame has it once a config value is declared; **N** = a new library entry or contract field; **H** = stays proteomics code behind a declared hook; **D** = deleted.

### Option column (`mod_prot_norm_ui.R:19-254` against `mod_omic_norm_panels.R:41-158`)

| Element | Proteomics | Frame | Fate |
|---|---|---|---|
| Section order | method, aesthetics, RUV, actions, export | the order of `norm$panels` | C |
| `norm_method` | four methods, default `cyclicloess`, test id `prot-norm-method` | `between_sample` with `methods`, `default`, `testids` | C |
| Method help text | "(default: Cyclic Loess)" | fixed text | N: optional `help` per section (decision 1) |
| Aesthetics | `h4` heading, same two selects | `h5` heading | F: takes the frame's heading level |
| RUV modes, grouping, automatic, manual, skip panels | same inputs, same defaults | `ruv_iii` section | C for `modes`, `default`, `percentage = 5`, `k = NULL`, `noun = "Proteins"` |
| RUV help texts (grouping, automatic x4, manual k, skip x2) | longer | fixed, lipid's | N: optional texts in `panels$ruv_iii` (decision 1) |
| Run button | "Run Normalization & RUV", no icon, test id `prot-norm-run` | "Run Normalization Pipeline", play icon | N: optional label; takes the frame's icon |
| Reset button | no icon | undo icon | F |
| Export button | `export_filtered_session`, "Export Filtered Data Session", test id `prot-norm-export-session`, its help text | `panels$export$button` (configurable id), fixed label and help | C for the id; N for label and help |
| S4 download | `export_s4_compatibility` download button and help | none | N: an `s4_download` section |
| `hr()` positions | before each section but the first | before each section | F: one extra rule above the method |

### QC tabset (`mod_prot_norm_ui.R:258-486` against `mod_omic_norm_panels.R:182-232`)

| Tab | Proteomics | Frame | Fate |
|---|---|---|---|
| RUV QC | `ruv_canonical_correlation_plot`, `ruv_optimization_summary`, `ruv_optimization_table`; one result | `ruv_qc_ui`: per assay `cancor_plot_<assay>`, `ruv_summary_<assay>`, `ruv_table_<assay>`, reading `normData$ruv_optimization_results[[assay]]` | C for heading and help; **H** for the content: proteomics' result has a different shape and a richer summary (`buildProtNormRuvOptimizationSummary`, with the mode and grouping) |
| PCA, Density, RLE, Correlation | twelve fixed outputs `<kind>_post_filtering` and siblings | `<kind>_<stage>_assay<n>` | F, with the image finding of section 6 |
| Post-Normalisation QC | `post_norm_filtering_summary`, `filtering_summary_text` | none | N: a library tab only proteomics declares (decision 2); its two renderers stay proteomics' (H) |
| Correlation Filtering | same title, slider, two buttons, `correlation_filter_summary`, `final_qc_plot`, plus `ruv_covariates_panel` under the summary | `omicNormCorrelationFilterTab` | F for the layout; N for a slot under the summary; H for the two outputs (section 3) |
| Log | none | `norm_log` | not declared (section 6) |

### Server (`mod_prot_norm_observer_helpers.R` against `mod_omic_norm.R:59-112`)

| # | Proteomics element | Frame today | Fate |
|---|---|---|---|
| 1 | Run: `runProtNormNormalizationWorkflow` | `norm$run$handler` (trx's shape) | H: the handler hook exists; its argument list must widen, or proteomics' handler closes over what it needs |
| 2 | Reset: `runProtNormResetObserver` | `norm$reset` | H: the hook exists |
| 3 | Correlation apply and skip: `runProtNormApplyCorrelationObserver`, `runProtNormSkipCorrelationObserver`; they write `correlation_filter_summary` themselves | the frame's own `applyOmicNormCorrelation`, no hook | N + H: a correlation handler pair, and the omic owns the two outputs |
| 4 | Export: `runProtNormExportObserver` | the frame's own `exportOmicNormSession` | N + H: an export handler |
| 5 | Tab entry: `runProtNormTabEntryWorkflow`, `ignoreInit = TRUE`, with F116-01's QC-status rule | `runOmicNormTabEntry`, `ignoreInit = FALSE`, class check only | N + H: an entry handler. The warning rule is behaviour and must not change |
| 6 | Aesthetic regeneration on colour or shape change | none | N + H |
| 7 | Restore: `observeEvent(workflowData$native_session_restored, ...)`, `:193-244` | none | N + H |
| 8 | DIA_LIMPA mode radio (`updateRadioButtons` once the route is known) | none | N + H: the route is known only at run time, so it cannot be static config |
| 9 | RUV covariates panel and `apply_ruv_covariates` observer | none | N + H |
| 10 | S4 download handler | none | N + H |
| 11 | `installWorkflowOmicsActions(session, workflowData, "proteomics", ...)` | the omics frame installs its own runner | H: part of proteomics' server hook |
| 12 | Design-driven choices | `omicNormDesignChoices` | **differs, see section 6** |
| 13 | Reactive state | `omicNormReactiveState` | N: proteomics' extra fields (`normalized_protein_obj`, `ruv_normalized_obj`, `best_k`, `control_genes_index`, `ruv_optimization_result`, `correlation_threshold`, `publication_paths`, and others) declared as initial state. The core reads them by name (`utils_workflow_dia_products.R`, `utils_workflow_dia_norm.R`, `utils_workflow_dia_sessions.R`), so no field is renamed |
| 14 | Plot aesthetics getter | `omicNormPlotAestheticsGetter` | F: the same rule |

Items 5 to 11 are proteomics-only registrations. The smallest contract is **one optional server hook** in `norm` (the omic registers its own observers and outputs with the frame's state and callbacks) plus **handler fields** for correlation, export and entry, beside the existing run and reset ones. W1 fixes the exact fields.

## 3. What is retired, and what stays

| Retired at W5 (about 1,150 lines) | Stays as proteomics' handlers (about 2,700 lines) |
|---|---|
| `mod_prot_norm_ui.R` body (486): becomes a wrapper | `mod_prot_norm_workflow_helpers.R` (786): the dispatchers and finishers of KIT-111, the entry rule |
| `mod_prot_norm_observer_helpers.R` shells (369): the registrations move to the hook | `mod_prot_norm_correlation_helpers.R` (441) |
| `mod_prot_norm_render_helpers.R`: the image renderer and registrar, `registerProtNormRenderOutputs`, the reactive state (about 300 of 396) | `mod_prot_norm_ruv_helpers.R` (356): the core calls `resolveProtNormRuvParameters` |
| | `mod_prot_norm_session_helpers.R` (478): export and reset |
| | `mod_prot_norm_qc_generation_helpers.R` (291), `mod_prot_norm_qc_support_helpers.R` (305) |

**Not done, with the reason:** running proteomics through the frame's own pipeline (`runOmicNormPipeline`, `applyOmicNormCorrelation`) instead of its handlers. The publication half writes proteomics' files and fields, and the journeys are exact on them. Unifying the handlers is a second, larger change with the journeys as its only safety net. It is not needed for this ticket's outcome.

## 4. Consumers

**Registry and exports:** `mod_omic_workflow.R:225`; `NAMESPACE:305-306` (`mod_prot_norm_server`, `mod_prot_norm_ui`); `man/normalizationAppletModule.Rd`. No shipped notebook or vignette calls the tab.

**Core callers of `*ProtNorm*` functions:** all but one are defined in `utils_workflow_dia_norm.R` already. The exception is `resolveProtNormRuvParameters` (`mod_prot_norm_ruv_helpers.R`), which stays.

**Tests** (24 files name a proteomics normalisation function). By weight: `test-prot-05b-norm-module-contracts.R` (233 references), `test-prot-05d-norm-module-characterization.R` (178), `test-prot-norm-core.R` (94), `test-prot-05c-norm-module-seams.R` (24), `test-workflow-saved-qc.R` (19), `test-module-ci-prot-norm.R` (16), `test-workflow-saved-exports.R` (13), `test-workflow-native-session.R` (12), `test-workflow-memory-mechanism-member.R` (10); fifteen more with one to seven. Most test the handlers, which stay. The functions W5 retires are named in one to three files each (the shells, the image renderer and registrar, `registerProtNormRenderOutputs`, `createProtNormReactiveState`); those blocks are the quarantine candidates.

**Entry lanes:** `helper-gui-entry.R:709-724` takes the defaults from `mod_prot_norm_ui` (the wrapper keeps working) and sets `export_filtered_session`, `ruv_mode`, `run_normalization`, `norm_qc_tabs`, `min_pearson_correlation_threshold`, `apply_correlation_filter`, `apply_ruv_covariates`. All are kept.

**Browser tests:** test ids `prot-norm-run`, `prot-norm-export-session` (`helper-e2e-browser.R`, `test-e2e-test-mode-integration.R`); also `prot-norm-method`, `prot-norm-ruv-mode` in the UI. All four fit `norm$testids`. No e2e test names an old plot output id, so W6 may have nothing to change.

**Module CI:** `helper-module-ci-prot-norm.R` reads `filtering_summary_text`, `final_qc_plot`; `test-module-ci-prot-norm.R` reads `correlation_filter_summary`. All three ids are kept.

**Harness** (`tools/profiling`, not edited here):
- Read at run time: the prefix `proteomics_workflow-normalization-`; the inputs `norm_method`, `ruv_mode`, the RUV settings, `run_normalization`, `norm_qc_tabs` = "Correlation Filtering", `min_pearson_correlation_threshold`, `apply_correlation_filter`, `apply_ruv_covariates`, `export_filtered_session`; the outputs `correlation_filter_summary` (text "Correlation filtering completed successfully!") and `ruv_covariates_summary`; the notifications (`workflow_kit_driver.R:194-296`, `workflow_checkpoint_dia_driver.R:1876-1897`). **All are kept, so no harness file needs to change for the journeys to run.**
- Labels only: the step table names the module `mod_prot_norm` (`workflow_kit.R:71-107`), and `kit-001/step-map.json:491-592` cites lines in `mod_prot_norm_workflow_helpers.R` and `mod_prot_norm_ruv_helpers.R`. Both files stay, so the module label stays true; the cited line numbers may shift. That is KIT-117's kind of correction, and nothing here forces it.
- `bindings-generator.R` binds no `mod_prot_norm_*` file.

## 5. Exit members (from `bindings.json`, the KIT-112 set)

| Member | Route | Corpus | Oracles | Use |
|---|---|---|---|---|
| `m` | DIA_LIMPA, covariates | `screen` | `limpa65_screen` | exact |
| `l` | protein (LFQ) | `lfq_fixture` | `lfq101_fixture` | exact |
| `t` | protein (TMT) | `tmt_fixture` | `tmt101_fixture` | exact |
| `q` | DIA | `full` | none | memory comparison |

Two gaps to settle when the runs are declared: the recent sets hold **no exact DIA (IQ) member** and **no DIA_LIMPA correction-mode member**. Either earlier members are re-bound for them, or the entry lanes (`prot_dia`, and `prot_dia_limpa`) are the evidence for those routes. The ticket's "both RUV modes" is not met by the current set as it stands.

## 6. Findings that change the plan

1. **A scientific default the frame does not have.** Proteomics preselects the RUV grouping variable as `batch` when the design has a batch column with values, else `group` (`mod_prot_norm_qc_support_helpers.R:157-170`; this is the TMT default). The frame always preselects `group` (`mod_omic_norm.R:159-176`). The join must keep proteomics' rule: `updateProtNormDesignDrivenChoices` stays, behind the server hook or a config field. A test must pin it before the switch.
2. **The Log tab.** The frame's log shows what its own pipeline writes through `addLog`. Proteomics' handlers write nothing to it, so the tab would be empty. The tabs are chosen by `norm$plot_tabs`, so proteomics simply does not declare it. The ticket said proteomics gains the Log tab: **this reverses that line.** Filling the log would mean editing every handler, which is new behaviour.
3. **Decision 3, the images: neither option as written.** KIT-111 D5's "images from recorded paths" is out of date: the display reads twelve fixed file names in `protein_qc_dir` (`pre_norm_pca.png`, `post_norm_pca.png`, `ruv_corrected_pca.png` and so on; `mod_prot_norm_render_helpers.R:20-93`). The frame builds `<stem>_<stage prefix>_<kind>.png`. Proteomics' stage prefixes are already the frame's defaults; only the stem is absent. **Answer: the frame accepts an empty stem (no stem, no underscore) and a QC directory key from `norm` (`protein_qc_dir`; `proteomicsDiaBinding` has no `qc_dir_role`).** No file is renamed.
4. **A visible addition:** the frame puts an assay label above each image row ("Assay: <name>"). Proteomics has one table, so it would show one label. Either it shows a fixed name, or the label row is omitted for a single unnamed assay (a small optional field). To decide at W1.
5. **Sizing.** With the handlers staying, the work is mostly contract, config and one registration hook. The estimate of 24 hours stands as an upper bound; the larger risk is the entry lanes' data captures at W4, where the stop rule applies.

## 7. Decision

W0 finds nothing that blocks the join or needs the harness changed. The plan proceeds with the corrections of section 6, once the user has seen findings 2 and 4 and the member gaps of section 5.
