# Stage 3b audit (KIT-072–080), 2026-11-01

Mode `audit_only` (user: "audit the ticketed work for completeness, rigour and reference accuracy"). Scope: `tickets/kit/KIT-072.md` … `KIT-080.md`, `tickets/kit/00-KIT-PLAN.md` v3.1 (the stage-3b rows, ledger and decision row), the README stage line; tree at `a7ac9d06`, clean. No ticket, plan, code or binding was changed by this audit; repairs need the user's go and are tracked in a separate remediation record.

Method: the read-only inventory (`audit_inventory.py`: graph, statuses, required sections, links), a mechanical check of every `file:line` and backticked identifier, then a line-by-line semantic check of every reference and factual premise (three read-only passes, each claim opened at its cited line, ±5 tolerance), the highest-impact claims re-checked directly, and a decision-design review against protocol v1.5.

## Summary

| Check | Result |
|---|---|
| Graph, statuses, required sections | 78 tickets, no cycle, no duplicate, no missing dependency; the 9 new tickets status-consistent with all required sections; selector returns KIT-072 |
| Links | 3 "missing links" are false positives (`[[field]]`, `[[assay]]`, `[[omic]]` R subscripts read as wiki links) |
| Files and line ranges | every cited existing file exists and every line is in range; 5 cited files are created by earlier stage-3b tickets (forward references; acceptable only in order) and 1 does not exist (`test-workflow-saved-table-codec*.R`) |
| Line-level accuracy | ~296 claims checked: ~258 confirmed, 18 wrong or misattributed, ~20 factual premises wrong or incomplete |
| Completeness | each ticket's touchpoint list misses sites (KIT-072, 075, 076, 077, 078 most) |
| Rigour | 4 findings block implementation or execution as written; 9 high |

Verdict: the series' shape and order hold, and most references are right, but **KIT-079 (ITSD), KIT-080 (operation pairs, control, metric) and KIT-078 (metric, gates) cannot meet their own acceptance as written**, and the shared "exactness set" is under-specified against protocol v1.5. The tickets need one revision pass before KIT-072 starts; KIT-072–074 need only reference and list corrections.

## Findings

Severity: **B** blocks implementation or execution as written; **H** the acceptance cannot support its claim or a significant touchpoint is missing; **M** incomplete list or wrong premise with a bounded effect; **L** reference or wording error.

### Cross-cutting

**S3B-01 (H) The exactness set's oracle source contradicts protocol v1.5 §7.1.** Every ticket re-runs `kit-journey-031l`/`032m` "against the bound oracles lipid31/metab32", which were recorded in process on `kit_031`/`kit_032`. §7.1: an in-process oracle "proves that the owner path reproduces the in-process path of the same kernels on the same build; it is not a cross-build claim". A later build judged against them is a cross-build comparison §7.1 does not define, and a member on a build that changes an in-process path (KIT-075 hydration, KIT-078 the worker) would be judged against an oracle that cannot see the change. Also undeclared: the per-ticket bindings (`0NNl`/`0NNm` single-journey keys, `kit_0NN` installs, `KIT_0NN_LIFECYCLE`), and `bindings-generator.R` is in no stage-3b `files_to_modify` except KIT-076/080. *Repair:* define the exactness set as (a) the recorder re-run on the ticket's build (the §7.1 oracle), (b) the member exact against (a), and (c) (a) equal to the stage-3 oracles state by state (the no-change claim, cross-build, as `kit-070/compare-oracles.R` does for proteomics); add the generator to each ticket's files. *Blocks:* execution of every ticket's gate.

**S3B-02 (H) Forward references.** KIT-074/075/076/078/079 cite `R/utils_workflow_omics_declarations.R` and `test-workflow-omics-declarations.R`, which KIT-073 creates. Valid only in order; the current locations are `utils_workflow_lipid_declarations.R:266` (core) and `:181-192` / metab `:164-175` (Pearson). *Repair:* note "created by KIT-073" where cited. *Blocks:* nothing if the order holds.

### KIT-072 (one path)

**S3B-03 (L) Mis-cited sites.** `readLipidDaSessionData` at `mod_lipid_da_session_helpers.R:340`: the branch is in `restoreLipidDaCurrentS4FromSession` (def `:296`, branch `:340`, double `saveState` `:347`). `handleLipidNormExportSession` (`:567`) and `runMetabNormExportSessionWorkflow` (`mod_metab_norm_session_helpers.R:740`) already refuse a double (KIT-071), so there is nothing to collapse. `test-metab-02h`, `test-metab-da-orchestration-direct-shared` use no state manager. `test-workflow-server-shared.R` uses an R6 stub from a mocked `createDIASessionOwner`, not a list double. `test-metab-02aa` and `test-metab-04` each match two files; only `-02aa-da-display-filter-characterization.R` and `-04-design-module-characterization.R` use doubles.

**S3B-04 (M) Omitted double paths.**
- Module sites: `mod_metab_norm_session_helpers.R:110` (`handleMetabNormApplyCorrelationOutcome`), `mod_metab_qc_duplicates_server_helpers.R:372` (`applyMetabDuplicateResolutionState(commit=)`), `mod_metab_qc_intensity.R:429` (`saveMetabIntensityFilterState`, `:228`), `runMetabNormRuvCorrectionStep` (`mod_metab_norm_pipeline_helpers.R` ~`:740-765`), the `!owned` sub-branches at `mod_lipid_qc_intensity_server_helpers.R:209` and `mod_lipid_norm_workflow_helpers.R:110,132,155,226,756`.
- The DA Open second path: `mod_lipid_da_session_helpers.R:476-477`, `mod_metab_da_observer_helpers.R:555-556`.
- Design import: `mod_metab_design_import_helpers.R:415`.
- Shared pass-throughs: `ingestWorkflowOmicsLegacySession` (`utils_workflow_dia_sessions.R:1327`) and `resolveWorkflowOmicsOwnedSession` (`:1402`).

**S3B-05 (M) Omitted tests.**
- Doubles: `test-lipid-00b` `:2720`, `test-lipid-00c` `:263`, `test-metab-00b` `:228`; `test-lipid-01e` `:198`, `test-lipid-01i` `:89`, `test-metab-01r` `:352,594,641`, `test-metab-01v` `:14`, `test-metab-qc-itsd-direct-shared` `:191`; `test-metab-05` `:428`; `test-metab-03c` `:39`; `test-module-ci-lipid-design.R:116,179,269`; `test-omics-audit-carry-overs.R:23`.
- The import module-CI lanes: `helper-module-ci-{lipid,metab}-import.R:21`.
- Tests asserting the pass-through that item 2 reverses (a recorded contract change, not in `files_to_modify`): `test-workflow-lipid-declarations.R:56-63`, `test-workflow-metab-declarations.R:57`, `test-workflow-omics-sessions.R:215-216,238-240`.

### KIT-073 (shared core)

**S3B-06 (M) The binding table cannot express the differences.**
- **Duplicate resolution:** lipid does a per-assay `lapply` over `resolveLipidDuplicateFeaturesByIntensityFn(assay_tibble, id_col, sample_cols)` with counts from `nrow`; metab makes one `resolveMetabDuplicateAssayData(assayList, metaboliteIdCol, logWarnFn)` call returning `statsList`. This needs a per-omic kernel adaptor, not a name and a result key.
- **Kernel names:** the call-site names are the local `…Fn` aliases (`resolveLipidDuplicateFeaturesByIntensityFn` `:62`, `runLipidPerAssayRuvOptimizationFn` `:152`, `runPerAssayRuvOptimizationFn` metab `:124`); only metab's intensity filter calls the bare name. The list also omits `extractLipidCtrlPerAssay` / `extractCtrlPerAssay` and `runLipidsDA` / `runMetabolitesDA`.
- **Parameter names:** the declared names differ (`lipids_*` / `metabolites_*`, `min_lipid_intensity_threshold` / `min_metabolite_intensity_threshold`).
- **Strings and separators:** product-type separators are mixed (`.`/`-`); checkpoint labels, inner DA schema (`*.da_references.v1`), design input keys and abort texts differ.
- **`@args` sections:** these live in `utils_workflow_operations.R:59-67`, not in the declaration files.
- **Replay branches:** the record has two pairs, the step replay (`:569-573`) and the DA replay (`:593-600`); the ticket names one.
- **Labels:** `lipidomicsStepLabels` / `metabolomicsStepLabels` are not in the alias list.
- **Range:** `:266-379` mixes the core (`:266`, `:353`, `:368`) with lipid wrappers and aliases.

### KIT-074 (neutral names)

**S3B-07 (M) The example list contradicts the inventory rule.**
- **Out of scope:** 14 of the 20 examples are called only by proteomics or kit internals, which the ticket's own rule excludes: `readDIAStateGuiMetadata`, `getDIAOwnedPayload`, `setDIAOwnedPayload`, `getDIAReactiveReference`, `resolveDIAOwnedValue`, `attachDIAStateSession`, `boundDIAProvenanceValues`, `takeDIAStepProvenance`, and the five `*DiaQc*` functions.
- **Already neutral:** `saveDIAOwnedValue` is already an alias of `saveWorkflowOwnedValue` (`utils_workflow_dia_state.R:1776`).
- **Outside the module globs:** `createDIASessionOwner` is called only from `mod_lipidomics.R:194` / `mod_metabolomics.R:194`.
- **Omitted:** `collectDIANativeReferences` (`utils_workflow_dia_sessions.R:284`), `publishDIAStateHead` and `withDIAStateScope` (both already neutral-aliased).

KIT-075, 076 and 077 will newly call some proteomics-only functions, and those then need neutral names. *Repair:* keep the rule, regenerate the examples from the tree, and add a "names a later stage-3b ticket starts calling" clause.

### KIT-075 (references)

**S3B-08 (H) Three premises of the mechanism are wrong.**
1. **`hydrate`:** changing the core default `hydrate = TRUE` (`utils_workflow_lipid_declarations.R:267`) does nothing.
   - Every wrapper declares its own `hydrate = TRUE` (lipid `:287,296,302,310,316,325,337`; metab `:242,251,257,265,271,280,292`).
   - Modules read `$object`: `mod_lipid_qc_intensity_server_helpers.R:163`, `mod_lipid_qc_duplicates_server_helpers.R:257`, `mod_lipid_norm_workflow_helpers.R:101,125,146,198,737`, `mod_metab_qc_intensity.R:410`, `mod_metab_qc_duplicates_server_helpers.R:416`, `mod_metab_norm_pipeline_helpers.R:341,417,502,604`, `mod_metab_norm_session_helpers.R:53`.
   - The DA wrapper re-attaches `theObject` at lipid `:437` / metab `:355`.
2. **`revertToState`:** returning a name breaks tests that use the returned object: `test-prot-02h-peptide-audit-shared.R:359`, `test-workflow-checkpoint-fixtures.R:90`, `test-workflow-checkpoints.R:311,337`, `test-workflow-saved-identity.R:187`. The in-memory `WorkflowState$revertToState` contract (`utils_workflow_state.R:136`) also returns the data. Proteomics GUI has no caller; it uses `selectDIAStateReference`.
3. **Owner lookup:** the gap is worse than "pass through". With proteomics selected in the same session, `savedReactiveVal` / `setDIAOwnedPayload` would save lipid/metab values into the **proteomics** owner's store (`getDIAActiveManager`, `:1613`). `getDIAOwnedInput` (`:1802`) is also implicit.

*Repair:* item 2 lists the wrappers and module `$object` readers. Item 6 either keeps `revertToState`'s return and adds a non-materialising variant for the GUI callers, or records the test-contract change. Item 1 adds `getDIAOwnedInput` and the cross-omic misroute test. *Blocks:* implementation.

**S3B-09 (H) The residency walk names `qcData`, which does not exist.** The real holders are:
- `workflow_data` (`mod_lipidomics.R:139`, `mod_metabolomics.R:139`);
- `normData`;
- `daData` (`mod_lipid_da.R:38`, `mod_metab_da_server_helpers.R:75`);
- the import `localData` (`mod_lipid_import_processing_helpers.R:236`, `mod_metab_import_server_helpers.R:298`);
- the summary state (`mod_lipid_summary_observer_helpers.R:11`, `mod_metab_summary_server_helpers.R:53`).

Also: `mod_lipidomics.R:155` / `mod_metabolomics.R:155` are dead slots, never assigned. The receipt probe is at `:1932-1937` (not `:1934-1939`); it counts only `is.data.frame || isS4` over `data_tbl`, `data_cln` and `ruv_normalised_for_da_analysis_obj`, and not `design_matrix`. *Blocks:* the acceptance test as written.

**S3B-10 (M) Omitted whole-object reads.**
- `getState` sites: `mod_lipid_da_analysis_helpers.R:34,37`, `mod_metab_da_observer_helpers.R:303`, `mod_lipid_norm_workflow_helpers.R:40,355,375`, `mod_lipid_norm_observer_helpers.R:260,323`, `mod_metab_norm_binding_helpers.R:104,124`, `mod_metab_norm_render_helpers.R:176`, `mod_metab_norm_pipeline_helpers.R:72`, `mod_metab_norm_session_helpers.R:528`, `mod_metab_qc_duplicates_server_helpers.R:335`.
- Integration RDS writes: `mod_lipid_summary_state_helpers.R:319`, `mod_metab_summary_observer_helpers.R:265`.
- `sessionData$current_s4_object` readers: `mod_lipid_da_session_helpers.R:55,112,304`, `mod_metab_da_observer_helpers.R:743`.
- The metab revert caller `:131,141`.

**S3B-11 (M) No terminal outcome for memory.** The ticket expects the lifecycle peak at or below the in-memory build, but says nothing about what happens if it is not. *Repair:* declare the observable before the run (the normalisation-step anonymous peak against kit_032 and against the in-memory build), with outcomes "reported; KIT-080 judges", plus a stop-and-report if it exceeds kit_032.

### KIT-076 (DA per contrast)

**S3B-12 (H) "Identical" holds bit-for-bit only for orthogonal designs.**
- **What is contrast-independent** (confirmed against limma 3.68.4): `lmFit`, `squeezeVar` with trend/robust on sigma/df/Amean, `treat`, per-contrast `qvalue` and BH. Nothing reads the all-contrast F or `F.p.value`; Glimma reads `fit.eb` only for coefficient names.
- **Empirically:**
  - `~ 0 + group`: single-contrast and all-contrast `topTable` / `topTreat` are `identical()`.
  - `~ 0 + group + batch` or a covariate: differences of about 4e-16 in t and 3e-16 in p. The cause is `contrasts.fit` dropping all-zero rows and taking `chol` of a different sub-matrix in the non-orthogonal path, plus dgemv vs dgemm.
- **Consequence:** "any difference is a stop" would stop on floating-point noise for any batch formula.

*Repair:* fix the rule before outcomes: `identical()` for orthogonal designs; for others, an absolute ≤ 1e-12 on t/p/q, identical significance calls and identical row sets, recorded as a deviation. *Blocks:* decision.

**S3B-13 (M) Undeclared behaviour changes.**
- **Failure isolation:** a failing contrast now fails alone, where before the whole assay was dropped silently. Declare this as a behaviour change. The metab twin `tryCatch` (`func_metab_da_orchestration.R:248-252`) is not cited.
- **Ordering and counts:** `significant_counts` sums all contrasts per assay (`:461-470`), and the long table is assay-major. The v2 assembly must reproduce both.
- **Product hash:** `fit.eb` sits inside the `da_results` product, so the product hash changes.

**S3B-14 (M) Omitted consumers.**
- Static volcanoes: `func_lipid_da_plotting.R:673`, `func_metab_da_volcano_static.R:32-50`.
- Exports calling the heatmap on `theObject`: `func_lipid_da_export.R:237`, metab `:239`.
- Download handlers: `mod_lipid_da.R:256`, `mod_metab_da_registration_helpers.R:126`.
- Selectors and tables: `mod_lipid_da_analysis_helpers.R:263-287`, `mod_metab_da_observer_helpers.R:273,675`, `mod_lipid_da_results_helpers.R:52,73,179`, `mod_metab_da_display_helpers.R:120,189,341`, `mod_lipid_da_display_helpers.R:7`.
- The v1 manifest builders: lipid `:206-224`, metab `:188-207`.
- Test fakes of `fit.eb`: `test-lipid-02-da-core-helpers.R:142,168`, `helper-module-ci-{lipid,metab}-da.R`.

No stage-5 code reads lipid/metab DA today. The templates read only exported files and the record.

**S3B-15 (L) Misattributions.**
- `applyProtDaContrastModel` is defined at `func_prot_da_model_stats.R:282-313`, with the fit half `fitProtDaBaseModel` at `:227`. `:480-499` is the call site; the extraction was commit `117c4d23`.
- `newWorkflowDiaBaseFitIdentity` (`utils_workflow_dia_da.R:704`) is proteomics-only, so KIT-074 gives it no neutral name. Its identity also hashes block, weights and the normalisation/RUV/software ids, so an omics fit identity is new code.

### KIT-077 (figures)

**S3B-16 (H) The proteomics figure products are protein-bound, not "extended by factors".**
- `newWorkflowDiaQcProductIdentity` aborts unless the object is `ProteinQuantitativeData` (`utils_workflow_dia_products.R:11`) and reads protein slots (`:17-20`).
- `saveWorkflowDiaQcProduct` requires schema `prot_norm_qc_plot_product_v1` (`:48`) and hard-codes `route = "iq"` (`:55,58,137`).
- `rebuildWorkflowDiaQcPlot` supports four types (`:71`).
- The filtering-QC and composite figures are `arrangeGrob` grids with no plot-data rebuild path.

*Repair:* an omics product identity and schema, plus a decision on how grids are stored (rendered asset only, no plot-data product). Declare it as new code, and declare that proteomics identities stay byte-identical.

**S3B-17 (M) The defects are real, but stated imprecisely.**
- **Metab:** the composite (`metabolomics_composite_QC_figure.*`, written by `savePlot` into `feature_qc_dir`) never reaches `QC_figures`. `:412` queues `composite_QC_figure.png`, which does not exist, and no other pattern matches. The report reads `qc_figures_dir/metabolomics_composite_QC_figure.png` (`metabolomics_report.rmd:69,1132`), so the composite is missing from the report.
- **Lipid:** the composite is copied by the unanchored pattern at `:533`. The inverted fallback (`:537-542`) only queues a missing file. The repair is to delete it; flipping the test would queue a duplicate.
- **Fixtures hide the metab defect:** `helper-module-ci-{metab,lipid}-summary.R:168-169,193-194` seed the unprefixed name, so the reproduction must not use them.

**S3B-18 (M) Omitted writers.**
- All-contrast multi-page PDFs: lipid `func_lipid_da_export.R:330,351`; metab `func_metab_da_results_output.R:332,353`.
- The volcano and heatmap PDFs beside the PNGs: `:218,271` / `:220,273`.
- The three composite files (pdf, png, rds).
- `updateMetaboliteFiltering` (`func_metab_qc.R:94`) as the metab filtering writer.

The Glimma interactive volcanos are never saved, unlike proteomics' `saveDIAInteractiveVolcanoAsset`. Decide whether they are in scope.

### KIT-078 (worker)

**S3B-19 (B) The decision metric cannot measure the objective.** The managed slice holds the controller, the worker and the monitor together (`utils_workflow_execution.R:259`), and the protocol's journey metrics are whole-slice. Moving the kernels into the worker therefore does not lower the whole-slice peak: a worker adds its own process baseline, about 0.85 GiB of packages, while the app keeps its own. As a result:
- "the app process holds no matrix" has no observable;
- KIT-080's per-step guard of ≤ 1.02 × control can fail on the normalisation step by design, and a valid unfavourable result would then be recorded against the stage's main mechanism.

*Repair (user decision):* bind a per-process observable before any run, for example the controller's anonymous peak per step from the monitor's per-unit trace if it exists, or a hydration probe proving the controller never materialises an assay during a managed step. Decide whether KIT-080 judges the worker on whole-slice not-worse, controller-process reduction, or both. If a new metric is needed, that is a protocol v1.6 amendment. *Blocks:* decision (KIT-078 acceptance, KIT-080).

**S3B-20 (H) The proteomics gate does not exercise the chain it changes.** The KIT-028 replay and the LFQ/TMT/DIA_LIMPA oracle records run in process, while KIT-078 rewrites the managed request/poll/accept/complete/publish chain that only managed proteomics sessions use. *Repair:* add a `gui`-mode managed proteomics member (the KIT-008 kind: DIA screen journey, named worker dispatches observed, populations reproduced) after the chain items and at closure.

**S3B-21 (H) The mechanism conflicts with existing contracts.**
- **Sections:** "Every worker step mapped to a config section" is impossible for log2 and between-sample (`utils_workflow_operations.R:60-62`: they write scalars, with no `@args` section), and for duplicates, Pearson and DA. The registry test `:226-245` needs a different rule.
- **Deadlines:** the registry test `:277-281` pins `worker_deadline` to NULL for every worker step; setting one is a pin change.
- **Inputs:** `validateWorkflowActionInput` (`:2824-2826`) requires `input_references == list(state=)`, which conflicts with sending ITSD ids, contrasts, pairs and the DA fit as references.
- **Results:** `acceptWorkflowActionResult` (`:1732-1733`) requires a component checkpoint, so the product-only declarations (`*_pearson_pairs`, `*_da_output`) would end `invalid_result`.
- **Export readiness:** `assertDIAExportReady` (`utils_workflow_saved_exports.R:91-119`) has a proteomics state list and gates retry on `publishDIAActionOutputs`.
- **Peptide-shaped completion:** `newWorkflowCompletionWorkflow` (`:547-554`), `validateWorkflowCompletionContext` (`:536-542`), `workflowPeptideProgressProduct` (`:559-563`), and the `readManagedPeptideResult` default (`:602`).

**S3B-22 (L) References and wording.**
- The registry test `:80-82` needs no change: `lipid` stays invalid.
- `:793` is a class branch, not a message, and needs no change.
- The driver's tab actions are at `workflow_kit_driver.R:536-603`, not `:499-538`.
- The omics members already launch under a managed application (`workflow_kit_coordinator.R:281`, `mode = "gui"`). The gap is dispatch, since no omics module calls `installWorkflowManagedStep`, so item 6 should be reworded.

### KIT-079 (partitions)

**S3B-23 (B) ITSD cannot be partitioned exactly without changing the kernel.** The ITSD kernel computes `average_norm_factor <- mean(...)` inside its body (`func_lipid_s4_normalization_methods.R:298`; metab `func_metab_s4_transformation_methods.R:453`) and takes no parameter for it. Run per sample batch, the unchanged kernel uses a batch-local mean, and a post-hoc rescale is unlikely to be hash-exact. That conflicts with both "kernels are never rewritten" and "exact against whole". The intensity filter is fine: `lipidIntensityFilteringHelper(min_lipid_intensity_threshold)` (`func_lipid_qc.R:40`) and `metaboliteIntensityFilteringHelper` (`func_metab_qc_filtering_helpers.R:13`) take the threshold. *Repair (user decision):* drop ITSD from KIT-079 (it stays whole_matrix in practice), or authorise an additive kernel parameter that accepts the pre-pass value, with the default behaviour byte-identical. *Blocks:* implementation.

**S3B-24 (M) The merge rules are under-specified.**
- **Axis:** "concatenate" names no axis. A complete-sample batch of a wide assay needs a column-bind (ids once, samples in design order); row and feature-group batches need a row-bind in source order.
- **Counts:** kernel `counts` are absent from the merges. `lipidomicsFeatureCounts` counts unique ids (`:28-32`) and duplicates still exist at the intensity filter, so per-batch sums would over-count. Counts must be recomputed after the merge.
- **Projection:** `partitionBySamples` is I/O-free, and the only wide projection applies one column set to every element (S3B-25).

*Repair:* state the axis and the count rule per step.

**S3B-25 (L)** `test-workflow-saved-table-codec*.R` does not exist. The codec and cursor tests are `test-workflow-saved-projection.R`, `test-perf-103.R` and `test-perf-103-audit-repairs.R`. Gaps 1–3 are confirmed real (`saved_table_codec.R:1866-1867`; `dia_state.R:440,467`, with `materializeWorkflowTable` aborting on a column a list element lacks, `:299-300`).

### KIT-080 (measurement)

**S3B-26 (B) The operation pairs cannot run, and the ticket forbids the fix.**
- **List assays:** `kitPopulationInputComponent` stops unless the assay component is a single entity id (`tools/profiling/workflow_kit.R:1159-1161`), and `verifyKitOperationMember` (`:1521-1535`) does the same. Both fail for a lists assay.
- **Pair slots:** `kitCases` allows at most 2 pairs per ticket key (`:356`), and operation stages take no suffix (`^kit-0[0-9]{2}$`; selector `:473`; protocol §4/§7.3). That is 4 operation members per ticket, not 16.
- **Ceiling:** operation members get a fixed 600 s (`:113`), not 1,800 s.
- **Decision rule:** `kitOperationDecision`'s `accepted` means not-worse (≤ 1.02); `reduction` (≤ 0.95 on `peak_bytes`) is reported, never required (`:392-399`). Item 4 would therefore switch defaults with no memory gain.
- **Scope:** KIT-080 lists no harness file, and its "one mechanism" allows only the default-variant change.

*Repair:* move the operation-member harness support into KIT-079 (or a tooling ticket before KIT-080): lists-assay input components, more pair slots or several keys, the operation ceiling from the bindings. Key the default switch on `reduction`. *Blocks:* execution.

**S3B-27 (B) The control and the corpus are outside protocol v1.5.**
- **Control:** §2 says "the control is the latest accepted KIT build" (`kit_036` in the bindings), and kit_032 has no outcomes row.
- **Corpus:** the ticket-journey-pair kind is defined "on the 300-sample corpus" with "the candidate's states verified against the KIT-005 oracles". §7.1's in-process oracle extends the oracle source but not the pair's corpus.

*Repair (user decision):* a protocol v1.6 note, or a recorded decision, naming kit_032 as the stage-3b control, the real lipid/metab corpora as pair corpora, and the oracle per member under §7.1. *Blocks:* execution.

**S3B-28 (M) Precondition availability.** The real lipid/metab corpora do not exist yet; KIT-066 says "the user supplies" them. The KC privacy rule permits aggregates and hashes only, and the ticket honours it. KIT-080 is blocked until the corpora exist. Record this as its blocked state, not as a defect.

**S3B-29 (L) Scheduling consequence.** KIT-066's new dependency on KIT-079 also delays its proteomics enrichment proof (the KC LFQ lane) and, through KIT-060, stage 5. The decision row records the dependency but not this consequence. *Repair:* one sentence in the decision row, or split KIT-066's proteomics item off if the user prefers.

## Evidence layers

| Layer | Established | Not established |
|---|---|---|
| `documentation_or_graph_integrity` | graph, statuses, sections, file and line existence, identifier existence | — |
| `implementation_conformance` (of the premises) | the cited code says what is listed above; limma per-contrast behaviour checked against source and empirically (reference BLAS) | behaviour under the optimised BLAS of the proof host |
| `runtime_behavior`, utility | nothing (no run was made) | every memory claim |

## Closeout

- `implementation_verified`: false (nothing implemented).
- `execution_complete`: true for the audit.
- `decision_supported`: false: S3B-19, 23, 26 and 27 need user decisions.
- `promotion_authorized`: false.
- `claim_scope`: the stage-3b ticket texts at `a7ac9d06`.

**Next bounded action.** The user decides S3B-19 (worker metric), S3B-23 (ITSD partition), S3B-26 (where the operation-member harness work goes) and S3B-27 (control and corpus). Then one revision pass of the tickets applies the reference, list and rule repairs, recorded in `STAGE3B-REMEDIATION-<date>.md`. KIT-072 can start after its own repairs (S3B-01, 03, 04, 05), which are independent of the four decisions.
