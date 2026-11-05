# KIT-096 audit: Transcriptomics GUI, part 1 (wiring, import, design, QC)

- **Mode:** `audit_only` (audit-evidence-contracts). No ticket, code or doc was edited; no test, install, journey or app run.
- **Target:** `tickets/kit/KIT-096.md` revision 1.0 (sha256 prefix `02dc7b47c3d2cd1c`), at HEAD `bbe1fae9`, branch `performance-v2`.
- **Read for context:** KIT-040, KIT-093, KIT-094, KIT-095, KIT-097, KIT-098, KIT-099; `00-KIT-PLAN.md` v3.8; the approved plan `~/.claude/plans/write-a-series-of-atomic-dragon.md`; `AGENTS.md`; `conventions/R.md`, `R-shiny.md`, `R-golem.md`, `R-bioinformatics.md`; the notebook `dev/distributable/plans/transcriptomics_workflow.rmd`.
- **Mechanical checks:** `next-openai-ticket.sh KIT-096` valid, no errors or warnings, status consistent (`pending`), `needs_planning: true`. `audit_inventory.py` exit 0: graph of 97 tickets, no cycles, no duplicates, no missing dependencies. Its one link finding (`00-KIT-PLAN.md` -> `omic`) is in the plan, not this ticket. The index title is the ticket title without "on the managed runner" (cosmetic).
- **Ticket state:** the ticket is planning-gated. A file-by-file plan is written when it starts. The findings below are what that plan has to settle. Some of them need a user decision first.

## Finding register

Severity: **blocking** means the ticket cannot be planned as written without a decision. **should-fix** means the ticket can be planned, but the gap would produce an unverifiable or wrong acceptance. **minor** means precision only.

### KIT096-F01: "the shared design builder" does not exist; reuse vs "other omics' modules unchanged" conflict

- **Severity / dimension:** blocking (needs a user decision) / implementation best practice, reference accuracy.
- **Evidence:**
  - There are three per-omic design builders and no shared one:
    - `R/mod_lipid_design_builder_{server,ui,helpers,state_helpers,action_helpers,display_helpers}.R`;
    - `R/mod_metab_design_builder_{server,ui}.R` (866 + 284 lines);
    - `R/mod_prot_design_builder_*.R`.
  - Nothing in `R/mod_metab_*` calls a lipid function (grep for `*Lipid*(` in `R/mod_metab_*.R`: none). The omics modules are copied per omic today. `R/mod_lipidomics.R` and `R/mod_metabolomics.R` are line-for-line twins (414 lines each, differing only in names).
  - The lipid builder is lipid-specific:
    - its UI text says "Lipidomics" (`R/mod_lipid_design_builder_ui.R:218`, `:246`);
    - its save path keeps `colMap$lipid_id_col` (`R/mod_lipid_design_builder_state_helpers.R:497-507`);
    - sample detection excludes MS-DIAL column patterns (`R/mod_lipid_design_validation.R:17-85`);
    - the results observer builds `createLipidomicsAssayData` and binds `"LIPIDOMICS"` (`R/mod_lipid_design_builder_helpers.R:29-160`, `:199-221`).
  - What the omics actually share in design is an input-id contract, not code. `tools/profiling/workflow_kit_driver.R` `runKitDesignActions` says "the lipidomics and metabolomics design builders share the proteomics builder's input ids under their own prefix".
  - The ticket says "design (the shared design builder reused)" and "shared UI helpers are reused, not copied", and also "Other omics' modules unchanged".
- **Impact:** "Reuse" has three readings, and each breaks a clause of the ticket:
  - (a) calling `mod_lipid_design_builder_*` unchanged, with `lipid_id_col` holding the gene id: the transcriptomics tab shows "Lipidomics" text, and the gene-id column is carried under a lipid name;
  - (b) parameterising the lipid builder (label, id-column name, class constructor): this edits a lipid module, against "Other omics' modules unchanged" unless "unchanged" is defined as behaviour plus lipid module-CI unchanged;
  - (c) copying it as `mod_trx_design_builder_*` (about 2-3k lines): against "reused, not copied".

  The planning step cannot pick one of these without a scope decision (AGENTS.md: "Request a scope decision before materially expanding the work").
- **Smallest repair:**
  - Replace "the shared design builder" with the concrete choice.
  - Recommended: (b), with the lipid builder's defaults unchanged and "unchanged" defined as the lipid and metab module-CI and the KIT-094 frozen dump unchanged.
  - Separately, list the genuinely shared parts to reuse as they are (see Verdicts, dimension 3).
- **Closure evidence:** the ticket revision names the option and the user's approval. The lipid and metab module-CI lanes pass unchanged after implementation.
- **Blocked state:** implementation (planning step).

### KIT096-F02: the module-CI infrastructure files that must change are outside `files_to_modify`

- **Severity / dimension:** should-fix / completeness.
- **Evidence:**
  - `tests/testthat/helper-module-ci-manifest.R:1`: `.MODULE_CI_VALID_OMICS <- c("all", "proteomics", "metabolomics", "lipidomics")`. Any transcriptomics scenario is refused at `:113-114`.
  - `tests/testthat/helper-module-ci-fixtures.R:17-23`: `module_ci_expected_import_formats()` has no transcriptomics entry, and packs are validated against it (`:102-113`, `:154-155`).
  - Scenario ids must match `^MCI-[0-9]{3}...` and `ticket_id` must match `^MCI-[0-9]{3}$` (`helper-module-ci-manifest.R` validator). The new scenarios therefore need new MCI ids (the last is MCI-025), not KIT-096.
  - `ci_lanes` in `tests/testdata/mci/manifest.json` is `foundation, proteomics-unit, metabolomics-unit, lipidomics-unit, browser-smoke, nightly, release`. There is no transcriptomics lane.
  - A fixture must exist at `tests/testdata/mci/<fixture.path>`. There is no `fx/transcriptomics/`, and `fixture_packs.json` and `generate-fixtures.R` are not listed.
  - `tests/testthat/test-module-ci-manifest.R:36-46` asserts the exact set of `import` x `unit-contract` scenarios. Adding a transcriptomics import scenario breaks it.
  - `files_to_modify` lists only `helper-module-ci-trx-*.R`, `test-module-ci-trx-*.R` and `tests/testdata/mci/manifest.json`.
- **Impact:** the test-first step cannot be written inside the declared file scope. The first scenario added breaks an existing test that is not in scope.
- **Smallest repair:**
  - Add to `files_to_modify`:
    - `helper-module-ci-manifest.R`, `helper-module-ci-fixtures.R`, `test-module-ci-manifest.R`;
    - `tests/testdata/mci/fixture_packs.json`, `tests/testdata/mci/fx/transcriptomics/`, `tests/testdata/mci/generate-fixtures.R`.
  - Name the MCI ids (MCI-026 onward) and the lane.
  - State that the fixtures come from the KIT-040 generator, not new data.
- **Closure evidence:** `test-module-ci-manifest.R`, `test-module-ci-fixture-integrity.R` and the new `test-module-ci-trx-*.R` pass. The lipid and metab module-CI scenario sets are unchanged.
- **Blocked state:** implementation.

### KIT096-F03: "No whole matrix held in the application" names no test, and its scope (held vs read; import and design) is undefined

- **Severity / dimension:** should-fix / completeness, scientific rigor (of the acceptance).
- **Evidence:**
  - The KIT-075 residency walk defaults to `classes = c("LipidomicsAssayData", "MetaboliteAssayData")` (`tests/testthat/helper-workflow-omics-residency.R:5`).
  - The KIT-082 probe is hard-coded to the two omics (`tests/testthat/test-workflow-omics-application-reads.R:5-21`).
  - Neither file, nor `test-workflow-omics-residency.R`, is in `files_to_modify`.
  - The lipid design step reads the whole tables in process for one action: `resolveWorkflowOwnedValue(held_data_cln, ...)`, then `createLipidomicsAssayData`, then `saveState` (`R/mod_lipid_design_builder_helpers.R:51-53`, `:114-145`). The import reads the file in process and then holds it by reference (`applyLipidImportResultToWorkflow`, `R/mod_lipid_import_processing_helpers.R:249-277`).
  - KIT-095 declares no worker import or design step (only the six compute steps, plus QC finalise and reset as store steps).
  - KIT-075's rule is "held" (no object kept on the bus). KIT-082's rule is "reads no assay outside the managed steps", and it covered QC, normalisation, DA and summary, not import or design.
- **Impact:** the acceptance cannot be verified as written:
  - Is an in-process read at import and design allowed (as in lipid and metab), or a breach?
  - Which probe proves "held"?

  An implementer could pass it with no test at all, or could expand scope by moving import and design onto the worker.
- **Smallest repair:** reword the acceptance as follows.
  - The KIT-075 residency walk, run with `TranscriptomicsAssayData` and the corpus sample columns, finds no object or assay table on the bus after import, design and QC, on both the managed and the in-process paths.
  - The KIT-082 probe is extended to the transcriptomics QC tab.
  - Import and design read the file or tables in process for one action, as lipid and metab do (stated as the accepted exception).
  - Add the three test files to `files_to_modify`.
- **Closure evidence:** the residency and probe tests fail on a deliberately holding variant (the positive control) and pass on the tabs.
- **Blocked state:** implementation (acceptance definition).

### KIT096-F04: what the import and design steps commit is not handed over by KIT-095 ("committing the KIT-095 states")

- **Severity / dimension:** should-fix / completeness (producer to consumer).
- **Evidence:**
  - KIT-095 lists the worker steps `trx_expression_filter` … `trx_ruviii_export` and the store steps "QC finalise and reset". It declares no initial or design state, no state name, and no input roles.
  - The lipid design step commits `lipid_raw_data_s4` through `saveState` after `bindWorkflowOwnerTypeAtDesign(…, "LIPIDOMICS", …)` and `recordLipidDesignProvenance` (`R/mod_lipid_design_builder_helpers.R:136-145`).
  - That provenance helper is the generic `recordWorkflowOmicsDesignProvenance(binding, …)` (`R/utils_workflow_omics_declarations.R:1011-1026`). It records each file in `source_dir/import_sources/` under `binding$design_input_prefix<n>`, plus `design_file`.
  - The owner refuses any workflow type outside `c("DIA","DIA_LIMPA","LFQ","TMT","LIPIDOMICS","METABOLOMICS")` (`R/utils_workflow_dia_state.R:1407-1414`; also `:1330`, and `workflowRoutesByType`, `R/utils_workflow_dia_sessions.R:817-818`). KIT-094 moves these lists to binding fields, and KIT-095 must add `TRANSCRIPTOMICS`.
  - The replay input check accepts only `^(lipid|metab)_assay` roles (`R/utils_workflow_analysis_record.R:465-466`, also KIT-094).
- **Impact:**
  - Import and design commit no "KIT-095 state", so the first acceptance criterion is partly unsatisfiable as worded.
  - If the GUI's design state (name, provenance, roles, type) differs from the one KIT-095's in-process journey creates, KIT-098's GUI journey is not comparable with the recorder. Its record may not validate or replay.
  - A separate annotation file copied to `import_sources/` would be recorded as a second assay input.
- **Smallest repair:**
  - KIT-096 states that the design step commits the initial `TranscriptomicsAssayData` under the state name, provenance (`recordWorkflowOmicsDesignProvenance(<trx binding>, …)`), workflow type `TRANSCRIPTOMICS` and input roles that KIT-095's journey uses.
  - KIT-095's hand-off names them.
  - The acceptance becomes: "import, design and QC commit the design state and the KIT-095 QC states (expression filter, QC finalise), with a record that validates".
- **Closure evidence:**
  - a test that the GUI design commit and KIT-095's journey setup give identical state metadata and provenance;
  - `validateWorkflowAnalysisRecord` valid after the QC tab.
- **Blocked state:** implementation.

### KIT096-F05: the import contract is undefined (formats, annotation, count policy surfaced)

- **Severity / dimension:** should-fix / completeness, scientific rigor.
- **Evidence:**
  - The ticket says only "import (counts and annotation)".
  - The notebook reads one file with columns `GeneID, GeneName, <samples>` and derives the annotation from it (`transcriptomics_workflow.rmd:116-152`). The KIT-040 corpus produces only that file.
  - The KIT-093 class has `annotation_id_column` (a column), no annotation table.
  - Common RNA-seq count tables are not addressed:
    - featureCounts adds `Chr/Start/End/Strand/Length` columns and a comment line;
    - salmon/RSEM "counts" are fractional estimates.
  - KIT-093's whole-number policy refuses fractional counts. The notebook rounds silently, but only at RUVs (`as.matrix(round(count_matrix_filtered))`, line 406). Its filter runs on unrounded values.
- **Impact:**
  - The scope of the import tab and its tests cannot be planned.
  - "Annotation" may grow into an unplanned file join.
  - A user with salmon output gets a refusal deep in a later step, unless the import surfaces the policy.
  - The ticket does not say that the GUI must not round silently.
- **Smallest repair:**
  - Define the formats: at minimum the notebook format, delimited, with the id and name columns chosen by the user. Any other format is listed or explicitly excluded.
  - Define annotation as the name column of the counts file (or, if a separate file, its join key and input role, routed to KIT-095).
  - State that the import applies the KIT-093 count policy and shows its classed refusal (NA, negative, fractional) by name, with no rounding option.
  - Tests use the KIT-040 invalid fixtures.
- **Closure evidence:** module-CI import scenarios accept the screen corpus and refuse each invalid fixture with the KIT-093 condition class shown to the user.
- **Blocked state:** implementation.

### KIT096-F06: reused lipid import and design defaults would change sample names and order relative to the notebook reference

- **Severity / dimension:** should-fix / scientific rigor.
- **Evidence:**
  - The lipid import sanitises sample names by default. The "Sanitize Sample Names" checkbox defaults to `value = TRUE` (`R/mod_lipid_import_ui_helpers.R:256-265`) and applies `janitor::make_clean_names` (`R/mod_lipid_import_processing_helpers.R:322-360`).
  - The lipid builder's fresh design orders runs by `gtools::mixedsort(sampleNames)` (`R/mod_lipid_design_builder_state_helpers.R:394-395`).
  - The notebook keeps the raw names and orders columns by the design file (`count_matrix[, design_matrix$Run]`, line 152).
  - KIT-098 requires "final states identical to the notebook reference".
- **Impact:** reusing these defaults gives logCPM, W and `newY` with renamed or reordered samples. The results are numerically the same science, but not `identical`. The exactness gate would then fail, or a tolerance would be introduced after the fact, which the contract forbids.
- **Smallest repair:**
  - State the sample contract: no name sanitisation for transcriptomics (or, if kept, the KIT-040 corpus names are already clean and this is asserted).
  - Sample order follows the design file or the counts file, as KIT-093's kernel and the KIT-040 reference define it.
  - Add a module-CI case: names and order after import and design are identical to the corpus design.
- **Closure evidence:** that test. Later, KIT-098's identity check.
- **Blocked state:** implementation, and KIT-098 execution.

### KIT096-F07: group level order is not captured by the design tab

- **Severity / dimension:** should-fix / scientific rigor.
- **Evidence:**
  - The notebook fixes `group = factor(group, levels = c("RPMI","Sera"))` (line 139).
  - It ranks control genes with `model.matrix(~ group)` and `glmQLFTest(fit_prelim, coef = 2)` (lines 232-236), and uses `X_bio <- design_prelim[, -1]` (line 241).
  - With two groups the p-values do not depend on the baseline. With more than two, `coef = 2` tests only one level against the baseline, so the ranking, the control set, k and W all depend on level order.
  - The design builders write `group` as character (`R/mod_lipid_design_builder_state_helpers.R:392-402`). No level order is recorded.
  - `filterByExpr(group = design$group)` also takes the group from the design.
- **Impact:** with three or more groups, the GUI's results depend on implicit alphabetical level order. Nothing records or controls it, and KIT-093's kernel has no stated rule for more than two groups.
- **Smallest repair:**
  - The design tab records the group level order (the reference level) in the design state.
  - The design refuses, or visibly warns about, more than two groups if KIT-093's kernel only supports the notebook's two-group ranking. That is a cross-ticket rule for KIT-093.
  - A module-CI case covers level order.
- **Closure evidence:** the design state carries the level order, and the kernel receives it unchanged (test).
- **Blocked state:** implementation (and a KIT-093 decision).

### KIT096-F08: the QC filter's "library sizes" pre-pass is TMM-scaled on the unfiltered matrix

- **Severity / dimension:** minor (routed: KIT-093 and KIT-095 own the kernel) / scientific rigor.
- **Evidence:**
  - The notebook calls `calcNormFactors(dge_raw, "TMM")` before `filterByExpr(dge_raw, group = …)` (lines 166-170).
  - `edgeR:::filterByExpr.DGEList` uses `lib.size <- y$samples$lib.size * y$samples$norm.factors` (read from the installed edgeR).
  - The notebook then drops zero-variance genes (lines 173-175).
  - KIT-093 says "library sizes as a declared global pre-pass".
- **Impact:** the QC tab's filter (and its partitioned variant) is exact only if the pre-pass is the effective (TMM-scaled) library size on the unfiltered matrix, followed by the zero-variance step. If the GUI passes raw column sums, it silently changes the mask.
- **Smallest repair:**
  - KIT-096 says that the QC tab calls the KIT-095 filter step with the design's group factor and no library-size argument of its own.
  - Name the zero-variance step as part of the filter result shown.
  - Route the effective-library-size wording to KIT-093.
- **Closure evidence:** the QC tab's mask is identical to the KIT-040 reference mask (module-CI or in-process).
- **Blocked state:** implementation.

### KIT096-F09: filtering metrics and QC figures are undefined for counts, and KIT-095 does not declare their worker product

- **Severity / dimension:** should-fix / completeness, producer to consumer.
- **Evidence:**
  - The lipid filtering metrics are intensity metrics: `n_lipids`, `missingness`, `sum_intensity_per_sample`, `cv_distribution`, `is_metrics` (`R/func_lipid_qc_filtering_helpers.R:53-110`). NA-missingness and internal-standard metrics have no meaning for counts.
  - Under a managed session the metrics come from the declaration (`binding$filtering_metrics`, `R/utils_workflow_omics_declarations.R:1119`, `:1150`, `:1205`), read by `readWorkflowOmicsFilteringMetrics` (`R/utils_workflow_omics_actions.R`) and drawn by `updateLipidFiltering` (`R/mod_lipid_qc_intensity_server_helpers.R:218-235`). KIT-095 does not list a filtering-metrics field for the transcriptomics binding.
  - "QC figures" is ambiguous. The notebook's figures are PCA and RLE for raw logCPM, TMM and RUVs (lines 492-498), and two of the three stages belong to KIT-097's normalisation tab.
- **Impact:** either the QC tab reads the matrix in the application to draw metrics (breaching F03), or it reuses intensity metrics that are wrong for counts. The split of figures between KIT-096 and KIT-097 is open.
- **Smallest repair:**
  - Define the QC tab's metrics: genes before, after `filterByExpr` and after zero-variance; library size and detected genes per sample.
  - Define its figures: the raw-logCPM PCA and RLE, and the filtering progress figure.
  - Require the metrics as the KIT-095 filter declaration's display product (`filtering_metrics`), and note that hand-off in KIT-095.
- **Closure evidence:** in the managed module-CI case the metrics come from the worker result. The figure renders from them with no assay read (the F03 probe).
- **Blocked state:** implementation.

### KIT096-F10: the technical-replicate column and replicate minimums that KIT-097 depends on are not stated

- **Severity / dimension:** minor / completeness (dependency correctness for KIT-097).
- **Evidence:**
  - KIT-097 offers RUVIII `M` from `technical_replicate`, and KIT-093's class has `technical_replicate_id`.
  - The builders derive `tech_rep_group = paste(group, replicates, sep = "_")` (`R/mod_lipid_design_builder_state_helpers.R:490-497`) and set `technical_replicate_id = "tech_rep_group"` (`R/mod_lipid_design_builder_helpers.R:153`). Technical replicates are assigned through `assign_tech_reps`.
  - RUVs `makeGroups(group)` (notebook line 410) and RUVIII `replicate.matrix` need at least two samples per group.
- **Impact:** KIT-097's alternative `M` can only be offered if KIT-096's design state carries a technical-replicate grouping. A design with a singleton group passes the builder but fails later in RUV.
- **Smallest repair:**
  - State that the design tab commits `technical_replicate_id` (the builder's `tech_rep_group`).
  - Apply KIT-093's validation at save, including at least two samples per group (or name where that is refused).
- **Closure evidence:** a module-CI design case with technical replicates. The singleton-group refusal is tested.
- **Blocked state:** implementation.

### KIT096-F11: KIT-098's driver depends on the builder's input ids and completion text; the ticket does not fix that contract

- **Severity / dimension:** minor / implementation best practice (consumer contract).
- **Evidence:** `tools/profiling/workflow_kit_driver.R` `runKitOmicsStep` reaches the design tab through:
  - the tab input `<route>_workflow-<route>_tabs`;
  - the prefix `<route>_workflow-design-builder-`;
  - the notification "Design matrix and contrasts saved successfully!";
  - `runKitDesignActions` with the ids `main_tabs`, `new_factor`, `add_factor`, `selected_runs`, `factor1_select`, `assign_metadata`.

  The lipid tabs carry `data-testid`s (`lipid-tab-import` …, `R/mod_lipidomics.R:63-117`).
- **Impact:** if the transcriptomics module departs from these ids, KIT-098 has to write a new design and import driver, which is unplanned work in a later ticket.
- **Smallest repair:**
  - State the ids: module id `transcriptomics_workflow`, tabs input `transcriptomics_tabs`, sub-module ids `import/design/qc`, and the builder id contract unchanged.
  - Name the `data-testid` prefix (for example `trx-tab-*`).
- **Closure evidence:** a module-CI or app-shell test that asserts the ids.
- **Blocked state:** none now (KIT-098 execution later).

### KIT096-F12: divergences from the convention suite are not recorded, and there is no size bound

- **Severity / dimension:** minor / implementation best practice.
- **Evidence:**
  - `conventions/R-shiny.md` requires bslib, `ExtendedTask` for long work, and `shinyFiles` for file selection.
  - The existing omics tabs use `shinydashboard`/`wellPanel`/`fluidRow` (`R/mod_lipid_qc.R`, `R/mod_lipidomics.R`), their own managed runner (`installWorkflowOmicsActions`), and a `fileInput` fallback in import.
  - `conventions/R.md:26-29` says to match the established idiom and flag the divergence.
  - The user asked to "match the stylistic aesthetics of other omics".
  - `R.md` asks for functions of about 50 lines or fewer. The lipid counterparts of these three tabs are 24 files, about 7.9k lines (import 2,215, design 3,280, QC 2,399).
- **Impact:** without a recorded divergence, a reviewer could demand bslib or `ExtendedTask` (a rewrite). Without a size bound, an 8 h estimate risks a large copied surface.
- **Smallest repair:**
  - One line recording the divergence: legacy shell and the managed runner kept, by user instruction and the convention scope rule.
  - A planned file and line budget in the planning step.
- **Closure evidence:** the plan and the implementation summary.
- **Blocked state:** none.

### KIT096-F13: reference accuracy, where the wiring is smaller than the ticket implies, and stale notes

- **Severity / dimension:** minor / reference accuracy.
- **Evidence (verified):**
  - `R/app_server.R` dispatches modules by name: `paste0("mod_", omic, "_ui")` (`:437`) and `paste0("mod_", omic, "_server")` (`:576`). The fallback "Module not yet implemented" is at `:475`, and the transcriptomics icon at `:407`.
  - The tile is already selectable (`R/app_ui.R:443-453`, with `toggleOmicSelection` at `:553`). Its "(finding m-2)" comment at `:443` becomes stale, and `:507` holds the hidden checkbox.
  - Defining `mod_transcriptomics_ui`/`_server` is therefore enough to wire the module. `R/app_server.R` probably needs no change.
  - `inst/workbooks/transcriptomics/` does not exist, so the watchdog logs a warning (`R/app_server.R:531-549`). This is harmless, but relevant to KIT-097's report template.
  - `setupDirectories` already has transcriptomics directories (`R/func_general_filemgmt_setupdirectories_helpers.R:25`, `:78-93`, `:216-219`).
  - The plan's "as `mod_lipidomics.R:180-196`" is accurate: bus at `:180`, owner at `:194`, runner at `:196`.
  - `omicsManagedCommit` is at `R/utils_workflow_omics_actions.R:172`.
  - `.MODULE_CI_VALID_OMICS` is at `helper-module-ci-manifest.R:1`.
  - `test-e2e-browser-harness.R:120` expects "Unsupported omic" for transcriptomics through `.E2E_VALID_OMICS` (`helper-e2e-browser.R:6`). That test does not break in KIT-096 (the e2e helper is untouched). It will break when KIT-097 adds its e2e test, and KIT-097's `files_to_modify` lists neither file.
  - The docs of `registerWorkflowOmicsBus`, `createWorkflowSessionOwner` and `installWorkflowOmicsActions` say "lipidomics or metabolomics", but the code is omic-generic.
- **Impact:** the planned edits to app_ui and app_server may be unnecessary churn. The e2e break is mis-scheduled into KIT-097 without file scope.
- **Smallest repair:**
  - Narrow the wiring item to "define `mod_transcriptomics_ui/_server`; remove the stale comment".
  - Drop `R/app_server.R`, unless the plan shows a need.
  - Note the e2e helper and test for KIT-097.
- **Closure evidence:** an app-shell test that selecting transcriptomics renders the module, not the fallback.
- **Blocked state:** none.

### KIT096-F14: session restore, adoption and help text are unassigned

- **Severity / dimension:** minor / completeness.
- **Evidence:**
  - Legacy adoption is lipid and metab only (`workflowOmicsSessionClass`, `R/utils_workflow_dia_sessions.R:1318`, `:1340`). This is correctly not needed for transcriptomics, which has no older build.
  - KIT-095 declares a session export, but neither KIT-096 nor KIT-097 names the GUI "Open" of a transcriptomics session.
  - The lipid tabs carry inline `helpText`. The ticket says nothing about help text, for example on the counts policy or the absence of a positional design.
- **Impact:** a gap in the GUI's completeness that no ticket owns.
- **Smallest repair:**
  - State: no legacy adoption for transcriptomics; native session Open in KIT-097.
  - The import and design help text states the whole-number policy and that the design is matched by sample name only.
- **Closure evidence:** the ticket text; a UI test of the help strings.
- **Blocked state:** none.

### KIT096-F15: no GUI-level test that the positional design fallback cannot return

- **Severity / dimension:** minor / scientific rigor.
- **Evidence:**
  - The notebook builds the design by position when the file is absent (`rep(c("RPMI","Sera"), each = 6)`, lines 125-134).
  - The reused builder has no positional path. A fresh design starts with NA groups (`R/mod_lipid_design_builder_state_helpers.R:391-402`), and an imported design is matched by `Run` name (`R/mod_lipid_design_import_helpers.R` around line 240, `intersect(imported_design$Run, all_cols)`).
  - The ticket repeats the contract but names no GUI test.
- **Impact:** low today. If a transcriptomics-specific import helper ever ports the notebook's load chunk, the fallback could return unseen.
- **Smallest repair:** module-CI design cases:
  - a design whose runs do not match the counts columns is refused, not mapped by position;
  - unassigned groups are refused at save;
  - the KIT-093 refusals (missing, extra, duplicate, unknown level, empty group) are shown by class.
- **Closure evidence:** those cases pass.
- **Blocked state:** implementation.

## Verdicts by dimension

1. **Completeness: not ready.**
   - The ticket names the three tabs and the wiring. It does not define:
     - the import formats or the annotation (F05);
     - the design state that is committed and its hand-off from KIT-095 (F04);
     - the filtering metrics and QC figures for counts (F09);
     - the technical-replicate output that KIT-097 needs (F10);
     - the module-CI infrastructure files the test-first step must touch (F02).
   - Its matrix-residency acceptance names no test (F03). Session Open and help text are unowned (F14).
   - Nothing breaks merely by removing the placeholder:
     - the app-shell tests stub their own modules;
     - the e2e "Unsupported omic" test breaks only when KIT-097 extends `.E2E_VALID_OMICS`.
2. **Scientific rigor: adequate in principle, under-specified in the GUI.**
   - The reused builder has no positional fallback (F15 asks only for a test).
   - Several GUI defaults would break `identical` against the notebook reference, or leave results controlled by implicit state: sample-name sanitisation and run order (F06), group level order (F07), and the effective library-size pre-pass the QC filter relies on (F08, routed).
   - Surfacing the counts policy at import, without silent rounding, is not stated (F05).
3. **Implementation best practice: blocked on what "reuse" means (F01).**
   - Genuinely shared today, to reuse as they are:
     - the runner: `installWorkflowOmicsActions`, `workflowOmicsManagedActions`, `omicsManagedCommit`/`omicsManagedProduct`, `workflowOmicsDisplay`;
     - owner and bus: `createWorkflowSessionOwner`, `registerWorkflowOmicsBus`, `bindWorkflowOmicsOwnerAtImport`/`bindWorkflowOwnerTypeAtDesign`, `setWorkflowOmicsBusTables`, `omicsBuilderTableVal`;
     - state reads: `readWorkflowOmicsMetadata`/`readWorkflowOmicsShell`, `recordWorkflowOmicsDesignProvenance`;
     - UI helpers: `render_workflow_stepper`, `testid`, `is_test_mode`.
   - The import, design and QC modules themselves are copied per omic (lipid and metab are twins), so "reuse" of those needs a decision between parameterising the lipid modules and copying them.
   - The divergences from the conventions should be recorded (F12). The input-id contract for KIT-098 should be fixed (F11). Test-first is stated, but its file scope is incomplete (F02).
4. **Codebase reference accuracy: mostly accurate, one false citation.**
   - Every cited file and line checked holds: `app_ui.R:443-453` and `:507`; `app_server.R:407`, `:437`, `:475`, `:576`; `mod_lipidomics.R:180-196`; `omicsManagedCommit`; the setup-directories transcriptomics entry; `.MODULE_CI_VALID_OMICS`; the e2e "Unsupported omic" test.
   - The exception is "the shared design builder", which does not exist (F01).
   - The wiring needs less than the ticket implies (F13). The hard-coded workflow-type allowlists the GUI depends on are correctly assigned to KIT-094 and KIT-095.

**Overall:** KIT-096 cannot be planned as written until the user decides what "reuse" of the import, design and QC modules means (F01). After that, the eight should-fix findings (F02 to F07, F09, F10) can be settled in the planning step, with two cross-ticket hand-offs to record in KIT-095 (F04, F09) and one in KIT-093 (F07, F08).

## Claim scope

This audit covers the static reading of the ticket and the code at `bbe1fae9`. It claims nothing about runtime, exactness or performance. No finding was closed.
