# KIT-102 W1: the omic GUI configuration contract (kit-101 design 4.1 to 4.4, values 4.11) and its
# validator. One plain list per omic with the groups `shell`, `binding` and `import`. It holds data,
# references to core functions (`utils_workflow_*`, `func_*`) and component ids, never a server
# or UI function. Every value is the one the GUI had at HEAD c1041c56 (before KIT-102), with the
# file:line it comes from there. KIT-102 W5 deleted the lipid and metab import helper files and the
# trx import section; in the provenance notes, `ui`, `processing`, `selection` and `output` are the
# lipid import's UI, processing, selection and output helper files at c1041c56 (for metab, `ui` is
# mod_metab_import_ui.R and `processing` and `display` its processing and display helper files) and
# `trx` is mod_transcriptomics.R at c1041c56. The values now render through the frame
# (mod_omic_import.R) and its source panels (mod_omic_import_feature_table.R,
# `protImportPanelUi` in mod_prot_import_ui.R). `import$preview` and `import$process` are the core's
# (`utils_workflow_omics_import.R`, `func_trx_counts.R`, `func_prot_import_*.R`); trx adds its count
# policy as `import$check`.
#
# KIT-103 W1 adds the `design` group (design 4.5) with the values the design tab and builder have at
# HEAD 3c657c14 (unchanged in R/ since 283f5eee), with their file:line there. KIT-103 W4 to W6
# deleted most of those design files (the lipid design helpers, mod_metab_design_builder_server.R,
# mod_metab_design_import_helpers.R, mod_metab_design_validation.R, the design section of
# mod_transcriptomics.R and mod_prot_design_builder_display_helpers.R), so the design group's
# file:line notes are provenance at 3c657c14; the values now render through the shared design tab
# and builder (mod_omic_design.R, mod_omic_design_builder.R). KIT-103 W7 wired the function fields
# to the core (`func_omics_design_shared.R`, `func_trx_design.R`, the declaration files,
# `utils_workflow_dia_state.R` and the proteomics annotation files) and made the ones every omic
# carries required (`validate`, `build_object`, `group_levels`; the feature-table import's
# `reader`, `record_design_file` and `check_files_first`), as KIT-102 W7 did for the import.
#
# KIT-104 W1 adds the `qc` group (design 4.6) with the values the QC tab and its sub-tabs have at
# HEAD b423b764 (R/ unchanged since f0e90318, the QC files since 321c2294), with their file:line
# there: the title, the tabset, the groups of sub-tabs (a component id each: the shared sub-tabs
# `omic_qc_*`, trx's filter plugin and proteomics' plugins, `mod_<component>_ui/server`), the start
# and completion, the tracking figure, the ITSD analysis, the finalise's per-omic differences, and
# the texts and notification ids the shared sub-tabs render per omic. KIT-104 W5 and W6 deleted
# seven of those files (the lipid intensity and duplicate helper files, mod_metab_qc_itsd.R and four
# metab helper files) and W7 the proteomics coordinator's body, so the qc group's file:line notes are
# provenance at b423b764; the values now render through the shared QC shell and sub-tabs
# (mod_omic_qc.R, mod_omic_qc_tabs.R). KIT-104 W8 made the
# finalise's function fields (`filtering_figure`, `finalize_alignment`), wired to the core by W2 to
# W6, required with the finalise sub-tab, which every feature omic has; `finalize_shell` (trx only)
# stays optional.
#
# KIT-105 W1 adds the `norm` group (design 4.7) with the values the lipid and metab normalisation
# tabs have at HEAD 9a98e765 (the lipid and metab norm files unchanged since 55a53d44), with their
# file:line there: the option sections in order with their options, the QC tabs, the stage grid, the
# RUV tab's texts, the managed run and correlation, the pre-QC on tab entry, the export's per-omic
# values, the ITSD selection table builder, the test ids and the per-omic display values of design
# D7 (the skip-correlation flag order, the final PCA's placeholder sizes). Proteomics (KIT-111) has
# none yet: its `norm` is NULL. The in-process run and correlation and the export's collector and
# saver are the core's (KIT-105 W2/W4) and required since KIT-105 W6. KIT-110 W3 adds trx's group
# (`omicGuiTrxNorm`, from mod_trx_norm.R and mod_transcriptomics.R at 38d9d25e): its option
# sections, its own run and reset handlers (the core's, utils_workflow_trx_norm.R) and the stage
# grid's fixed assay, file stem, stage prefixes and registered rows. Since then the validator
# requires a field only with the tab or panel that uses it (KIT-110 D8): the ITSD table with the
# `itsd` tab; the correlation, skip-correlation flags, final PCA and managed-correlation notice with
# the `correlation_filter` tab; the export with the `export` panel.
#
# KIT-106 W1 adds the `da` group (design 4.8, the ticket's D1 and D8) for all four omics, with the
# values the DA tabs have at the KIT-106 start commit b75561de (R/ unchanged since 46ec9e01), with
# their file:line there: the tab and status key, the test ids, the model (the formula input, or trx's
# model shown), the threshold and fold-change inputs, the table's columns, the summary, the nouns, the
# results tabs, the assay and table-filter rules, the views (static volcano, heatmap, table), the run
# and its source, the session load and its assay check, the refusal order (`preflight_first`), the
# texts, the download and proteomics' plugins (init, load, run, archive). Function fields are the
# core's; those still in an R/mod_*.R file at the start commit (lipid's preflight, trx's model text,
# static volcano, source read, contrast entries and export report) are wired when W2 moves them.
# `testids`, `summary`, `tabs` and `table_filter` are required with D8's always-required fields, as
# every omic carries them.
#
# KIT-107 W1 adds the `enrichment` group (design 4.9, the ticket's D1 and D9): proteomics' enrichment
# tab at the KIT-107 start commit 8e7294b6 (R/ unchanged since b9820d9e; the ten enrichment module
# files since 3880ddb0), with each value's file:line there: the tab, its status key and value, the
# declaration and its run, the DA source, the identifier column, the cutoffs, the correction, the
# excluded IEA flag, the taxon and its prefill, the result services, the organism filter's toggle and
# the completed state's name. Lipid, metab and trx carry `enrichment = NULL` (no enrichment tab). The
# function fields whose functions are still in R/mod_prot_enrich_*.R at the start commit
# (`organism`, `annotation`, `kernel`, and the organism filter's `mapping` and `apply`) are optional
# and absent until the W2 merge wires them; W4 makes them required. The declared parameters and the
# background are the declaration's, not repeated here; labels and help texts stay in the markup (D1).
#
# KIT-108 W1 adds the `summary` group (design 4.10, the ticket's D1, D4 and D9) for all four omics: the Summary &
# Export tab at the KIT-108 start commit 5be34ee5 (R/ identical to a19a6bf2 and to kit_107's source 8b92d7ac), with
# each value's file:line there: the module id and test id prefix, the label placeholder and description, the final
# state, the record's owner class, the publication copy (its gate, its required artifacts and the `project_dirs`
# global's lookup), the integration export (its file name, managed operation, in-process writer and trx's gate on the
# RUVIII export), the report (its template directory, template, label and template status), the export steps (trx's
# RUVIII export, D4) and the texts that differ per omic. The function fields whose functions are still in an R/mod_*.R
# file at the start commit (proteomics' final state, copy gate, integration file name and writer, report template and
# status; trx's final state, integration writer and the export step's choices, availability, parameters and status
# text; every omic's session-state payload) are optional and absent until the W3 (lipid, metab, trx) and W4
# (proteomics) merges wire them; W6 makes them required.

# ------------------------------------------------------------------------------------------------
# The entry point.
# ------------------------------------------------------------------------------------------------

#' The shared GUI layer's configuration of one omic, validated.
#' @param omic_type One of `proteomics`, `lipidomics`, `metabolomics` or `transcriptomics`.
#' @return A plain list with the groups `shell`, `binding`, `import`, `design`, `qc`, `norm`, `da`,
#'   `enrichment` and `summary`.
#' @noRd
omicGuiConfig <- function(omic_type) {
    builders <- list(proteomics = omicGuiConfigProteomics, lipidomics = omicGuiConfigLipidomics,
        metabolomics = omicGuiConfigMetabolomics, transcriptomics = omicGuiConfigTranscriptomics)
    if (!is.character(omic_type) || length(omic_type) != 1L || !omic_type %in% names(builders)) {
        cli::cli_abort(
            "{.arg omic_type} must be one of {.val {names(builders)}}, not {.val {omic_type}}.",
            class = "omic_gui_config_invalid")
    }
    validateOmicGuiConfig(builders[[omic_type]]())
}

# ------------------------------------------------------------------------------------------------
# The validator.
# ------------------------------------------------------------------------------------------------

#' Check an omic GUI configuration: each required field present with its type, the runner and the
#' other choices, every reference a function, the tab and status keys unique.
#' @param config A configuration, as `omicGuiConfig()` builds it.
#' @return `config`, invisibly. Otherwise an `omic_gui_config_invalid` error naming the field.
#' @noRd
validateOmicGuiConfig <- function(config) {
    if (!is.list(config)) omicGuiConfigAbort("config", "must be a list")
    schema <- omicGuiConfigSchema()
    groups <- c(shell = "named_list", binding = "named_list", import = "named_list",
        design = "named_list", qc = "named_list", norm = "named_list_or_null", da = "named_list",
        enrichment = "named_list_or_null", summary = "named_list")
    checkOmicGuiFields(config, groups)
    validateOmicGuiShell(config$shell, schema)
    validateOmicGuiBinding(config$binding, config$shell$omic_type, schema)
    validateOmicGuiImport(config$import, schema)
    validateOmicGuiDesign(config$design, config, schema)
    validateOmicGuiQc(config$qc, config$shell$routes, schema)
    if (!is.null(config$norm)) validateOmicGuiNorm(config$norm, schema)
    validateOmicGuiDa(config$da, schema)
    if (!is.null(config$enrichment)) validateOmicGuiEnrichment(config$enrichment, schema)
    validateOmicGuiSummary(config$summary, config$shell, schema)
    invisible(config)
}

#' The contract's fields and their types (design 4.2 to 4.5), by group. The type codes are
#' `omicGuiFieldHasType`'s.
#' @noRd
omicGuiConfigSchema <- function() {
    list(
        shell = c(omic_type = "chr1", prefix = "chr1", testid = "chr1", label = "chr1_or_null",
            label_lower = "chr1", header_icon = "chr1_or_null", tabset_id = "chr1",
            tab_style = "chr1", routes = "chr", tabs = "list", bus = "named_list",
            bus_registration = "chr1"),
        shell_optional = c(chain_rules = "named_list", extras = "chr"),
        tab = c(status_key = "chr1", value = "chr1", module_id = "chr1", component = "chr1",
            title = "chr1", icon = "chr1", testid = "chr1", step_name = "chr1", step_icon = "chr1"),
        binding = c(binding = "list", runner = "chr1", owner = "ref", read_state = "ref",
            read_metadata = "ref"),
        dia_binding = c(omic_type = "chr1", routes = "chr", route_for_format = "ref", owner = "ref",
            read_state = "ref", read_metadata = "ref", import_metadata = "refs"),
        import = c(source_panel = "chr1", title = "chr1", formats = "named_chr",
            default_format = "chr1", detect = "ref", help = "chr", success_text = "chr1",
            preview = "ref", process = "ref"),
        import_optional = c(route_choice = "named_chr", default_route = "chr1", aux_inputs = "chr0",
            testids = "named_chr", validate_mapping = "ref", check = "ref"),
        feature_table = c(format_heading = "chr1", assay_name = "chr1", id_label = "chr1",
            id_custom_label = "chr1", id_custom_placeholder = "chr1", annotation_label = "chr1",
            annotation_custom_placeholder = "chr1", second_assay = "lgl1", id_mapping_key = "chr1",
            aux_inputs = "chr0", validate_mapping = "ref", inputs = "named_chr",
            outputs = "named_chr", testids = "named_chr", status_modes = "named_list",
            format_labels = "named_chr", validation_noun = "chr1", assay_names = "named_list",
            notification_id = "chr1", refusal_prefix = "chr1", summary_file = "chr1",
            count_key = "chr1"),
        second_assay = c(second_assay_columns = "lgl1", readers = "refs"),
        assay_names = c(rule = "chr1", path_pattern = "chr1", messages = "named_chr"),
        # KIT-103: `run_order` is required and has no shared default.
        design = c(module_id = "chr1", help = "chr1", formula_default = "chr1", run_order = "ref",
            run_order_rule = "chr1", record_provenance = "ref", state_name = "by_route",
            state_description = "by_route", data_adapter = "chr1",
            design_import = "named_list_or_null", commit = "chr1", builder_ui = "chr1",
            record_excluded_runs = "lgl1", label_lower = "chr1", multi_assay_banner = "lgl1",
            no_data_text = "chr1", preview_text = "chr1", preview_assays = "lgl1",
            testids = "named_chr0", builder_texts = "named_chr", group_levels = "ref_or_null",
            validate = "ref_or_null", build_object = "ref"),
        design_optional = c(annotation_step = "refs"),
        design_wide = c(id_mapping_key = "chr1"),
        design_shared_builder = c(multi_assay_note = "chr1", matrix_heading = "chr1"),
        builder_texts = c(save_success = "chr1", no_assignment = "chr1",
            contrast_group_prefix = "chr1", contrast_as_is = "chr1", replicate_suffix = "chr1"),
        design_import = c(plugin = "chr1", import_dir_path = "chr1"),
        design_import_feature_table = c(id_candidates = "named_list", qc_baseline = "ref",
            preflight = "lgl1", preflight_contrasts = "lgl1", state_description = "chr1", reader = "ref",
            record_design_file = "lgl1", check_files_first = "lgl1"),
        id_candidates = c(id = "chr", annotation = "chr"),
        dia_design = c(initial_state = "ref", rename = "ref", subset = "ref", commit = "chr1"),
        # KIT-104 W1 (design 4.6). `complete_state` and `itsd` are NULL where the omic has none.
        qc = c(title = "chr1", tabset_id = "by_route", groups = "list", auto_init = "lgl1",
            complete_state = "route_map_or_null", itsd = "named_list_or_null", texts = "named_chr0",
            notification_ids = "named_chr0"),
        # KIT-104 W6 (F104-04): the finalise's shell read, `function(binding, manager)`, for an omic
        # whose class cannot be rebuilt as the generic shell (trx); optional.
        qc_read = c(finalize_shell = "ref"),
        qc_group = c(id = "chr1_or_null", title = "chr1_or_null", component = "chr1_or_null",
            tabs = "list"),
        qc_tab = c(id = "chr1", title = "chr1", component = "chr1", routes = "chr",
            needs_paths = "lgl1"),
        # Required with the shared finalise sub-tab; its function fields since KIT-104 W8, once their
        # functions were in the core (every feature omic carries them).
        qc_finalize = c(finalize_zero_rows = "lgl1", finalize_opens_norm = "lgl1",
            finalize_per_assay = "lgl1", filtering_figure = "ref", finalize_alignment = "ref"),
        qc_itsd = c(name_columns = "chr", metrics = "ref", quant = "ref", tip = "chr1"),
        # The texts and the notification ids each shared sub-tab renders, by its component.
        qc_texts = list(
            omic_qc_intensity = c("intensity_heading", "intensity_help", "intensity_proportion_help",
                "intensity_original_label", "intensity_working", "intensity_result_heading",
                "intensity_total", "intensity_success", "intensity_error"),
            omic_qc_duplicates = c("duplicates_help", "duplicates_none"),
            omic_qc_itsd = "itsd_pattern_help",
            omic_qc_finalize = c("finalize_total_label", "finalize_id_label", "finalize_total",
                "finalize_no_assays", "finalize_no_design")),
        qc_notifications = c(omic_qc_intensity = "intensity", omic_qc_duplicates = "duplicates",
            omic_qc_itsd = "itsd"),
        # KIT-105 W1 (design 4.7), checked when the omic has the group. KIT-110 D8: these always; the
        # others with the tab or panel that uses them (`norm_by_tab`, `norm_by_panel`), or by type when present.
        norm = c(module_id = "chr1", panels = "named_list", plot_tabs = "chr", stage_labels = "chr",
            assay_slots = "int1", ruv_heading = "chr1", ruv_help = "chr1", testids = "named_chr"),
        # KIT-120 W1: the fields of the frame's own actions, required unless the omic registers its own server
        # (`norm$server`): the run and the pre-QC on entry always, the others with their tab or panel.
        norm_frame = c(run = "named_list", pre_qc_on_entry = "lgl1"),
        norm_frame_by_tab = list(
            correlation_filter = c(correlation = "named_list", skip_correlation_flags = "chr1", final_pca = "named_list",
                # KIT-105 W5 (D7): whether a failed managed correlation is also shown as an error notification.
                managed_correlation_notify = "lgl1")),
        norm_frame_by_panel = list(export = c(export = "named_list")),
        # KIT-120 W1: the Post-Normalisation QC tab's heading and help (`norm_post_norm_qc`).
        norm_by_tab = list(itsd = c(itsd_table = "ref"), post_norm_qc = c(post_norm_qc = "named_list")),
        norm_post_norm_qc = c(heading = "chr1", help = "chr1"),
        # KIT-110 D3, D4, D7: an omic's own reset handler; the stage grid's fixed assay names, file stem, stage file
        # prefixes (`norm_stage_prefixes`) and registered assay rows.
        # KIT-120 W1: an omic's own server registrations, by name as the DA plugins are (a module function in an
        # R/mod_*.R file: `server(input, output, session, context)`, its actions and its outputs, in place of the frame's); the QC directory's key in the experiment paths (else the binding's
        # `qc_dir_role`); whether the stage grid shows an assay label row (shown unless FALSE) and names its files
        # `<stage prefix>_<kind>.png`, with no assay stem; the id of a UI output under the correlation tab's summary.
        # KIT-122 (audit F5): whether the RUV QC tab shows one result's plot, summary and table as static outputs (the
        # library's `omicNormRuvSingleResult`) in place of the per-assay `ruv_qc_ui`.
        norm_optional = c(reset = "ref", assay_names = "chr", image_stem = "chr1", stage_prefixes = "named_chr",
            output_assay_slots = "int1", server = "chr1", qc_dir = "chr1", assay_labels = "lgl1", image_stemless = "lgl1",
            correlation_slot = "chr1", ruv_single_result = "lgl1",
            # KIT-129 W2: the RUV design presets (`validateOmicGuiRuvDesign`), checked when the group has them.
            ruv_design = "named_list"),
        # KIT-129 W2 (ticket sections 1 and 2): the RUV design presets. `enabled` is the gate of the panel
        # section; `tiers` the tier labels; `k_rules` the rules for k with the texts shown beside the choice;
        # `tables` one preset table per route family, the first the one rendered, `route_tables` the workflow
        # types that take another; `conditions` the conditions the omic's simulations cover, shown once.
        norm_ruv_design = c(enabled = "lgl1", testid = "chr1", tiers = "named_chr", k_rules = "named_list",
            tables = "named_list"),
        norm_ruv_design_optional = c(route_tables = "named_chr", conditions = "chr1"),
        norm_ruv_tiers = c("recommended", "valid_with_conditions", "legacy", "none"),
        norm_ruv_k_rules = c("num_sv_be", "sweep_bounded", "user"),
        norm_ruv_k_rule = c(label = "chr1", reads = "chr1", pros = "chr1", cons = "chr1", diagnostic = "chr1"),
        # A table: the design selected when the tab opens, the recommended one, its options; the design selected
        # when the design matrix declares technical replicates, by type.
        norm_ruv_table = c(selected = "chr1", recommended = "chr1", options = "list"),
        # KIT-129 fourth pass item 2: the table's features in words, its texts (reasons, labels and refusals of the
        # run; one wording table, `omicGuiRuvTexts`) and the k rules the run path does not take on this route, with
        # the reason.
        norm_ruv_table_optional = c(selected_with_technical_replicates = "chr1", features = "chr1", texts = "named_chr",
            k_rules_not_offered = "named_chr"),
        # Every table offers each of these once (the user's rule: every design in every omic); a route may add
        # the optional ones.
        norm_ruv_designs = c("groups_all_features", "technical_replicates", "supplied_controls", "first_pass_nulls",
            "legacy", "none"),
        norm_ruv_designs_optional = "legacy_correction",
        # An option: `status` says whether its figure is measured, pending (a simulation still running) or
        # unavailable (with the reason). A measured option carries `norm_ruv_figure`; no other carries a number.
        norm_ruv_option = c(id = "chr1", tier = "chr1", label = "chr1", description = "chr1", tooltip = "chr1",
            note = "chr1", status = "chr1", source = "chr1", test_form = "chr1", k_rules = "chr0"),
        norm_ruv_option_optional = c(replicate_structure = "chr1", control_rule = "chr1", fdp = "num1",
            fdp_where = "chr1", power = "num1", single_scenario = "lgl1", k_rule = "chr1", k_rules_pending = "chr",
            k_rule_figures = "named_list",
            unavailable = "chr1",
            # KIT-129 fourth pass item 2: a precondition on the design (`technical_replicates` or
            # `no_technical_replicates`; the reason is the table's text of that name), and a reason the run path does
            # not take the option on this route (the option keeps its figure).
            needs = "chr1", not_offered = "chr1"),
        norm_ruv_status = c("measured", "pending", "unavailable"),
        norm_ruv_figure = c("fdp", "fdp_where", "power"),
        norm_stage_prefixes = c("pre_norm", "post_norm", "ruv_corrected"),
        # `run` and `correlation`: the managed function and its in-process twin (the core's; D6).
        norm_action = c(managed = "ref", in_process = "ref"),
        # KIT-110 D4: an omic's own run handler and the settings it takes from the tab's inputs.
        norm_run_optional = c(settings = "ref", handler = "ref"),
        # `export`: the per-omic values (D7, P-18), and the core's collector and saver.
        norm_export = c(prefix = "chr1", counts_field = "chr1", samples = "chr1", progress = "named_list",
            summary_time = "chr1", collect = "ref", save = "ref"),
        # KIT-105 W5 (D7): the export's progress display: the gathering step's text, the closing step's text (none
        # when empty) and whether the saver steps the display through its files.
        norm_export_progress = c(gather = "chr1", complete = "chr0", steps = "lgl1"),
        norm_testids = c(method = "chr1", ruv_mode = "chr1", run = "chr1", export = "chr1", reset = "chr1"),
        # KIT-110 D8: the panels that need a test id (`run` always; `reset` optional).
        norm_testid_panels = list(method = "between_sample", ruv_mode = c("ruv_iii", "ruv_grid"), export = "export"),
        # The managed final PCA's placeholders (D7: "Preparing ..." and "Error: ..." text sizes).
        norm_final_pca = c(preparing_size = "num1", error_size = "num1"),
        # The shared panel library's option sections, in their order, and each one's options.
        norm_panels = list(aesthetics = character(), itsd = c(aggregation = "named_chr"),
            log2 = c(offset = "num1"), between_sample = c(methods = "named_chr", default = "chr1"),
            # KIT-110 D2: trx's TMM/logCPM and RUV grid.
            tmm_logcpm = character(),
            # KIT-120 W1: the manual k may be unset (proteomics fills it from the canonical correlation).
            ruv_iii = c(modes = "named_chr", default = "chr1", percentage = "num1", k = "num1_or_null", noun = "chr1"),
            ruv_grid = c(modes = "named_chr", default = "chr1", percentage_min = "num1", percentage_max = "num1",
                percentage_step = "num1", k_penalty_weight = "num1", percentage = "num1", k = "num1"),
            actions = character(), export = c(button = "chr1"),
            # KIT-120 W1: a download button (proteomics' S4 compatibility file; its handler is the omic's server's).
            s4_download = c(button = "chr1", label = "chr1", help = "chr1")),
        # KIT-120 W1: a section's own texts in place of the library's, by type when present; the RUV-III section's help
        # texts by name (`norm_ruv_help`).
        norm_panel_optional = list(between_sample = c(help = "chr1"), ruv_iii = c(help = "named_chr"),
            actions = c(run_label = "chr1"), export = c(label = "chr1", help = "chr1")),
        norm_ruv_help = c("grouping", "automatic", "range", "metric", "penalty", "adaptive", "manual", "manual_k", "skip",
            "skip_alternative"),
        # The QC tabs a frame can show, in their order.
        norm_plot_tabs = c("itsd", "ruv", "pca", "density", "rle", "correlation", "post_norm_qc", "correlation_filter", "log"),
        # KIT-106 W1 (design 4.8, D1, D8): the DA tab. These always; `views$heatmap` with the heatmap tab;
        # `session_load` with the load section (its test id) unless the `load` plugin renders it; `run` unless the
        # `run` plugin is present; `source` unless the `init` plugin is present; the rest by type when present.
        da = c(module_id = "chr1", status_key = "chr1", testids = "named_chr", model = "named_list",
            threshold = "named_list", lfc = "named_list", columns = "named_list", summary = "named_list",
            tabs = "named_chr", assays = "chr1", table_filter = "chr1", views = "named_list", texts = "named_chr",
            download = "named_list"),
        da_optional = c(status_keys = "chr", nouns = "named_chr", status_text = "ref", preflight = "ref",
            preflight_first = "lgl1", contrast_entries = "ref", contrast_entry = "ref", plugins = "named_chr"),
        # `editable`: the formula input (its label, default and help); otherwise the model shown (its heading,
        # output and help), with its text and realised formula by type when present.
        da_model = c(editable = "lgl1", help = "chr"),
        da_model_editable = c(label = "chr1", default = "chr1"),
        da_model_shown = c(heading = "chr1", output = "chr1"),
        da_model_optional = c(text = "ref", formula = "ref"),
        da_input = c(label = "chr1", value = "num1", min = "num1", max = "num1", step = "num1"),
        da_lfc = c(tested = "lgl1"),
        # The table's columns: an omic that names its identifier by the object's slot (`id_slot`) has no `id`.
        da_columns = c(fold_change = "chr1", pvalue = "chr1", significance = "chr1"),
        da_columns_optional = c(id = "chr1", label = "chr1", id_slot = "chr1", significance_label = "chr1"),
        da_summary = c(noun = "chr1", significance_label = "chr1"),
        da_tabs = c("volcano", "heatmap", "table"),
        da_texts = c(status_waiting = "chr1", contrasts_none = "chr1"),
        da_texts_optional = c(status_ready = "chr1"),
        da_static_volcano = c(output = "chr1", show_labels = "lgl1", n_labels = "num1"),
        da_heatmap = c(kernel = "ref", labels_argument = "chr1", tree_cut = "named_list", save_prefix = "chr1",
            save_note = "lgl1"),
        da_tree_cut = c(default = "chr1", n_clusters = "num1", cut_height = "num1", min_cluster_size = "num1"),
        da_table = c(significance = "chr1", max_rows = "num1"),
        da_download = c(stem = "chr1", extension = "chr1"),
        da_run = c(managed = "ref", in_process = "ref"),
        da_source = c(class = "chr1", no_data = "chr1"),
        da_session_load = c(file = "chr1", assay_rule = "chr1", assay_message = "chr1"),
        da_plugins = c("init", "load", "run", "archive"),
        # KIT-107 W1 (design 4.9, D1, D9): the enrichment tab, checked when the omic has the group. The function
        # function groups `organism`, `annotation` and `kernel` (the core's, wired at the W2 merge) are required since
        # KIT-107 W4 and checked whole; `organism_filter` is NULL or its toggle's bus key (`mixed_species_from`) with
        # its `mapping` and `apply`.
        enrichment = c(module_id = "chr1", status_key = "chr1", tab_value = "chr1", operation = "ref", run = "ref",
            da_source = "named_list", id_column = "named_list", cutoffs = "named_list", correction = "named_list",
            exclude_iea = "lgl1", taxon = "named_list", services = "chr", organism_filter = "named_list_or_null",
            state_name = "chr1"),
        enrichment_optional = c(organism = "named_list", annotation = "named_list", kernel = "named_list"),
        enrichment_da_source = c(manifest = "ref", components = "chr", read = "ref"),
        enrichment_id_column = c(slot = "chr1", default = "chr1", gprofiler_override = "chr1"),
        enrichment_cutoffs = c("up", "down", "q"),
        enrichment_cutoff = c(value = "num1", min = "num1", max = "num1", step = "num1"),
        enrichment_correction = c(choices = "named_chr", selected = "chr1"),
        enrichment_taxon = c(default = "chr1", prefill_from = "chr1"),
        enrichment_services = c("gprofiler2", "clusterprofiler", "stringdb"),
        enrichment_organism = c(supported = "ref", method = "ref"),
        enrichment_annotation = c(resolve = "ref", match = "ref"),
        enrichment_kernel = c(select = "ref", test_double = "ref"),
        enrichment_organism_filter = c(mixed_species_from = "chr1"),
        enrichment_organism_filter_optional = c(mapping = "ref", apply = "ref"),
        # KIT-108 W1 (design 4.10, D1, D4, D9): the Summary & Export tab, every omic. The function fields W3 and W4 wired
        # to the core are required since KIT-108 W6 (the session state, the final state's context or states, the
        # integration file name, the report template, an export step's functions). `export_steps` is a list, possibly
        # empty, of steps; `texts` the per-omic texts.
        summary = c(module_id = "chr1", testid = "chr1", placeholder = "chr1", description = "chr1",
            final_state = "named_list", copy = "named_list", integration = "named_list", report = "named_list",
            texts = "named_chr0", session_state = "ref"),
        summary_optional = c(record_class = "chr1_or_null"),
        # The final state: the head's by its `context` function, or the first of `states` with its `read_ids` reader.
        summary_final_state = c(context = "ref", states = "chr", read_ids = "ref"),
        # `required`: NULL (failures not classified) or the required artifacts' display names, checked by the validator.
        summary_copy = c(inherits = "lgl1"),
        summary_copy_optional = c(gate = "ref_or_null"),
        # `operation` and `requires` are present, NULL or a string; `filename` a string or a function; `write` the
        # in-process writer (NULL: the final object's RDS). KIT-136: `requires` may also be a function of the session's
        # bus giving the export's name or NULL (transcriptomics needs no RUVIII export when RUV was not applied).
        summary_integration = c(requires = "chr1_or_ref_or_null"),
        summary_integration_optional = c(write = "ref_or_null", missing_message = "chr1"),
        # `template`: a file name, or the route's resolver; a file name needs the `label` of the one-template status.
        # KIT-108 W1's captures: the report flow's per-omic order (`render_after_failed_template`: metab renders after
        # a failed template retrieval; `checking_progress_always`: proteomics shows the template check's progress when
        # the template is already in the project; `record_after_success`: proteomics exports the record bundle after
        # the report's success notice and outputs, the others before).
        summary_report = c(dir = "chr1", label = "chr1_or_null", render_after_failed_template = "lgl1",
            checking_progress_always = "lgl1", record_after_success = "lgl1"),
        summary_report_optional = c(status = "ref_or_null"),
        summary_export_step = c(id = "chr1", title = "chr1", help = "chr", choice = "named_list", button = "named_list",
            status = "named_list", figure = "named_list", operation = "chr1", accept = "ref", run_reference = "ref",
            publications = "list", error = "named_list"),
        summary_export_step_optional = c(parameters = "ref"),
        summary_export_choice = c(output = "chr1", input = "chr1", label = "chr1", selected = "chr1",
            unavailable = "chr"),
        summary_export_choice_optional = c(choices = "ref", available = "ref"),
        summary_export_button = c(input = "chr1", label = "chr1", testid = "chr1", icon = "chr1", class = "chr1"),
        summary_export_status = c(output = "chr1", idle = "chr1"),
        summary_export_status_optional = c(text = "ref"),
        summary_export_figure = c(output = "chr1", file = "chr1", dir = "chr1", height = "chr1", image = "named_chr",
            alt = "chr1", alt_idle = "chr1"),
        summary_export_publication = c(operation = "chr1", parameters = "named_list"),
        summary_export_error = c(prefix = "chr1", type = "chr1", duration = "num1")
    )
}

#' @noRd
validateOmicGuiShell <- function(shell, schema) {
    checkOmicGuiFields(shell, schema$shell, "shell")
    checkOmicGuiFields(shell, schema$shell_optional, "shell", required = FALSE)
    checkOmicGuiChoice(shell$tab_style, "shell$tab_style", c("icon_argument", "icon_in_title"))
    checkOmicGuiChoice(shell$bus_registration, "shell$bus_registration", c("omics", "userData"))
    for (i in seq_along(shell$tabs)) {
        checkOmicGuiFields(shell$tabs[[i]], schema$tab, sprintf("shell$tabs[[%d]]", i))
    }
    for (key in c("status_key", "value", "module_id", "testid")) {
        keys <- vapply(shell$tabs, `[[`, character(1), key)
        if (anyDuplicated(keys)) {
            omicGuiConfigAbort(paste0("shell$tabs$", key),
                sprintf("must be unique across the tabs (repeated: %s)", keys[anyDuplicated(keys)]))
        }
    }
    # The initial `tab_status` is derived from the tabs; the bus names its place.
    if (!"tab_status" %in% names(shell$bus)) {
        omicGuiConfigAbort("shell$bus", "must name `tab_status`")
    }
    invisible(shell)
}

#' @noRd
validateOmicGuiBinding <- function(binding, omic_type, schema) {
    checkOmicGuiFields(binding, schema$binding, "binding")
    checkOmicGuiChoice(binding$runner, "binding$runner", c("omics_actions", "dia_managed_step"))
    if (identical(binding$runner, "dia_managed_step")) {
        checkOmicGuiFields(binding$binding, schema$dia_binding, "binding$binding")
    }
    if (!identical(binding$binding$omic_type, omic_type)) {
        omicGuiConfigAbort("binding$binding", sprintf("must be the %s binding", omic_type))
    }
    invisible(binding)
}

#' @noRd
validateOmicGuiImport <- function(import, schema) {
    checkOmicGuiFields(import, schema$import, "import")
    checkOmicGuiFields(import, schema$import_optional, "import", required = FALSE)
    checkOmicGuiChoice(import$source_panel, "import$source_panel", c("feature_table", "proteomics"))
    if (!identical(import$source_panel, "feature_table")) return(invisible(import))
    checkOmicGuiFields(import, schema$feature_table, "import")
    checkOmicGuiFields(import$assay_names, schema$assay_names, "import$assay_names")
    checkOmicGuiChoice(import$assay_names$rule, "import$assay_names$rule",
        c("require_pair", "skip_incomplete"))
    if (import$second_assay) checkOmicGuiFields(import, schema$second_assay, "import")
    invisible(import)
}

#' The design group: its fields and choices; a per-route state names every route; the wide-table
#' adapter needs the column-mapping key, the shared builder UI its texts, and the `dia_long` adapter
#' the DIA adaptor's design entries with the same commit hook.
#' @noRd
validateOmicGuiDesign <- function(design, config, schema) {
    checkOmicGuiFields(design, schema$design, "design")
    checkOmicGuiFields(design, schema$design_optional, "design", required = FALSE)
    checkOmicGuiChoice(design$run_order_rule, "design$run_order_rule", c("mixedsort", "input"))
    checkOmicGuiChoice(design$data_adapter, "design$data_adapter", c("wide_assay_list", "dia_long"))
    checkOmicGuiChoice(design$commit, "design$commit", c("shared", "prot_hook"))
    checkOmicGuiChoice(design$builder_ui, "design$builder_ui", c("shared", "proteomics"))
    checkOmicGuiFields(design$builder_texts, schema$builder_texts, "design$builder_texts")
    for (field in c("state_name", "state_description")) {
        value <- design[[field]]
        if (length(value) > 1L && !setequal(names(value), config$shell$routes)) {
            omicGuiConfigAbort(paste0("design$", field), "must name each of the omic's routes")
        }
    }
    if (identical(design$builder_ui, "shared")) {
        checkOmicGuiFields(design, schema$design_shared_builder, "design")
    }
    if (identical(design$data_adapter, "wide_assay_list")) {
        checkOmicGuiFields(design, schema$design_wide, "design")
    } else {
        adaptor <- config$binding$binding$design
        checkOmicGuiFields(adaptor, schema$dia_design, "binding$binding$design")
        if (!identical(design$commit, adaptor$commit)) {
            omicGuiConfigAbort("design$commit", "must be the DIA adaptor's commit hook")
        }
    }
    if (!is.null(design$design_import)) validateOmicGuiDesignImport(design$design_import, schema)
    invisible(design)
}

#' The import-existing-design plugin: the feature table's reader values, or proteomics' own panel.
#' @noRd
validateOmicGuiDesignImport <- function(design_import, schema) {
    path <- "design$design_import"
    checkOmicGuiFields(design_import, schema$design_import, path)
    checkOmicGuiChoice(design_import$plugin, paste0(path, "$plugin"),
        c("feature_table", "proteomics"))
    checkOmicGuiChoice(design_import$import_dir_path, paste0(path, "$import_dir_path"),
        c("req", "empty"))
    if (!identical(design_import$plugin, "feature_table")) return(invisible(design_import))
    checkOmicGuiFields(design_import, schema$design_import_feature_table, path)
    checkOmicGuiFields(design_import$id_candidates, schema$id_candidates,
        paste0(path, "$id_candidates"))
    invisible(design_import)
}

#' The qc group (design 4.6): its fields; per-route values naming the omic's routes; the groups of
#' sub-tabs; and what the shared sub-tabs on them need (the finalise fields, the ITSD analysis, the
#' texts and the notification ids).
#' @param routes The omic's routes (`shell$routes`).
#' @noRd
validateOmicGuiQc <- function(qc, routes, schema) {
    checkOmicGuiFields(qc, schema$qc, "qc")
    checkOmicGuiFields(qc, schema$qc_read, "qc", required = FALSE)
    if (length(qc$tabset_id) > 1L && !setequal(names(qc$tabset_id), routes)) {
        omicGuiConfigAbort("qc$tabset_id", "must name each of the omic's routes")
    }
    if (!all(names(qc$complete_state) %in% routes)) {
        omicGuiConfigAbort("qc$complete_state", "must name only the omic's routes")
    }
    components <- validateOmicGuiQcGroups(qc$groups, routes, schema)
    if ("omic_qc_finalize" %in% components) {
        checkOmicGuiFields(qc, schema$qc_finalize, "qc")
        if (qc$finalize_per_assay) checkOmicGuiFields(qc$texts, c(finalize_count_column = "chr1"), "qc$texts")
    }
    if ("omic_qc_itsd" %in% components) {
        if (is.null(qc$itsd)) omicGuiConfigAbort("qc$itsd", "is required by the ITSD sub-tab")
        checkOmicGuiFields(qc$itsd, schema$qc_itsd, "qc$itsd")
    }
    for (component in intersect(components, names(schema$qc_texts))) {
        keys <- schema$qc_texts[[component]]
        checkOmicGuiFields(qc$texts, stats::setNames(rep("chr1", length(keys)), keys), "qc$texts")
    }
    for (component in intersect(components, names(schema$qc_notifications))) {
        checkOmicGuiFields(qc$notification_ids,
            stats::setNames("chr1", schema$qc_notifications[[component]]), "qc$notification_ids")
    }
    invisible(qc)
}

#' The QC groups: one unnamed group of flat tabs (lipid, metab, trx), or named groups rendered by
#' their group plugins (proteomics); tab ids unique within a group, tab routes the omic's.
#' @return The components the groups' tabs name.
#' @noRd
validateOmicGuiQcGroups <- function(groups, routes, schema) {
    for (i in seq_along(groups)) {
        path <- sprintf("qc$groups[[%d]]", i)
        checkOmicGuiFields(groups[[i]], schema$qc_group, path)
        if (!is.null(groups[[i]]$id) && (is.null(groups[[i]]$title) || is.null(groups[[i]]$component))) {
            omicGuiConfigAbort(path, "is a named group, so it needs a title and a component")
        }
        tabs <- groups[[i]]$tabs
        for (j in seq_along(tabs)) {
            checkOmicGuiFields(tabs[[j]], schema$qc_tab, sprintf("%s$tabs[[%d]]", path, j))
            if (!all(tabs[[j]]$routes %in% routes)) {
                omicGuiConfigAbort(sprintf("%s$tabs[[%d]]$routes", path, j), "must name only the omic's routes")
            }
        }
        ids <- vapply(tabs, `[[`, character(1), "id")
        if (anyDuplicated(ids)) {
            omicGuiConfigAbort(paste0(path, "$tabs$id"),
                sprintf("must be unique within the group (repeated: %s)", ids[anyDuplicated(ids)]))
        }
    }
    unnamed <- vapply(groups, function(group) is.null(group$id), logical(1))
    if (any(unnamed) && length(groups) > 1L) {
        omicGuiConfigAbort("qc$groups", "may hold an unnamed group only as its one group")
    }
    ids <- unlist(lapply(groups, `[[`, "id"))
    if (anyDuplicated(ids)) {
        omicGuiConfigAbort("qc$groups$id", sprintf("must be unique (repeated: %s)", ids[anyDuplicated(ids)]))
    }
    unique(unlist(lapply(groups, function(group) vapply(group$tabs, `[[`, character(1), "component"))))
}

#' The norm group (design 4.7): its fields; the option sections from the shared panel library, in
#' its order, each with its options; the QC tabs from the frame's set; three stage labels; the
#' managed run and its in-process twin; the correlation, the export's values and choices and the D7
#' display values with the tab or panel that uses them (KIT-110 D8); the test ids; and the optional
#' fields by type (an omic's own run and reset handlers, the stage grid's fixed values).
#' @noRd
validateOmicGuiNorm <- function(norm, schema) {
    checkOmicGuiFields(norm, schema$norm, "norm")
    sections <- names(schema$norm_panels)
    panels <- names(norm$panels)
    if (length(panels) == 0L || !all(panels %in% sections) || is.unsorted(match(panels, sections))) {
        omicGuiConfigAbort("norm$panels", sprintf("must name option sections of the panel library in its order (%s)",
            paste(sections, collapse = ", ")))
    }
    for (panel in panels) {
        checkOmicGuiFields(norm$panels[[panel]], schema$norm_panels[[panel]], paste0("norm$panels$", panel))
        checkOmicGuiFields(norm$panels[[panel]], schema$norm_panel_optional[[panel]], paste0("norm$panels$", panel),
            required = FALSE)
    }
    if (!all(names(norm$panels$ruv_iii$help) %in% schema$norm_ruv_help)) {
        omicGuiConfigAbort("norm$panels$ruv_iii$help", sprintf("must name help texts of the RUV-III section (%s)",
            paste(schema$norm_ruv_help, collapse = ", ")))
    }
    for (panel in intersect(c("between_sample", "ruv_iii", "ruv_grid"), panels)) {
        choices <- norm$panels[[panel]][[if (identical(panel, "between_sample")) "methods" else "modes"]]
        checkOmicGuiChoice(norm$panels[[panel]]$default, paste0("norm$panels$", panel, "$default"), choices)
    }
    tabs <- norm$plot_tabs
    if (anyDuplicated(tabs) || !all(tabs %in% schema$norm_plot_tabs) || is.unsorted(match(tabs, schema$norm_plot_tabs))) {
        omicGuiConfigAbort("norm$plot_tabs", sprintf("must name the frame's QC tabs once each, in their order (%s)",
            paste(schema$norm_plot_tabs, collapse = ", ")))
    }
    for (tab in intersect(names(schema$norm_by_tab), tabs)) checkOmicGuiFields(norm, schema$norm_by_tab[[tab]], "norm")
    if ("post_norm_qc" %in% tabs) checkOmicGuiFields(norm$post_norm_qc, schema$norm_post_norm_qc, "norm$post_norm_qc")
    checkOmicGuiFields(norm, schema$norm_optional, "norm", required = FALSE)
    if (!is.null(norm$ruv_design)) validateOmicGuiRuvDesign(norm$ruv_design, schema)
    # KIT-120 W1: an omic with its own server registers its actions itself; the frame's action fields are then unused.
    framed <- is.null(norm$server)
    # (audit F4) The frame looks nothing up for an omic with its own server, so the group names the stage grid's
    # QC directory and assay itself.
    if (!framed) checkOmicGuiFields(norm, c(qc_dir = "chr1", assay_names = "chr"), "norm")
    if (framed) {
        checkOmicGuiFields(norm, schema$norm_frame, "norm")
        for (tab in intersect(names(schema$norm_frame_by_tab), tabs)) checkOmicGuiFields(norm, schema$norm_frame_by_tab[[tab]], "norm")
        for (panel in intersect(names(schema$norm_frame_by_panel), panels)) {
            checkOmicGuiFields(norm, schema$norm_frame_by_panel[[panel]], "norm")
        }
    }
    if (!is.null(norm$stage_prefixes) && !setequal(names(norm$stage_prefixes), schema$norm_stage_prefixes)) {
        omicGuiConfigAbort("norm$stage_prefixes", sprintf("must name exactly the stages %s",
            paste(schema$norm_stage_prefixes, collapse = ", ")))
    }
    if (length(norm$stage_labels) != 3L) omicGuiConfigAbort("norm$stage_labels", "must hold the three stages' captions")
    if (framed) {
        checkOmicGuiFields(norm$run, schema$norm_action, "norm$run")
        checkOmicGuiFields(norm$run, schema$norm_run_optional, "norm$run", required = FALSE)
        # KIT-115 (A110-04): the frame calls an omic's own handler with `run$settings(input)`.
        if (!is.null(norm$run$handler) && is.null(norm$run$settings)) {
            omicGuiConfigAbort("norm$run$settings", "is required with norm$run$handler")
        }
    }
    if (framed && "correlation_filter" %in% tabs) {
        checkOmicGuiFields(norm$correlation, schema$norm_action, "norm$correlation")
        checkOmicGuiChoice(norm$skip_correlation_flags, "norm$skip_correlation_flags", c("before_commit", "after_commit"))
        checkOmicGuiFields(norm$final_pca, schema$norm_final_pca, "norm$final_pca")
    }
    if (framed && "export" %in% panels) {
        checkOmicGuiFields(norm$export, schema$norm_export, "norm$export")
        checkOmicGuiChoice(norm$export$samples, "norm$export$samples", c("design", "columns"))
        checkOmicGuiFields(norm$export$progress, schema$norm_export_progress, "norm$export$progress")
        checkOmicGuiChoice(norm$export$summary_time, "norm$export$summary_time", c("session", "written"))
    }
    needed <- c("run", names(Filter(function(owners) any(owners %in% panels), schema$norm_testid_panels)))
    checkOmicGuiFields(norm$testids, schema$norm_testids[needed], "norm$testids")
    checkOmicGuiFields(norm$testids, schema$norm_testids[setdiff(names(schema$norm_testids), needed)], "norm$testids",
        required = FALSE)
    invisible(norm)
}

#' The RUV design presets (KIT-129 W2; ticket sections 1 and 2): the gate, the four tiers' labels, the three k
#' rules with their texts, and each table.
#' @param design The norm group's `ruv_design`.
#' @param schema `omicGuiConfigSchema()`.
#' @return `design`, invisibly. Otherwise an `omic_gui_config_invalid` error naming the field.
#' @noRd
validateOmicGuiRuvDesign <- function(design, schema) {
    path <- "norm$ruv_design"
    checkOmicGuiFields(design, schema$norm_ruv_design, path)
    checkOmicGuiFields(design, schema$norm_ruv_design_optional, path, required = FALSE)
    if (!setequal(names(design$tiers), schema$norm_ruv_tiers)) {
        omicGuiConfigAbort(paste0(path, "$tiers"), sprintf("must name exactly the tiers %s",
            paste(schema$norm_ruv_tiers, collapse = ", ")))
    }
    if (!identical(names(design$k_rules), schema$norm_ruv_k_rules)) {
        omicGuiConfigAbort(paste0(path, "$k_rules"), sprintf("must name the k rules %s, in this order",
            paste(schema$norm_ruv_k_rules, collapse = ", ")))
    }
    for (rule in names(design$k_rules)) {
        checkOmicGuiFields(design$k_rules[[rule]], schema$norm_ruv_k_rule, paste0(path, "$k_rules$", rule))
    }
    if (length(design$tables) == 0L) omicGuiConfigAbort(paste0(path, "$tables"), "must hold at least one table")
    if (!all(design$route_tables %in% names(design$tables))) {
        omicGuiConfigAbort(paste0(path, "$route_tables"), sprintf("must name tables of the group (%s)",
            paste(names(design$tables), collapse = ", ")))
    }
    for (name in names(design$tables)) {
        validateOmicGuiRuvTable(design$tables[[name]], paste0(path, "$tables$", name), schema)
    }
    invisible(design)
}

#' One preset table: every design offered once, the selected ones among the available options, and the one
#' option of the recommended tier named as the recommended design.
#' @param table The table.
#' @param path Its path in the configuration, for the refusal.
#' @param schema `omicGuiConfigSchema()`.
#' @return `table`, invisibly.
#' @noRd
validateOmicGuiRuvTable <- function(table, path, schema) {
    checkOmicGuiFields(table, schema$norm_ruv_table, path)
    checkOmicGuiFields(table, schema$norm_ruv_table_optional, path, required = FALSE)
    for (i in seq_along(table$options)) {
        validateOmicGuiRuvOption(table$options[[i]], sprintf("%s$options[[%d]]", path, i), schema)
    }
    ids <- vapply(table$options, `[[`, character(1), "id")
    allowed <- c(schema$norm_ruv_designs, schema$norm_ruv_designs_optional)
    if (anyDuplicated(ids) || !all(schema$norm_ruv_designs %in% ids) || !all(ids %in% allowed)) {
        omicGuiConfigAbort(paste0(path, "$options"), sprintf("must offer each design once (%s; optionally %s)",
            paste(schema$norm_ruv_designs, collapse = ", "), paste(schema$norm_ruv_designs_optional, collapse = ", ")))
    }
    statuses <- vapply(table$options, `[[`, character(1), "status")
    selectable <- ids[statuses != "unavailable"]
    # KIT-129 fourth pass item 2: a precondition's reason is one of the table's texts (the reason when it is not met,
    # `omicNormRuvUnavailableReason`); a rule left out is a k rule.
    reasons <- c(technical_replicates = "no_technical_replicates", no_technical_replicates = "technical_replicates_declared")
    for (i in seq_along(table$options)) {
        needs <- table$options[[i]]$needs
        if (!is.null(needs) && !reasons[[needs]] %in% names(table$texts)) {
            omicGuiConfigAbort(sprintf("%s$options[[%d]]$needs", path, i), "must name one of the table's texts")
        }
    }
    if (!all(names(table$k_rules_not_offered) %in% schema$norm_ruv_k_rules)) {
        omicGuiConfigAbort(paste0(path, "$k_rules_not_offered"), "must name k rules")
    }
    for (field in intersect(c("selected", "recommended", "selected_with_technical_replicates"), names(table))) {
        checkOmicGuiChoice(table[[field]], paste0(path, "$", field), selectable)
    }
    tiers <- vapply(table$options, `[[`, character(1), "tier")
    if (!identical(ids[tiers == "recommended"], table$recommended)) {
        omicGuiConfigAbort(paste0(path, "$recommended"), "must be the one option of the recommended tier")
    }
    invisible(table)
}

#' One preset: its fields; the figure with a measured status and with no other; the reason with an unavailable
#' one; the k rules it offers, the one selected, the ones without a measured figure and the figures measured with
#' another rule than the selected one; never the separation
#' rule with all features as controls (the two canonical-correlation curves coincide).
#' @param option The option.
#' @param path Its path in the configuration, for the refusal.
#' @param schema `omicGuiConfigSchema()`.
#' @return `option`, invisibly.
#' @noRd
validateOmicGuiRuvOption <- function(option, path, schema) {
    at <- \(field) paste0(path, "$", field)
    checkOmicGuiFields(option, schema$norm_ruv_option, path)
    checkOmicGuiFields(option, schema$norm_ruv_option_optional, path, required = FALSE)
    checkOmicGuiChoice(option$tier, at("tier"), schema$norm_ruv_tiers)
    checkOmicGuiChoice(option$status, at("status"), schema$norm_ruv_status)
    measured <- identical(option$status, "measured")
    for (field in schema$norm_ruv_figure) {
        if (measured && !field %in% names(option)) omicGuiConfigAbort(at(field), "is missing for a measured option")
        if (!measured && field %in% names(option)) {
            omicGuiConfigAbort(at(field), "must be absent unless the option's status is measured")
        }
    }
    for (field in intersect(c("fdp", "power"), names(option))) {
        if (option[[field]] < 0 || option[[field]] > 1) omicGuiConfigAbort(at(field), "must be a proportion (0 to 1)")
    }
    if (identical(option$status, "unavailable") != ("unavailable" %in% names(option))) {
        omicGuiConfigAbort(at("unavailable"), "must give the reason of an unavailable option, and of no other")
    }
    if (!is.null(option$needs)) checkOmicGuiChoice(option$needs, at("needs"), c("technical_replicates", "no_technical_replicates"))
    rules <- option$k_rules
    if (!all(rules %in% schema$norm_ruv_k_rules) || anyDuplicated(rules)) {
        omicGuiConfigAbort(at("k_rules"), sprintf("must name k rules once each (%s)",
            paste(schema$norm_ruv_k_rules, collapse = ", ")))
    }
    if (identical(option$control_rule, "all_features") && "sweep_bounded" %in% rules) {
        omicGuiConfigAbort(at("k_rules"), "must not offer the separation rule with all features as controls")
    }
    if ((length(rules) > 0L) != ("k_rule" %in% names(option)) || !all(option[["k_rule"]] %in% rules)) {
        omicGuiConfigAbort(at("k_rule"), "must be one of the option's k rules, and absent when it offers none")
    }
    if (!all(option$k_rules_pending %in% setdiff(rules, option[["k_rule"]]))) {
        omicGuiConfigAbort(at("k_rules_pending"), "must name k rules the option offers, other than the selected one")
    }
    # The figure measured with another rule than the selected one: a rule the option offers, not a pending one.
    others <- setdiff(rules, c(option[["k_rule"]], option$k_rules_pending, "user"))
    if (!all(names(option$k_rule_figures) %in% others)) {
        omicGuiConfigAbort(at("k_rule_figures"),
            "must be named by k rules the option offers, other than the selected, a pending one and the user's number")
    }
    for (rule in names(option$k_rule_figures)) {
        figure <- option$k_rule_figures[[rule]]
        checkOmicGuiFields(figure, c(fdp = "num1", fdp_where = "chr1", power = "num1", source = "chr1"),
            paste0(at("k_rule_figures"), "$", rule))
        if (any(c(figure$fdp, figure$power) < 0 | c(figure$fdp, figure$power) > 1)) {
            omicGuiConfigAbort(paste0(at("k_rule_figures"), "$", rule), "must hold proportions (0 to 1)")
        }
    }
    invisible(option)
}

#' The da group (KIT-106 D1, D8): its fields; the model, editable or shown; the threshold and fold-change
#' inputs; the table's columns, the summary, the results tabs from the frame's set in its order, the assay
#' and table-filter rules, the views, texts and download; the heatmap view with the heatmap tab; the
#' session load with the load section unless a plugin renders it; the run and the source unless a plugin
#' replaces them; the optional fields by type.
#' @noRd
validateOmicGuiDa <- function(da, schema) {
    checkOmicGuiFields(da, schema$da, "da")
    checkOmicGuiFields(da, schema$da_optional, "da", required = FALSE)
    checkOmicGuiFields(da$model, schema$da_model, "da$model")
    checkOmicGuiFields(da$model, if (da$model$editable) schema$da_model_editable else schema$da_model_shown, "da$model")
    # KIT-106 W6: a shown model's text and realised formula are required (trx's, wired at the W2 merge).
    checkOmicGuiFields(da$model, schema$da_model_optional, "da$model", required = !isTRUE(da$model$editable))
    checkOmicGuiFields(da$threshold, schema$da_input, "da$threshold")
    checkOmicGuiFields(da$lfc, c(schema$da_input, schema$da_lfc), "da$lfc")
    checkOmicGuiFields(da$columns, schema$da_columns, "da$columns")
    checkOmicGuiFields(da$columns, schema$da_columns_optional, "da$columns", required = FALSE)
    if (is.null(da$columns[["id"]]) == is.null(da$columns[["id_slot"]])) {
        omicGuiConfigAbort("da$columns", "must name the identifier by its column (`id`) or by the object's slot (`id_slot`)")
    }
    checkOmicGuiFields(da$summary, schema$da_summary, "da$summary")
    tabs <- names(da$tabs)
    if (anyDuplicated(tabs) || !all(c("volcano", "table") %in% tabs) || !all(tabs %in% schema$da_tabs) ||
        is.unsorted(match(tabs, schema$da_tabs))) {
        omicGuiConfigAbort("da$tabs", sprintf("must name the volcano and table tabs (and the heatmap's) once each, in their order (%s)",
            paste(schema$da_tabs, collapse = ", ")))
    }
    checkOmicGuiChoice(da$assays, "da$assays", c("multi", "single", "none"))
    checkOmicGuiChoice(da$table_filter, "da$table_filter", c("all_only", "all_or_combined", "store_page"))
    checkOmicGuiFields(da$views, c(static_volcano = "named_list", table = "named_list"), "da$views")
    checkOmicGuiFields(da$views$static_volcano, schema$da_static_volcano, "da$views$static_volcano")
    # KIT-106 W6: every omic's static volcano names its kernel (trx's wired at the W2 merge).
    checkOmicGuiFields(da$views$static_volcano, c(plot = "ref"), "da$views$static_volcano")
    checkOmicGuiFields(da$views$table, schema$da_table, "da$views$table")
    if ("heatmap" %in% tabs) {
        checkOmicGuiFields(da$views, c(heatmap = "named_list"), "da$views")
        checkOmicGuiFields(da$views$heatmap, schema$da_heatmap, "da$views$heatmap")
        checkOmicGuiFields(da$views$heatmap$tree_cut, schema$da_tree_cut, "da$views$heatmap$tree_cut")
        checkOmicGuiChoice(da$views$heatmap$tree_cut$default, "da$views$heatmap$tree_cut$default",
            c("none", "k_clusters", "height_cutoff", "dynamic"))
    }
    checkOmicGuiFields(da$texts, schema$da_texts, "da$texts")
    checkOmicGuiFields(da$texts, schema$da_texts_optional, "da$texts", required = FALSE)
    checkOmicGuiFields(da$download, schema$da_download, "da$download")
    plugins <- names(da$plugins)
    if (!all(plugins %in% schema$da_plugins)) {
        omicGuiConfigAbort("da$plugins", sprintf("must name only the plugins %s", paste(schema$da_plugins, collapse = ", ")))
    }
    if ("load_session" %in% names(da$testids) && !"load" %in% plugins) {
        checkOmicGuiFields(da, c(session_load = "named_list"), "da")
        checkOmicGuiFields(da$session_load, schema$da_session_load, "da$session_load")
        checkOmicGuiFields(da$session_load, c(assay_omic_type = "chr1", assay_check_first = "lgl1"), "da$session_load",
            required = FALSE)
        checkOmicGuiChoice(da$session_load$assay_rule, "da$session_load$assay_rule", c("setequal", "identical"))
    }
    if (!"run" %in% plugins) {
        checkOmicGuiFields(da, c(run = "named_list"), "da")
        checkOmicGuiFields(da$run, schema$da_run, "da$run")
        checkOmicGuiFields(da$run, c(export_report = "ref"), "da$run", required = FALSE)
    }
    if (!"init" %in% plugins) {
        checkOmicGuiFields(da, c(source = "named_list"), "da")
        checkOmicGuiFields(da$source, schema$da_source, "da$source")
        # KIT-106 W6: the source read is required with the frame's own init (trx's wired at the W2 merge).
        checkOmicGuiFields(da$source, c(read_metadata = "ref"), "da$source")
    }
    invisible(da)
}

#' The enrichment group (KIT-107 D1, D9): its fields and types; the DA source, the identifier column; the three
#' cutoffs (`up`, `down`, `q`), each with its four numbers; the correction's choice among its choices; the taxon;
#' the result services from the tab's set (`gprofiler2`, `clusterprofiler`, `stringdb`); the organism filter, NULL or
#' its fields; and the optional function groups, whole when present.
#' @noRd
validateOmicGuiEnrichment <- function(enrichment, schema) {
    checkOmicGuiFields(enrichment, schema$enrichment, "enrichment")
    # KIT-107 W4: the function groups are the core's (wired at the W2 merge) and required.
    checkOmicGuiFields(enrichment, schema$enrichment_optional, "enrichment")
    checkOmicGuiFields(enrichment$da_source, schema$enrichment_da_source, "enrichment$da_source")
    checkOmicGuiFields(enrichment$id_column, schema$enrichment_id_column, "enrichment$id_column")
    if (!setequal(names(enrichment$cutoffs), schema$enrichment_cutoffs) || length(enrichment$cutoffs) != 3L) {
        omicGuiConfigAbort("enrichment$cutoffs", sprintf("must name exactly the cutoffs %s",
            paste(schema$enrichment_cutoffs, collapse = ", ")))
    }
    for (cutoff in schema$enrichment_cutoffs) {
        checkOmicGuiFields(enrichment$cutoffs[[cutoff]], schema$enrichment_cutoff, paste0("enrichment$cutoffs$", cutoff))
    }
    checkOmicGuiFields(enrichment$correction, schema$enrichment_correction, "enrichment$correction")
    checkOmicGuiChoice(enrichment$correction$selected, "enrichment$correction$selected",
        unname(enrichment$correction$choices))
    checkOmicGuiFields(enrichment$taxon, schema$enrichment_taxon, "enrichment$taxon")
    services <- enrichment$services
    if (anyDuplicated(services) || !all(services %in% schema$enrichment_services)) {
        omicGuiConfigAbort("enrichment$services", sprintf("must name result services among %s, once each",
            paste(schema$enrichment_services, collapse = ", ")))
    }
    if (!is.null(enrichment$organism_filter)) {
        checkOmicGuiFields(enrichment$organism_filter, schema$enrichment_organism_filter, "enrichment$organism_filter")
        checkOmicGuiFields(enrichment$organism_filter, schema$enrichment_organism_filter_optional,
            "enrichment$organism_filter")
    }
    for (group in names(schema$enrichment_optional)) {
        if (!is.null(enrichment[[group]])) {
            checkOmicGuiFields(enrichment[[group]], schema[[paste0("enrichment_", group)]], paste0("enrichment$", group))
        }
    }
    invisible(enrichment)
}

#' The summary group (KIT-108 D1, D4, D9): its fields and types; the module id is the shell's summary tab's (when the
#' shell has one) and the test id the shell's; the final state is named by its `context`, or by `states` with their `read_ids` reader, not both; the copy's
#' required artifacts are NULL or names; the integration's operation is NULL or a `*_integration_export` operation, its
#' gate NULL or an export with its message, its file name a string or a function; a report template given as a file
#' name has the label of the one-template status; each export step's panel strings and references.
#' @noRd
validateOmicGuiSummary <- function(summary, shell, schema) {
    checkOmicGuiFields(summary, schema$summary, "summary")
    checkOmicGuiFields(summary, schema$summary_optional, "summary", required = FALSE)
    # A shell configured without its summary tab (a shell test's reduced tabs) has no module id to match.
    tab <- Filter(function(tab) identical(tab$status_key, "session_summary"), shell$tabs)
    if (length(tab) == 1L && !identical(summary$module_id, tab[[1L]]$module_id)) {
        omicGuiConfigAbort("summary$module_id", "must be the module id of the shell's summary tab")
    }
    if (!identical(summary$testid, shell$testid)) omicGuiConfigAbort("summary$testid", "must be the shell's test id")
    final <- summary$final_state
    checkOmicGuiFields(final, schema$summary_final_state, "summary$final_state", required = FALSE)
    if (!is.null(final$context) && !is.null(final$states)) {
        omicGuiConfigAbort("summary$final_state", "must name the final state by `context` or by `states`, not both")
    }
    if (!is.null(final$states)) checkOmicGuiFields(final, c(read_ids = "ref"), "summary$final_state")
    if (is.null(final$context) && is.null(final$states)) {
        omicGuiConfigAbort("summary$final_state", "must name the final state by `context` or by `states`")
    }
    copy <- summary$copy
    checkOmicGuiFields(copy, schema$summary_copy, "summary$copy")
    checkOmicGuiFields(copy, schema$summary_copy_optional, "summary$copy", required = FALSE)
    if (!"required" %in% names(copy)) omicGuiConfigAbort("summary$copy$required", "is missing")
    if (!is.null(copy$required) && !omicGuiFieldHasType(copy$required, "chr")) {
        omicGuiConfigAbort("summary$copy$required", "must be NULL or the required artifacts' names")
    }
    integration <- summary$integration
    checkOmicGuiFields(integration, schema$summary_integration, "summary$integration")
    checkOmicGuiFields(integration, schema$summary_integration_optional, "summary$integration", required = FALSE)
    if (!"operation" %in% names(integration)) omicGuiConfigAbort("summary$integration$operation", "is missing")
    operation <- integration$operation
    if (!is.null(operation) && !(omicGuiFieldHasType(operation, "chr1") && endsWith(operation, "_integration_export"))) {
        omicGuiConfigAbort("summary$integration$operation", "must be NULL or an `*_integration_export` operation")
    }
    if (!is.null(integration$requires)) {
        checkOmicGuiFields(integration, c(missing_message = "chr1"), "summary$integration")
    }
    if (!"filename" %in% names(integration)) omicGuiConfigAbort("summary$integration$filename", "is missing")
    if ("filename" %in% names(integration) && !(omicGuiFieldHasType(integration$filename, "chr1") ||
        is.function(integration$filename))) {
        omicGuiConfigAbort("summary$integration$filename", "must be a file name or a function")
    }
    report <- summary$report
    checkOmicGuiFields(report, schema$summary_report, "summary$report")
    checkOmicGuiFields(report, schema$summary_report_optional, "summary$report", required = FALSE)
    if (!"template" %in% names(report)) omicGuiConfigAbort("summary$report$template", "is missing")
    if ("template" %in% names(report)) {
        if (omicGuiFieldHasType(report$template, "chr1")) {
            checkOmicGuiFields(report, c(label = "chr1"), "summary$report")
        } else if (!is.function(report$template)) {
            omicGuiConfigAbort("summary$report$template", "must be a file name or a function")
        }
    }
    steps <- summary$export_steps
    if (!"export_steps" %in% names(summary)) omicGuiConfigAbort("summary$export_steps", "is missing")
    if (!is.list(steps) || is.object(steps)) omicGuiConfigAbort("summary$export_steps", "must be a list of steps")
    for (i in seq_along(steps)) {
        path <- sprintf("summary$export_steps[[%d]]", i)
        step <- steps[[i]]
        checkOmicGuiFields(step, schema$summary_export_step, path)
        checkOmicGuiFields(step, schema$summary_export_step_optional, path)
        checkOmicGuiFields(step$choice, schema$summary_export_choice, paste0(path, "$choice"))
        checkOmicGuiFields(step$choice, schema$summary_export_choice_optional, paste0(path, "$choice"))
        checkOmicGuiFields(step$button, schema$summary_export_button, paste0(path, "$button"))
        checkOmicGuiFields(step$status, schema$summary_export_status, paste0(path, "$status"))
        checkOmicGuiFields(step$status, schema$summary_export_status_optional, paste0(path, "$status"))
        checkOmicGuiFields(step$figure, schema$summary_export_figure, paste0(path, "$figure"))
        for (j in seq_along(step$publications)) {
            checkOmicGuiFields(step$publications[[j]], schema$summary_export_publication,
                sprintf("%s$publications[[%d]]", path, j))
        }
        checkOmicGuiFields(step$error, schema$summary_export_error, paste0(path, "$error"))
    }
    invisible(summary)
}

#' Check the fields a schema names: present (when required) and of their type.
#' @param values The list holding the fields.
#' @param schema Named character vector: field name to type code.
#' @param path The fields' path in the configuration, for the refusal; NULL at the top level.
#' @param required Whether an absent field is refused.
#' @noRd
checkOmicGuiFields <- function(values, schema, path = NULL, required = TRUE) {
    for (field in names(schema)) {
        where <- if (is.null(path)) field else paste0(path, "$", field)
        if (!field %in% names(values)) {
            if (required) omicGuiConfigAbort(where, "is missing")
            next
        }
        if (!omicGuiFieldHasType(values[[field]], schema[[field]])) {
            omicGuiConfigAbort(where, sprintf("must be of type %s", schema[[field]]))
        }
    }
    invisible(values)
}

#' Whether a value has a schema type: `chr1` one non-empty string; `chr1_or_null`; `chr` a
#' non-empty character vector; `chr0` a character vector; `named_chr`; `named_chr0` a possibly
#' empty named character vector; `by_route` one non-empty string, or non-empty strings named by
#' route; `route_map_or_null` NULL, or non-empty strings named by route (a route not named has no
#' value); `lgl1`; `num1` one number; `num1_or_null`; `int1` one positive whole number; `ref` a function;
#' `ref_or_null`; `refs` a named list of functions; `list` a non-empty list; `named_list`;
#' `named_list_or_null`.
#' @noRd
omicGuiFieldHasType <- function(value, type) {
    is_chr <- is.character(value) && !anyNA(value)
    named <- length(value) == 0L ||
        (!is.null(names(value)) && all(nzchar(names(value))) && !anyDuplicated(names(value)))
    is_chr1 <- is_chr && length(value) == 1L && nzchar(value)
    is_num1 <- is.numeric(value) && length(value) == 1L && !is.na(value)
    isTRUE(switch(type,
        num1 = is_num1,
        num1_or_null = is.null(value) || is_num1,
        int1 = is_num1 && value >= 1 && value == round(value),
        chr1 = is_chr1,
        chr1_or_null = is.null(value) || is_chr1,
        chr1_or_ref_or_null = is.null(value) || is_chr1 || is.function(value),
        chr = is_chr && length(value) >= 1L,
        chr0 = is_chr,
        named_chr = is_chr && length(value) >= 1L && named,
        named_chr0 = is_chr && named,
        by_route = is_chr && length(value) >= 1L && all(nzchar(value)) &&
            (length(value) == 1L || named),
        route_map_or_null = is.null(value) ||
            (is_chr && length(value) >= 1L && all(nzchar(value)) && !is.null(names(value)) && named),
        lgl1 = is.logical(value) && length(value) == 1L && !is.na(value),
        ref = is.function(value),
        ref_or_null = is.null(value) || is.function(value),
        refs = is.list(value) && length(value) >= 1L && named &&
            all(vapply(value, is.function, logical(1))),
        list = is.list(value) && length(value) >= 1L,
        named_list = is.list(value) && named,
        named_list_or_null = is.null(value) || (is.list(value) && named)))
}

#' @noRd
checkOmicGuiChoice <- function(value, field, choices) {
    if (!value %in% choices) {
        omicGuiConfigAbort(field,
            sprintf("must be one of %s, not %s", paste(choices, collapse = ", "), value))
    }
    invisible(value)
}

#' @noRd
omicGuiConfigAbort <- function(field, problem) {
    cli::cli_abort("omic GUI configuration field {.field {field}} {problem}.",
        class = "omic_gui_config_invalid")
}

# ------------------------------------------------------------------------------------------------
# Shared pieces: tabs, bus, binding and the feature-table panel's status displays.
# ------------------------------------------------------------------------------------------------

#' One workflow tab. Today the stepper icon is the tab icon for every omic.
#' @param component The tab's component id: today's public module is `mod_<component>_ui/server`.
#' @noRd
omicGuiTab <- function(status_key, value, module_id, component, title, icon, testid, step_name) {
    list(status_key = status_key, value = value, module_id = module_id, component = component,
        title = title, icon = icon, testid = testid, step_name = step_name, step_icon = icon)
}

#' The six tabs of lipid, metab and trx (mod_lipidomics.R:59-116, stepper :271-276; the same lines
#' in mod_metabolomics.R; mod_transcriptomics.R:220-225, stepper :265-269). Titles keep their
#' leading space: the shell writes them after the icon (`tab_style = "icon_in_title"`).
#' @noRd
omicGuiFeatureTabs <- function(prefix, da_title) {
    tab <- function(status_key, value, component, title, icon, step_name) {
        omicGuiTab(status_key, value, value, paste0(prefix, "_", component), title, icon,
            paste0(prefix, "-tab-", value), step_name)
    }
    list(
        tab("setup_import", "import", "import", " Setup & Import", "file-import", "Import"),
        tab("design_matrix", "design", "design", " Design Matrix", "th", "Design"),
        tab("quality_control", "qc", "qc", " Quality Control", "check-double", "QC"),
        tab("normalization", "norm", "norm", " Normalization", "balance-scale", "Normalize"),
        tab("differential_analysis", "de", "da", da_title, "chart-line", "DE"),
        tab("session_summary", "summary", "summary", " Summary & Export", "download", "Summary")
    )
}

#' The initial bus of lipid, metab and trx, in today's order (mod_lipidomics.R:139-176,
#' mod_metabolomics.R:139-176, mod_transcriptomics.R:232-235). `tab_status` is derived from the
#' tabs.
#' @param assay_slot The omic's assay slot (`lipid_assay_data`, `metabolite_assay_data`); NULL
#'   for trx.
#' @noRd
omicGuiFeatureBus <- function(data_type, assay_slot = NULL) {
    bus <- list(data_tbl = NULL, data_format = NULL, data_type = data_type, column_mapping = NULL,
        config_list = list(), state_manager = NULL)
    if (!is.null(assay_slot)) bus[assay_slot] <- list(NULL)
    c(bus, list(contrasts = list(), contrasts_tbl = NULL, tab_status = NULL,
        state_update_trigger = NULL, processing_log = list()))
}

#' The shell of lipid, metab and trx: the header `paste0(" ", label, " Workflow")` after
#' `header_icon`, the pill tabset `<omic>_tabs` (mod_lipidomics.R:43-55;
#' mod_transcriptomics.R:217-219) and the bus registered with `registerWorkflowOmicsBus`
#' (mod_lipidomics.R:180, mod_transcriptomics.R:246).
#' @noRd
omicGuiFeatureShell <- function(omic_type, prefix, label, header_icon, route, da_title, bus) {
    list(omic_type = omic_type, prefix = prefix, testid = prefix, label = label,
        label_lower = tolower(label), header_icon = header_icon,
        tabset_id = paste0(omic_type, "_tabs"), tab_style = "icon_in_title", routes = route,
        tabs = omicGuiFeatureTabs(prefix, da_title), bus = bus, bus_registration = "omics")
}

#' The omics runner's binding group (`installWorkflowOmicsActions`, mod_lipidomics.R:196; the owner
#' `createWorkflowSessionOwner`, :194).
#' @noRd
omicGuiFeatureBinding <- function(binding, read_state, read_metadata) {
    list(binding = binding, runner = "omics_actions", owner = createWorkflowSessionOwner,
        read_state = read_state, read_metadata = read_metadata)
}

#' The feature-table panel's column status displays, the same for lipid, metab and trx: the
#' builder's success and empty modes, and whether a custom name matches case-insensitively (lipid
#' `output`:384, :398-399, :413-415, :429-431; metab `display`:63-195; now
#' mod_omic_import_feature_table.R).
#' @noRd
omicGuiFeatureStatusModes <- function() {
    list(id = c(success = "unique_ids"),
        annotation = c(success = "ok", empty = "optional"),
        id_custom = c(success = "found_unique_ids", empty = "prompt", case = "insensitive"),
        annotation_custom = c(success = "found", empty = "optional", case = "insensitive"))
}

#' The design tab of lipid, metab and trx (at 3c657c14): one frame (mod_lipid_design_ui.R; the same
#' in mod_metab_design.R:43-111) and one builder whose UI differs only in the omic word
#' (mod_lipid_design_builder_ui.R:218, :246; mod_metab_design_builder_ui.R:245, :273), with the
#' same texts in their servers. Trx rendered the lipid frame and builder under its profile. Since
#' KIT-103 W4/W5 they are the shared frame and builder (mod_omic_design.R, mod_omic_design_builder.R).
#' @param builder_word The omic word in the builder UI (trx showed lipid's until P-12, 2026-11-11).
#' @noRd
omicGuiFeatureDesign <- function(prefix, label_lower, builder_word, id_mapping_key) {
    list(module_id = "design",
        # mod_lipid_design_validation.R:5; metab inline (mod_metab_design.R:235-238)
        formula_default = "~ 0 + group", label_lower = label_lower,
        data_adapter = "wide_assay_list", id_mapping_key = id_mapping_key, commit = "shared",
        builder_ui = "shared",
        # mod_lipid_design_ui.R:85, :63, :73-74 (mod_metab_design.R:105, :83, :93-94)
        no_data_text = paste(" Please complete the 'Import' step first. The builder will appear",
            "once data is available."),
        preview_text = "This section shows the design matrix and contrasts saved to the workflow.",
        preview_assays = TRUE,
        # mod_lipid_design_server.R:47, :50, :57 (mod_metab_design.R:430, :441); trx has no import
        testids = c(import_dir = paste0(prefix, "-design-import-dir"),
            import_confirm = paste0(prefix, "-design-import-confirm")),
        multi_assay_note = paste0("<b>Multi-Assay Note:</b> This ", builder_word,
            " workflow supports multiple assays "),
        matrix_heading = paste0("Current Design Matrix (All Samples) - ",
            paste0(toupper(substr(builder_word, 1L, 1L)), substring(builder_word, 2L))),
        # mod_lipid_design_builder_state_helpers.R:621, :626; mod_lipid_design_builder_display_helpers.R
        # :41, :43, :132 (mod_metab_design_builder_server.R:783, :859, :493, :495)
        builder_texts = c(save_success = "Design saved successfully!",
            no_assignment = "No samples have been assigned to groups.",
            contrast_group_prefix = paste("Note: Contrasts will use 'group' prefix based on",
                "formula: ~ 0 + group"),
            contrast_as_is = "Note: Contrasts will use group names as-is",
            replicate_suffix = "selected samples:"))
}

#' The import-existing-design values of lipid and metab (mod_lipid_design_import_helpers.R;
#' mod_metab_design_import_helpers.R at 3c657c14): the column-mapping inference when
#' `column_mapping.json` is absent, the QC baseline, whether the import runs the DA preflight
#' (`preflight`: both since P-11, an invalid imported design is refused before anything is written),
#' whether the preflight's normalised contrasts replace the file's (`preflight_contrasts`: lipid, as
#' before; metab keeps the file's, so no valid metab import changes), and the state's description. The output `import_dir_path` waits for a path (`req`). The reader is
#' the core's (`readOmicsDesignImport`, KIT-103 W2). `record_design_file`: whether the provenance
#' records the imported design file (lipid) or the owner's source directory (metab, P-14).
#' `check_files_first`: whether the files are checked before the modal closes (metab, P-11).
#' @noRd
omicGuiFeatureDesignImport <- function(annotation_candidates, qc_baseline, preflight, class,
                                       record_design_file, check_files_first, preflight_contrasts) {
    list(plugin = "feature_table", import_dir_path = "req",
        id_candidates = list(id = c("Peak ID", "Alignment ID", "Alignment.ID", "peak id",
            "alignment id"), annotation = annotation_candidates),
        qc_baseline = qc_baseline, preflight = preflight, preflight_contrasts = preflight_contrasts,
        state_description = paste(class, "S4 object created from imported design"),
        reader = readOmicsDesignImport, record_design_file = record_design_file,
        check_files_first = check_files_first)
}

#' One QC sub-tab: its module id, title and component (`mod_<component>_ui/server`), the routes it
#' shows on, and whether its server takes the experiment paths (proteomics' protein specs).
#' @noRd
omicGuiQcTab <- function(id, title, component, routes, needs_paths = FALSE) {
    list(id = id, title = title, component = component, routes = routes, needs_paths = needs_paths)
}

#' The QC sub-tabs of lipid and metab: one unnamed group of the four shared sub-tabs, under today's
#' module ids and titles (mod_lipid_qc.R:58-84, :120-125, the titles at mod_lipid_qc_intensity.R:41,
#' mod_lipid_qc_duplicates.R:42, mod_lipid_qc_itsd.R:41, mod_lipid_qc_s4.R:43; the same in
#' mod_metab_qc.R:53-80, :116-121 and the metab sub-tab UIs).
#' @noRd
omicGuiFeatureQcGroups <- function(route) {
    list(list(id = NULL, title = NULL, component = NULL, tabs = list(
        omicGuiQcTab("intensity", "Intensity Filter", "omic_qc_intensity", route),
        omicGuiQcTab("duplicates", "Duplicate Resolution", "omic_qc_duplicates", route),
        omicGuiQcTab("itsd", "Internal Standards", "omic_qc_itsd", route),
        omicGuiQcTab("s4_finalize", "Finalize QC", "omic_qc_finalize", route))))
}

#' The texts of lipid's and metab's QC sub-tabs, which differ only in the feature noun and the class
#' (lipid at b423b764: mod_lipid_qc_intensity.R:46, :47, :70; mod_lipid_qc_intensity_server_helpers.R
#' :22, :150, :202, :219, :235, :263; mod_lipid_qc_duplicates.R:48;
#' mod_lipid_qc_duplicates_server_helpers.R:26; mod_lipid_qc_itsd.R:57; mod_lipid_qc_s4.R:191, :203,
#' :251, :325, :366, :372. Metab: mod_metab_qc_intensity.R:46, :47, :70, :131, :364, :256, :265,
#' :310, :429; mod_metab_qc_duplicates_ui.R:24; mod_metab_qc_duplicates_ui_helpers.R:57;
#' mod_metab_qc_itsd_ui.R:25; mod_metab_qc_s4_server_helpers.R:60, :72, :194, :292, :338, :344).
#' `%d` stands for today's counts, in today's order.
#' @param noun The feature noun (`lipid`, `metabolite`).
#' @noRd
omicGuiFeatureQcTexts <- function(noun, class) {
    capitalised <- paste0(toupper(substr(noun, 1L, 1L)), substring(noun, 2L))
    plural <- paste0(noun, "s")
    c(intensity_heading = paste(capitalised, "Intensity Filter Parameters"),
        intensity_help = sprintf(paste("Filter %s based on intensity thresholds. %ss with low intensities",
            "across many samples will be removed."), plural, capitalised),
        intensity_proportion_help = sprintf(paste("Remove %s where this proportion of samples are below the",
            "intensity threshold (default: 0.5)"), plural),
        intensity_original_label = sprintf("Original %s:", plural),
        intensity_working = sprintf("Applying %s intensity filter...", noun),
        intensity_result_heading = paste(capitalised, "Intensity Filter Applied Successfully"),
        intensity_total = sprintf("Total: %%d -> %%d %s (removed %%d)", plural),
        intensity_success = sprintf("Intensity filter applied: %%d %s retained", plural),
        intensity_error = sprintf("Error applying %s intensity filter:", noun),
        duplicates_help = sprintf(paste("Identify and resolve duplicate %s IDs within each assay. Duplicates",
            "are resolved by keeping the feature with the highest average intensity across samples."), noun),
        duplicates_none = sprintf("All %s IDs are unique. No resolution needed.", noun),
        itsd_pattern_help = sprintf(paste("Regular expression to identify internal standards in %s IDs.",
            "Leave empty to use pattern from S4 object."), noun),
        finalize_total_label = sprintf("Total %ss:", capitalised),
        finalize_id_label = sprintf("%s ID Column:", capitalised),
        finalize_count_column = paste0(capitalised, "s"),
        finalize_total = sprintf("Total %s retained: %%d", plural),
        finalize_no_assays = sprintf("Current %s object has no assay data to finalize.", class),
        finalize_no_design = sprintf("Current %s object has no design matrix rows to finalize.", class))
}

#' The normalisation tab of lipid and metab at 9a98e765 (design 4.7): one layout, written twice
#' (lipid generated by buildLipidNormOptionsControlPanel and buildLipidNormQcTabsetPanel,
#' mod_lipid_norm_ui_helpers.R; metab by hand, mod_metab_norm_ui.R:25-658), differing in the RUV
#' noun, the test ids and the per-omic values of design D7. The option sections, in order: the
#' aesthetics (:8-25), ITSD with its aggregation (:29-52, default median), log2 with its offset
#' (:57-65), the between-sample methods (:70-85), RUV-III with its modes, percentage and k
#' (:90-204), the actions (run, reset; :225-246) and the export button (:251-261). The QC tabs, the
#' stage captions, the assay rows and the RUV tab's texts: buildLipidNormQcTabsetPanel (:285-307;
#' metab :290-658). The managed run and correlation: the core's (utils_workflow_omics_actions.R:254,
#' :347), called by both tabs (mod_lipid_norm_workflow_helpers.R:150, :738;
#' mod_metab_norm_binding_helpers.R:628, mod_metab_norm_session_helpers.R:19). Both start the
#' pre-normalisation QC when the tab is opened (mod_lipid_norm_workflow_helpers.R:871-893;
#' mod_metab_norm_entry_helpers.R:84-96).
#' @param prefix The test-id prefix (`lipid`, `metab`).
#' @param noun The RUV-III control noun (`Lipids`, `Metabolites`).
#' @param export The export's `prefix`, `counts_field` and `samples` rule.
#' @param itsd_table The ITSD selection table builder, `function(assay_data, id_col, annotation_cols)`.
#' @param skip_correlation_flags When skip-correlation sets the tab's flags: before or after the
#'   store commit.
#' @param final_pca The managed final PCA's placeholder sizes.
#' @param managed_correlation_notify Whether a failed managed correlation is also notified as an error (KIT-105 W5).
#' @noRd
omicGuiFeatureNorm <- function(prefix, noun, export, itsd_table, skip_correlation_flags, final_pca, managed_correlation_notify) {
    list(module_id = "norm",
        panels = list(aesthetics = list(),
            itsd = list(aggregation = c("Median (recommended)" = "median", "Mean" = "mean", "Sum" = "sum")),
            # KIT-172: plain log2 (the offset is no longer offered; a record's offset replays).
            log2 = list(offset = 0),
            between_sample = list(methods = c("Cyclic Loess" = "cyclicloess", "Quantile" = "quantile",
                "Scale (Median Absolute Values)" = "scale", "None" = "none"), default = "cyclicloess"),
            ruv_iii = list(modes = c("Automatic (recommended)" = "automatic", "Manual (advanced users)" = "manual",
                "Skip RUV" = "skip"), default = "automatic", percentage = 5, k = 2, noun = noun),
            actions = list(), export = list(button = "export_session")),
        plot_tabs = c("itsd", "ruv", "pca", "density", "rle", "correlation", "correlation_filter", "log"),
        stage_labels = c("Post-Filtering", "Post-Normalization", "RUV-Corrected"), assay_slots = 2L,
        ruv_heading = "RUV Parameter Optimization Results",
        ruv_help = paste("Canonical correlation plots and optimization summary per assay. Each assay is optimized",
            "independently."),
        # The in-process twins and the export's collector and saver: the core's (KIT-105 W4; required since W6).
        run = list(managed = runWorkflowOmicsManagedNormalisation, in_process = runWorkflowOmicsInProcessNormalisation),
        correlation = list(managed = runWorkflowOmicsManagedCorrelation, in_process = runWorkflowOmicsInProcessCorrelation),
        pre_qc_on_entry = TRUE,
        export = c(export, list(collect = collectWorkflowOmicsNormSession, save = saveWorkflowOmicsNormSessionArtifacts)),
        itsd_table = itsd_table,
        testids = c(method = paste0(prefix, "-norm-method"), ruv_mode = paste0(prefix, "-norm-ruv-mode"),
            run = paste0(prefix, "-norm-run"), export = paste0(prefix, "-norm-export-session")),
        skip_correlation_flags = skip_correlation_flags, final_pca = final_pca,
        managed_correlation_notify = managed_correlation_notify,
        ruv_design = omicGuiRuvDesign(prefix, list(features = omicGuiRuvTableFeature(tolower(noun)))))
}

# ------------------------------------------------------------------------------------------------
# KIT-129 W2: the RUV design presets (ticket sections 1 and 2; the user's gate of 2026-11-18). One
# catalogue of designs (tier, label, description, tooltip, what the design does, its k rules) and one table of
# measured figures and notes per route family. The panel section (mod_omic_norm_panels.R) reads them;
# the study parameters, the record and the report templates are to read the same fields.
#
# The figures are KIT-124's: the mean false discovery proportion at a nominal 0.05 in the worst core
# cell with true effects (`fdp`, with the cell in words in `fdp_where`), and the power in the centre
# cell (`power`; 6 per group, 20% changed, one factor); `single_scenario` marks a figure measured in one cell. Each row's
# `source` gives the simulator, arm and cells; it is not shown. To update a figure when a simulation
# finishes, change the row's numbers and `source`, set `status = "measured"`, and drop the rule from
# `k_rules_pending` when a k rule's own figure arrives. The row's figure is the one measured with the pre-selected
# k rule (`k_rule`: the estimated number wherever it is defined, the user's decision M15 of 2026-11-18); a figure
# measured with the other rule is in `k_rule_figures`, by rule, with its own `source`. Nothing is computed from
# these numbers.
# ------------------------------------------------------------------------------------------------

#' Whether the RUV design section is shown (the gate of KIT-129 W2). Off until the switch commit: the
#' normalisation tab is then today's, byte for byte. `options(MultiScholaR.ruv_design_panel = TRUE)` shows it.
#' @return One logical.
#' @noRd
omicGuiRuvDesignEnabled <- function() {
    isTRUE(getOption("MultiScholaR.ruv_design_panel", FALSE))
}

#' The log2 step's texts of lipidomics and metabolomics (KIT-172; proposed wording, implemented as proposed, the owner
#' edits before release, VAL-41): the section's sentence, the recorded-offset sentence of a project whose record used an
#' offset (`%s` the offset), the step-log line of the values set to missing (`%s` the counts per assay), the
#' interactive volcano's sentence where missing cells are drawn, and the warning shown when `scale` runs on log2 values
#' below 0 (`%s` the counts per assay).
#' @noRd
omicGuiLog2Texts <- function() {
    c(help = "log2 of each value; zeros are read as not detected at import",
        recorded_offset = paste("This project's record used log2(x + %s) for this step; its results replay with that offset.",
            "A new run uses plain log2."),
        nonpositive_log = "Log2 step: values set to missing (0 or below; not finite): %s",
        volcano_missing = paste("Samples with no value for a feature are drawn at 0 in this plot's sample panel; the",
            "interactive plot cannot draw missing values."),
        scale_warning = paste("Scale normalisation rescales each sample so that the median of its absolute log2 values is the",
            "same in every sample. Because it divides log2 values instead of shifting them, its result depends on the unit of",
            "your table (for example \u00b5M against nM). In this data, %s values reaching the log2 step are below 1, so",
            "their log2 values are negative, and taking absolute values treats a log2 value of \u22122 like +2. Cyclic Loess",
            "and Quantile do not depend on the unit."))
}

#' The proteomics label-free and TMT log2 texts (KIT-177; proposed wording, implemented as proposed, the owner edits
#' before release, VAL-41): the reopen sentence of a design saved by an earlier version, the notification and step-log
#' lines of the values set to missing at the design step (`%d` the counts), the refusal of a negative or infinite
#' quantity (`%s` the run and the values), and the warning shown when `scale` runs on log2 values below 0 (`%d` the count).
#' @noRd
omicGuiProtLog2Texts <- function() {
    c(reopen = paste("This design was saved by an earlier version, which added 1 to every value before taking log2.",
            "Your saved results are unchanged. If you save the design again, values are log2-transformed without adding 1,",
            "and later results will differ slightly, most for small values."),
        notice_zero = "%d values of 0",
        notice_nan = "%d values that were not a number",
        notice_both = "%d values of 0 and %d values that were not a number",
        notice_tail = "were read as not measured (missing). They are not used in the analysis.",
        log_zero = "Values equal to 0 set to missing: %d",
        log_nan = "Values that were not a number set to missing: %d",
        refusal = paste("Run '%s' contains negative or infinite values (for example %s). Quantities must be zero or",
            "positive; correct the file and import it again."),
        scale_warning = paste("Scale normalisation divides each sample by the median of its absolute log2 values. In this",
            "data, %d values are below 1, so their log2 values are negative, and taking absolute values treats a log2",
            "value of \u22122 like +2. The result then depends on the unit of your table. Cyclic Loess and Quantile do not",
            "have this problem."),
        scale_log = "Scale normalisation: log2 values below 0 in the input: %d")
}

#' An omic's RUV design group: the gate, the tier labels, the k rules and its tables.
#' @param prefix The omic's test-id prefix.
#' @param tables The preset tables, named; the first is the one the tab opens with.
#' @param route_tables The workflow types that take another table, named by type; NULL for none.
#' @param conditions The conditions the omic's simulations cover, shown under the options; NULL for none.
#' @return The norm group's `ruv_design`.
#' @noRd
omicGuiRuvDesign <- function(prefix, tables, route_tables = NULL, conditions = NULL) {
    design <- list(enabled = omicGuiRuvDesignEnabled(), testid = paste0(prefix, "-norm-ruv-design"),
        tiers = c(recommended = "Recommended", valid_with_conditions = "Valid with conditions",
            legacy = "Legacy", none = "None"),
        k_rules = omicGuiRuvKRules(), tables = tables)
    design$route_tables <- route_tables
    design$conditions <- conditions
    design
}

#' The three rules for the number of factors (k), with the texts shown behind the choice's explanation (ticket
#' section 2; kit-124/design-review-20261118.md section 5: what each reads, for and against, the measured counts).
#' The ids are the record's `ruv_k_rule` values. Wording: docs/qa/workflow-execution/kit/kit-129/panel-ui/
#' wording-review.md (plain words first, the method's name second).
#' @return A named list of `label`, `reads`, `pros`, `cons` and `diagnostic`.
#' @noRd
omicGuiRuvKRules <- function() {
    list(
        num_sv_be = list(label = "Estimated number",
            reads = paste("The software estimates how many unwanted factors the data contain. It first sets the",
                "differences between the groups aside, then counts the patterns (principal components) in what is",
                "left that explain more of the variation than they do in 20 random shuffles (permutations) of the",
                "data. Method: sva::num.sv, method \"be\"."),
            pros = paste("It counts the factors directly, needs no control features and gives the same answer every",
                "time (fixed random seed). It is the pre-selected rule of every approach that estimates k, and the",
                "rule their main figures were measured with. In simulation it returned on average 1.00 where there was truly 1 factor, 1.94 where there",
                "were 2 and 3.00 where there were 3."),
            cons = paste("It comes from surrogate variable analysis (the sva package), not from the RUV publications.",
                "Because the group differences are set aside first, a factor that is tied to the groups can be",
                "under-counted. It is a rough estimate with only 3 samples per group. It needs the optional R",
                "package sva."),
            diagnostic = paste("For each pattern, the share of the remaining variation it explains in your data",
                "beside the shuffled data, with the estimated number marked (RUV QC tab).")),
        sweep_bounded = list(label = "Separation rule",
            reads = paste("Uses the canonical-correlation plot, which shows for each k how strongly the first k",
                "patterns in the data are related to the groups: one curve for all features, one for the control",
                "features only. The rule takes the k at which the two curves are furthest apart, with an upper",
                "limit on k in small experiments."),
            pros = paste("It is this package's established diagnostic, links k to the controls you chose and needs",
                "no extra package. It is the rule the figures of the approach of earlier versions were measured",
                "with."),
            cons = paste("In simulation it returned on average 2.9 where there was truly 1 factor, and 1.55 where",
                "there were 3. It shows how free the controls are of group differences, not how many factors",
                "there are. It cannot be used when all features are the controls."),
            diagnostic = "The canonical-correlation plot and its table (RUV QC tab)."),
        user = list(label = "Your own number",
            reads = "You enter the number of factors.",
            pros = paste("The RUV publications ask for judgement guided by diagnostic plots rather than one",
                "automatic rule. Every approach accepts a number."),
            cons = paste("The figures shown for each approach were measured with the estimated number or the",
                "separation rule, not with your number."),
            diagnostic = paste("The estimated number's plot and the canonical-correlation plot. When all features",
                "are the controls, the canonical-correlation plot is drawn with controls chosen from the data, and",
                "says so.")))
}

#' The designs' common part (ticket section 1's table): tier, label (what is done, in plain words), description
#' (the visible second line: how, with the method's name), tooltip (the exact statement, extra detail only),
#' replicate structure, control rule, test form (the record's values) and k rules (section 2: the estimated number
#' for every design with a k, and pre-selected there (decision M15); the separation rule where a separate control
#' set exists; `k_rules_pending` are the combinations without a measured figure of their own: none since the new
#' arms were scored). The two legacy designs' descriptions end with what their test reads (`omicGuiRuvTable`).
#' @param features The omic's features in words (`proteins`, `genes`, `lipids`, `metabolites`).
#' @param method The RUV method's name by kind of design: `replicates` and `controls`.
#' @return A named list of options without their figures.
#' @noRd
omicGuiRuvDesignCatalogue <- function(features, method) {
    in_model <- paste("The estimated factors are included as covariates in the differential model, which is fitted",
        "to the uncorrected values; the corrected values are used for plots and export only.")
    earlier <- sprintf(paste("The %s that differ least between the groups in a first test serve as controls; an",
        "automatic search chooses how many of them, and the number of factors."), features)
    list(
        groups_all_features = list(id = "groups_all_features", tier = "recommended",
            label = sprintf("Estimate from all %s, within your groups", features),
            description = sprintf(paste("Samples of the same group are treated as replicates of each other, and all",
                "%s are used to estimate the unwanted variation (%s)."), features, method[["replicates"]]),
            tooltip = sprintf("%s with the experimental groups as replicate sets and all %s as negative controls. %s",
                method[["replicates"]], features, in_model),
            replicate_structure = "group", control_rule = "all_features", test_form = "factors_in_model",
            k_rules = c("num_sv_be", "user"), k_rule = "num_sv_be"),
        technical_replicates = list(id = "technical_replicates", tier = "valid_with_conditions",
            label = sprintf("Estimate from all %s, within your technical replicates", features),
            description = sprintf(paste("The technical replicates declared in the experimental design are compared",
                "with each other, and all %s are used to estimate the unwanted variation (%s)."), features,
                method[["replicates"]]),
            tooltip = sprintf(paste("%s with the technical replicates declared in the design as replicate sets and all",
                "%s as negative controls. %s"), method[["replicates"]], features, in_model),
            replicate_structure = "technical_replicates", control_rule = "all_features",
            test_form = "factors_in_model", k_rules = c("num_sv_be", "user"), k_rule = "num_sv_be"),
        supplied_controls = list(id = "supplied_controls", tier = "valid_with_conditions",
            label = sprintf("Estimate from your own list of control %s", features),
            description = sprintf(paste("You supply %s that you know do not differ between the groups (negative",
                "controls); the unwanted variation is estimated from them (%s)."), features, method[["controls"]]),
            tooltip = sprintf(paste("%s with a list of %s you know not to change between the groups as negative",
                "controls. %s"), method[["controls"]], features, in_model),
            replicate_structure = "group", control_rule = "supplied_list", test_form = "factors_in_model",
            k_rules = c("num_sv_be", "sweep_bounded", "user"), k_rule = "num_sv_be"),
        first_pass_nulls = list(id = "first_pass_nulls", tier = "valid_with_conditions",
            label = sprintf("Estimate from control %s chosen from the data", features),
            description = sprintf(paste("A first test of the groups is run; the 25%% most significant %s are set",
                "aside and all the others serve as controls (%s)."), features, method[["controls"]]),
            tooltip = sprintf(paste("%s with negative controls chosen from the data: a first test of the groups is",
                "run and the 25%% most significant %s are left out. %s"), method[["controls"]], features, in_model),
            replicate_structure = "group", control_rule = "first_pass_drop_top_25", test_form = "factors_in_model",
            k_rules = c("num_sv_be", "sweep_bounded", "user"), k_rule = "num_sv_be"),
        legacy = list(id = "legacy", tier = "legacy", label = "Approach of earlier versions", description = earlier,
            tooltip = sprintf(paste("The approach of earlier versions: the least significant %s of a first test are",
                "the negative controls, and a search over the share of controls (a sweep) chooses how many and k."),
                features),
            replicate_structure = "group", control_rule = "sweep_least_significant", k_rules = character()),
        legacy_correction = list(id = "legacy_correction", tier = "legacy",
            label = "Approach of earlier versions, values corrected before the test", description = earlier,
            tooltip = sprintf(paste("The approach of earlier versions in its correction mode: the least significant",
                "%s of a first test are the negative controls, the quantified matrix is corrected and the model",
                "carries no factors."), features),
            replicate_structure = "group", control_rule = "sweep_least_significant",
            test_form = "corrected_then_tested", k_rules = character()),
        none = list(id = "none", tier = "none", label = "No removal of unwanted variation",
            description = "The differential test is run on the normalised values as they are.",
            tooltip = "No RUV: the differential test reads the normalised values.",
            test_form = "none", k_rules = character()))
}

#' One preset table: the catalogue's designs in the rows' order, each with its measured figure and note. A legacy
#' design's description ends with what its test reads in this table (its `test_form`).
#' @param features,method `omicGuiRuvDesignCatalogue()`'s arguments.
#' @param rows The omic's rows, named by design: `status`, the figure, `note`, `source`, and any field of the
#'   catalogue the omic words differently.
#' @param selected_with_technical_replicates The design selected when the design declares technical replicates;
#'   NULL where the omic has none.
#' @param k_rules_not_offered The k rules the run path does not take on this route, named, with the reason; NULL for
#'   none.
#' @return The table. The selected design is the recommended one (decision M2; KIT-129 fourth pass item 2): the
#'   section is shown only with its gate on, so today's tab is unchanged.
#' @noRd
omicGuiRuvTable <- function(features, method, rows, selected_with_technical_replicates = "technical_replicates",
    k_rules_not_offered = NULL) {
    catalogue <- omicGuiRuvDesignCatalogue(features, method)
    reads <- c(corrected_then_tested = "The values are corrected first, and the differential test is run on the corrected values.",
        factors_in_model = paste("The estimated factors are included as covariates in the differential test, which",
            "is run on the uncorrected values."))
    table <- list(selected = "groups_all_features", recommended = "groups_all_features",
        options = unname(lapply(names(rows), \(id) {
            option <- utils::modifyList(catalogue[[id]], rows[[id]])
            if (identical(option$tier, "legacy")) {
                option$description <- paste(option$description, reads[[option$test_form]],
                    "Its settings appear below when it is selected.")
            }
            option
        })))
    table$selected_with_technical_replicates <- selected_with_technical_replicates
    table$features <- features
    table$texts <- omicGuiRuvTexts(features)
    table$k_rules_not_offered <- k_rules_not_offered
    table
}

#' The RUV design section's texts beside the preset tables (KIT-129 fourth pass): the reasons an option or a k rule is
#' not available, the labels of the inputs the run reads, and the run's refusals of the panel's choices. One wording
#' table, under the user's review with the presets (docs/qa/workflow-execution/kit/kit-129/panel-ui/wording-review.md).
#' @param features The omic's features in words (`proteins`, `genes`, `lipids`, `metabolites`).
#' @return A named character vector.
#' @noRd
omicGuiRuvTexts <- function(features) {
    feature <- sub("s$", "", features)
    c(no_technical_replicates = "the experimental design declares no technical replicates.",
        technical_replicates_declared = paste("the experimental design declares technical replicates, so the approach",
            "within your technical replicates is used instead."),
        not_available = "This approach is not available:",
        k_rule_not_available = "This rule for the number of factors is not available here:",
        control_list = sprintf("Your control %s: one identifier per line, or separated by commas or spaces", features),
        seed = "Random seed for the estimated number:",
        k_missing = "Enter the number of factors (k) as a whole number of 1 or more.",
        seed_invalid = "Enter the random seed for the estimated number as a whole number.",
        list_missing = sprintf("Enter at least one control %s identifier.", feature),
        # KIT-129 fourth pass item 3: what the run did, under the options.
        result_heading = "What the last run did",
        result_design = "Approach used:",
        result_k = "Number of factors (k):",
        result_rule_num_sv_be = "estimated number",
        result_rule_sweep_bounded = "separation rule",
        result_rule_user = "your own number",
        result_k_estimate = "estimated number for comparison:",
        result_k_estimate_none = "the estimated number could not be computed:",
        result_sets_group = "Replicate sets: your groups.",
        result_sets_technical_replicates = "Replicate sets: your technical replicates (the measurements of each biological sample).",
        result_not_applied = "Unwanted variation was not removed:",
        result_seed = "Random seed for the estimated number:",
        result_assay = "Assay",
        # KIT-129 fifth pass item 3: a supplied list. %d: the identifiers listed, then those found in the data.
        result_list = "Your list: %d identifiers, of which %d were found in your data.",
        result_controls = "Controls used (measured in every sample):",
        separation_plot_heading = "Number of factors (k): the separation rule",
        # KIT-129 sixth pass (decisions M35, M36; drafts, the user sets the final wording): the check that an estimated
        # unwanted factor is associated with the groups, as a detection limit with the approach's measured consequences.
        # detectable: the largest association (0 to 1), the detection limit, then the approach's false discovery
        # proportion in simulation with no association, at 0.3 and at 0.6, and without correction at 0.6.
        association_detectable = paste("An estimated unwanted factor is associated with your groups (%s; values above %s",
            "are detectable for this approach). In simulation, this approach's false discovery proportion at a 5%% target",
            "was %s with no association, %s at 0.3 and %s at 0.6. Without correction it was %s at 0.6."),
        # not detectable: the largest association, then the detection limit.
        association_not_detectable = paste("No association between the estimated unwanted factors and your groups was",
            "detectable (largest %s; detection limit %s)."),
        # no limit: the largest association.
        association_no_limit = paste("Largest association of an estimated unwanted factor with your groups: %s. No detection",
            "limit has been established for this approach with this rule for the number of factors."),
        # M35: a typed k equal to the estimated number is checked against the estimated number's limit.
        association_typed_k_equals = paste("Your number of factors equals the estimated number, so the detection limit of the",
            "estimated number applies."),
        association_not_measured = "not measured",
        # DIA with limpa: the largest association.
        association_dia_limpa = paste("Largest association of an estimated unwanted factor with your groups: %s. In DIA with",
            "limpa no detection limit has been established: in simulation a candidate limit was exceeded in 86%% of data",
            "sets at an association of 0.6 and in only 37%% at 0.3, so only strong associations would be detectable."),
        association_legacy = paste("No detection limit exists for the approach of earlier versions: the association of its",
            "factors with your groups does not show this condition, so it is not checked."),
        # KIT-129 fifth pass item 5: the RUV QC tab's column heading of the check.
        association_column = "Unwanted factors and your groups")
}

#' The new design's fields in words, for the saved study parameters, the report's methods paragraph and the record's
#' workbook (KIT-129 sixth pass items 4 and 6; drafts, the user sets the final wording; listed in
#' docs/qa/workflow-execution/kit/kit-129/panel-ui/wording-review.md). No value holds "=", which the reports' parser of
#' `study_parameters.txt` reads as a separator.
#' @return A named character vector.
#' @noRd
omicGuiRuvRecordTexts <- function() {
    c(section = "Unwanted variation (RUV)",
        approach = "Approach",
        design_groups_all_features = "all features as controls, the experimental groups as replicates (recommended)",
        design_technical_replicates = "all features as controls, the declared technical replicates as replicates",
        design_supplied_controls = "a supplied list of control features, the experimental groups as replicates",
        design_first_pass_nulls = "controls chosen from the data: the features outside the 25% most significant of a first test",
        design_none = "no correction for unwanted variation",
        replicates = "Replicates",
        # %s: the design's column
        replicates_group = "the experimental groups (column %s)",
        replicates_technical_replicates = "the declared technical replicates (column %s)",
        k_rule = "Rule for the number of factors (k)",
        rule_num_sv_be = "estimated number (sva num.sv, method be)",
        rule_sweep_bounded = "separation rule (canonical correlations of the controls against all features)",
        rule_user = "chosen by the analyst",
        k_typed = "Number of factors typed",
        k = "Number of factors applied",
        k_estimate = "Estimated number of factors",
        # %s: the estimated number, its random seed and its permutations
        k_estimate_value = "%s (random seed %s, %s permutations)",
        controls = "Controls used (measured in every sample)",
        list = "Supplied control list",
        # %s: identifiers listed, found in the data, the list's sha256
        list_value = "%s identifiers, %s found in the data (sha256 %s)",
        first_pass = "First test for the controls",
        first_limma_f_test_groups = "limma F test of the groups",
        first_dpcde_f_test_groups = "limpa dpcDE F test of the groups",
        first_limma_smallest_pairwise_p = "limma, the smallest p-value over the pairwise comparisons of the groups",
        first_edger_qlf_f_test_groups = "edgeR quasi-likelihood F test of the groups",
        se = "Standard errors of the comparisons",
        se_exact_per_comparison = "exact for each comparison (the comparison refitted as a model coefficient, limma contrastAsCoef)",
        se_contrasts_fit = "limma contrasts.fit (exact when no factor is in the model)",
        se_edger_contrast = "edgeR, each comparison tested as a contrast of one fit",
        test = "Test",
        test_limma_treat_trend_robust = "limma treat with the eBayes trend and robust settings",
        test_limpa_dpcde = "limpa dpcDE with limma eBayes",
        test_limma_da_settings = "limma with the differential analysis tab's eBayes settings",
        test_edger_qlf = "edgeR quasi-likelihood F test",
        association = "Association of the factors with the groups",
        # %s: the largest association, then the outcome's words
        association_value = "largest %s; %s",
        outcome_detectable = "detectable (above the detection limit of %s for this approach)",
        outcome_not_detectable = "not detectable (detection limit %s)",
        outcome_no_limit = "no detection limit established for this approach with this rule for the number of factors",
        not_applied = "Not applied",
        # %s: the assay's name
        assay = "assay %s",
        methods = "Methods text",
        # %s: the method, the replicate sets, the controls, k, the k rule, the standard errors, the test
        methods_applied = paste("Unwanted variation was estimated with %s, using %s as replicate sets and %s as negative",
            "controls. %s factor(s) were used (%s) and included as covariates in the differential model fitted to the",
            "uncorrected values; the corrected matrix was used for figures and export only. Standard errors: %s. Test: %s."),
        # RUVg uses no replicate sets. %s: the method, the controls, k, the k rule, the standard errors, the test
        methods_applied_ruvg = paste("Unwanted variation was estimated with %s, using %s as negative controls. %s factor(s)",
            "were used (%s) and included as covariates in the differential model fitted to the uncorrected values; the",
            "corrected matrix was used for figures and export only. Standard errors: %s. Test: %s."),
        methods_none = "No correction for unwanted variation was applied.",
        # %s: the method, the reason
        methods_not_applied = "%s could not be applied (%s); no correction was made.",
        # %s: the largest association, the detection limit
        methods_association = paste("At least one estimated factor was associated with the experimental groups (largest",
            "association %s, above the detection limit of %s for this approach)."),
        method_ruviii = "RUV-III (Molania et al., 2019)",
        method_ruvs = "RUVs (Risso et al., 2014)",
        method_ruvg = "RUVg (Risso et al., 2014)",
        sets_group = "the experimental groups",
        sets_technical_replicates = "the declared technical replicates",
        # %s: the features in words (proteins, genes, lipids, metabolites)
        controls_all = "all %s",
        # %s: the number of listed controls found, the features
        controls_supplied = "a supplied list of %s %s",
        controls_first_pass = "the %s outside the 25%% most significant of a first test")
}

#' The scenarios of KIT-124's simulations in words, for `fdp_where`: the condition only (the column heading says
#' it is the worst scenario tested; a figure measured in one scenario carries `single_scenario`).
#' @return A named character vector by scenario.
#' @noRd
omicGuiRuvCells <- function() {
    c(c1 = "12 samples per group", c2 = "10% of features truly different", c3 = "no unwanted variation present",
        c5 = "unwanted variation tied to the groups (correlation 0.3)",
        c6 = "4 samples per group, each measured twice (technical duplicates)", c7 = "three groups")
}

#' The note of a supplied control list, with the false discovery proportion measured when one control in five
#' truly changes (KIT-124's arm 9).
#' @param features The omic's features in words.
#' @param contaminated That proportion.
#' @return One string.
#' @noRd
omicGuiRuvSuppliedNote <- function(features, contaminated) {
    sprintf(paste("Use when you have %s known not to differ between the groups. Only as good as the list: when one",
        "listed feature in five truly differed, false discoveries rose to %.1f%%."), features, 100 * contaminated)
}

#' The caution that leads the note of controls chosen from the data (ruling R3, the user, 2026-11-19), with the omic's
#' own false discovery proportion under association with the groups (`ruvGroupAssociationConsequences()`). Its last
#' sentence follows `ruvGroupAssociationThresholds()`: it points at the check's verdict only where the check has a
#' detection limit for this approach in this omic.
#' @param omic The check's omic key (`proteomics` for classic, `dia_limpa`, `transcriptomics`, `lipidomics`,
#'   `metabolomics`).
#' @param more How far above the target (`far more`, `more`).
#' @param measured The omic's measured figure and where it was measured.
#' @return One string.
#' @noRd
omicGuiRuvFirstPassCaution <- function(omic, more, measured) {
    thresholds <- ruvGroupAssociationThresholds()
    has_limit <- any(thresholds$omic == omic & thresholds$ruv_design == "first_pass_nulls")
    check <- if (has_limit) {
        "If the check below reports a detectable association, choose the recommended approach."
    } else {
        paste("The check below shows how strongly the estimated unwanted factors are associated with your groups, but",
            "gives no verdict for this approach. Unless you know the unwanted variation does not differ between your",
            "groups, choose the recommended approach.")
    }
    paste(sprintf(paste("Use this approach only when the unwanted variation is unlikely to differ between your groups.",
        "When it did in simulation, this approach reported %s false discoveries than the 5%% target (%s)."), more, measured),
        check)
}

#' Classic proteomics (DIA with IQ, LFQ, TMT): ticket section 1a, simulator S1.
#' @return The table.
#' @noRd
omicGuiRuvTableClassic <- function() {
    cell <- omicGuiRuvCells()
    omicGuiRuvTable("proteins", c(replicates = "RUV-III", controls = "RUV-III"), list(
        groups_all_features = list(status = "measured", fdp = 0.067, fdp_where = cell[["c1"]], power = 0.92,
            note = paste("The general-purpose choice. With 12 samples per group it reported slightly too many false",
                "discoveries: 11.7% of simulated data sets had more than 10% false discoveries (acceptance limit",
                "5.7%). 5.5% false discoveries in the reference scenario."),
            source = paste("KIT-129 1a; KIT-124 S1 arm 32 (scored/exact-gha; the shipped arm 34 identical, confirmation round",
                "Part 1): worst cell C1 0.067 (more than 10% in 0.117, limit 0.057), C0 0.055, power C0 0.92")),
        technical_replicates = list(status = "measured", fdp = 0.055, fdp_where = cell[["c6"]], power = 0.90,
            single_scenario = TRUE, needs = "technical_replicates",
            note = "Use when the experimental design declares technical replicates. Measured in one scenario only.",
            source = "KIT-129 1a; KIT-124 S1 arm 38 (scored/confirm-gha; identical to arm 33): C6 0.055, power C6 0.90"),
        supplied_controls = list(status = "measured", fdp = 0.054, fdp_where = cell[["c5"]], power = 0.91,
            note = omicGuiRuvSuppliedNote("proteins", 0.121),
            k_rule_figures = list(sweep_bounded = list(fdp = 0.060, fdp_where = cell[["c6"]], power = 0.91,
                source = "KIT-124 S1 arm 41 (scored/confirm-gha): worst cell C6 0.060, power C0 0.91")),
            source = "KIT-129 1a; KIT-124 S1 arm 35 (scored/confirm-gha): worst cell C5 0.054, power C0 0.91; arm 36: C1 0.121"),
        first_pass_nulls = list(status = "measured", fdp = 0.079, fdp_where = cell[["c5"]], power = 0.92,
            note = paste(omicGuiRuvFirstPassCaution("proteomics", "more", paste("7.9% at a correlation of 0.3, with the estimated",
                "number of factors; 28.0% at a correlation of 0.6, with the separation rule for the number of factors")),
                "Not when the unwanted variation is tied to the groups: it fails there with either rule",
                "for the number of factors (7.9% false discoveries, as shown). At most 5.3% in the other scenarios",
                "tested."),
            k_rule_figures = list(sweep_bounded = list(fdp = 0.079, fdp_where = cell[["c5"]], power = 0.91,
                source = "KIT-124 S1 arm 40 (scored/confirm-gha): worst cell C5 0.079, at most 0.056 elsewhere, power C0 0.91")),
            source = "KIT-129 1a; KIT-124 S1 arm 39 (scored/confirm-gha): worst cell C5 0.079, at most 0.053 elsewhere, power C0 0.92"),
        legacy = list(status = "measured", fdp = 0.385, fdp_where = cell[["c5"]], power = 0.92,
            test_form = "corrected_then_tested",
            note = paste("For reproducing analyses made with earlier versions. False discoveries ranged from 3% to",
                "39% across the scenarios tested, against a target of 5%; 9.8% in the reference scenario."),
            source = paste("KIT-129 1a; KIT-124 S1 arm 3 (scored/confirm-gha; C6 and C8 depend on the runner machine):",
                "worst cell C5 0.385, C0 0.098, power C0 0.92")),
        none = list(status = "measured", fdp = 0.052, fdp_where = cell[["c3"]], power = 0.49,
            note = paste("Only when little unwanted variation is expected. False discoveries stay near the 5%",
                "target, but about half the power is lost when unwanted variation is present."),
            source = "KIT-129 1a; KIT-124 S1 arm 37 (scored/confirm-gha; identical to arm 27): worst cell C3 0.052, power C0 0.49")))
}

#' DIA with limpa: ticket section 1b, simulator S2 (run on the `fastlimpadpc` engine).
#' @return The table.
#' @noRd
omicGuiRuvTableDiaLimpa <- function() {
    cell <- omicGuiRuvCells()
    omicGuiRuvTable("proteins", c(replicates = "RUV-III", controls = "RUV-III"), list(
        groups_all_features = list(status = "measured", fdp = 0.056, fdp_where = cell[["c6"]], power = 0.91,
            # The user's answer U2: on a design with technical replicates the groups are selectable in classic
            # proteomics only.
            needs = "no_technical_replicates",
            # Decision M42 (the user, 2026-11-21): the note states the result at a correlation of 0.6.
            note = paste("The general-purpose choice. In every scenario tested up to a correlation of 0.3 with the groups",
                "it did as well as an analysis that was given the true unwanted factor (5.3% false discoveries with",
                "unwanted variation tied to the groups, correlation 0.3). With unwanted variation tied to the groups at",
                "a correlation of 0.6 it had 6.6% false discoveries, against 4.7% for an analysis given the true",
                "unwanted factor, and 9.4% of simulated data sets had more than 10% false discoveries, against 0.6%.",
                "Measured with the fastlimpadpc engine; expected to hold for the limpa engine, not measured there."),
            source = paste("KIT-129 1b; KIT-124 S2 arm 32 (scored/exact-gha; the shipped arm 34 identical, confirmation round",
                "Part 1): worst cell C6 0.056, C5 0.053, power C0 0.91")),
        technical_replicates = list(status = "measured", fdp = 0.056, fdp_where = cell[["c6"]], power = 0.77,
            single_scenario = TRUE, needs = "technical_replicates",
            note = paste("Use when the experimental design declares technical replicates. Measured in one scenario",
                "only. In DIA with limpa, when technical replicates are declared, this and the recommended approach",
                "are the same calculation."),
            source = "KIT-129 1b; KIT-124 S2 arm 32 in C6 (scored/exact-gha; the shipped arm 34 identical): C6 0.056, power C6 0.77"),
        supplied_controls = list(status = "measured", fdp = 0.056, fdp_where = cell[["c6"]], power = 0.91,
            note = omicGuiRuvSuppliedNote("proteins", 0.090),
            k_rule_figures = list(sweep_bounded = list(fdp = 0.057, fdp_where = cell[["c5"]], power = 0.90,
                source = "KIT-124 S2 arm 41 (scored/confirm-gha): worst cell C5 0.057, power C0 0.90")),
            source = "KIT-129 1b; KIT-124 S2 arm 35 (scored/confirm-gha): worst cell C6 0.056, power C0 0.91; arm 36: C1 0.090"),
        first_pass_nulls = list(status = "measured", fdp = 0.147, fdp_where = cell[["c5"]], power = 0.92,
            note = paste(omicGuiRuvFirstPassCaution("dia_limpa", "far more", paste("14.7% at a correlation of 0.3, with the",
                "estimated number of factors; 30.2% at a correlation of 0.6, with the separation rule for the number of factors")),
                "Fails in DIA with limpa, with either rule for the number of factors: too many false",
                "discoveries in every scenario tested that had unwanted variation (7.2% to 14.7%; 9.2% in the",
                "reference scenario). Offered for completeness."),
            k_rule_figures = list(sweep_bounded = list(fdp = 0.151, fdp_where = cell[["c5"]], power = 0.92,
                source = "KIT-124 S2 arm 40 (scored/confirm-gha): worst cell C5 0.151, C0 0.095, power C0 0.92")),
            source = "KIT-129 1b; KIT-124 S2 arm 39 (scored/confirm-gha): worst cell C5 0.147, C1 0.072, C0 0.092, power C0 0.92"),
        legacy = list(status = "measured", fdp = 0.435, fdp_where = cell[["c5"]], power = 0.90,
            label = "Approach of earlier versions, factors in the test", test_form = "factors_in_model",
            note = paste("For reproducing analyses made with earlier versions. It chooses its controls as the limpa",
                "workflow script does. False discoveries ranged from 5% to 44% across the scenarios tested; 29.0%",
                "in the reference scenario."),
            source = "KIT-129 1b; KIT-124 S2 arm 4 (scored/confirm-gha): worst cell C5 0.435, C0 0.290, power C0 0.90"),
        legacy_correction = list(status = "measured", fdp = 0.466, fdp_where = cell[["c5"]], power = 0.91,
            note = paste("For reproducing analyses made with earlier versions. False discoveries ranged from 8% to",
                "47% across the scenarios tested."),
            source = "KIT-129 1b; KIT-124 S2 arm 3 (scored/confirm-gha): worst cell C5 0.466, power C0 0.91"),
        none = list(status = "measured", fdp = 0.191, fdp_where = cell[["c5"]], power = 0.56,
            note = paste("In DIA with limpa, false discoveries do not stay at the 5% target when unwanted variation",
                "is present: 15% of simulated data sets had more than 10% false discoveries; 6.5% false discoveries",
                "in the reference scenario."),
            source = "KIT-129 1b; KIT-124 S2 arm 1 (the shipped arm 37 identical, confirmation round Part 1): worst cell C5 0.191, C0 0.065, power C0 0.56")))
}

#' Transcriptomics: ticket section 1c, simulator S3. Technical replicates are summed into samples at the design
#' step, so that design is shown disabled with the reason (the user's answer U3).
#' @return The table.
#' @noRd
omicGuiRuvTableTrx <- function() {
    cell <- omicGuiRuvCells()
    omicGuiRuvTable("genes", c(replicates = "RUVSeq's RUVs", controls = "RUVSeq's RUVg"), list(
        groups_all_features = list(status = "measured", fdp = 0.055, fdp_where = cell[["c2"]], power = 0.77,
            note = paste("The general-purpose choice. In every scenario tested it did as well as an analysis that",
                "was given the true unwanted factor (5.2% false discoveries with unwanted variation tied to the",
                "groups at a correlation of 0.6)."),
            source = paste("KIT-129 1c; KIT-124 S3 arm 31 (scored/exact-gha; the shipped arm 34 identical, confirmation round",
                "Part 1): worst cell C2 0.055, power C0 0.77; arm 34 stress S3 0.052 (scored/confirm-gha)")),
        technical_replicates = list(status = "unavailable",
            unavailable = paste("the technical replicates of a biological sample are added together into one set",
                "of counts before the analysis (edgeR User's Guide, section 4.9.1), so none are left for RUV to",
                "compare."),
            note = "Not offered in transcriptomics.",
            source = "KIT-129 1c, answer U3; KIT-128 (the sum at the design step)"),
        supplied_controls = list(status = "measured", fdp = 0.056, fdp_where = cell[["c2"]], power = 0.76,
            note = paste(omicGuiRuvSuppliedNote("genes", 0.076),
                "The separation rule keeps false discoveries on target but costs power (61% against 76%)."),
            k_rule_figures = list(sweep_bounded = list(fdp = 0.053, fdp_where = cell[["c2"]], power = 0.61,
                source = "KIT-124 S3 arm 41 (scored/confirm-gha; identical to arm 22): worst cell C2 0.053, power C0 0.61")),
            source = "KIT-129 1c; KIT-124 S3 arm 35 (scored/confirm-gha): worst cell C2 0.056, power C0 0.76; arm 36: C1 0.076"),
        first_pass_nulls = list(status = "measured", fdp = 0.071, fdp_where = cell[["c5"]], power = 0.77,
            note = paste(omicGuiRuvFirstPassCaution("transcriptomics", "more", paste("6.4% at a correlation of 0.3 and 9.6% at a",
                "correlation of 0.6, with the separation rule for the number of factors")),
                "Not when the unwanted variation is tied to the groups: it fails there (7.1% false",
                "discoveries, as shown, and 7.1% of unchanged genes with p below 0.05, more than the acceptance",
                "criterion allows). At most 5.5% in the other scenarios tested, with as much power as an analysis",
                "that was given the true unwanted factor. With the separation rule it finds few of the real",
                "differences when the unwanted variation is tied to the groups (power 25%) or with three groups",
                "(5%)."),
            k_rule_figures = list(sweep_bounded = list(fdp = 0.064, fdp_where = cell[["c5"]], power = 0.60,
                source = paste("KIT-124 S3 arm 40 (scored/confirm-gha; identical to arm 24): worst cell C5 0.064, power C0 0.60,",
                    "C5 0.25, C7 0.049"))),
            source = paste("KIT-129 1c; KIT-124 S3 arm 39 (scored/confirm-gha): worst cell C5 0.071",
                "(type I 0.071, rule 1 fail), at most 0.055 elsewhere, power C0 0.77 (arm 2: 0.77)")),
        legacy = list(status = "measured", fdp = 0.299, fdp_where = cell[["c5"]], power = 0.73,
            test_form = "factors_in_model",
            note = paste("For reproducing analyses made with earlier versions. It fails when the unwanted variation",
                "is tied to the groups: 30% false discoveries at a correlation of 0.3, 72% at 0.6. At most 5.4% in",
                "the other scenarios tested."),
            source = paste("KIT-129 1c; KIT-124 S3 arm 3 (scored/confirm-gha): worst cell C5 0.299, stress S3 0.724, at most",
                "0.054 elsewhere, power C0 0.73")),
        none = list(status = "measured", fdp = 0.044, fdp_where = cell[["c3"]], power = 0.005,
            note = paste("Only when little unwanted variation is expected. False discoveries stay near the 5%",
                "target, but almost no real differences are found when unwanted variation is strong."),
            source = "KIT-129 1c; KIT-124 S3 arm 1 (the shipped arm 37 identical, confirmation round Part 1): worst cell C3 0.044, power C0 0.005")),
        selected_with_technical_replicates = NULL)
}

#' Lipidomics and metabolomics: ticket section 1d, simulator S4 (one table for both, scored over two assays).
#' Every figure is the shipped method's, from KIT-124's confirmation round (scored/confirm-gha; KIT-129 final figures,
#' 2026-11-20).
#' @param features `lipids` or `metabolites`.
#' @return The table.
#' @noRd
omicGuiRuvTableFeature <- function(features) {
    cell <- omicGuiRuvCells()
    omicGuiRuvTable(features, c(replicates = "RUV-III", controls = "RUV-III"), list(
        groups_all_features = list(status = "measured", fdp = 0.085, fdp_where = cell[["c5"]], power = 0.47,
            needs = "no_technical_replicates",
            note = sprintf(paste("The general-purpose choice, with limits: no approach met every acceptance",
                "criterion for %s. Too many false discoveries when the unwanted variation is strongly tied to the",
                "groups (15%% at a correlation of 0.6) and with 3 samples per group (11%%); 7.2%% in the reference",
                "scenario. If you know the batches, put batch in the statistical model: 18%% false discoveries when",
                "an unbalanced batch was left for RUV to find."), features),
            source = paste("KIT-129 1d; KIT-124 S4 arm 34 (scored/confirm-gha): worst cell C5 0.085, C0 0.072, stress S3 0.148,",
                "S1 0.111, S8 0.180, power C0 0.47")),
        technical_replicates = list(status = "measured", fdp = 0.063, fdp_where = cell[["c6"]], power = 0.34,
            single_scenario = TRUE, needs = "technical_replicates",
            note = "Use when the experimental design declares technical replicates. Measured in one scenario only.",
            source = paste("KIT-129 1d; KIT-124 S4 arm 38 (scored/confirm-gha; the cell is reported, not core): C6 0.063,",
                "power C6 0.34")),
        supplied_controls = list(status = "measured", fdp = 0.079, fdp_where = cell[["c5"]], power = 0.40,
            label = sprintf("Your own list of control %s, for example internal standards that were not used for normalisation",
                features),
            note = paste("Check that the standards do not differ between the groups: when 2 of 10 did, false",
                "discoveries rose to 40%. The separation rule gives more false discoveries and less power here (as",
                "shown). With 5 standards per assay, the separation rule's number of factors was not below the number",
                "of standards in about a quarter of simulated data sets (2,327 of 9,000), which were then analysed",
                "without RUV."),
            k_rule_figures = list(sweep_bounded = list(fdp = 0.095, fdp_where = cell[["c5"]], power = 0.26,
                source = paste("KIT-124 S4 arm 41 (scored/confirm-gha; 2,327 of 9,000 rows refused by M32 and scored",
                    "without RUV): worst cell C5 0.095, C0 0.081, power C0 0.26"))),
            source = "KIT-129 1d; KIT-124 S4 arm 35 (standards; scored/confirm-gha): worst cell C5 0.079, power C0 0.40; arm 36: C2 0.399"),
        first_pass_nulls = list(status = "measured", fdp = 0.203, fdp_where = cell[["c5"]], power = 0.48,
            note = paste(omicGuiRuvFirstPassCaution(c(lipids = "lipidomics", metabolites = "metabolomics")[[features]],
                "far more", paste("20.3% at a correlation of 0.3, with the estimated number of factors; 56.2% at a correlation",
                    "of 0.6, with the separation rule for the number of factors")),
                "Not when the unwanted variation is tied to the groups: it fails there (20.3% false",
                "discoveries, as shown). With three groups it missed two acceptance criteria (8.4% false",
                "discoveries). 5.2% to 7.3% in the other scenarios tested. With the separation rule it reports too",
                "many false discoveries in most scenarios tested (10.8% in the reference scenario)."),
            k_rule_figures = list(sweep_bounded = list(fdp = 0.261, fdp_where = cell[["c5"]], power = 0.53,
                source = "KIT-124 S4 arm 40 (scored/confirm-gha): worst cell C5 0.261, C0 0.108, power C0 0.53")),
            source = paste("KIT-129 1d; KIT-124 S4 arm 39 (scored/confirm-gha; C4 void): worst cell C5 0.203,",
                "C7 0.084 (rules 2, 3 fail), C0 0.063, C1 0.052, C2 0.073, C3 0.058, power C0 0.48")),
        legacy = list(status = "measured", fdp = 0.557, fdp_where = cell[["c2"]], power = 0.42,
            test_form = "corrected_then_tested",
            note = paste("Measured in earlier versions, when this design took its control features from the top of",
                "the table instead of choosing them by comparing the groups; in the simulations the rows were in random",
                "order, so the controls were a random set of features, some with real differences. More than 90% of",
                "those simulated data sets had more than 10% false discoveries (47.0% in the reference scenario). These",
                "figures do not describe the current selection, which has not yet been measured."),
            source = paste("KIT-129 1d; KIT-124 S4 arm 3 (scored/confirm-gha): worst cell C2 0.557, C0 0.470, power C0",
                "0.42; controls positional (random rows) before KIT-190")),
        none = list(status = "measured", fdp = 0.059, fdp_where = cell[["c7"]], power = 0.003,
            note = paste("Only when little unwanted variation is expected. False discoveries stay near the 5%",
                "target, but almost no real differences are found when the signal drifts during the run."),
            source = "KIT-129 1d; KIT-124 S4 arm 37 (scored/confirm-gha; identical to arm 1): worst cell C7 0.059, power C0 0.003")))
}

#' The DA tab of lipid and metab at the KIT-106 start commit (b75561de; the DA files unchanged since
#' 35f855bd for lipid and 0aee1526 for metab): one layout written twice (mod_lipid_da_ui.R:9-444,
#' mod_metab_da.R:38-523), differing in the test-id prefix, the nouns and the per-omic values of design
#' D7. Shared: the tab and status key (mod_lipid_da_analysis_helpers.R:272; mod_metab_da_observer_helpers.R:671),
#' the formula input and its help (mod_lipid_da_ui.R:24-30), the threshold and fold-change inputs
#' (:35-50), the tabs (:120, :177, :373), the static volcano's labels (mod_lipid_da_display_helpers.R:
#' 137-138; mod_metab_da_display_helpers.R:436-437), the heatmap's tree-cut defaults (mod_lipid_da_ui.R:
#' 310-331; mod_metab_da.R:333-366, KIT-101 F-03) and save note (mod_lipid_da.R:74-87;
#' mod_metab_da_registration_helpers.R:40), the table's defaults (mod_lipid_da_ui.R:399-420), the status and
#' contrasts texts (mod_lipid_da_display_helpers.R:29-40; mod_metab_da_display_helpers.R:278-329) and the
#' managed DA (utils_workflow_omics_actions.R:521, called at mod_lipid_da_analysis_helpers.R:198,
#' mod_metab_da_observer_helpers.R:722).
#' @param prefix The test-id prefix (`lipid`, `metab`).
#' @param values The per-omic values: `columns`, `summary`, `nouns`, `table_filter`, `static_volcano`
#'   (the plot), `heatmap` (its kernel, its label argument and save prefix), `run_in_process`, `source`,
#'   `session_load`, `preflight_first`, `download` (the file stem).
#' @noRd
omicGuiFeatureDa <- function(prefix, values) {
    list(module_id = "de", status_key = "differential_analysis",
        testids = c(load_session = paste0(prefix, "-da-load-session"), run = paste0(prefix, "-da-run"),
            download = paste0(prefix, "-da-download"), results_tabs = paste0(prefix, "-da-results-tabs")),
        model = list(editable = TRUE, label = "Model Formula:", default = "~ 0 + group",
            help = "Formula fetched from S4 object @args if available"),
        threshold = list(label = "Q-value Threshold:", value = 0.05, min = 0.001, max = 0.2, step = 0.005),
        lfc = list(label = "Log Fold-Change Cutoff:", value = 0, min = 0, max = 2, step = 0.1, tested = TRUE),
        columns = values$columns, summary = values$summary, nouns = values$nouns,
        tabs = c(volcano = "Volcano Plot", heatmap = "Heatmap", table = "DA Results Table"),
        assays = "multi", table_filter = values$table_filter,
        views = list(
            static_volcano = list(output = "volcano_static", plot = values$static_volcano, show_labels = TRUE, n_labels = 15),
            heatmap = list(kernel = values$heatmap$kernel, labels_argument = values$heatmap$labels_argument,
                tree_cut = list(default = "none", n_clusters = 4, cut_height = 10, min_cluster_size = 5),
                save_prefix = values$heatmap$save_prefix, save_note = TRUE),
            table = list(significance = "significant", max_rows = 1000)),
        run = list(managed = runWorkflowOmicsManagedDa, in_process = values$run_in_process),
        source = values$source, session_load = values$session_load,
        preflight_first = values$preflight_first,
        texts = c(status_ready = "Data loaded. Ready to run analysis.",
            status_waiting = "Waiting for data.\nClick 'Load Filtered Session' to begin.",
            contrasts_none = "No contrasts defined.\nLoad a filtered session or define contrasts in the design tab."),
        download = list(stem = values$download, extension = "csv"))
}

#' The lipid assay-name rule and messages (lipid `processing`:443-468; now
#' `validateOmicsImportAssayNames()`): a second name needs a second file and the reverse. Trx runs
#' the same rule (`processing`:652).
#' @noRd
omicGuiLipidAssayNames <- function() {
    list(rule = "require_pair", path_pattern = "[/\\\\]",
        messages = c(required = "Primary lipidomics assay name is required",
            path = "Assay names must not contain path separators",
            pair = "Second lipidomics assay requires both a file and an assay name",
            duplicate = "Duplicate lipidomics assay names are not allowed"))
}

#' The lipid panel's column ids (lipid `ui`:167, :188, :176, :196); trx renders them too.
#' @noRd
omicGuiLipidColumnIds <- function() {
    list(inputs = c(id_col = "lipid_id_col", id_col_custom = "lipid_id_col_custom"),
        outputs = c(id_status = "lipid_id_status", id_status_custom = "lipid_id_status_custom"))
}

# ------------------------------------------------------------------------------------------------
# Proteomics.
# ------------------------------------------------------------------------------------------------

#' The proteomics configuration: the DIA binding adaptor (utils_workflow_prot_declarations.R), the
#' owner `createDIASessionOwner` (mod_proteomics.R:266) and the DIA reader `readDIAStateMetadata`
#' (utils_workflow_dia_state.R:237).
#' @noRd
omicGuiConfigProteomics <- function() {
    list(shell = omicGuiProteomicsShell(),
        binding = list(binding = proteomicsDiaBinding, runner = "dia_managed_step",
            owner = createDIASessionOwner, read_state = readDIAStateMetadata,
            read_metadata = readDIAStateMetadata),
        import = omicGuiProteomicsImport(),
        design = omicGuiProteomicsDesign(),
        qc = omicGuiProteomicsQc(),
        # KIT-120: proteomics' tab on the shared frame (P-32 closed).
        norm = omicGuiProteomicsNorm(),
        da = omicGuiProteomicsDa(),
        enrichment = omicGuiProteomicsEnrichment(),
        summary = omicGuiProteomicsSummary())
}

#' Proteomics' normalisation tab on the shared frame (KIT-120; mod_prot_norm_ui.R and
#' mod_prot_norm_observer_helpers.R at 8ca68e93, mapped in kit-120/sites.md). The sections are the library's in its
#' order, with proteomics' own inputs, defaults, help texts, labels and test ids; the manual k starts unset (it is
#' filled from the canonical correlation). The tabs are the frame's without ITSD and the log, with the
#' Post-Normalisation QC tab. The stage grid shows the protein table's twelve figures under their existing file names
#' (`pre_norm_pca.png` ...) in `protein_qc_dir`, without an assay label. The actions and the other outputs are
#' proteomics' own (`protNormFrameServer`), for all five routes; the DIA_LIMPA route's covariates panel is rendered by
#' that server under the correlation summary.
#' @noRd
omicGuiProteomicsNorm <- function() {
    list(module_id = "normalization",
        panels = list(aesthetics = list(),
            between_sample = list(methods = c("Cyclic Loess" = "cyclicloess", "Quantile" = "quantile",
                "Scale (Median Absolute Values)" = "scale", "None" = "none"), default = "cyclicloess",
                help = "Method for between-sample normalization (default: Cyclic Loess)"),
            ruv_iii = list(modes = c("Automatic (recommended)" = "automatic", "Manual (advanced users)" = "manual",
                "Skip RUV (troublesome datasets)" = "skip"), default = "automatic", percentage = 5, k = NULL, noun = "Proteins",
                help = c(
                    grouping = paste("Experimental factor to preserve during batch correction. RUV-III will remove unwanted",
                        "variation while preserving differences between groups defined by this variable. For TMT workflows",
                        "with batches, select 'batch' to remove batch effects."),
                    automatic = paste("These parameters control the automatic optimization of RUV-III negative control",
                        "selection. The function tests different percentages and finds the optimal balance between",
                        "separation quality and k value."),
                    range = paste("Range of percentages to test (1-50%). Start with 1-20% for most datasets. Larger datasets",
                        "may benefit from higher ranges."),
                    metric = paste("Method for evaluating separation between All and Control groups in canonical correlation",
                        "plots. Max difference is most robust for typical datasets."),
                    penalty = paste("Controls how strongly high k values are penalized (0.1 = lenient, 0.9 = strict). Higher",
                        "values favor lower k to preserve biological signal."),
                    adaptive = "Automatically adjusts penalty thresholds based on dataset size. Recommended for most analyses.",
                    manual_k = "Will be auto-populated from canonical correlation analysis",
                    skip = paste("Batch correction will not be applied. Use this option for datasets where RUV-III is",
                        "problematic or unnecessary."),
                    skip_alternative = paste("Alternative approach: Consider including batch as a covariate in your",
                        "differential expression model if batch effects are present."))),
            actions = list(run_label = "Run Normalization & RUV"),
            export = list(button = "export_filtered_session", label = "Export Filtered Data Session",
                help = "Save filtered data and contrasts for later DE analysis"),
            s4_download = list(button = "export_s4_compatibility", label = "Export S4 compatibility RDS",
                help = "Optional export for older tools; requires a guarded session.")),
        plot_tabs = c("ruv", "pca", "density", "rle", "correlation", "post_norm_qc", "correlation_filter"),
        stage_labels = c("Post-Filtering", "Post-Normalization", "RUV-Corrected"), assay_slots = 1L,
        ruv_heading = "RUV Parameter Optimization Results",
        ruv_help = paste("Canonical correlation plot and optimization summary for RUV-III parameter selection. This plot",
            "shows the separation between 'All' and 'Control' protein groups across different k values."),
        ruv_single_result = TRUE,
        post_norm_qc = list(heading = "Filtering Progression Summary",
            help = paste("Summary of protein filtering steps through normalization and RUV correction. This step removes",
                "proteins with excessive missing values after RUV processing.")),
        correlation_slot = "ruv_covariates_panel",
        server = "protNormFrameServer",
        assay_names = "proteins", image_stemless = TRUE, qc_dir = "protein_qc_dir", assay_labels = FALSE,
        testids = c(method = "prot-norm-method", ruv_mode = "prot-norm-ruv-mode", run = "prot-norm-run",
            export = "prot-norm-export-session"),
        # KIT-129 W2: the classic routes' table, and DIA with limpa's once the workflow type is known.
        ruv_design = omicGuiRuvDesign("prot", list(classic = omicGuiRuvTableClassic(), dia_limpa = omicGuiRuvTableDiaLimpa()),
            route_tables = c(DIA_LIMPA = "dia_limpa"),
            # KIT-129 1 ("Conditions of use still open"): until the stress cells of these two simulators are scored.
            conditions = "The simulations covered up to 25% missing values and 6 or more samples per group."))
}

#' The proteomics Summary & Export tab at 5be34ee5 (mod_prot_summary_*.R; KIT-108 D1): the module id and test ids (the
#' shell's tab, `omicGuiProteomicsShell`; mod_prot_summary_ui.R:33, :46, :64, :72, :106), the placeholder and
#' description (:25, :27), the record and its bundle only on a `DIAStateAdapter` owner
#' (mod_prot_summary_session_helpers.R:59; mod_prot_summary_report_helpers.R:545), the copy's failures unclassified and
#' its `project_dirs` global looked up with `inherits = FALSE` (mod_prot_summary_publication_helpers.R:89-153, :105), no
#' managed integration operation and no gate (mod_prot_summary_session_helpers.R:82-98), the package's template directory
#' (mod_prot_summary_report_helpers.R:254), no one-template label (the route's status lists four templates, :29-61), no
#' export step, the native integration's notification (mod_prot_summary_session_helpers.R:97-98), and the report flow's
#' order (W1's captures): the progress "Detecting workflow type..." and the template check's progress always
#' (mod_prot_summary_report_helpers.R:512, :520), a failed retrieval stopping the report (:528-530), the record bundle
#' after the success notice and outputs (:534-550; `activateProtSummaryRenderedReport`, :389-408). Wired by KIT-108 W4
#' (utils_workflow_dia_summary.R; the integration writer returns the save's notice, :97-98, and the session-state
#' payload is `protSummarySessionState`, today's export in the core's shape): `final_state$context`
#' (`resolveProtSummaryFinalS4State`, mod_prot_summary_report_helpers.R:1-27), `copy$gate` (today's
#' `assertDIAExportReady(getDIAActiveManager())`, mod_prot_summary_publication_helpers.R:104), `integration$filename`
#' and `$write` (native `%s_%s_final.msdia.tar` with the archive's refresh and copy, else `%s_%s_final_s4.RDS` and
#' `saveRDS`; mod_prot_summary_session_helpers.R:82-93), `report$template` (`resolveProtSummaryReportTemplate`,
#' mod_prot_summary_report_helpers.R:131-248), `report$status` (`buildProtSummaryTemplateStatus`, :29-61) and
#' `session_state` (`prepareProtSummarySessionStateExport`, mod_prot_summary_session_helpers.R:181-215).
#' @noRd
omicGuiProteomicsSummary <- function() {
    list(module_id = "session_summary", testid = "prot", placeholder = "e.g., my_proteomics_analysis",
        description = "Full protein analysis workflow with config parameters",
        # KIT-108 W4: the route's plugins (utils_workflow_dia_summary.R); the integration writer returns its notice
        # (KIT-061: "Saved Integration input", the contract's folder `integration/proteomics/`).
        final_state = list(context = resolveProtSummaryFinalS4State),
        record_class = "DIAStateAdapter",
        copy = list(gate = protSummaryCopyGate, required = NULL, inherits = FALSE),
        integration = list(filename = protSummaryIntegrationFilename, operation = NULL, write = protSummaryIntegrationWrite,
            requires = NULL),
        report = list(dir = "proteomics", template = resolveProtSummaryReportTemplate, label = NULL,
            status = buildProtSummaryTemplateStatus, render_after_failed_template = FALSE, checking_progress_always = TRUE,
            record_after_success = TRUE),
        session_state = protSummarySessionState,
        export_steps = list(),
        # KIT-061: the integration input is not written while the owner holds no final matrix.
        texts = c(progress_detect = "Detecting workflow type...",
            integration_not_written = "not written (no final proteomics matrix yet)"))
}

#' The proteomics enrichment tab at 8e7294b6 (mod_prot_enrich*.R; KIT-107 D1): the tab, its status key and value
#' (the shell's tab, utils_omic_gui_config.R:1080-1081 at the start commit, module id = status key; the run's
#' status write, mod_prot_enrich_state_helpers.R:264-276; the tab observer's test,
#' mod_prot_enrich_registration_helpers.R:336), the
#' declaration and its run (`diaEnrichmentOperation`, `runDIAReferenceEnrichment`; utils_workflow_dia_da.R:539-559,
#' :572-601), the DA source (the manifest test, utils_workflow_dia_da.R:115; the source state's components the run
#' hydrates, mod_prot_enrich_analysis_helpers.R:533-534; the DA tables' reader, func_general_s4_objects.R:164-196), the
#' identifier column (the object's slot, its default and g:Profiler's gene-name override,
#' mod_prot_enrich_analysis_helpers.R:1-21), the cutoffs (mod_prot_enrich.R:53-60, :63-70, :73-80), the correction
#' (:144-152), the excluded IEA flag (not an input; mod_prot_enrich_analysis_helpers.R:41), the taxon and its prefill
#' from the bus (mod_prot_enrich.R:83-88; mod_prot_enrich_registration_helpers.R:81-108), the result tabs
#' (mod_prot_enrich.R:161-376; STRING-DB the placeholder), the organism filter's toggle and the bus key that switches
#' it on (mod_prot_enrich.R:93-99; mod_prot_enrich_registration_helpers.R:110-153) and the completed state
#' (mod_prot_enrich_state_helpers.R:227). The W2 merge wires the function fields whose functions are in
#' R/mod_prot_enrich_*.R at the start commit: `organism` (`buildProtEnrichSupportedOrganisms`,
#' mod_prot_enrich_observer_helpers.R:164-181; `resolveProtEnrichAnalysisMethod`, mod_prot_enrich_execution_helpers.R:
#' 206-236), `annotation` (`resolveProtEnrichUniprotAnnotations`, `resolveProtEnrichAnnotationMatching`,
#' mod_prot_enrich_input_helpers.R:599-717, :719-774), `organism_filter$mapping` and `$apply`
#' (`resolveProtEnrichOrganismMapping`, `applyProtEnrichOrganismFilter`, :1-180, :182-274) and `kernel`
#' (`resolveProtEnrichProcessEnrichmentsFn`, `runProtEnrichDeterministicProcessEnrichments`,
#' mod_prot_enrich_execution_helpers.R:196-204, :126-194).
#' @noRd
omicGuiProteomicsEnrichment <- function() {
    cutoff <- function(value, min, max, step) list(value = value, min = min, max = max, step = step)
    list(module_id = "enrichment_analysis", status_key = "enrichment_analysis", tab_value = "enrichment",
        operation = diaEnrichmentOperation, run = runDIAReferenceEnrichment,
        da_source = list(manifest = isDIAResultManifest, components = c("design_matrix", "protein_id_table"),
            read = createDAResultsForEnrichment),
        id_column = list(slot = "protein_id_column", default = "uniprot_acc", gprofiler_override = "gene_name"),
        cutoffs = list(up = cutoff(0, 0, 5, 0.1), down = cutoff(0, 0, 5, 0.1), q = cutoff(0.05, 0.001, 0.2, 0.005)),
        correction = list(choices = c(`g:SCS (Conservative)` = "gSCS", `FDR (Less Conservative)` = "fdr"),
            selected = "gSCS"),
        exclude_iea = FALSE,
        taxon = list(default = "9606", prefill_from = "taxon_id"),
        services = c("gprofiler2", "clusterprofiler", "stringdb"),
        organism_filter = list(mixed_species_from = "mixed_species_analysis", mapping = resolveProtEnrichOrganismMapping,
            apply = applyProtEnrichOrganismFilter),
        # KIT-107 W2 merge: the core's functions (utils_workflow_enrichment.R, utils_workflow_dia_enrichment.R).
        organism = list(supported = buildProtEnrichSupportedOrganisms, method = resolveProtEnrichAnalysisMethod),
        annotation = list(resolve = resolveProtEnrichUniprotAnnotations, match = resolveProtEnrichAnnotationMatching),
        kernel = list(select = resolveProtEnrichProcessEnrichmentsFn, test_double = runProtEnrichDeterministicProcessEnrichments),
        state_name = "enrichment_completed")
}

#' The proteomics DA tab at b75561de (mod_prot_da_*.R; KIT-106 D5): the tab and status key (mod_proteomics.R, the
#' shell's tab), the keys the run plugin also writes (utils_workflow_dia_da.R:451-452), the test ids
#' (mod_prot_da_ui.R:33, :98, :112, :125), the formula input and its help (:45-52), the threshold and fold-change
#' inputs (:57-71), the table's columns (the object's `protein_id_column`; mod_prot_da_render_table.R:28, :61-71) and
#' its summary (:103-149), the nouns (mod_prot_da_ui.R:183, :227, :254, :377; mod_prot_da_render_heatmap.R:117), the
#' tabs (mod_prot_da_ui.R:138, :208, :436), no assay selectors, the table paged from the store
#' (mod_prot_da_render_table.R:12-28), the static volcano with its label inputs' defaults (mod_prot_da_render_volcano.R:
#' 32-44), the heatmap (generateProtDAHeatmap, which routes a DA manifest to generateDIAReferenceHeatmap,
#' func_prot_da_heatmap.R:56-57; mod_prot_da_render_heatmap.R:31-99), its tree-cut defaults (mod_prot_da_ui.R:306-362)
#' and save prefix (mod_prot_da_render_heatmap.R:152), no save note (P-25), the status and contrasts texts (the
#' completed status, `diaDaStatusText`, is the plugins' own: KIT-109 removed the unread field), the
#' download (mod_prot_da_render_table.R:151-160). The init, load, run and archive steps are plugins (D5): each is the
#' handler today's server calls (mod_prot_da_server.R:36-45) or the archive chooser (mod_prot_da_ui.R:80-88).
#' @noRd
omicGuiProteomicsDa <- function() {
    # KIT-109 (A106-05): `status_keys` and `status_text` removed; the plugin server read neither.
    list(module_id = "differential_expression", status_key = "differential_expression",
        testids = c(load_session = "prot-da-load-session", run = "prot-da-run", download = "prot-da-download",
            results_tabs = "prot-da-results-tabs"),
        model = list(editable = TRUE, label = "Model Formula:", default = "~ 0 + group",
            help = "Formula fetched from S4 object @args"),
        threshold = list(label = "Q-value Threshold:", value = 0.05, min = 0.001, max = 0.2, step = 0.005),
        lfc = list(label = "Log Fold-Change Cutoff:", value = 0, min = 0, max = 2, step = 0.1, tested = TRUE),
        columns = list(id_slot = "protein_id_column", fold_change = "log2FC", pvalue = "raw_pvalue",
            significance = "fdr_qvalue"),
        summary = list(noun = "genes", significance_label = "q"),
        nouns = c(top_n = "Top N Genes:", scaling = "Row (gene) scaling", labels = "Show Gene Labels",
            label_top = "Label Top Proteins", cluster = "proteins"),
        tabs = c(volcano = "Volcano Plot", heatmap = "Heatmap", table = "DA Results Table"),
        assays = "none", table_filter = "store_page",
        views = list(
            static_volcano = list(output = "volcano_plot_static", plot = generateProtDAVolcanoStatic, show_labels = TRUE,
                n_labels = 10),
            heatmap = list(kernel = generateProtDAHeatmap, labels_argument = "show_gene_names",
                tree_cut = list(default = "k_clusters", n_clusters = 4, cut_height = 0.5, min_cluster_size = 3),
                save_prefix = "prot_", save_note = FALSE),
            table = list(significance = "significant", max_rows = 1000)),
        texts = c(status_waiting = "[WAITING] Waiting for analysis...\nClick 'Run DA Analysis' to start",
            contrasts_none = "No contrasts available.\nComplete normalization and\ncorrelation filtering first."),
        download = list(stem = "proteomics_da_results_", extension = "csv"),
        plugins = c(init = "da_server_init_handlers", load = "da_server_load_session_handler",
            run = "da_server_run_analysis_handler", archive = "prot_da_session_archive"))
}

#' The proteomics QC tab (mod_prot_qc.R, mod_prot_qc_peptide.R, mod_prot_qc_protein.R at b423b764):
#' the peptide and protein groups, rendered by their group plugins `mod_prot_qc_{peptide,protein}`,
#' each tab one of the 17 sub-tab plugins with the routes it shows on. The tabs per route are
#' today's spec tables (`getProtQcPeptideModuleSpecs`, :28-84; `getProtQcProteinModuleSpecs`,
#' :27-87): a tab's component is its `uiFn`/`serverFn` stem; `needs_paths` is the protein specs'
#' `needsExperimentPaths` (the peptide runner passes the paths by name when it has them, :96-103).
#' TMT and LFQ show the protein group only (mod_prot_qc.R:98-133). The sub-tabs start on the
#' design's trigger only (:150-201); QC completes on `protein_replicate_filtered` (:204-216), which
#' the DIA_LIMPA route never reaches (P-16, fixed by KIT-104 W7).
#' @noRd
omicGuiProteomicsQc <- function() {
    peptide <- function(id, title, stem, routes) {
        omicGuiQcTab(id, title, paste0("prot_qc_peptide_", stem), routes)
    }
    protein <- function(id, title, stem, routes) {
        omicGuiQcTab(id, title, paste0("prot_qc_protein_", stem), routes, needs_paths = TRUE)
    }
    dia <- c("DIA", "DIA_LIMPA")
    common <- c("DIA", "TMT", "LFQ")
    list(title = "Quality Control & Filtering", # mod_prot_qc.R:43
        # mod_prot_qc.R:121, :130 (TMT, LFQ), :138 (DIA, DIA_LIMPA)
        tabset_id = c(DIA = "qc_tabs_lfq", DIA_LIMPA = "qc_tabs_lfq", TMT = "qc_tabs_tmt",
            LFQ = "qc_tabs_tmt"),
        groups = list(
            # mod_prot_qc.R:140; mod_prot_qc_peptide.R:32-38 (DIA_LIMPA), :40-83 (DIA)
            list(id = "peptide_qc", title = "Peptide QC", component = "prot_qc_peptide", tabs = list(
                peptide("qvalue_filter", "Q-Value Filter", "qvalue", dia),
                peptide("protein_peptide_filter", "Protein Evidence", "protein", "DIA"),
                peptide("rollup", "Precursor Rollup", "rollup", dia),
                peptide("intensity_filter", "Intensity Filter", "intensity", "DIA"),
                peptide("sample_filter", "Sample Quality", "sample", "DIA"),
                peptide("replicate_filter", "Replicate Filter", "replicate", "DIA"),
                peptide("imputation", "Imputation", "impute", "DIA"),
                peptide("log2_matrix", "Log2 Matrix", "log2", "DIA_LIMPA"),
                peptide("accession_cleanup", "Peptide Accession Cleanup", "accession", "DIA_LIMPA"))),
            # mod_prot_qc.R:123, :131, :141; mod_prot_qc_protein.R:27-35 (DIA), :37-68, :73-79
            # (DIA_LIMPA)
            list(id = "protein_qc", title = "Protein QC", component = "prot_qc_protein", tabs = list(
                protein("rollup", "IQ Protein Rollup", "rollup", "DIA"),
                protein("cleanup", "Accession Cleanup", "cleanup", common),
                protein("intensity_filter", "Protein Intensity Filter", "intensity", common),
                protein("duplicate_removal", "Duplicate Removal", "dedup", common),
                protein("replicate_filter", "Protein Replicate Filter", "replicate", common),
                protein("dpc_quant", "DPC Quantification", "dpc_quant", "DIA_LIMPA"),
                protein("normalisation", "Between-Sample Normalisation", "normalisation", "DIA_LIMPA"),
                protein("peptide_count", "Peptide Count Filter", "peptide_count", "DIA_LIMPA")))),
        auto_init = FALSE,
        # P-16 (KIT-104 W7, the lead's decision 2026-11-08): DIA_LIMPA completes at its last QC
        # state (mod_prot_qc_protein_peptide_count.R), so normalisation opens as on the other routes.
        complete_state = c(DIA = "protein_replicate_filtered", DIA_LIMPA = "peptide_count_filtered",
            TMT = "protein_replicate_filtered", LFQ = "protein_replicate_filtered"),
        itsd = NULL, texts = stats::setNames(character(), character()),
        notification_ids = stats::setNames(character(), character()))
}

#' The proteomics shell (mod_proteomics.R): no header; the icon as `tabPanel`'s argument
#' (:117-184, stepper :380-386); the bus (:205-245) on `session$userData` (:248); the
#' peptide-checkpoint panel (:101-113, :275-304). TMT and LFQ open the QC tab after the design, as
#' DIA does (VAL-26, KIT-180).
#' @noRd
omicGuiProteomicsShell <- function() {
    tab <- function(status_key, value, component, title, icon, testid, step_name) {
        omicGuiTab(status_key, value, status_key, component, title, icon, testid, step_name)
    }
    bus <- list(data_tbl = NULL, fasta_file_path = NULL, config_list = NULL, taxon_id = NULL,
        organism_name = NULL, design_matrix = NULL, data_cln = NULL, contrasts_tbl = NULL,
        state_manager = NULL, peptide_data = NULL, protein_log2_quant = NULL, protein_data = NULL,
        ruv_normalised_for_da_analysis_obj = NULL, da_analysis_results_list = NULL,
        uniprot_dat_cln = NULL, enrichment_results = NULL, tab_status = NULL,
        state_update_trigger = NULL, processing_log = list())
    list(omic_type = "proteomics", prefix = "prot", testid = "prot", label = NULL,
        label_lower = "proteomics", header_icon = NULL, tabset_id = "workflow_tabs",
        tab_style = "icon_argument", routes = c("DIA", "DIA_LIMPA", "TMT", "LFQ"),
        tabs = list(
            tab("setup_import", "setup", "prot_import", "Setup & Import", "upload",
                "prot-tab-setup", "Import"),
            tab("design_matrix", "design", "prot_design", "Design Matrix", "table",
                "prot-tab-design", "Design"),
            tab("quality_control", "qc", "prot_qc", "Quality Control", "chart-line",
                "prot-tab-qc", "QC"),
            tab("normalization", "normalization", "prot_norm", "Normalization", "balance-scale",
                "prot-tab-normalization", "Normalize"),
            tab("differential_expression", "da", "prot_da", "Differential Expression", "chart-bar",
                "prot-tab-da", "DE"),
            tab("enrichment_analysis", "enrichment", "prot_enrich", "Enrichment Analysis",
                "network-wired", "prot-tab-enrichment", "Enrich"),
            tab("session_summary", "session_summary", "prot_summary", "Session Summary & Report",
                "file-export", "prot-tab-session-summary", "Summary")),
        bus = bus, bus_registration = "userData",
        extras = "peptide_checkpoint")
}

#' The proteomics source panel's data: the title (the frame, mod_omic_import.R), and the format
#' override and DIA-NN route choices, the step headings, the auxiliary inputs and the test ids that
#' `protImportPanelUi` (mod_prot_import_ui.R) renders; the success text is the one
#' `completeProtImportSuccessState` shows (mod_prot_import_processing_helpers.R).
#' @noRd
omicGuiProteomicsImport <- function() {
    list(source_panel = "proteomics", title = "Setup & Import Data",
        formats = c("Auto-detect" = "auto", "DIA-NN" = "diann", "Spectronaut DIA" = "spectronaut",
            "FragPipe LFQ" = "fragpipe", "MaxQuant LFQ" = "maxquant",
            "Proteome Discoverer TMT" = "pd_tmt"),
        default_format = "auto",
        route_choice = c("DIA, iq MaxLFQ (optional replicate-mean imputation)" = "DIA",
            "DIA with imputation leveraging limpa (Li et al. 2025)" = "DIA_LIMPA"),
        default_route = "DIA",
        # func_prot_import_detection.R:14
        detect = detectProteomicsFormat,
        # The header read and the reading by format (func_prot_import_detection.R,
        # func_prot_import_readers.R); the source panel calls them in its own order.
        preview = readProtImportHeaders, process = importProteomicsData,
        help = c(search = "Step 1: Import Searched Data File", fasta = "Step 2: Import FASTA File",
            organism = "Step 3: Organism Information"),
        aux_inputs = c("fasta", "organism", "mapping_files"),
        testids = c(search_results_file = "prot-import-search-results-file",
            fasta_file = "prot-import-fasta-file", process = "prot-import-process"),
        success_text = "Data import successful!")
}

#' The proteomics design tab (mod_prot_design.R) over the `dia_long` adapter (the DIA adaptor's
#' `design` entries), with its own builder UI (`mod_prot_design_builder_ui`, back-reference #4),
#' its own import panel and its own commit (`runProtDesignBuilderSaveFlow`). `build_object`
#' (`buildProtDesignStateCheckpoint`) and `annotation_step` (the UniProt and FASTA steps) are the
#' core's (KIT-103 W3, utils_workflow_dia_state.R and the proteomics annotation files). File:line
#' notes are at 3c657c14; mod_prot_design_builder_display_helpers.R was deleted by KIT-103 W6.
#' @noRd
omicGuiProteomicsDesign <- function() {
    list(module_id = "design_matrix",
        # mod_prot_design.R:82
        help = paste("Use the tools below to define your experimental groups and the contrasts for",
            "differential analysis. Alternatively, import an existing design from a previous",
            "analysis."),
        # No save-time preflight reads it (D7); the builder's contrast note names it
        # (mod_prot_design_builder_display_helpers.R:400-406 at 3c657c14).
        formula_default = "~ 0 + group",
        # mod_prot_design_builder_state_helpers.R:348
        run_order = gtools::mixedsort, run_order_rule = "mixedsort", group_levels = NULL,
        validate = NULL,
        # utils_workflow_dia_state.R:742, called in buildProtDesignStateCheckpoint
        # (mod_prot_design_state_helpers.R:15, :98)
        record_provenance = setDIAStepProvenance,
        # KIT-103 W3: the design state's checkpoint (utils_workflow_dia_state.R) and the annotation
        # step (func_prot_annotation_uniprot_full_helpers.R, func_prot_annotation_fasta_processing_helpers.R).
        build_object = buildProtDesignStateCheckpoint,
        annotation_step = list(uniprot = retrieveProtDesignUniprotAnnotations,
            fasta = processProtDesignImportedFasta),
        # mod_prot_design_state_helpers.R:19-21, :90-91, by route (`workflowType`)
        state_name = c(DIA = "raw_data_s4", DIA_LIMPA = "raw_data_s4", TMT = "protein_s4_initial",
            LFQ = "protein_s4_initial"),
        state_description = c(DIA = "Initial Peptide S4 object created after design matrix.",
            DIA_LIMPA = "Initial Peptide S4 object created after design matrix.",
            TMT = "Initial Protein S4 object created from TMT data after design matrix.",
            LFQ = "Initial Protein S4 object created from LFQ data after design matrix."),
        data_adapter = "dia_long",
        # mod_prot_design_import_helpers.R:564-573: in test mode `import_dir_path` renders "" until
        # a path is typed.
        design_import = list(plugin = "proteomics", import_dir_path = "empty"),
        commit = "prot_hook", builder_ui = "proteomics", record_excluded_runs = FALSE,
        label_lower = "proteomics", multi_assay_banner = FALSE,
        # mod_prot_design.R:117, :98, :101-107
        no_data_text = paste(" Please complete the 'Setup & Import' step first. The builder will",
            "appear here once data is available."),
        preview_text = paste("This section shows the design matrix and contrasts that have been",
            "saved to the workflow."),
        preview_assays = FALSE,
        # mod_prot_design.R:78; mod_prot_design_builder_action_helpers.R:7
        testids = c(import_existing = "prot-design-import-existing",
            replicate_start = "prot-design-replicate-start"),
        # mod_prot_design_builder_state_helpers.R:109, :122; mod_prot_design_builder_display_helpers.R
        # :403, :405; mod_prot_design_builder_action_helpers.R:4 (all at 3c657c14)
        builder_texts = c(save_success = "Design saved successfully. You can close this builder.",
            no_assignment = paste("No samples have been assigned to groups. Please assign metadata",
                "before saving."),
            contrast_group_prefix = paste0("Note: Contrasts will use 'group' prefix (e.g., ",
                "groupGA_Control-groupGA_Elevated)\nbased on current formula: ~ 0 + group"),
            contrast_as_is = "Note: Contrasts will use group names as-is (e.g., GA_Control-GA_Elevated)",
            replicate_suffix = "selected runs:"))
}

# ------------------------------------------------------------------------------------------------
# Lipidomics.
# ------------------------------------------------------------------------------------------------

#' The lipidomics configuration (mod_lipidomics.R; the readers utils_workflow_lipid_declarations.R
#' :120, :123).
#' @noRd
omicGuiConfigLipidomics <- function() {
    list(shell = omicGuiFeatureShell("lipidomics", "lipid", "Lipidomics", "flask", "LIPIDOMICS",
            " Differential Analysis", omicGuiFeatureBus("lipid", "lipid_assay_data")),
        binding = omicGuiFeatureBinding(lipidomicsDeclarationBinding, readLipidWorkflowState,
            readLipidWorkflowMetadata),
        import = omicGuiLipidImport(),
        design = omicGuiLipidDesign(),
        qc = omicGuiLipidQc(),
        norm = omicGuiLipidNorm(),
        da = omicGuiLipidDa(),
        # KIT-107 D9: no enrichment tab.
        enrichment = NULL,
        summary = omicGuiLipidSummary())
}

#' The publication copy of lipid, metab and trx at 5be34ee5 (KIT-108 D1): no gate; the three required artifacts whose
#' failure stops the copy (`classifyLipidSummaryCopyFailures`, mod_lipid_summary_state_helpers.R:211-227, called at
#' mod_lipid_summary_workflow_helpers.R:268-277; metab's mod_metab_summary_observer_helpers.R:147-163, :617-626); the
#' `project_dirs` global looked up with `exists`' default `inherits = TRUE` (mod_lipid_summary_workflow_helpers.R:255;
#' mod_metab_summary_observer_helpers.R:603).
#' @noRd
omicGuiFeatureSummaryCopy <- function() {
    list(gate = NULL, required = c("Contrasts Table", "Design Matrix", "Study Parameters"), inherits = TRUE)
}

#' The lipid Summary & Export tab at 5be34ee5 (mod_lipid_summary*.R; KIT-108 D1): the module id and test ids (the
#' shell's tab; mod_lipid_summary.R:41, :61-134), the placeholder and description (:53, :55), the final state: the first
#' of the six data states, head first, read as ids on a managed session (`readLipidWorkflowIds`,
#' utils_workflow_lipid_declarations.R:130), else `getState` (`collectLipidSummaryWorkflowArgsContext`,
#' mod_lipid_summary_state_helpers.R:312-377, states :322-329, :332, :336-337), the record whenever there is a bus
#' (mod_lipid_summary_workflow_helpers.R:67, :479), the copy (`omicGuiFeatureSummaryCopy`), the integration file and its
#' managed operation, else the final object's RDS in process (:89, :124-132), the template, its directory and its
#' label (:359, :387; mod_lipid_summary_state_helpers.R:7-24), the report flow's order (the template check's progress
#' only when the template is missing, :407-408; a failed retrieval stops the report, :444-446; the record bundle before
#' the success notice, :477-494), no export step and no text of its own. Wired at the W3 merge: `session_state`
#' (`prepareLipidSummarySessionState`, mod_lipid_summary_state_helpers.R:133-209).
#' @noRd
omicGuiLipidSummary <- function() {
    list(module_id = "summary", testid = "lipid", placeholder = "e.g., my_lipidomics_analysis",
        description = "Full lipidomics analysis workflow with normalization and DA",
        final_state = list(states = c("lipid_correlation_filtered", "lipid_norm_complete", "lipid_ruv_corrected",
            "lipid_normalized", "loaded_for_de", "lipid_qc_complete"), read_ids = readLipidWorkflowIds),
        record_class = NULL,
        copy = omicGuiFeatureSummaryCopy(),
        # KIT-135: the in-process writer writes the RDS as before, then the declared integration input beside it.
        integration = list(filename = "%s_%s_final_s4.RDS", operation = "lipid_integration_export", write = lipidSummaryIntegrationWrite,
            requires = NULL),
        report = list(dir = "lipidomics", template = "lipidomics_report.rmd", label = "Lipidomics", status = NULL,
            render_after_failed_template = FALSE, checking_progress_always = FALSE, record_after_success = FALSE),
        export_steps = list(),
        texts = character(),
        session_state = prepareLipidSummarySessionState)
}

#' The lipid DA tab at b75561de (the shared values are `omicGuiFeatureDa`'s). The columns
#' (mod_lipid_da_results_helpers.R:106), the summary's noun and label (:6-7), the nouns (mod_lipid_da_ui.R:195,
#' :222, :301; mod_lipid_da_display_helpers.R:239), the table filter where only "All" means every assay
#' (mod_lipid_da_results_helpers.R:20, :117), the static volcano and heatmap kernels (mod_lipid_da_display_helpers.R:
#' 121-140, :158, :173) and the heatmap's save prefix (mod_lipid_da_export_helpers.R:74), the in-process reference DA
#' (mod_lipid_da_analysis_helpers.R:218), the source and its no-data text (:34-36), the session file and assay check:
#' a set comparison, for a lipidomics session only (mod_lipid_da_session_helpers.R:56-64, :227), the preflight after
#' the running notice (mod_lipid_da_analysis_helpers.R:139-157), the download's stem (mod_lipid_da_results_helpers.R:
#' 173-193). The preflight (`validateLipidDaAnalysisPreflight`, mod_lipid_da_analysis_helpers.R:84) moved
#' to the core at KIT-106 W2 (func_omics_design_shared.R) and is wired.
#' @noRd
omicGuiLipidDa <- function() {
    da <- omicGuiFeatureDa("lipid", list(
        columns = list(id = "lipid_id", label = "lipid_name", fold_change = "logFC", pvalue = "raw_pvalue",
            significance = "fdr_qvalue"),
        # KIT-172: the count of rows not tested because a compared group has no observed value.
        summary = list(noun = "lipids", significance_label = "Q", not_tested = TRUE),
        nouns = c(top_n = "Top N Lipids:", scaling = "Row (lipid) scaling", labels = "Show Lipid Labels", cluster = "lipids"),
        table_filter = "all_only",
        static_volcano = generateLipidDAVolcanoStatic,
        heatmap = list(kernel = generateLipidDAHeatmap, labels_argument = "show_lipid_names", save_prefix = "lipid_"),
        run_in_process = runLipidReferenceDaAnalysis,
        source = list(read_metadata = readLipidWorkflowMetadata, class = "LipidomicsAssayData",
            no_data = "No lipidomics data loaded. Please load a filtered session first."),
        session_load = list(file = "lipid_filtered_session_data_latest.rds", assay_rule = "setequal",
            assay_omic_type = "lipidomics",
            assay_message = "Loaded lipid session assay_names do not match S4 lipid_data assays: session [%s], S4 [%s]"),
        preflight_first = FALSE,
        download = "lipidomics_da_results_"))
    # KIT-106 W2 moved the preflight to the core (func_omics_design_shared.R); wired at the W2 merge.
    da$preflight <- validateLipidDaAnalysisPreflight
    da
}

#' The lipid normalisation tab (mod_lipid_norm_*.R at 9a98e765; the shared values are
#' `omicGuiFeatureNorm`'s). The RUV noun (mod_lipid_norm_ui_helpers.R:187, :194) and the test ids
#' (:84, :103, :234, :259). The export: the `lipid_` files (mod_lipid_norm_workflow_helpers.R:412-416,
#' :576), the `lipid_counts` field (:382) and the per-assay sample count of design samples among the
#' columns (:336-346; P-18). The ITSD selection table: buildLipidItsdSelectionTable
#' (func_lipid_norm_support_helpers.R:19), as the tab calls it (mod_lipid_norm_observer_helpers.R:181,
#' mod_lipid_norm_workflow_helpers.R:212). D7: skip-correlation sets the flags before the store commit
#' (mod_lipid_norm_workflow_helpers.R:633-639); the managed final PCA's placeholders are sizes 5 and 4
#' (mod_lipid_norm_support_helpers.R:350-351). KIT-105 W5 (D7, at f8d4c5fc): the export's progress shows the
#' gathering and the completion steps only, and the summary carries the session's export time
#' (mod_lipid_norm_workflow_helpers.R:469, :568-577); a failed managed correlation is logged, not notified (:738-743).
#' @noRd
omicGuiLipidNorm <- function() {
    omicGuiFeatureNorm("lipid", "Lipids",
        export = list(prefix = "lipid_", counts_field = "lipid_counts", samples = "design",
            progress = list(gather = "Gathering data and publishing outputs...", complete = "Export complete", steps = FALSE),
            summary_time = "session"),
        itsd_table = function(assay_data, id_col, annotation_cols) {
            buildLipidItsdSelectionTable(assay_data = assay_data, lipid_id_col = id_col, annotation_cols = annotation_cols)
        },
        skip_correlation_flags = "before_commit", final_pca = list(preparing_size = 5, error_size = 4),
        managed_correlation_notify = FALSE)
}

#' The lipid QC tab (mod_lipid_qc*.R at b423b764). `itsd` holds the inline analysis's entries
#' (mod_lipid_qc_itsd.R): the name columns tried after the id and annotation columns (:257-263), the
#' metrics kernel with its id-column argument (:282-287), the quant getter (:305) and the refusal's
#' tip (:337).
#' @noRd
omicGuiLipidQc <- function() {
    list(title = "Lipid Quality Control & Filtering", # mod_lipid_qc.R:45
        tabset_id = "lipid_qc_tabs", # mod_lipid_qc.R:120
        groups = omicGuiFeatureQcGroups("LIPIDOMICS"),
        # Started on the trigger or when the state's class is the binding's (mod_lipid_qc.R:129-159).
        auto_init = TRUE, complete_state = NULL,
        # mod_lipid_qc_intensity_server_helpers.R:142; mod_lipid_qc_s4.R:450 (func_lipid_qc.R:134)
        filtering_figure = updateLipidFiltering,
        itsd = list(name_columns = c("Lipid name", "Name", "lipid"),
            metrics = function(assay_data, is_pattern, id_col, sample_id_col) {
                getLipidInternalStandardMetrics(assay_data = assay_data, is_pattern = is_pattern,
                    lipid_id_col = id_col, sample_id_col = sample_id_col)
            },
            quant = getLipidQuantData,
            tip = "Tip: Check if your data has a column with lipid names (e.g., 'Lipid name')"),
        # Any zero-row assay is refused (mod_lipid_qc_s4.R:365); normalisation opens (:420-436); a
        # per-assay table (:75-79).
        finalize_zero_rows = TRUE, finalize_opens_norm = TRUE, finalize_per_assay = TRUE,
        # D11 (KIT-104 W5, the shared finalise's first consumer): the shared alignment with lipid's
        # arguments (mod_lipid_qc_s4.R:374-394 at b423b764).
        finalize_alignment = function(shell) omicsQcFinalizeAlignment(shell, "lipidomics"),
        texts = omicGuiFeatureQcTexts("lipid", "LipidomicsAssayData"),
        # mod_lipid_qc_intensity_server_helpers.R:151; mod_lipid_qc_duplicates_server_helpers.R:373,
        # :415; mod_lipid_qc_itsd.R:223
        notification_ids = c(intensity = "lipid_intensity_filter_working",
            duplicates = "lipid_dup_resolve_working", itsd = "is_analysis_working"))
}

#' The lipid design tab (mod_lipid_design_*.R at 3c657c14). `validate` is the core's save preflight
#' (`validateOmicsDesignSave`, as `validateLipidDesignDaPreflight` was called) and `build_object` the
#' core's object (`lipidDesignObject`, the inline `createLipidomicsAssayData` call), KIT-103 W2.
#' @noRd
omicGuiLipidDesign <- function() {
    c(omicGuiFeatureDesign("lipid", "lipidomics", "lipidomics", "lipid_id_col"), list(
        # mod_lipid_design_ui.R:39
        help = paste("Use the tools below to define your experimental groups and contrasts for",
            "differential analysis. All assays (Pos/Neg modes) share the same design matrix."),
        # mod_lipid_design_builder_server.R:43; mod_lipid_design_builder_state_helpers.R:298, :395
        run_order = gtools::mixedsort, run_order_rule = "mixedsort", group_levels = NULL,
        # utils_workflow_lipid_declarations.R:211; mod_lipid_design_builder_helpers.R:175-183
        record_provenance = recordLipidDesignProvenance, state_name = "lipid_raw_data_s4",
        state_description = "Initial LipidomicsAssayData S4 object created after design matrix",
        validate = function(results, columnMapping) {
            validateOmicsDesignSave(results, columnMapping, "lipidomics")
        },
        build_object = lipidDesignObject,
        # mod_lipid_design_ui.R:42-49; the builder's removed runs are not returned
        # (mod_lipid_design_builder_server.R:180-183)
        record_excluded_runs = FALSE, multi_assay_banner = TRUE,
        # mod_lipid_design_import_helpers.R:230, :388 (func_lipid_qc.R:134), :293-305, :383
        design_import = omicGuiFeatureDesignImport(c("Lipid name", "Lipid.name", "Name", "name"),
            updateLipidFiltering, preflight = TRUE, class = "LipidomicsAssayData",
            record_design_file = TRUE, check_files_first = FALSE, preflight_contrasts = TRUE)))
}

#' The lipid panel (the provenance abbreviations are the file header's).
#' @noRd
omicGuiLipidImport <- function() {
    c(omicGuiLipidColumnIds(), list(source_panel = "feature_table",
        # ui:359, :87
        title = "Lipidomics Data Import", format_heading = "Step 1: Select Vendor Format",
        # ui:77-84, :93
        formats = c("MS-DIAL" = "msdial", "Progenesis QI" = "progenesis", "XCMS" = "xcms",
            "Compound Discoverer" = "compound_discoverer", "LipidSearch" = "lipidsearch",
            "Other/Custom" = "custom"),
        default_format = "msdial",
        # selection:23
        detect = detectLipidomicsFormat,
        # The core's preview and import (utils_workflow_omics_import.R).
        preview = previewOmicsImport, process = importOmicsTables,
        # ui:100, :105, :168, :189, :191, :210, :233
        help = "Import one or more assay files (e.g., LCMS_Pos, LCMS_Neg, GCMS)",
        assay_name = "LCMS_Pos", id_label = "Lipid ID Column", id_custom_label = "Lipid ID Column",
        id_custom_placeholder = "e.g., Compound_ID", annotation_label = "Annotation Column",
        annotation_custom_placeholder = "e.g., Lipid_Name",
        # ui:113-125. The reader per detected format (processing:395-406; the preview,
        # selection:39-43); the second assay's reader is given the id and annotation columns
        # (processing:408-425).
        second_assay = TRUE,
        readers = list(default = importLipidMSDIALData, lipidsearch = importLipidSearchData),
        second_assay_columns = TRUE,
        # processing:470-472; ui:255-285 (in order); output:5
        id_mapping_key = "lipid_id_col",
        aux_inputs = c("exclude_norm", "sanitize_names", "is_pattern"),
        validate_mapping = validateLipidColumnMapping,
        # ui:96, :110-111, :123-124, :171, :213, :271, :283, :336
        testids = c(vendor_format = "lipid-import-vendor-format",
            assay1_name = "lipid-import-assay1-name", assay1_file = "lipid-import-assay1-file",
            assay2_name = "lipid-import-assay2-name", assay2_file = "lipid-import-assay2-file",
            id_col = "lipid-import-lipid-id-col", annotation_col = "lipid-import-annotation-col",
            sanitize_names = "lipid-import-sanitize", is_pattern = "lipid-import-is-pattern",
            process = "lipid-import-process"),
        status_modes = omicGuiFeatureStatusModes(),
        # processing:5-13, :153
        format_labels = c(msdial = "MS-DIAL", lipidsearch = "LipidSearch",
            progenesis = "Progenesis QI", xcms = "XCMS",
            compound_discoverer = "Compound Discoverer", trx_counts = "Count table"),
        validation_noun = "Lipids",
        assay_names = omicGuiLipidAssayNames(),
        # processing:646, :683, :561, :296, :721
        notification_id = "lipid_import_working", refusal_prefix = "Invalid lipidomics import",
        summary_file = "lipidomics_import_summary.tsv", count_key = "n_lipids",
        success_text = "Data imported successfully!"))
}

# ------------------------------------------------------------------------------------------------
# Metabolomics.
# ------------------------------------------------------------------------------------------------

#' The metabolomics configuration (mod_metabolomics.R; the readers
#' utils_workflow_metab_declarations.R:111, :114).
#' @noRd
omicGuiConfigMetabolomics <- function() {
    list(shell = omicGuiFeatureShell("metabolomics", "metab", "Metabolomics", "flask",
            "METABOLOMICS", " Differential Analysis",
            omicGuiFeatureBus("metabolite", "metabolite_assay_data")),
        binding = omicGuiFeatureBinding(metabolomicsDeclarationBinding, readMetabWorkflowState,
            readMetabWorkflowMetadata),
        import = omicGuiMetabImport(),
        design = omicGuiMetabDesign(),
        qc = omicGuiMetabQc(),
        norm = omicGuiMetabNorm(),
        da = omicGuiMetabDa(),
        # KIT-107 D9: no enrichment tab.
        enrichment = NULL,
        summary = omicGuiMetabSummary())
}

#' The metab Summary & Export tab at 5be34ee5 (mod_metab_summary*.R; KIT-108 D1): the module id and test ids (the
#' shell's tab; mod_metab_summary.R:57, :70, :88, :96, :130), the placeholder and description (:49, :51), the final
#' state: the first of the six data states, head first, read as ids on a managed session (`readMetabWorkflowIds`,
#' utils_workflow_metab_declarations.R:124), else `getState` (mod_metab_summary_observer_helpers.R:310-328), the record
#' whenever there is a bus (:395, :885), the copy (`omicGuiFeatureSummaryCopy`), the integration file and its managed
#' operation, else the final object's RDS in process (:420-425, :459-467), the template, its directory and its label
#' (:757, :784; mod_metab_summary_server_helpers.R:1-27), the report flow's order (the template check's progress only
#' when the template is missing, :769-773; after a failed retrieval the report still renders, :819-846; the record bundle
#' before the success notice, :883-900), no export step and no text of its own. Wired at the W3 merge:
#' `session_state` (`prepareMetabSummarySessionState`, mod_metab_summary_observer_helpers.R:70-145).
#' @noRd
omicGuiMetabSummary <- function() {
    list(module_id = "summary", testid = "metab", placeholder = "e.g., my_metabolomics_analysis",
        description = "Full metabolomics analysis workflow with normalization and DA",
        final_state = list(states = c("metab_correlation_filtered", "metab_norm_complete", "metab_ruv_corrected",
            "metab_normalized", "loaded_for_de", "metab_qc_complete"), read_ids = readMetabWorkflowIds),
        record_class = NULL,
        copy = omicGuiFeatureSummaryCopy(),
        # KIT-135: the in-process writer writes the RDS as before, then the declared integration input beside it.
        integration = list(filename = "%s_%s_final_s4.RDS", operation = "metab_integration_export", write = metabSummaryIntegrationWrite,
            requires = NULL),
        report = list(dir = "metabolomics", template = "metabolomics_report.rmd", label = "Metabolomics", status = NULL,
            render_after_failed_template = TRUE, checking_progress_always = FALSE, record_after_success = FALSE),
        export_steps = list(),
        texts = character(),
        session_state = prepareMetabSummarySessionState)
}

#' The metab DA tab at b75561de (the shared values are `omicGuiFeatureDa`'s). The columns
#' (mod_metab_da_display_helpers.R:131-143), the summary's noun and label (:93-99), the nouns (mod_metab_da.R:222,
#' :249, :328; mod_metab_da_display_helpers.R:224), the table filter where "All" and "Combined" mean every assay and a
#' missing `comparison` or `friendly_name` column matches nothing (:13-29), the static volcano and heatmap kernels
#' (:423-440, :480, :491). In mod_metab_da_observer_helpers.R: the heatmap's save prefix (:37), the in-process
#' reference DA (:739), the source and its no-data text (:303, :309, :337), the session file and assay check, in
#' order (:448, :789-803), the preflight before the running notice (:306-382, :384-424; `validateMetabDesignDaPreflight`,
#' func_omics_design_shared.R:788) and the download's stem (:264).
#' @noRd
omicGuiMetabDa <- function() {
    da <- omicGuiFeatureDa("metab", list(
        columns = list(id = "metabolite_id", label = "metabolite_name", fold_change = "logFC", pvalue = "raw_pvalue",
            significance = "fdr_qvalue"),
        # KIT-172: the count of rows not tested because a compared group has no observed value.
        summary = list(noun = "metabolites", significance_label = "Q", not_tested = TRUE),
        nouns = c(top_n = "Top N Metabolites:", scaling = "Row (metabolite) scaling", labels = "Show Metabolite Labels",
            cluster = "metabolites"),
        table_filter = "all_or_combined",
        static_volcano = generateMetabDAVolcanoStatic,
        heatmap = list(kernel = generateMetabDAHeatmap, labels_argument = "show_metabolite_names", save_prefix = "metab_"),
        run_in_process = runMetabReferenceDaAnalysis,
        source = list(read_metadata = readMetabWorkflowMetadata, class = "MetaboliteAssayData",
            no_data = "No metabolomics data loaded. Please load a filtered session first."),
        session_load = list(file = "metab_filtered_session_data_latest.rds", assay_rule = "identical",
            assay_message = "Session assay_names do not match current_s4_object assays: session=%s; object=%s",
            # KIT-106 W2 (D7, from W1's captures): metab checks the assays before the session is adopted; lipid after.
            assay_check_first = TRUE),
        preflight_first = TRUE,
        download = "metabolomics_da_results_"))
    # KIT-106 W4: the frame calls the preflight with the state read for the run, the contrasts and the formula
    # (`omicDaPreflight`, mod_omic_da.R); metab's validator takes the state's design and assays, and is not called without
    # a formula (mod_metab_da_observer_helpers.R:351-362 at fe97aac1). It is looked up by name when it runs.
    da$preflight <- function(currentS4, contrastsTbl, formulaString) {
        if (is.null(formulaString)) return(list(valid = TRUE))
        validateMetabDesignDaPreflight(designMatrix = currentS4@design_matrix, assayList = currentS4@metabolite_data,
            contrastsTbl = contrastsTbl, formulaString = formulaString, sampleIdCol = currentS4@sample_id,
            groupCol = currentS4@group_id, requireContrasts = TRUE)
    }
    da
}

#' The metab normalisation tab (mod_metab_norm*.R at 9a98e765; the shared values are
#' `omicGuiFeatureNorm`'s, written by hand in mod_metab_norm_ui.R:25-658). The RUV noun
#' (mod_metab_norm_ui.R:209, :216) and the test ids (:106, :125, :256, :281). The export: the
#' `metab_` files (mod_metab_norm_session_helpers.R:181, :576, :757), the `metabolite_counts` field
#' (:538) and the per-assay sample count of the columns less the id and annotation columns (:480-497;
#' P-18). The ITSD selection table: buildItsdSelectionTable (func_metab_norm_support_helpers.R:111),
#' as the tab calls it (resolveMetabNormManualItsdFeatureIds, mod_metab_norm_pipeline_helpers.R:177).
#' D7: skip-correlation sets the flags after the store commit (mod_metab_norm_session_helpers.R:349-361,
#' :427-436); the managed final PCA's placeholders are both size 5 (mod_metab_norm_render_helpers.R:420-422).
#' KIT-105 W5 (D7, at 037bc246): the export's progress shows the gathering step and then each of the saver's files,
#' with no completion step, and the summary carries the time it is written (mod_metab_norm_session_helpers.R:566, :580,
#' :672, :720-749); a failed managed correlation is also notified as an error (:15-21, :91-94).
#' @noRd
omicGuiMetabNorm <- function() {
    omicGuiFeatureNorm("metab", "Metabolites",
        export = list(prefix = "metab_", counts_field = "metabolite_counts", samples = "design",
            progress = list(gather = "Gathering data...", complete = character(0), steps = TRUE), summary_time = "written"),
        itsd_table = function(assay_data, id_col, annotation_cols) {
            buildItsdSelectionTable(assay_data = assay_data, metabolite_id_col = id_col, annotation_cols = annotation_cols)
        },
        skip_correlation_flags = "after_commit", final_pca = list(preparing_size = 5, error_size = 5),
        managed_correlation_notify = TRUE)
}

#' The metab QC tab (mod_metab_qc*.R at b423b764). `itsd` holds `analyzeMetabQcItsdData`'s entries
#' (mod_metab_qc_itsd.R): the name columns (`metabQcItsdColumns`, :36-44), the metrics kernel with
#' its id-column argument (:102-107), the quant getter (:125) and the refusal's tip (:157).
#' @noRd
omicGuiMetabQc <- function() {
    list(title = "Metabolite Quality Control & Filtering", # mod_metab_qc.R:44
        tabset_id = "metab_qc_tabs", # mod_metab_qc.R:116
        groups = omicGuiFeatureQcGroups("METABOLOMICS"),
        # Started on the trigger or when the state's class is the binding's (mod_metab_qc.R:125-167).
        auto_init = TRUE, complete_state = NULL,
        # mod_metab_qc_intensity.R:206; mod_metab_qc_s4_server_helpers.R:398 (func_metab_qc.R:268)
        filtering_figure = updateMetaboliteFiltering,
        itsd = list(name_columns = c("Metabolite name", "Name", "metabolite"),
            metrics = function(assay_data, is_pattern, id_col, sample_id_col) {
                getInternalStandardMetrics(assay_data = assay_data, is_pattern = is_pattern,
                    metabolite_id_col = id_col, sample_id_col = sample_id_col)
            },
            quant = getMetaboliteQuantData,
            tip = "Tip: Check if your data has a column with metabolite names (e.g., 'Metabolite name')"),
        # Only an empty assay list is refused (mod_metab_qc_s4_server_helpers.R:335-339); only
        # `quality_control` is set (:368-378); a per-assay table (mod_metab_qc_s4.R:73-77).
        finalize_zero_rows = FALSE, finalize_opens_norm = FALSE, finalize_per_assay = TRUE,
        # D11, the parity check (KIT-104 W6): the shared alignment with lipid's arguments; today's
        # metab finalise checked the class, the assay list and the design rows only
        # (mod_metab_qc_s4_server_helpers.R:320-350 at 827f2328).
        finalize_alignment = function(shell) omicsQcFinalizeAlignment(shell, "metabolomics"),
        texts = omicGuiFeatureQcTexts("metabolite", "MetaboliteAssayData"),
        # mod_metab_qc_intensity.R:365; mod_metab_qc_duplicates_server_helpers.R:418;
        # mod_metab_qc_itsd.R:402
        notification_ids = c(intensity = "metab_intensity_filter_working",
            duplicates = "metab_dup_resolve_working", itsd = "is_analysis_working"))
}

#' The metab design tab (mod_metab_design*.R at 3c657c14). `validate` is the core's save preflight
#' (`validateOmicsDesignSave`, as `validateMetabDesignDaPreflight` was called) and `build_object` the
#' core's object (`metabDesignObject`, the inline `createMetaboliteAssayData` call), KIT-103 W2.
#' @noRd
omicGuiMetabDesign <- function() {
    c(omicGuiFeatureDesign("metab", "metabolomics", "metabolomics", "metabolite_id_col"), list(
        # mod_metab_design.R:59
        help = paste("Use the tools below to define your experimental groups and contrasts for",
            "differential analysis. All assays (Pos/Neg modes) share the same design matrix."),
        # mod_metab_design_builder_server.R:237
        run_order = gtools::mixedsort, run_order_rule = "mixedsort", group_levels = NULL,
        # utils_workflow_metab_declarations.R:193; mod_metab_design.R:341-351
        record_provenance = recordMetabDesignProvenance, state_name = "metab_raw_data_s4",
        state_description = "Initial MetaboliteAssayData S4 object created after design matrix",
        validate = function(results, columnMapping) {
            validateOmicsDesignSave(results, columnMapping, "metabolomics")
        },
        build_object = metabDesignObject,
        # mod_metab_design.R:62-68
        record_excluded_runs = FALSE, multi_assay_banner = TRUE,
        # mod_metab_design_import_helpers.R:250, :433 (func_metab_qc.R:268), :424; the import ran
        # no DA preflight until P-11 (GUI audit R4): it now refuses an invalid design as lipid's does,
        # and keeps the file's contrasts
        design_import = omicGuiFeatureDesignImport(c("Metabolite name", "Metabolite.name", "Name",
            "name"), updateMetaboliteFiltering, preflight = TRUE, class = "MetaboliteAssayData",
            record_design_file = TRUE, check_files_first = TRUE, preflight_contrasts = FALSE)))
}

#' The metab panel (the provenance abbreviations are the file header's). Its ids, labels and
#' messages stay metab's; its layout is the shared panel's (KIT-102 contract).
#' @noRd
omicGuiMetabImport <- function() {
    list(source_panel = "feature_table",
        # ui:26, :31
        title = "Metabolomics Data Import", format_heading = "Step 1: Select Vendor Format",
        # ui:36-43
        formats = c("MS-DIAL" = "msdial", "Progenesis QI" = "progenesis", "XCMS" = "xcms",
            "Compound Discoverer" = "compound_discoverer", "Other/Custom" = "custom"),
        default_format = "msdial",
        # processing:323
        detect = detectMetabolomicsFormat,
        # The core's preview and import (utils_workflow_omics_import.R).
        preview = previewOmicsImport, process = importOmicsTables,
        # ui:50, :62, :176, :195, :197, :215, :236
        help = "Import one or more assay files (e.g., LCMS_Pos, LCMS_Neg, GCMS)",
        assay_name = "LCMS_Pos", id_label = "Metabolite ID Column",
        id_custom_label = "Metabolite ID Column", id_custom_placeholder = "e.g., Compound_ID",
        annotation_label = "Annotation Column",
        annotation_custom_placeholder = "e.g., Metabolite_Name",
        # ui:95-137. Every format is read with the MS-DIAL reader (processing:19, :36; the preview,
        # :324-325); the second assay's reader is given the file only.
        second_assay = TRUE, readers = list(default = importMetabMSDIALData),
        second_assay_columns = FALSE,
        # processing:69; ui:259-273 (in order); display:231
        id_mapping_key = "metabolite_id_col", aux_inputs = c("is_pattern", "sanitize_names"),
        validate_mapping = validateMetabColumnMapping,
        # ui:175, :194, :182, :202
        inputs = c(id_col = "metabolite_id_col", id_col_custom = "metabolite_id_col_custom"),
        outputs = c(id_status = "metabolite_id_status",
            id_status_custom = "metabolite_id_status_custom"),
        # ui:46, :64, :70, :108, :114, :318
        testids = c(vendor_format = "metab-import-vendor-format",
            assay1_name = "metab-import-assay1-name", assay1_file = "metab-import-assay1-file",
            assay2_name = "metab-import-assay2-name", assay2_file = "metab-import-assay2-file",
            process = "metab-import-process"),
        status_modes = omicGuiFeatureStatusModes(),
        # display:46-52, :264
        format_labels = c(msdial = "MS-DIAL", progenesis = "Progenesis QI", xcms = "XCMS",
            compound_discoverer = "Compound Discoverer"),
        validation_noun = "Metabolites",
        # A second assay counts only with both a file and a name (processing:24-29); the path rule
        # runs in the validation (:137-138) and the duplicate message lists the names (:29).
        assay_names = list(rule = "skip_incomplete", path_pattern = "[/\\\\\r\n\t]",
            messages = c(required = "Metabolomics import requires non-empty assay names.",
                path = "Assay names must not contain path separators or control characters.",
                duplicate = "Duplicate metabolomics assay names are not allowed: ")),
        # processing:557, :588, :270, :79, :435
        notification_id = "metab_import_working", refusal_prefix = "Invalid metabolomics import",
        summary_file = "metabolomics_import_summary.tsv", count_key = "n_metabolites",
        success_text = "Data imported successfully!")
}

# ------------------------------------------------------------------------------------------------
# Transcriptomics.
# ------------------------------------------------------------------------------------------------

#' The transcriptomics configuration: the transcriptomics profile's fields at c1041c56
#' (mod_transcriptomics.R:32-49; KIT-102 W5 removed the import-only ones) and the lipid panel's
#' values it rendered under that profile (design 4.12).
#' @noRd
omicGuiConfigTranscriptomics <- function() {
    list(shell = omicGuiFeatureShell("transcriptomics", "trx", "Transcriptomics", "chart-line",
            "TRANSCRIPTOMICS", " Differential Expression", omicGuiFeatureBus("transcriptomics")),
        # The metadata reader as the profile's (:48), with the reference the other readers take.
        binding = omicGuiFeatureBinding(transcriptomicsDeclarationBinding, readTrxWorkflowState,
            function(manager, reference = NULL) {
                readWorkflowOmicsMetadata(transcriptomicsDeclarationBinding, manager, reference)
            }),
        import = omicGuiTrxImport(),
        design = omicGuiTrxDesign(),
        qc = omicGuiTrxQc(),
        norm = omicGuiTrxNorm(),
        da = omicGuiTrxDa(),
        # KIT-107 D9: no enrichment tab.
        enrichment = NULL,
        summary = omicGuiTrxSummary())
}

#' The trx Summary & Export tab at 5be34ee5: the lipid tab under the profile's summary entries (KIT-108 D1, D4;
#' the profile's normalisation, DA and summary entries, mod_transcriptomics.R:45-57; the profile, :39-41): the module id and test ids (the
#' shell's tab; the profile's `testid`, :39; mod_lipid_summary.R:41), the placeholder and description (:52-53), the record
#' whenever there is a bus, the copy (`omicGuiFeatureSummaryCopy`), the integration file and its managed operation
#' (:54; mod_lipid_summary_workflow_helpers.R:89, :124-127) gated on the accepted RUVIII export with its message
#' (mod_transcriptomics.R:56; `trxSummaryIntegrationReady`, mod_trx_summary.R:143; mod_lipid_summary_workflow_helpers.R:
#' 119-122), the template, its directory and its label (mod_transcriptomics.R:39; mod_lipid_summary_state_helpers.R:17,
#' :20), lipid's report flow (mod_lipid_summary_workflow_helpers.R:352-534), the RUVIII export step
#' (`omicGuiTrxRuviiiExportStep`) and the save's text before that export
#' (mod_lipid_summary_workflow_helpers.R:121). Wired at the W3 merge: `final_state$context` (`trxSummaryContext`,
#' mod_trx_summary.R:121-126), `integration$write` (`trxSummaryIntegrationWrite`, :136-139) and `session_state` (the
#' lipid payload with `trxSummaryStateObject`, :130-132; mod_lipid_summary_state_helpers.R:158, :168, :172).
#' @noRd
omicGuiTrxSummary <- function() {
    list(module_id = "summary", testid = "trx", placeholder = "e.g., my_transcriptomics_analysis",
        description = "Full transcriptomics analysis workflow with normalization, DE and the RUVIII integration export",
        final_state = list(context = trxSummaryContext),
        record_class = NULL,
        copy = omicGuiFeatureSummaryCopy(),
        integration = list(filename = "transcriptomics_normalized_data.tsv", operation = "trx_integration_export",
            write = trxSummaryIntegrationWrite, requires = trxSummaryIntegrationRequires,
            missing_message = "The integration export needs the RUVIII export: run it in the panel above, then save again."),
        report = list(dir = "transcriptomics", template = "transcriptomics_report.rmd", label = "Transcriptomics",
            status = NULL, render_after_failed_template = FALSE, checking_progress_always = FALSE,
            record_after_success = FALSE),
        export_steps = list(omicGuiTrxRuviiiExportStep()),
        texts = c(integration_not_written = "not written (no integration export yet)"),
        session_state = trxSummarySessionState)
}

#' Trx's RUVIII export at 5be34ee5 as an export step (KIT-108 D4; mod_trx_summary.R:23-117): the panel
#' (`trxSummaryExportUi`, :43-54: title, help text, the choice's output, the button, the status and the figure), the
#' choice (the radio input, its label and default, the unavailable text; `trxSummaryServer`, :99-104), the status's idle
#' text (:110), the figure (`ruviii_before_and_after.png` in `publication_graphs_dir`, :111-115), the managed product
#' and its acceptance on the head (`acceptTrxRuviiiExport`, utils_workflow_trx_declarations.R:411; mod_trx_summary.R:
#' 76-83), the in-process run (`runTrxReferenceRuviiiExport`, utils_workflow_trx_declarations.R:815; mod_trx_summary.R:
#' 88), the publications (:73) and the error notification (:68-71). Wired at the W3 merge: `choice$choices`
#' (`trxRuviiiMChoices`, :25-29), `choice$available` (`trxHeadHasTechnicalReplicates`, :33-40), `parameters`
#' (`trxRuviiiExportParameters`, new: :76-78) and `status$text` (`trxRuviiiStatusText`, :58-61).
#' @noRd
omicGuiTrxRuviiiExportStep <- function() {
    list(id = "ruviii_export", title = "RUVIII Integration Export",
        help = c("RUVIII (ruv::RUVIII) on the TMM logCPM, with the RUVs inferred negative controls and k, for the",
            "integration export only: the DA results are not changed. M is the replicate structure: the group (the notebook),",
            "or the design's technical replicates when it has them."),
        choice = list(output = "ruviii_m_ui", input = "ruviii_m_source", label = "Replicate structure M:",
            selected = "group", unavailable = c("Technical replicates: not available; the design has no technical",
                "replicates (trx_ruviii_replicates_invalid)."), choices = trxRuviiiMChoices,
            available = trxHeadHasTechnicalReplicates),
        button = list(input = "run_ruviii_export", label = "Run RUVIII Export", testid = "trx-summary-ruviii-export",
            icon = "play", class = "btn-primary"),
        status = list(output = "ruviii_status", idle = "Not run yet.", text = trxRuviiiStatusText),
        figure = list(output = "ruviii_before_after", file = "ruviii_before_and_after.png", dir = "publication_graphs_dir",
            height = "500px", image = c(content_type = "image/png", width = "100%", height = "auto"),
            alt = "RUVIII before and after", alt_idle = "RUVIII export not run yet"),
        operation = "trx_ruviii_export", accept = acceptTrxRuviiiExport, run_reference = runTrxReferenceRuviiiExport,
        parameters = trxRuviiiExportParameters,
        publications = list(list(operation = "trx_qc_figures", parameters = list(stage = "ruviii")),
            list(operation = "trx_composite_figure", parameters = list(figure = "ruviii_before_and_after"))),
        error = list(prefix = "RUVIII export error:", type = "error", duration = 15))
}

#' The trx DA tab at b75561de: lipid's tab under the profile's DA entries (mod_transcriptomics.R:49-53, mod_trx_da.R).
#' No load section (mod_lipid_da_ui.R:60) and no heatmap (:176); the test ids by the profile's prefix (:11). The model
#' shown in place of the formula input, its text and the realised formula of the managed settings (mod_trx_da.R:24-39;
#' mod_lipid_da_analysis_helpers.R:197); the threshold and display-only fold-change labels (mod_transcriptomics.R:49-50;
#' the binding's `da_fixed`, utils_workflow_trx_declarations.R:722); the columns (mod_transcriptomics.R:51); the
#' summary (mod_trx_da.R:157); one assay, on which the volcano opens (:123); the static volcano (:74-78); the in-process
#' reference DA (mod_transcriptomics.R:53); the source's class and no-data text (mod_trx_da.R:129-130); the texts
#' (:51-54, and lipid's contrasts text); the download's stem (:167). Wired at the KIT-106 W2 merge, after W2 moved them
#' from mod_trx_da.R to utils_workflow_trx_da.R: the model text (`trxDaModelText`), the static volcano (`trxDaStaticVolcano`), the source read (`trxDaReadMetadata`), the
#' contrast entries (`trxDaContrastEntries`, `trxDaContrastEntry`) and the export report (`trxDaExport`).
#' @noRd
omicGuiTrxDa <- function() {
    list(module_id = "de", status_key = "differential_analysis",
        testids = c(run = "trx-da-run", download = "trx-da-download", results_tabs = "trx-da-results-tabs"),
        model = list(editable = FALSE, heading = "Model", output = "formula_display",
            help = c("edgeR quasi-likelihood: estimateDisp and glmQLFit(robust = TRUE), then glmQLFTest per contrast on the",
                "one fit; Benjamini-Hochberg FDR. The formula is realised from the design: ~ 0 + group, plus the RUV factors",
                "W_1..W_k as covariates when the normalisation added them. No fold-change test is applied."),
            formula = trxRealisedFormula, text = trxDaModelText),
        threshold = list(label = "FDR (BH) threshold:", value = 0.05, min = 0.001, max = 0.2, step = 0.005),
        lfc = list(label = "Fold-change line (display only; no test):", value = 0, min = 0, max = 2, step = 0.1,
            tested = FALSE),
        columns = list(id = "GeneID", label = "GeneName", fold_change = "logFC", pvalue = "PValue", significance = "FDR",
            significance_label = "FDR (BH)"),
        summary = list(noun = "genes", significance_label = "FDR (BH)"),
        tabs = c(volcano = "Volcano Plot", table = "DA Results Table"),
        assays = "single", table_filter = "all_only",
        views = list(static_volcano = list(output = "volcano_static", plot = trxDaStaticVolcano, show_labels = TRUE,
            n_labels = 15),
            table = list(significance = "significant", max_rows = 1000)),
        contrast_entries = trxDaContrastEntries, contrast_entry = trxDaContrastEntry,
        run = list(managed = runWorkflowOmicsManagedDa, in_process = runTrxReferenceDaAnalysis,
            export_report = trxDaExport),
        source = list(class = "TranscriptomicsAssayData", read_metadata = trxDaReadMetadata,
            no_data = "Run the normalisation first: the differential analysis runs on the normalised data."),
        preflight_first = FALSE,
        texts = c(status_ready = "Ready. The DA runs on the normalised (RUVs) state.",
            status_waiting = "Ready. The DA runs on the normalised (RUVs) state.",
            contrasts_none = "No contrasts defined.\nLoad a filtered session or define contrasts in the design tab."),
        download = list(stem = "transcriptomics_da_results_", extension = "csv"))
}

#' The trx normalisation tab (KIT-110): `trxNormOptionsPanel`, `mod_trx_norm_server` and the profile's
#' `norm_*` fields at 38d9d25e (mod_trx_norm.R:25-53, :215-235; mod_transcriptomics.R:48-52), as values. The options:
#' TMM and logCPM (mod_trx_norm.R:40-41); the RUV grid, its modes and default (:43-44) and its inputs' values
#' (:28-31, :35-36); the actions with the reset's test id (:48-52). The profile's tabs, captions and RUV texts, one
#' assay row (mod_transcriptomics.R:48-52; mod_lipid_norm_ui.R:26 at 38d9d25e). The run and reset are trx's own
#' handlers (D4; the core's, utils_workflow_trx_norm.R, from mod_trx_norm.R:57-64, :130-187, :191-207). The stage grid:
#' the one count table (mod_trx_norm.R:219), the file stem `trx` and the stage prefixes (:222-223), and the
#' registrars' default of two assay rows (:221-224, P-21). No pre-QC on tab entry: the shell calls trx's server
#' without `selected_tab` (mod_omic_workflow.R:214). No ITSD, correlation-filter tab or export panel (D7).
#' @noRd
omicGuiTrxNorm <- function() {
    list(module_id = "norm",
        panels = list(tmm_logcpm = list(),
            ruv_grid = list(modes = c("Automatic (recommended)" = "automatic", "Manual (advanced users)" = "manual"),
                default = "automatic", percentage_min = 5, percentage_max = 30, percentage_step = 5, k_penalty_weight = 0.05,
                percentage = 10, k = 2),
            actions = list()),
        plot_tabs = c("ruv", "pca", "rle", "log"),
        stage_labels = c("Raw logCPM", "TMM logCPM", "RUVs-corrected logCPM (display only)"), assay_slots = 1L,
        ruv_heading = "RUVs Control Selection (canonical correlation)",
        ruv_help = paste("The CCA panels and the grid of percentages: the controls are inferred negative controls (the least",
            "significant genes in a preliminary group GLM)."),
        run = list(managed = runTrxNormManaged, in_process = runTrxNormInProcess, settings = trxNormGridParameters,
            handler = runTrxNormalisation),
        reset = resetTrxNormalisation,
        pre_qc_on_entry = FALSE,
        assay_names = "counts", image_stem = "trx", stage_prefixes = c(pre_norm = "raw", post_norm = "tmm", ruv_corrected = "ruvs"),
        output_assay_slots = 2L,
        testids = c(ruv_mode = "trx-norm-ruv-mode", run = "trx-norm-run", reset = "trx-norm-reset"),
        ruv_design = omicGuiRuvDesign("trx", list(genes = omicGuiRuvTableTrx())))
}

#' The trx QC tab: at b423b764 the lipid QC orchestrator and finalise under the transcriptomics profile's QC
#' fields (mod_transcriptomics.R:38-41; mod_trx_qc_filter.R): the expression-filter plugin and the
#' shared finalise with one table. Since KIT-104 W6 the profile has no QC field: the tab is the
#' shared shell under this configuration, with the figure, alignment and shell read below.
#' @noRd
omicGuiTrxQc <- function() {
    list(title = "Transcriptomics Quality Control & Filtering", # mod_transcriptomics.R:38
        tabset_id = "trx_qc_tabs", # mod_trx_qc_filter.R:179
        # mod_trx_qc_filter.R:179-180, the titles at :77 and mod_lipid_qc_s4.R:43
        groups = list(list(id = NULL, title = NULL, component = NULL, tabs = list(
            omicGuiQcTab("filter", "Expression Filter", "trx_qc_filter", "TRANSCRIPTOMICS"),
            omicGuiQcTab("s4_finalize", "Finalize QC", "omic_qc_finalize", "TRANSCRIPTOMICS")))),
        # The lipid orchestrator's class check under the profile (mod_lipid_qc.R:151-158).
        auto_init = TRUE, complete_state = NULL, itsd = NULL,
        # The finalise's tracking figure from the worker's metrics, or from the state in process, as
        # trx's finalise drew it (trxQcFinalizeCollaborators, mod_trx_qc_filter.R:234-237 at 827f2328).
        filtering_figure = function(theObject = NULL, metrics = NULL, ...) {
            trxQcFilteringFigure(NULL, metrics %||% trxFilteringMetrics(theObject))
        },
        # D10, the parity check (KIT-104 W6): trx's own design check (func_trx_design.R); today's trx
        # finalise checked the class, the row counts and the design rows only (mod_trx_qc_filter.R:224-229).
        finalize_alignment = trxDesignAlignment,
        # F104-04 (lead decision, option 1): the class's validity refuses a count table with no rows,
        # so the finalise reads a structural shell: the stored slots and design, the count table as its
        # zero-row prototype (`trxReadTemplate`, which runs no validity) and the stored row counts.
        finalize_shell = function(binding, manager) {
            if (!inherits(manager, "DIAStateAdapter")) {
                cli::cli_abort(paste(binding$step_scope, "needs the session's state owner"), class = "saved_entity_scope")
            }
            reference <- getWorkflowStateReference(manager)
            if (is.null(reference)) return(NULL)
            withWorkflowStateScope(manager, function() {
                shell <- trxReadTemplate(manager, reference)
                tables <- readWorkflowStateMetadata(manager, reference)$data$components[[binding$data_slot]]
                attr(shell, "assay_rows") <- vapply(tables, function(table) as.integer(table$metadata$shape[[1L]]), integer(1))
                shell
            })
        },
        # Any zero-row assay is refused (mod_trx_qc_filter.R:226); lipid's tab status opens
        # normalisation (mod_lipid_qc_s4.R:420-436); one table, no per-assay table (:75).
        finalize_zero_rows = TRUE, finalize_opens_norm = TRUE, finalize_per_assay = FALSE,
        # mod_trx_qc_filter.R:202-203, :212, :227 (one refusal for both checks)
        texts = c(finalize_total_label = "Total Genes:", finalize_id_label = "Gene ID Column:",
            finalize_total = "Total genes retained: %d",
            finalize_no_assays = "Current TranscriptomicsAssayData object has no genes or no design rows to finalize.",
            finalize_no_design = "Current TranscriptomicsAssayData object has no genes or no design rows to finalize."),
        notification_ids = stats::setNames(character(), character()))
}

#' The trx design tab: at 3c657c14 the lipid frame and builder under the transcriptomics profile's design
#' fields (mod_transcriptomics.R:36-41; design 4.12; KIT-103 W5 removed them from the profile).
#' `group_levels` (`trxDesignGroupLevels`), `validate` (`trxDesignCheck`) and `build_object`
#' (`trxDesignObject`) are func_trx_design.R's (KIT-103 W2).
#' @noRd
omicGuiTrxDesign <- function() {
    design <- c(omicGuiFeatureDesign("trx", "transcriptomics", "transcriptomics", "lipid_id_col"), list(
        # mod_transcriptomics.R:36-38
        help = paste("Samples are matched to the count columns by sample name only, in the counts",
            "file's order; there is no positional mapping. Each group needs at least two samples.",
            "The order in which you add the factor levels is the group level order the analysis",
            "uses."),
        # mod_transcriptomics.R:39-41
        run_order = identity, run_order_rule = "input",
        record_provenance = recordTrxDesignProvenance, state_name = "trx_raw_data_s4",
        state_description = "Initial TranscriptomicsAssayData S4 object created after design matrix",
        group_levels = trxDesignGroupLevels, validate = trxDesignCheck, build_object = trxDesignObject,
        # The removed runs on the result (mod_lipid_design_builder_server.R:180-183); no banner and
        # no import button under a profile (mod_lipid_design_ui.R:34, :42).
        record_excluded_runs = TRUE, multi_assay_banner = FALSE, design_import = NULL))
    design$testids <- character()
    design
}

#' The trx panel: at c1041c56 the lipid panel under the trx profile (`ui` with the profile, the
#' lipid `processing` called from trx:139); now the shared feature-table panel.
#' @noRd
omicGuiTrxImport <- function() {
    c(omicGuiLipidColumnIds(), list(source_panel = "feature_table",
        # trx:35; ui:87
        title = "Transcriptomics Data Import", format_heading = "Step 1: Format",
        # trx:34-35
        formats = c("Count table (gene ID, gene name, sample counts)" = "trx_counts"),
        default_format = "trx_counts",
        # One format, reported with full confidence (trx:89; now `previewTrxCountsImport()`).
        detect = function(headers, filename = NULL) list(format = "trx_counts", confidence = 1),
        # The count table's preview and the count policy (func_trx_counts.R); the core's import.
        preview = previewTrxCountsImport, process = importOmicsTables, check = checkTrxImportCounts,
        # trx:37-39, :36; ui:189, :191, :233 (the lipid texts, unchanged under the profile)
        help = paste("One delimited count table: a gene ID column, a gene name column and one",
            "column per sample, named exactly as in the design. Counts must be whole numbers: a",
            "missing, negative or fractional count is refused, with no rounding. Sample names are",
            "kept as they are in the file."),
        assay_name = "counts", id_label = "Gene ID Column", id_custom_label = "Lipid ID Column",
        id_custom_placeholder = "e.g., Compound_ID", annotation_label = "Gene Name Column",
        annotation_custom_placeholder = "e.g., Lipid_Name",
        # ui:113; trxDesignObject reads `lipid_id_col` (trx:195); ui:255
        second_assay = FALSE, id_mapping_key = "lipid_id_col", aux_inputs = character(),
        # output:5; trx:117
        validate_mapping = validateLipidColumnMapping,
        # ui:96, :110-111, :171, :213, :336
        testids = c(vendor_format = "trx-import-format", assay1_name = "trx-import-assay1-name",
            assay1_file = "trx-import-assay1-file", id_col = "trx-import-id-col",
            annotation_col = "trx-import-annotation-col", process = "trx-import-process"),
        status_modes = omicGuiFeatureStatusModes(),
        # processing:11, :153
        format_labels = c(trx_counts = "Count table"), validation_noun = "Lipids",
        assay_names = omicGuiLipidAssayNames(),
        # processing:646; trx:119, :144; processing:296, :721
        notification_id = "lipid_import_working", refusal_prefix = "Invalid transcriptomics import",
        summary_file = "transcriptomics_import_summary.tsv", count_key = "n_lipids",
        success_text = "Data imported successfully!"))
}
