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
            preflight = "lgl1", state_description = "chr1", reader = "ref",
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
            assay_slots = "int1", ruv_heading = "chr1", ruv_help = "chr1", run = "named_list",
            pre_qc_on_entry = "lgl1", testids = "named_chr"),
        norm_by_tab = list(itsd = c(itsd_table = "ref"),
            correlation_filter = c(correlation = "named_list", skip_correlation_flags = "chr1", final_pca = "named_list",
                # KIT-105 W5 (D7): whether a failed managed correlation is also shown as an error notification.
                managed_correlation_notify = "lgl1")),
        norm_by_panel = list(export = c(export = "named_list")),
        # KIT-110 D3, D4, D7: an omic's own reset handler; the stage grid's fixed assay names, file stem, stage file
        # prefixes (`norm_stage_prefixes`) and registered assay rows.
        norm_optional = c(reset = "ref", assay_names = "chr", image_stem = "chr1", stage_prefixes = "named_chr",
            output_assay_slots = "int1"),
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
            ruv_iii = c(modes = "named_chr", default = "chr1", percentage = "num1", k = "num1", noun = "chr1"),
            ruv_grid = c(modes = "named_chr", default = "chr1", percentage_min = "num1", percentage_max = "num1",
                percentage_step = "num1", k_penalty_weight = "num1", percentage = "num1", k = "num1"),
            actions = character(), export = c(button = "chr1")),
        # The QC tabs a frame can show, in their order.
        norm_plot_tabs = c("itsd", "ruv", "pca", "density", "rle", "correlation", "correlation_filter", "log"),
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
        # KIT-108 W1 (design 4.10, D1, D4, D9): the Summary & Export tab, every omic. The function fields still in an
        # R/mod_*.R file at the start commit are optional (`*_optional`), checked by type when present, until the W3 and
        # W4 merges wire them. `export_steps` is a list, possibly empty, of steps; `texts` the per-omic texts.
        summary = c(module_id = "chr1", testid = "chr1", placeholder = "chr1", description = "chr1",
            final_state = "named_list", copy = "named_list", integration = "named_list", report = "named_list",
            texts = "named_chr0"),
        summary_optional = c(record_class = "chr1_or_null", session_state = "ref"),
        # The final state: the head's by its `context` function, or the first of `states` with its `read_ids` reader.
        summary_final_state = c(context = "ref", states = "chr", read_ids = "ref"),
        # `required`: NULL (failures not classified) or the required artifacts' display names, checked by the validator.
        summary_copy = c(inherits = "lgl1"),
        summary_copy_optional = c(gate = "ref_or_null"),
        # `operation` and `requires` are present, NULL or a string; `filename` a string or a function; `write` the
        # in-process writer (NULL: the final object's RDS).
        summary_integration = c(requires = "chr1_or_null"),
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
    for (panel in intersect(names(schema$norm_by_panel), panels)) checkOmicGuiFields(norm, schema$norm_by_panel[[panel]], "norm")
    checkOmicGuiFields(norm, schema$norm_optional, "norm", required = FALSE)
    if (!is.null(norm$stage_prefixes) && !setequal(names(norm$stage_prefixes), schema$norm_stage_prefixes)) {
        omicGuiConfigAbort("norm$stage_prefixes", sprintf("must name exactly the stages %s",
            paste(schema$norm_stage_prefixes, collapse = ", ")))
    }
    if (length(norm$stage_labels) != 3L) omicGuiConfigAbort("norm$stage_labels", "must hold the three stages' captions")
    checkOmicGuiFields(norm$run, schema$norm_action, "norm$run")
    checkOmicGuiFields(norm$run, schema$norm_run_optional, "norm$run", required = FALSE)
    if ("correlation_filter" %in% tabs) {
        checkOmicGuiFields(norm$correlation, schema$norm_action, "norm$correlation")
        checkOmicGuiChoice(norm$skip_correlation_flags, "norm$skip_correlation_flags", c("before_commit", "after_commit"))
        checkOmicGuiFields(norm$final_pca, schema$norm_final_pca, "norm$final_pca")
    }
    if ("export" %in% panels) {
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
    if ("filename" %in% names(integration) && !(omicGuiFieldHasType(integration$filename, "chr1") ||
        is.function(integration$filename))) {
        omicGuiConfigAbort("summary$integration$filename", "must be a file name or a function")
    }
    report <- summary$report
    checkOmicGuiFields(report, schema$summary_report, "summary$report")
    checkOmicGuiFields(report, schema$summary_report_optional, "summary$report", required = FALSE)
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
        checkOmicGuiFields(step, schema$summary_export_step_optional, path, required = FALSE)
        checkOmicGuiFields(step$choice, schema$summary_export_choice, paste0(path, "$choice"))
        checkOmicGuiFields(step$choice, schema$summary_export_choice_optional, paste0(path, "$choice"), required = FALSE)
        checkOmicGuiFields(step$button, schema$summary_export_button, paste0(path, "$button"))
        checkOmicGuiFields(step$status, schema$summary_export_status, paste0(path, "$status"))
        checkOmicGuiFields(step$status, schema$summary_export_status_optional, paste0(path, "$status"), required = FALSE)
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
#' value); `lgl1`; `num1` one number; `int1` one positive whole number; `ref` a function;
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
        int1 = is_num1 && value >= 1 && value == round(value),
        chr1 = is_chr1,
        chr1_or_null = is.null(value) || is_chr1,
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
#' @param builder_word The omic word in the builder UI (trx: lipid's).
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
#' `column_mapping.json` is absent, the QC baseline, whether the import runs the DA preflight, and
#' the state's description. The output `import_dir_path` waits for a path (`req`). The reader is
#' the core's (`readOmicsDesignImport`, KIT-103 W2). `record_design_file`: whether the provenance
#' records the imported design file (lipid) or the owner's source directory (metab, P-14).
#' `check_files_first`: whether the files are checked before the modal closes (metab, P-11).
#' @noRd
omicGuiFeatureDesignImport <- function(annotation_candidates, qc_baseline, preflight, class,
                                       record_design_file, check_files_first) {
    list(plugin = "feature_table", import_dir_path = "req",
        id_candidates = list(id = c("Peak ID", "Alignment ID", "Alignment.ID", "peak id",
            "alignment id"), annotation = annotation_candidates),
        qc_baseline = qc_baseline, preflight = preflight,
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
            log2 = list(offset = 1),
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
        managed_correlation_notify = managed_correlation_notify)
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
        # KIT-105 D3 (P-20): proteomics keeps its own normalisation frame until KIT-111.
        norm = NULL,
        da = omicGuiProteomicsDa(),
        enrichment = omicGuiProteomicsEnrichment(),
        summary = omicGuiProteomicsSummary())
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
#' after the success notice and outputs (:534-550; `activateProtSummaryRenderedReport`, :389-408). Wired at the W4
#' merge (their functions are in R/mod_prot_summary_*.R at the start commit): `final_state$context`
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
        final_state = list(),
        record_class = "DIAStateAdapter",
        copy = list(required = NULL, inherits = FALSE),
        integration = list(operation = NULL, requires = NULL),
        report = list(dir = "proteomics", label = NULL, render_after_failed_template = FALSE,
            checking_progress_always = TRUE, record_after_success = TRUE),
        export_steps = list(),
        texts = c(integration_saved = "Saved Integration artifact", progress_detect = "Detecting workflow type..."))
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
#' and save prefix (mod_prot_da_render_heatmap.R:152), no save note (P-25), the status and contrasts texts and the
#' completed status (mod_prot_da_handlers_init.R:320-335; `diaDaStatusText`, utils_workflow_dia_da.R:441), the
#' download (mod_prot_da_render_table.R:151-160). The init, load, run and archive steps are plugins (D5): each is the
#' handler today's server calls (mod_prot_da_server.R:36-45) or the archive chooser (mod_prot_da_ui.R:80-88).
#' @noRd
omicGuiProteomicsDa <- function() {
    list(module_id = "differential_expression", status_key = "differential_expression",
        status_keys = c("differential_abundance", "enrichment_analysis"),
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
        status_text = diaDaStatusText,
        texts = c(status_waiting = "[WAITING] Waiting for analysis...\nClick 'Run DA Analysis' to start",
            contrasts_none = "No contrasts available.\nComplete normalization and\ncorrelation filtering first."),
        download = list(stem = "DA_results_", extension = "zip"),
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
#' (:117-184, stepper :380-386); the bus (:205-245) on `session$userData` (:248); the TMT/LFQ QC
#' bypass (:406-416); the peptide-checkpoint panel (:101-113, :275-304).
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
        chain_rules = list(qc_bypass = c("TMT", "LFQ")), bus = bus, bus_registration = "userData",
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
        integration = list(filename = "%s_%s_final_s4.RDS", operation = "lipid_integration_export", write = NULL,
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
        summary = list(noun = "lipids", significance_label = "Q"),
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
            record_design_file = TRUE, check_files_first = FALSE)))
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
        integration = list(filename = "%s_%s_final_s4.RDS", operation = "metab_integration_export", write = NULL,
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
        summary = list(noun = "metabolites", significance_label = "Q"),
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
        export = list(prefix = "metab_", counts_field = "metabolite_counts", samples = "columns",
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
        # mod_metab_design_import_helpers.R:250, :433 (func_metab_qc.R:268), :424; the import runs
        # no DA preflight
        design_import = omicGuiFeatureDesignImport(c("Metabolite name", "Metabolite.name", "Name",
            "name"), updateMetaboliteFiltering, preflight = FALSE, class = "MetaboliteAssayData",
            record_design_file = FALSE, check_files_first = TRUE)))
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

#' The transcriptomics configuration: `trxOmicProfile()`'s fields at c1041c56
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
#' `trxNormDaSummaryProfile`, mod_transcriptomics.R:45-57; `trxOmicProfile`, :39-41): the module id and test ids (the
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
            write = trxSummaryIntegrationWrite, requires = "trx_ruviii_export",
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
                "one fit; Benjamini-Hochberg FDR. The formula is realised from the design: ~ 0 + group plus the RUVs factors",
                "W_1..W_k as covariates. No fold-change test is applied."),
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
            no_data = "Run the normalisation first: the DA runs on the normalised (RUVs) state."),
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
        testids = c(ruv_mode = "trx-norm-ruv-mode", run = "trx-norm-run", reset = "trx-norm-reset"))
}

#' The trx QC tab: at b423b764 the lipid QC orchestrator and finalise under `trxOmicProfile()`'s QC
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

#' The trx design tab: at 3c657c14 the lipid frame and builder under `trxOmicProfile()`'s design
#' fields (mod_transcriptomics.R:36-41; design 4.12; KIT-103 W5 removed them from the profile).
#' `group_levels` (`trxDesignGroupLevels`), `validate` (`trxDesignCheck`) and `build_object`
#' (`trxDesignObject`) are func_trx_design.R's (KIT-103 W2).
#' @noRd
omicGuiTrxDesign <- function() {
    design <- c(omicGuiFeatureDesign("trx", "transcriptomics", "lipidomics", "lipid_id_col"), list(
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
