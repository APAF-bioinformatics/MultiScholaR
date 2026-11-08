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
# the texts and notification ids the shared sub-tabs render per omic. Fields whose functions are
# still in R/mod_*.R (trx's `filtering_figure`; every `finalize_alignment`, which W2 writes) are
# absent and optional until KIT-104 W8 wires them.

# ------------------------------------------------------------------------------------------------
# The entry point.
# ------------------------------------------------------------------------------------------------

#' The shared GUI layer's configuration of one omic, validated.
#' @param omic_type One of `proteomics`, `lipidomics`, `metabolomics` or `transcriptomics`.
#' @return A plain list with the groups `shell`, `binding`, `import`, `design` and `qc`.
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
        design = "named_list", qc = "named_list")
    checkOmicGuiFields(config, groups)
    validateOmicGuiShell(config$shell, schema)
    validateOmicGuiBinding(config$binding, config$shell$omic_type, schema)
    validateOmicGuiImport(config$import, schema)
    validateOmicGuiDesign(config$design, config, schema)
    validateOmicGuiQc(config$qc, config$shell$routes, schema)
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
        # Required at KIT-104 W8, once their functions are in the core (absent until then).
        qc_optional = c(filtering_figure = "ref", finalize_alignment = "ref"),
        # KIT-104 W6 (F104-04): the finalise's shell read, `function(binding, manager)`, for an omic
        # whose class cannot be rebuilt as the generic shell (trx); optional.
        qc_read = c(finalize_shell = "ref"),
        qc_group = c(id = "chr1_or_null", title = "chr1_or_null", component = "chr1_or_null",
            tabs = "list"),
        qc_tab = c(id = "chr1", title = "chr1", component = "chr1", routes = "chr",
            needs_paths = "lgl1"),
        # Required with the shared finalise sub-tab.
        qc_finalize = c(finalize_zero_rows = "lgl1", finalize_opens_norm = "lgl1",
            finalize_per_assay = "lgl1"),
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
            omic_qc_itsd = "itsd")
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
    checkOmicGuiFields(qc, schema$qc_optional, "qc", required = FALSE)
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
#' value); `lgl1`; `ref` a function; `ref_or_null`; `refs` a named list of functions; `list` a
#' non-empty list; `named_list`; `named_list_or_null`.
#' @noRd
omicGuiFieldHasType <- function(value, type) {
    is_chr <- is.character(value) && !anyNA(value)
    named <- length(value) == 0L ||
        (!is.null(names(value)) && all(nzchar(names(value))) && !anyDuplicated(names(value)))
    is_chr1 <- is_chr && length(value) == 1L && nzchar(value)
    isTRUE(switch(type,
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
        qc = omicGuiProteomicsQc())
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
        qc = omicGuiLipidQc())
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
        qc = omicGuiMetabQc())
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
        qc = omicGuiTrxQc())
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
