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

# ------------------------------------------------------------------------------------------------
# The entry point.
# ------------------------------------------------------------------------------------------------

#' The shared GUI layer's configuration of one omic, validated.
#' @param omic_type One of `proteomics`, `lipidomics`, `metabolomics` or `transcriptomics`.
#' @return A plain list with the groups `shell`, `binding` and `import`.
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
    groups <- c(shell = "named_list", binding = "named_list", import = "named_list")
    checkOmicGuiFields(config, groups)
    validateOmicGuiShell(config$shell, schema)
    validateOmicGuiBinding(config$binding, config$shell$omic_type, schema)
    validateOmicGuiImport(config$import, schema)
    invisible(config)
}

#' The contract's fields and their types (design 4.2 to 4.4), by group. The type codes are
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
        assay_names = c(rule = "chr1", path_pattern = "chr1", messages = "named_chr")
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
#' non-empty character vector; `chr0` a character vector; `named_chr`; `lgl1`; `ref` a function;
#' `refs` a named list of functions; `list` a non-empty list; `named_list`.
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
        lgl1 = is.logical(value) && length(value) == 1L && !is.na(value),
        ref = is.function(value),
        refs = is.list(value) && length(value) >= 1L && named &&
            all(vapply(value, is.function, logical(1))),
        list = is.list(value) && length(value) >= 1L,
        named_list = is.list(value) && named))
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
        import = omicGuiProteomicsImport())
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
        import = omicGuiLipidImport())
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
        import = omicGuiMetabImport())
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
        import = omicGuiTrxImport())
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
