# KIT-102 W1: the omic GUI configuration contract (kit-101 design 4.2 to 4.4 and 4.11), its
# validator and the DIA binding adaptor. The field lists below are written from the ticket and the
# design, independently of the implementation's own schema, so that a field the validator forgets
# is caught here.

guiConfigOmics <- c("proteomics", "lipidomics", "metabolomics", "transcriptomics")

guiConfigRequired <- list(
    shell = c("omic_type", "prefix", "testid", "label", "label_lower", "header_icon", "tabset_id",
        "tab_style", "routes", "tabs", "bus", "bus_registration"),
    binding = c("binding", "runner", "owner", "read_state", "read_metadata"),
    import = c("source_panel", "title", "formats", "default_format", "detect", "help",
        "success_text", "preview", "process"),
    # KIT-103 W1 (design 4.5 and the ticket's rendering fields). The function fields still in
    # R/mod_*.R at W1 (validate, build_object, trx group_levels, annotation_step, the import reader)
    # are optional until W7.
    design = c("module_id", "help", "formula_default", "run_order", "run_order_rule",
        "record_provenance", "state_name", "state_description", "data_adapter", "design_import",
        "commit", "builder_ui", "record_excluded_runs", "label_lower", "multi_assay_banner",
        "no_data_text", "preview_text", "preview_assays", "testids", "builder_texts")
)

guiConfigBuilderTexts <- c("save_success", "no_assignment", "contrast_group_prefix",
    "contrast_as_is", "replicate_suffix")

# The feature-table panel's fields (lipid, metab, trx).
guiConfigFeatureTable <- c("format_heading", "assay_name", "id_label", "id_custom_label",
    "id_custom_placeholder", "annotation_label", "annotation_custom_placeholder", "second_assay",
    "id_mapping_key", "aux_inputs", "validate_mapping", "inputs", "outputs", "testids",
    "status_modes", "format_labels", "validation_noun", "assay_names", "notification_id",
    "refusal_prefix", "summary_file", "count_key")

guiConfigTabFields <- c("status_key", "value", "module_id", "component", "title", "icon", "testid",
    "step_name", "step_icon")

expectGuiConfigRefused <- function(config, field) {
    expect_error(validateOmicGuiConfig(config), regexp = field, fixed = TRUE,
        class = "omic_gui_config_invalid")
}

# Every function reachable through the configuration's lists.
guiConfigFunctions <- function(x) {
    if (is.function(x)) return(list(x))
    if (!is.list(x)) return(list())
    unlist(lapply(x, guiConfigFunctions), recursive = FALSE)
}

test_that("each omic's configuration validates, with the shell, binding, import and design groups", {
    for (omic in guiConfigOmics) {
        config <- omicGuiConfig(omic)
        expect_identical(names(config), c("shell", "binding", "import", "design"), info = omic)
        expect_identical(config$shell$omic_type, omic, info = omic)
        expect_identical(config$binding$binding$omic_type, omic, info = omic)
        expect_identical(validateOmicGuiConfig(config), config, info = omic)
        for (group in names(guiConfigRequired)) {
            expect_true(all(guiConfigRequired[[group]] %in% names(config[[group]])), info = omic)
        }
    }
    expect_error(omicGuiConfig("genomics"), class = "omic_gui_config_invalid")
})

test_that("the shell keys, runners and panels are today's per omic", {
    prot <- omicGuiConfig("proteomics")
    lipid <- omicGuiConfig("lipidomics")
    metab <- omicGuiConfig("metabolomics")
    expect_identical(prot$binding$runner, "dia_managed_step")
    expect_identical(prot$import$source_panel, "proteomics")
    expect_null(prot$shell$label)
    expect_identical(prot$shell$tabset_id, "workflow_tabs")
    expect_identical(vapply(prot$shell$tabs, `[[`, "", "status_key"), c("setup_import",
        "design_matrix", "quality_control", "normalization", "differential_expression",
        "enrichment_analysis", "session_summary"))
    expect_identical(prot$shell$routes, proteomicsDiaBinding$routes)
    for (config in list(lipid, metab, omicGuiConfig("transcriptomics"))) {
        expect_identical(config$binding$runner, "omics_actions")
        expect_identical(config$import$source_panel, "feature_table")
        expect_identical(vapply(config$shell$tabs, `[[`, "", "status_key"), c("setup_import",
            "design_matrix", "quality_control", "normalization", "differential_analysis",
            "session_summary"))
    }
    expect_identical(lipid$import$id_mapping_key, "lipid_id_col")
    expect_identical(metab$import$id_mapping_key, "metabolite_id_col")
    expect_identical(metab$import$count_key, "n_metabolites")
    expect_identical(vapply(metab$shell$tabs, `[[`, "", "testid")[[1L]], "metab-tab-import")
})

test_that("removing a required field, or giving it the wrong type, is refused by class", {
    for (omic in c("proteomics", "lipidomics")) {
        config <- omicGuiConfig(omic)
        for (group in names(guiConfigRequired)) {
            for (field in guiConfigRequired[[group]]) {
                path <- paste0(group, "$", field)
                removed <- config
                removed[[group]][[field]] <- NULL
                expectGuiConfigRefused(removed, path)
                mistyped <- config
                mistyped[[group]][field] <- list(42L)
                expectGuiConfigRefused(mistyped, path)
            }
        }
    }
    lipid <- omicGuiConfig("lipidomics")
    for (field in guiConfigFeatureTable) {
        removed <- lipid
        removed$import[[field]] <- NULL
        expectGuiConfigRefused(removed, paste0("import$", field))
        mistyped <- lipid
        mistyped$import[field] <- list(42L)
        expectGuiConfigRefused(mistyped, paste0("import$", field))
    }
    for (field in guiConfigTabFields) {
        removed <- lipid
        removed$shell$tabs[[2L]][[field]] <- NULL
        expectGuiConfigRefused(removed, paste0("shell$tabs[[2]]$", field))
    }
})

test_that("the runner, the references and the tab and status keys are checked", {
    lipid <- omicGuiConfig("lipidomics")
    bad <- lipid
    bad$binding$runner <- "worker"
    expectGuiConfigRefused(bad, "binding$runner")
    bad <- lipid
    bad$import$detect <- "detectLipidomicsFormat"
    expectGuiConfigRefused(bad, "import$detect")
    bad <- lipid
    bad$binding$owner <- list(createWorkflowSessionOwner)
    expectGuiConfigRefused(bad, "binding$owner")
    bad <- lipid
    bad$shell$tabs[[3L]]$status_key <- "design_matrix"
    expectGuiConfigRefused(bad, "status_key")
    bad <- lipid
    bad$shell$tabs[[4L]]$value <- "import"
    expectGuiConfigRefused(bad, "value")
    bad <- lipid
    bad$binding$binding <- metabolomicsDeclarationBinding
    expectGuiConfigRefused(bad, "binding$binding")
    prot <- omicGuiConfig("proteomics")
    for (field in c("omic_type", "routes", "route_for_format", "owner", "read_state",
        "read_metadata", "import_metadata")) {
        removed <- prot
        removed$binding$binding[[field]] <- NULL
        expectGuiConfigRefused(removed, paste0("binding$binding$", field))
    }
})

test_that("no configuration function is defined in an R/mod_*.R file (srcref, under load_all)", {
    if (is.null(utils::getSrcref(omicGuiConfig))) {
        testthat::skip("source references are kept under load_all only")
    }
    expect_identical(utils::getSrcFilename(omicGuiConfig), "utils_omic_gui_config.R")
    for (omic in guiConfigOmics) {
        functions <- guiConfigFunctions(omicGuiConfig(omic))
        files <- unlist(lapply(functions, utils::getSrcFilename))
        expect_gt(length(files), 0L)
        expect_false(any(startsWith(files, "mod_")),
            info = paste(omic, paste(unique(files[startsWith(files, "mod_")]), collapse = ", ")))
        # KIT-103 W1: the design group, and the DIA adaptor's design entries, by name.
        config <- omicGuiConfig(omic)
        design <- c(guiConfigFunctions(config$design), guiConfigFunctions(config$binding$binding$design))
        expect_gt(length(design), 0L)
        design_files <- unlist(lapply(design, utils::getSrcFilename))
        expect_false(any(startsWith(design_files, "mod_")), info = paste(omic, "design"))
    }
})

test_that("the trx shell and import values equal trxOmicProfile()'s fields today (design 4.12)", {
    profile <- trxOmicProfile()
    config <- omicGuiConfig("transcriptomics")
    for (field in c("omic_type", "prefix", "testid", "label", "label_lower")) {
        expect_identical(config$shell[[field]], profile[[field]], info = field)
    }
    expect_identical(config$shell$routes, profile$workflow_type)
    # KIT-102 W5 removed the profile's import-only fields; their values (mod_transcriptomics.R:34-39
    # at c1041c56) are pinned here as literals, unchanged.
    expect_identical(config$import$title, "Transcriptomics Data Import")
    expect_identical(config$import$formats,
        c("Count table (gene ID, gene name, sample counts)" = "trx_counts"))
    expect_identical(config$import$default_format, "trx_counts")
    expect_identical(config$import$help, paste("One delimited count table: a gene ID column, a",
        "gene name column and one column per sample, named exactly as in the design. Counts must",
        "be whole numbers: a missing, negative or fractional count is refused, with no rounding.",
        "Sample names are kept as they are in the file."))
    expected <- c(assay_name = "counts", id_label = "Gene ID Column",
        annotation_label = "Gene Name Column")
    for (field in names(expected)) {
        expect_identical(config$import[[field]], expected[[field]], info = field)
    }
    expect_identical(config$binding$binding, profile$binding)
    expect_identical(config$binding$read_state, profile$read_state)
    expect_identical(config$binding$binding$class, profile$class)
})

test_that("the DIA binding adaptor answers with today's core functions", {
    expect_identical(proteomicsDiaBinding$omic_type, "proteomics")
    expect_identical(proteomicsDiaBinding$routes, c("DIA", "DIA_LIMPA", "TMT", "LFQ"))
    route <- proteomicsDiaBinding$route_for_format
    for (format in c("diann", "pd_tmt", "maxquant", "fragpipe", "spectronaut", "unknown")) {
        expect_identical(route(format), workflowTypeForImportFormat(format), info = format)
    }
    expect_identical(route("diann", dia_route = "DIA_LIMPA"), "DIA_LIMPA")
    expect_identical(route("spectronaut", data_type = "peptide"), NA_character_)
    testthat::local_mocked_bindings(
        recordDIAImportMetadata = function(manager, table, result, reference) {
            list(table, reference)
        },
        getDIAImportMetadata = function(manager) "metadata",
        readDIAStateMetadata = function(manager, reference = NULL) c("state", reference)
    )
    metadata <- proteomicsDiaBinding$import_metadata
    expect_identical(metadata$record(NULL, "table", NULL, "reference"), list("table", "reference"))
    expect_identical(metadata$read(NULL), "metadata")
    expect_identical(proteomicsDiaBinding$read_state(NULL, "r"), c("state", "r"))
    expect_identical(proteomicsDiaBinding$read_metadata(NULL), "state")
})

# ------------------------------------------------------------------------------------------------
# KIT-103 W1: the design group (kit-101 design 4.5) and the DIA adaptor's design entries.
# ------------------------------------------------------------------------------------------------

test_that("the design group holds today's values per omic", {
    configs <- lapply(stats::setNames(guiConfigOmics, guiConfigOmics), omicGuiConfig)
    design <- lapply(configs, `[[`, "design")
    # The tab's module id, as the shell's (the harness DOM ids <omic>_workflow-design-builder-* and
    # proteomics_workflow-design_matrix-builder-*).
    for (omic in guiConfigOmics) {
        tab <- Filter(function(tab) identical(tab$status_key, "design_matrix"), configs[[omic]]$shell$tabs)
        expect_identical(design[[omic]]$module_id, tab[[1L]]$module_id, info = omic)
        expect_identical(design[[omic]]$label_lower, configs[[omic]]$shell$label_lower, info = omic)
        expect_named(design[[omic]]$builder_texts, guiConfigBuilderTexts, ignore.order = TRUE)
    }
    expect_identical(vapply(design, `[[`, "", "module_id"), c(proteomics = "design_matrix",
        lipidomics = "design", metabolomics = "design", transcriptomics = "design"))
    # Run order: mixedsort for a fresh lipid, metab and prot design, input order for trx; no default.
    for (omic in c("proteomics", "lipidomics", "metabolomics")) {
        expect_identical(design[[omic]]$run_order, gtools::mixedsort, info = omic)
        expect_identical(design[[omic]]$run_order_rule, "mixedsort", info = omic)
        expect_identical(design[[omic]]$run_order(c("S10", "S2", "S1")), c("S1", "S2", "S10"))
    }
    expect_identical(design$transcriptomics$run_order, identity)
    expect_identical(design$transcriptomics$run_order_rule, "input")
    expect_identical(design$transcriptomics$run_order(c("S10", "S2", "S1")), c("S10", "S2", "S1"))
    # Adapters, commits and builder UIs.
    expect_identical(vapply(design, `[[`, "", "data_adapter"), c(proteomics = "dia_long",
        lipidomics = "wide_assay_list", metabolomics = "wide_assay_list",
        transcriptomics = "wide_assay_list"))
    expect_identical(vapply(design, `[[`, "", "commit"), c(proteomics = "prot_hook",
        lipidomics = "shared", metabolomics = "shared", transcriptomics = "shared"))
    expect_identical(vapply(design, `[[`, "", "builder_ui"), c(proteomics = "proteomics",
        lipidomics = "shared", metabolomics = "shared", transcriptomics = "shared"))
    expect_identical(vapply(design, `[[`, TRUE, "record_excluded_runs"), c(proteomics = FALSE,
        lipidomics = FALSE, metabolomics = FALSE, transcriptomics = TRUE))
    # States, descriptions and provenance.
    expect_identical(design$lipidomics$state_name, "lipid_raw_data_s4")
    expect_identical(design$lipidomics$state_description,
        "Initial LipidomicsAssayData S4 object created after design matrix")
    expect_identical(design$metabolomics$state_name, "metab_raw_data_s4")
    expect_identical(design$metabolomics$state_description,
        "Initial MetaboliteAssayData S4 object created after design matrix")
    expect_identical(design$proteomics$state_name, c(DIA = "raw_data_s4", DIA_LIMPA = "raw_data_s4",
        TMT = "protein_s4_initial", LFQ = "protein_s4_initial"))
    expect_identical(design$proteomics$state_description, c(
        DIA = "Initial Peptide S4 object created after design matrix.",
        DIA_LIMPA = "Initial Peptide S4 object created after design matrix.",
        TMT = "Initial Protein S4 object created from TMT data after design matrix.",
        LFQ = "Initial Protein S4 object created from LFQ data after design matrix."))
    expect_identical(design$lipidomics$record_provenance, recordLipidDesignProvenance)
    expect_identical(design$metabolomics$record_provenance, recordMetabDesignProvenance)
    expect_identical(design$transcriptomics$record_provenance, recordTrxDesignProvenance)
    expect_identical(design$proteomics$record_provenance, setDIAStepProvenance)
    # Proteomics has no save-time preflight (D7) and no group-level rule; nor have lipid and metab.
    expect_true("validate" %in% names(design$proteomics))
    expect_null(design$proteomics$validate)
    for (omic in c("proteomics", "lipidomics", "metabolomics")) {
        expect_true("group_levels" %in% names(design[[omic]]), info = omic)
        expect_null(design[[omic]]$group_levels, info = omic)
    }
    # The feature omics' column-mapping key is the import's.
    for (omic in c("lipidomics", "metabolomics", "transcriptomics")) {
        expect_identical(design[[omic]]$id_mapping_key, configs[[omic]]$import$id_mapping_key)
    }
    # Frame and builder texts: trx renders the lipid frame and builder under its profile.
    expect_identical(vapply(design, `[[`, TRUE, "multi_assay_banner"), c(proteomics = FALSE,
        lipidomics = TRUE, metabolomics = TRUE, transcriptomics = FALSE))
    expect_identical(design$metabolomics$matrix_heading,
        "Current Design Matrix (All Samples) - Metabolomics")
    expect_identical(design$transcriptomics$matrix_heading,
        "Current Design Matrix (All Samples) - Lipidomics")
    expect_identical(design$transcriptomics$multi_assay_note,
        "<b>Multi-Assay Note:</b> This lipidomics workflow supports multiple assays ")
    expect_identical(design$proteomics$no_data_text, paste(" Please complete the 'Setup & Import'",
        "step first. The builder will appear here once data is available."))
    expect_identical(design$lipidomics$no_data_text,
        " Please complete the 'Import' step first. The builder will appear once data is available.")
    expect_identical(vapply(design, `[[`, TRUE, "preview_assays"), c(proteomics = FALSE,
        lipidomics = TRUE, metabolomics = TRUE, transcriptomics = TRUE))
    expect_identical(design$lipidomics$builder_texts[["save_success"]], "Design saved successfully!")
    expect_identical(design$proteomics$builder_texts[["save_success"]],
        "Design saved successfully. You can close this builder.")
    expect_identical(design$proteomics$builder_texts[["no_assignment"]],
        "No samples have been assigned to groups. Please assign metadata before saving.")
    expect_identical(design$proteomics$builder_texts[["contrast_as_is"]],
        "Note: Contrasts will use group names as-is (e.g., GA_Control-GA_Elevated)")
    expect_identical(design$metabolomics$builder_texts[["contrast_group_prefix"]],
        "Note: Contrasts will use 'group' prefix based on formula: ~ 0 + group")
    expect_identical(design$proteomics$builder_texts[["replicate_suffix"]], "selected runs:")
    expect_identical(design$lipidomics$builder_texts[["replicate_suffix"]], "selected samples:")
    expect_identical(design$lipidomics$testids, c(import_dir = "lipid-design-import-dir",
        import_confirm = "lipid-design-import-confirm"))
    expect_identical(design$proteomics$testids, c(import_existing = "prot-design-import-existing",
        replicate_start = "prot-design-replicate-start"))
    expect_length(design$transcriptomics$testids, 0L)
    # The import-existing-design plugin: none for trx; proteomics' own; the feature table's reader
    # values for lipid and metab.
    expect_null(design$transcriptomics$design_import)
    expect_true("design_import" %in% names(design$transcriptomics))
    expect_identical(design$proteomics$design_import$plugin, "proteomics")
    expect_identical(design$proteomics$design_import$import_dir_path, "empty")
    expect_identical(design$lipidomics$design_import$import_dir_path, "req")
    expect_identical(design$lipidomics$design_import$qc_baseline, updateLipidFiltering)
    expect_identical(design$metabolomics$design_import$qc_baseline, updateMetaboliteFiltering)
    expect_true(design$lipidomics$design_import$preflight)
    expect_false(design$metabolomics$design_import$preflight)
    expect_identical(design$metabolomics$design_import$id_candidates$annotation,
        c("Metabolite name", "Metabolite.name", "Name", "name"))
    expect_identical(design$lipidomics$design_import$state_description,
        "LipidomicsAssayData S4 object created from imported design")
})

test_that("run_order is required with no default; each required design field is refused by class", {
    for (omic in guiConfigOmics) {
        config <- omicGuiConfig(omic)
        removed <- config
        removed$design$run_order <- NULL
        expectGuiConfigRefused(removed, "design$run_order")
        mistyped <- config
        mistyped$design$run_order <- "mixedsort"
        expectGuiConfigRefused(mistyped, "design$run_order")
    }
    lipid <- omicGuiConfig("lipidomics")
    prot <- omicGuiConfig("proteomics")
    # The choices.
    for (field in c("run_order_rule", "data_adapter", "commit", "builder_ui")) {
        bad <- lipid
        bad$design[[field]] <- "other"
        expectGuiConfigRefused(bad, paste0("design$", field))
    }
    # The function fields still optional at W1 are checked when present.
    for (field in c("validate", "build_object", "group_levels", "annotation_step")) {
        bad <- lipid
        bad$design[[field]] <- "a function name"
        expectGuiConfigRefused(bad, paste0("design$", field))
    }
    # The builder's texts.
    for (field in guiConfigBuilderTexts) {
        bad <- lipid
        bad$design$builder_texts <- bad$design$builder_texts[names(bad$design$builder_texts) != field]
        expectGuiConfigRefused(bad, paste0("design$builder_texts$", field))
    }
    # Required by the wide-table adapter and the shared builder UI.
    for (field in c("id_mapping_key", "multi_assay_note", "matrix_heading")) {
        removed <- lipid
        removed$design[[field]] <- NULL
        expectGuiConfigRefused(removed, paste0("design$", field))
    }
    # A per-route state name names every route.
    bad <- prot
    bad$design$state_name <- bad$design$state_name[c("DIA", "TMT")]
    expectGuiConfigRefused(bad, "design$state_name")
    # The import-existing-design plugin.
    bad <- lipid
    bad$design$design_import$plugin <- "folder"
    expectGuiConfigRefused(bad, "design$design_import$plugin")
    bad <- lipid
    bad$design$design_import$import_dir_path <- "blank"
    expectGuiConfigRefused(bad, "design$design_import$import_dir_path")
    for (field in c("id_candidates", "qc_baseline", "preflight", "state_description")) {
        removed <- lipid
        removed$design$design_import[[field]] <- NULL
        expectGuiConfigRefused(removed, paste0("design$design_import$", field))
    }
    bad <- lipid
    bad$design$design_import$reader <- "a function name"
    expectGuiConfigRefused(bad, "design$design_import$reader")
    removed <- lipid
    removed$design$design_import$id_candidates$annotation <- NULL
    expectGuiConfigRefused(removed, "design$design_import$id_candidates$annotation")
    # The DIA adapter needs the adaptor's design entries, and the adaptor's commit hook.
    for (field in c("initial_state", "rename", "subset", "commit")) {
        removed <- prot
        removed$binding$binding$design[[field]] <- NULL
        expectGuiConfigRefused(removed, paste0("binding$binding$design$", field))
    }
    bad <- prot
    bad$design$commit <- "shared"
    expectGuiConfigRefused(bad, "design$commit")
})

test_that("the trx design values equal trxOmicProfile()'s eight design fields today (design 4.12)", {
    # KIT-103 W5 removed the eight fields from the profile; their values at 283f5eee are frozen here
    # (mod_transcriptomics.R:36-41 then; group_levels, validate and build_object are the
    # func_trx_design.R copies, W2).
    design <- omicGuiConfig("transcriptomics")$design
    frozen <- list(help = paste("Samples are matched to the count columns by sample name only, in the counts file's",
            "order; there is no positional mapping. Each group needs at least two samples. The order in which you add",
            "the factor levels is the group level order the analysis uses."),
        run_order = identity, group_levels = trxDesignGroupLevels, validate = trxDesignCheck,
        build_object = trxDesignObject, record_provenance = recordTrxDesignProvenance, state_name = "trx_raw_data_s4",
        state_description = "Initial TranscriptomicsAssayData S4 object created after design matrix")
    # group_levels, validate and build_object join the configuration when W7 wires them, and then they
    # are compared too.
    deferred <- c("group_levels", "validate", "build_object")
    for (field in names(frozen)) {
        if (field %in% deferred && !field %in% names(design)) next
        expect_identical(design[[field]], frozen[[field]], info = field)
    }
    expect_false(any(c("design_help", "run_order", "group_levels", "validate_design", "build_object",
        "record_provenance", "state_name", "state_description") %in% names(trxOmicProfile())))
    expect_true(design$record_excluded_runs)
    expect_identical(design$label_lower, "transcriptomics")
    expect_identical(trxOmicProfile()$omic_type, "transcriptomics")
})

test_that("the DIA binding adaptor's design entries call today's core functions", {
    entries <- proteomicsDiaBinding$design
    expect_identical(entries$commit, "prot_hook")
    testthat::local_mocked_bindings(
        transformDIADesignTable = function(manager, reference, samples = NULL, rename = NULL,
            numeric_columns = character(), strip_parser = FALSE) {
            list(manager = manager, reference = reference, samples = samples, rename = rename,
                numeric_columns = numeric_columns, strip_parser = strip_parser)
        },
        buildDIAReferenceDesignState = function(manager, reference, config_list, column_mapping) {
            list("reference", manager, reference, config_list, column_mapping)
        },
        buildProtDesignInitialState = function(dataTbl, configList, columnMapping) {
            list("table", dataTbl, configList, columnMapping)
        }
    )
    reference <- list(schema = 1L, entity_id = strrep("a", 64L), payload = list())
    table <- data.frame(Run = c("S2", "S1"))
    expect_identical(entries$initial_state(NULL, reference, "config", "mapping"),
        list("reference", NULL, reference, "config", "mapping"))
    expect_identical(entries$initial_state(NULL, table, "config", "mapping"),
        list("table", table, "config", "mapping"))
    renamed <- entries$rename("manager", reference, c(S1 = "A1"))
    expect_identical(renamed[c("manager", "reference", "samples", "rename", "strip_parser")],
        list(manager = "manager", reference = reference, samples = NULL, rename = c(S1 = "A1"),
            strip_parser = FALSE))
    subset <- entries$subset("manager", reference, c("S1", "S2"), strip_parser = TRUE)
    expect_identical(subset[c("samples", "rename", "strip_parser")],
        list(samples = c("S1", "S2"), rename = NULL, strip_parser = TRUE))
    expect_identical(entries$subset("manager", reference, "S1")$strip_parser, FALSE)
})
