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
    # KIT-103 W1 (design 4.5 and the ticket's rendering fields). KIT-103 W7 wired the core's
    # function fields; every omic carries validate, build_object and group_levels (NULL where the
    # omic has none). Only proteomics has an annotation step, so it stays optional.
    design = c("module_id", "help", "formula_default", "run_order", "run_order_rule",
        "record_provenance", "state_name", "state_description", "data_adapter", "design_import",
        "commit", "builder_ui", "record_excluded_runs", "label_lower", "multi_assay_banner",
        "no_data_text", "preview_text", "preview_assays", "testids", "builder_texts",
        "group_levels", "validate", "build_object"),
    # KIT-104 W1 (design 4.6 and the ticket's rendering fields). `complete_state` and `itsd` are
    # NULL where the omic has none; the finalise fields, required with a finalise tab, are checked by
    # their own blocks below (KIT-104 W8 made its function fields required there).
    qc = c("title", "tabset_id", "groups", "auto_init", "complete_state", "itsd", "texts",
        "notification_ids"),
    # KIT-106 W1 (D8, and the fields every omic carries); the conditional and optional ones are checked
    # by the da blocks below.
    da = c("module_id", "status_key", "testids", "model", "threshold", "lfc", "columns", "summary", "tabs",
        "assays", "table_filter", "views", "texts", "download")
)

# KIT-107 W1 (D1, D9): the enrichment group's required fields (proteomics only; NULL for the other omics). The
# function fields still in R/mod_*.R at the start commit (`organism`, `annotation`, `kernel`) are optional until the
# W2 merge wires them; the enrichment blocks below check them.
guiConfigEnrichmentRequired <- c("module_id", "status_key", "tab_value", "operation", "run", "da_source", "id_column",
    "cutoffs", "correction", "exclude_iea", "taxon", "services", "organism_filter", "state_name")

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

test_that("each omic's configuration validates, with the shell, binding, import, design, qc, norm, da and enrichment groups", {
    for (omic in guiConfigOmics) {
        config <- omicGuiConfig(omic)
        expect_identical(names(config), c("shell", "binding", "import", "design", "qc", "norm", "da", "enrichment"),
            info = omic)
        expect_identical(config$shell$omic_type, omic, info = omic)
        expect_identical(config$binding$binding$omic_type, omic, info = omic)
        expect_identical(validateOmicGuiConfig(config), config, info = omic)
        for (group in names(guiConfigRequired)) {
            expect_true(all(guiConfigRequired[[group]] %in% names(config[[group]])), info = omic)
        }
        # KIT-107 W1 (D1, D9): proteomics alone has the enrichment group; the other omics carry NULL.
        if (identical(omic, "proteomics")) {
            expect_true(all(guiConfigEnrichmentRequired %in% names(config$enrichment)), info = omic)
        } else {
            expect_true("enrichment" %in% names(config), info = omic)
            expect_null(config$enrichment, info = omic)
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
        # KIT-104 W1: the qc group. Since W8 every feature omic carries its tracking figure and its
        # finalise alignment (lipid and metab also their ITSD kernels).
        qc_files <- as.character(unlist(lapply(guiConfigFunctions(config$qc), utils::getSrcFilename)))
        if (omic != "proteomics") expect_gt(length(qc_files), 0L)
        expect_false(any(startsWith(qc_files, "mod_")), info = paste(omic, "qc"))
        # KIT-105 W1: the norm group (lipid and metab: the managed run and correlation and the ITSD
        # selection table; the in-process twins and the export's collector and saver from W6).
        norm_files <- as.character(unlist(lapply(guiConfigFunctions(config$norm), utils::getSrcFilename)))
        # KIT-110 W3: trx's too (its run and reset handlers, utils_workflow_trx_norm.R).
        if (omic %in% c("lipidomics", "metabolomics", "transcriptomics")) expect_gt(length(norm_files), 0L)
        expect_false(any(startsWith(norm_files, "mod_")), info = paste(omic, "norm"))
        # KIT-106 W1 (D8): the da group (the managed and in-process DA, the plots and the readers; a function still in an
        # R/mod_*.R file is wired only when it moves to the core).
        da_files <- as.character(unlist(lapply(guiConfigFunctions(config$da), utils::getSrcFilename)))
        expect_gt(length(da_files), 0L)
        expect_false(any(startsWith(da_files, "mod_")), info = paste(omic, "da"))
        # KIT-107 W1 (D9): the enrichment group (proteomics: the declaration, its run and the DA source's readers; a
        # function still in an R/mod_*.R file is wired only when W2 moves it to the core).
        enrichment_files <- as.character(unlist(lapply(guiConfigFunctions(config$enrichment), utils::getSrcFilename)))
        if (identical(omic, "proteomics")) expect_gt(length(enrichment_files), 0L)
        expect_false(any(startsWith(enrichment_files, "mod_")), info = paste(omic, "enrichment"))
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
    # KIT-104 W8 removed the profile's `read_state` and `class` (mod_transcriptomics.R:39-40 at
    # 691aaa71); their values are pinned here as literals, unchanged.
    expect_identical(config$binding$read_state, readTrxWorkflowState)
    expect_identical(config$binding$binding$class, "TranscriptomicsAssayData")
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
    # KIT-103 W7: the function fields are the core's (W2, W3), wired.
    expect_identical(design$lipidomics$build_object, lipidDesignObject)
    expect_identical(design$metabolomics$build_object, metabDesignObject)
    expect_identical(design$transcriptomics$build_object, trxDesignObject)
    expect_identical(design$proteomics$build_object, buildProtDesignStateCheckpoint)
    expect_identical(design$transcriptomics$group_levels, trxDesignGroupLevels)
    expect_identical(design$transcriptomics$validate, trxDesignCheck)
    testthat::local_mocked_bindings(validateOmicsDesignSave = function(results, columnMapping, omicType) {
        list(results = results, columnMapping = columnMapping, omicType = omicType)
    })
    for (omic in c("lipidomics", "metabolomics")) {
        expect_true(is.function(design[[omic]]$validate), info = omic)
        if (is.function(design[[omic]]$validate)) {
            expect_identical(design[[omic]]$validate("results", "mapping"),
                list(results = "results", columnMapping = "mapping", omicType = omic), info = omic)
        }
        expect_identical(design[[omic]]$design_import$reader, readOmicsDesignImport, info = omic)
    }
    expect_true(isTRUE(design$lipidomics$design_import$record_design_file))
    expect_identical(design$metabolomics$design_import$record_design_file, FALSE)
    # Metab checks the import's files before the modal closes (P-11); lipid after.
    expect_identical(design$lipidomics$design_import$check_files_first, FALSE)
    expect_true(isTRUE(design$metabolomics$design_import$check_files_first))
    expect_identical(design$proteomics$annotation_step, list(uniprot = retrieveProtDesignUniprotAnnotations,
        fasta = processProtDesignImportedFasta))
    for (omic in c("lipidomics", "metabolomics", "transcriptomics")) {
        expect_false("annotation_step" %in% names(design[[omic]]), info = omic)
    }
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
    # The function fields are checked by type; the optional annotation step when present.
    for (field in c("validate", "build_object", "group_levels", "annotation_step")) {
        bad <- lipid
        bad$design[[field]] <- "a function name"
        expectGuiConfigRefused(bad, paste0("design$", field))
    }
    bad <- prot
    bad$design$annotation_step$fasta <- "processProtDesignImportedFasta"
    expectGuiConfigRefused(bad, "design$annotation_step")
    bad <- prot
    bad$design$annotation_step <- unname(bad$design$annotation_step)
    expectGuiConfigRefused(bad, "design$annotation_step")
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
    for (field in c("id_candidates", "qc_baseline", "preflight", "state_description", "reader",
        "record_design_file", "check_files_first")) {
        removed <- lipid
        removed$design$design_import[[field]] <- NULL
        expectGuiConfigRefused(removed, paste0("design$design_import$", field))
    }
    bad <- lipid
    bad$design$design_import$reader <- "a function name"
    expectGuiConfigRefused(bad, "design$design_import$reader")
    for (field in c("record_design_file", "check_files_first")) {
        bad <- lipid
        bad$design$design_import[[field]] <- NA
        expectGuiConfigRefused(bad, paste0("design$design_import$", field))
    }
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
    # func_trx_design.R functions, W2, wired by W7).
    design <- omicGuiConfig("transcriptomics")$design
    frozen <- list(help = paste("Samples are matched to the count columns by sample name only, in the counts file's",
            "order; there is no positional mapping. Each group needs at least two samples. The order in which you add",
            "the factor levels is the group level order the analysis uses."),
        run_order = identity, group_levels = trxDesignGroupLevels, validate = trxDesignCheck,
        build_object = trxDesignObject, record_provenance = recordTrxDesignProvenance, state_name = "trx_raw_data_s4",
        state_description = "Initial TranscriptomicsAssayData S4 object created after design matrix")
    for (field in names(frozen)) {
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

# ------------------------------------------------------------------------------------------------
# KIT-104 W1: the qc group (kit-101 design 4.6) and what the shared QC sub-tabs render per omic.
# ------------------------------------------------------------------------------------------------

guiConfigQcGroupFields <- c("id", "title", "component", "tabs")
guiConfigQcTabFields <- c("id", "title", "component", "routes", "needs_paths")
guiConfigQcFinalize <- c("finalize_zero_rows", "finalize_opens_norm", "finalize_per_assay")
# KIT-104 W8: the finalise's function fields, required with a finalise tab (every feature omic).
guiConfigQcFinalizeFunctions <- c("filtering_figure", "finalize_alignment")
guiConfigQcItsd <- c("name_columns", "metrics", "quant", "tip")
# The texts and the notification ids each shared sub-tab renders, by its component.
guiConfigQcTexts <- list(
    omic_qc_intensity = c("intensity_heading", "intensity_help", "intensity_proportion_help",
        "intensity_original_label", "intensity_working", "intensity_result_heading", "intensity_total",
        "intensity_success", "intensity_error"),
    omic_qc_duplicates = c("duplicates_help", "duplicates_none"),
    omic_qc_itsd = "itsd_pattern_help",
    omic_qc_finalize = c("finalize_total_label", "finalize_id_label", "finalize_total",
        "finalize_no_assays", "finalize_no_design"))
guiConfigQcNotifications <- c(omic_qc_intensity = "intensity", omic_qc_duplicates = "duplicates",
    omic_qc_itsd = "itsd")

# A group's tabs on a route.
guiConfigQcTabs <- function(group, route) Filter(function(tab) route %in% tab$routes, group$tabs)

# A group's tabs on a route in the shape of today's proteomics spec tables (getProtQc*ModuleSpecs).
guiConfigQcSpecs <- function(group, route) {
    lapply(guiConfigQcTabs(group, route), function(tab) {
        c(list(id = tab$id, title = tab$title, uiFn = paste0("mod_", tab$component, "_ui"),
            serverFn = paste0("mod_", tab$component, "_server")),
            if (isTRUE(tab$needs_paths)) list(needsExperimentPaths = TRUE))
    })
}

test_that("the qc group holds today's values per omic", {
    configs <- lapply(stats::setNames(guiConfigOmics, guiConfigOmics), omicGuiConfig)
    qc <- lapply(configs, `[[`, "qc")
    # mod_prot_qc.R:43; mod_lipid_qc.R:45; mod_metab_qc.R:44; mod_transcriptomics.R:38 (qc_title).
    expect_identical(vapply(qc, `[[`, "", "title"), c(proteomics = "Quality Control & Filtering",
        lipidomics = "Lipid Quality Control & Filtering",
        metabolomics = "Metabolite Quality Control & Filtering",
        transcriptomics = "Transcriptomics Quality Control & Filtering"))
    # mod_lipid_qc.R:120; mod_metab_qc.R:116; mod_trx_qc_filter.R:179; mod_prot_qc.R:121, :130, :138.
    expect_identical(qc$lipidomics$tabset_id, "lipid_qc_tabs")
    expect_identical(qc$metabolomics$tabset_id, "metab_qc_tabs")
    expect_identical(qc$transcriptomics$tabset_id, "trx_qc_tabs")
    expect_identical(qc$proteomics$tabset_id, c(DIA = "qc_tabs_lfq", DIA_LIMPA = "qc_tabs_lfq",
        TMT = "qc_tabs_tmt", LFQ = "qc_tabs_tmt"))
    expect_identical(vapply(qc, `[[`, TRUE, "auto_init"), c(proteomics = FALSE, lipidomics = TRUE,
        metabolomics = TRUE, transcriptomics = TRUE))
    # QC completes on protein_replicate_filtered (mod_prot_qc.R:204-216 at 827f2328), and on the
    # DIA_LIMPA route at its last QC state, peptide_count_filtered (P-16, fixed by KIT-104 W7); the
    # other omics' finalise tab completes QC.
    expect_identical(qc$proteomics$complete_state, c(DIA = "protein_replicate_filtered",
        DIA_LIMPA = "peptide_count_filtered", TMT = "protein_replicate_filtered",
        LFQ = "protein_replicate_filtered"))
    for (omic in c("lipidomics", "metabolomics", "transcriptomics")) {
        expect_true("complete_state" %in% names(qc[[omic]]), info = omic)
        expect_null(qc[[omic]]$complete_state, info = omic)
    }
    # One unnamed group of flat tabs for lipid, metab and trx; the sub-tab module ids, titles and
    # components; the omic's one route.
    feature_tabs <- list(
        lipidomics = c(intensity = "Intensity Filter", duplicates = "Duplicate Resolution",
            itsd = "Internal Standards", s4_finalize = "Finalize QC"),
        metabolomics = c(intensity = "Intensity Filter", duplicates = "Duplicate Resolution",
            itsd = "Internal Standards", s4_finalize = "Finalize QC"),
        transcriptomics = c(filter = "Expression Filter", s4_finalize = "Finalize QC"))
    components <- list(lipidomics = c("omic_qc_intensity", "omic_qc_duplicates", "omic_qc_itsd",
        "omic_qc_finalize"), transcriptomics = c("trx_qc_filter", "omic_qc_finalize"))
    components$metabolomics <- components$lipidomics
    for (omic in names(feature_tabs)) {
        groups <- qc[[omic]]$groups
        expect_length(groups, 1L)
        expect_null(groups[[1L]]$id)
        expect_null(groups[[1L]]$title)
        expect_null(groups[[1L]]$component)
        tabs <- groups[[1L]]$tabs
        expect_identical(vapply(tabs, `[[`, "", "id"), names(feature_tabs[[omic]]), info = omic)
        expect_identical(vapply(tabs, `[[`, "", "title"), unname(feature_tabs[[omic]]), info = omic)
        expect_identical(vapply(tabs, `[[`, "", "component"), components[[omic]], info = omic)
        for (tab in tabs) {
            expect_identical(tab$routes, configs[[omic]]$shell$routes, info = omic)
            expect_false(tab$needs_paths, info = omic)
        }
    }
    # Proteomics: the peptide and protein groups, rendered by their group plugins.
    groups <- qc$proteomics$groups
    expect_identical(vapply(groups, `[[`, "", "id"), c("peptide_qc", "protein_qc"))
    expect_identical(vapply(groups, `[[`, "", "title"), c("Peptide QC", "Protein QC"))
    expect_identical(vapply(groups, `[[`, "", "component"), c("prot_qc_peptide", "prot_qc_protein"))
    # The finalise: lipid refuses any zero-row assay and opens normalisation; metab neither; trx as
    # lipid, with one table (design 4.6; D2, D3).
    finalize <- function(omic) unlist(qc[[omic]][guiConfigQcFinalize])
    expect_identical(finalize("lipidomics"), c(finalize_zero_rows = TRUE, finalize_opens_norm = TRUE,
        finalize_per_assay = TRUE))
    expect_identical(finalize("metabolomics"), c(finalize_zero_rows = FALSE,
        finalize_opens_norm = FALSE, finalize_per_assay = TRUE))
    expect_identical(finalize("transcriptomics"), c(finalize_zero_rows = TRUE,
        finalize_opens_norm = TRUE, finalize_per_assay = FALSE))
    expect_false(any(c(guiConfigQcFinalize, "filtering_figure", "finalize_alignment") %in%
        names(qc$proteomics)))
    # The tracking figure: the core's (trx's arrives with W2's move of trxQcFilteringFigure).
    expect_identical(qc$lipidomics$filtering_figure, updateLipidFiltering)
    expect_identical(qc$metabolomics$filtering_figure, updateMetaboliteFiltering)
    # The ITSD analysis: lipid and metab only.
    expect_identical(qc$lipidomics$itsd$name_columns, c("Lipid name", "Name", "lipid"))
    expect_identical(qc$metabolomics$itsd$name_columns, c("Metabolite name", "Name", "metabolite"))
    expect_identical(qc$lipidomics$itsd$quant, getLipidQuantData)
    expect_identical(qc$metabolomics$itsd$quant, getMetaboliteQuantData)
    expect_identical(qc$lipidomics$itsd$tip,
        "Tip: Check if your data has a column with lipid names (e.g., 'Lipid name')")
    expect_identical(qc$metabolomics$itsd$tip,
        "Tip: Check if your data has a column with metabolite names (e.g., 'Metabolite name')")
    for (omic in c("proteomics", "transcriptomics")) {
        expect_true("itsd" %in% names(qc[[omic]]), info = omic)
        expect_null(qc[[omic]]$itsd, info = omic)
    }
    # The metrics entry calls the omic's kernel with its id-column argument.
    testthat::local_mocked_bindings(
        getLipidInternalStandardMetrics = function(assay_data, is_pattern, lipid_id_col, sample_id_col,
            sample_columns = NULL) list("lipid", assay_data, is_pattern, lipid_id_col, sample_id_col),
        getInternalStandardMetrics = function(assay_data, is_pattern, metabolite_id_col, sample_id_col,
            sample_columns = NULL) list("metab", assay_data, is_pattern, metabolite_id_col, sample_id_col))
    expect_identical(qc$lipidomics$itsd$metrics("assay", "^IS_", "Name", "Run"),
        list("lipid", "assay", "^IS_", "Name", "Run"))
    expect_identical(qc$metabolomics$itsd$metrics("assay", "^IS_", "Name", "Run"),
        list("metab", "assay", "^IS_", "Name", "Run"))
    # The texts and notification ids: every key its sub-tabs need, and today's values.
    for (omic in guiConfigOmics) {
        used <- unlist(lapply(qc[[omic]]$groups, function(group) vapply(group$tabs, `[[`, "", "component")))
        expect_true(all(unlist(guiConfigQcTexts[intersect(used, names(guiConfigQcTexts))]) %in%
            names(qc[[omic]]$texts)), info = omic)
        expect_setequal(names(qc[[omic]]$notification_ids),
            unname(guiConfigQcNotifications[intersect(used, names(guiConfigQcNotifications))]))
    }
    expect_length(qc$proteomics$texts, 0L)
    expect_identical(qc$lipidomics$notification_ids, c(intensity = "lipid_intensity_filter_working",
        duplicates = "lipid_dup_resolve_working", itsd = "is_analysis_working"))
    expect_identical(qc$metabolomics$notification_ids, c(intensity = "metab_intensity_filter_working",
        duplicates = "metab_dup_resolve_working", itsd = "is_analysis_working"))
    expect_identical(qc$lipidomics$texts[["intensity_total"]], "Total: %d -> %d lipids (removed %d)")
    expect_identical(qc$metabolomics$texts[["intensity_success"]],
        "Intensity filter applied: %d metabolites retained")
    expect_identical(qc$metabolomics$texts[["duplicates_none"]],
        "All metabolite IDs are unique. No resolution needed.")
    expect_identical(qc$lipidomics$texts[["finalize_no_assays"]],
        "Current LipidomicsAssayData object has no assay data to finalize.")
    expect_identical(qc$transcriptomics$texts[["finalize_total"]], "Total genes retained: %d")
    expect_identical(qc$transcriptomics$texts[["finalize_id_label"]], "Gene ID Column:")
    # Trx has one refusal for both checks (mod_trx_qc_filter.R:227).
    expect_identical(qc$transcriptomics$texts[["finalize_no_assays"]],
        "Current TranscriptomicsAssayData object has no genes or no design rows to finalize.")
    expect_identical(qc$transcriptomics$texts[["finalize_no_design"]],
        qc$transcriptomics$texts[["finalize_no_assays"]])
})

test_that("proteomics' QC groups per route are today's spec tables (getProtQc{Peptide,Protein}ModuleSpecs)", {
    groups <- omicGuiConfig("proteomics")$qc$groups
    names(groups) <- vapply(groups, `[[`, "", "id")
    for (route in c("DIA", "DIA_LIMPA")) {
        expect_identical(guiConfigQcSpecs(groups$peptide_qc, route), getProtQcPeptideModuleSpecs(route),
            info = route)
    }
    # TMT and LFQ render the protein group only (mod_prot_qc.R:98-133).
    for (route in c("TMT", "LFQ")) expect_length(guiConfigQcTabs(groups$peptide_qc, route), 0L)
    for (route in c("DIA", "DIA_LIMPA", "TMT", "LFQ")) {
        expect_identical(guiConfigQcSpecs(groups$protein_qc, route),
            getProtQcProteinModuleSpecs(workflowType = route), info = route)
    }
    # Each plugin's module exists under its stem; so do trx's filter and the group plugins.
    stems <- c(unlist(lapply(groups, function(group) vapply(group$tabs, `[[`, "", "component"))),
        "trx_qc_filter", vapply(groups, `[[`, "", "component"))
    for (stem in stems) {
        for (side in c("_ui", "_server")) {
            expect_true(exists(paste0("mod_", stem, side), mode = "function"), info = stem)
        }
    }
})

test_that("each qc field, required or conditional, is refused by class when missing or mistyped", {
    lipid <- omicGuiConfig("lipidomics")
    prot <- omicGuiConfig("proteomics")
    trx <- omicGuiConfig("transcriptomics")
    # The finalise fields are required with a finalise tab; proteomics has none.
    for (field in guiConfigQcFinalize) {
        removed <- lipid
        removed$qc[[field]] <- NULL
        expectGuiConfigRefused(removed, paste0("qc$", field))
        mistyped <- lipid
        mistyped$qc[[field]] <- NA
        expectGuiConfigRefused(mistyped, paste0("qc$", field))
    }
    # KIT-104 W8: so are its function fields, for every feature omic.
    for (omic in c("lipidomics", "metabolomics", "transcriptomics")) {
        config <- omicGuiConfig(omic)
        for (field in guiConfigQcFinalizeFunctions) {
            expect_true(is.function(config$qc[[field]]), info = paste(omic, field))
            removed <- config
            removed$qc[[field]] <- NULL
            expectGuiConfigRefused(removed, paste0("qc$", field))
        }
    }
    # The texts and notification ids each tab's component renders.
    for (component in names(guiConfigQcTexts)) {
        for (key in guiConfigQcTexts[[component]]) {
            removed <- lipid
            removed$qc$texts <- removed$qc$texts[names(removed$qc$texts) != key]
            expectGuiConfigRefused(removed, paste0("qc$texts$", key))
        }
    }
    for (key in guiConfigQcNotifications) {
        removed <- lipid
        removed$qc$notification_ids <- removed$qc$notification_ids[names(removed$qc$notification_ids) != key]
        expectGuiConfigRefused(removed, paste0("qc$notification_ids$", key))
    }
    removed <- trx
    removed$qc$texts <- removed$qc$texts[names(removed$qc$texts) != "finalize_total"]
    expectGuiConfigRefused(removed, "qc$texts$finalize_total")
    # The ITSD tab needs the analysis entries.
    bad <- lipid
    bad$qc["itsd"] <- list(NULL)
    expectGuiConfigRefused(bad, "qc$itsd")
    for (field in guiConfigQcItsd) {
        removed <- lipid
        removed$qc$itsd[[field]] <- NULL
        expectGuiConfigRefused(removed, paste0("qc$itsd$", field))
        mistyped <- lipid
        mistyped$qc$itsd[field] <- list(42L)
        expectGuiConfigRefused(mistyped, paste0("qc$itsd$", field))
    }
    # The groups and their tabs.
    for (field in guiConfigQcGroupFields) {
        removed <- prot
        removed$qc$groups[[2L]][[field]] <- NULL
        expectGuiConfigRefused(removed, paste0("qc$groups[[2]]$", field))
    }
    for (field in guiConfigQcTabFields) {
        removed <- lipid
        removed$qc$groups[[1L]]$tabs[[3L]][[field]] <- NULL
        expectGuiConfigRefused(removed, paste0("qc$groups[[1]]$tabs[[3]]$", field))
        mistyped <- lipid
        mistyped$qc$groups[[1L]]$tabs[[3L]][field] <- list(42L)
        expectGuiConfigRefused(mistyped, paste0("qc$groups[[1]]$tabs[[3]]$", field))
    }
    bad <- lipid
    bad$qc$groups[[1L]]$tabs[[2L]]$id <- "intensity"
    expectGuiConfigRefused(bad, "qc$groups[[1]]$tabs$id")
    bad <- lipid
    bad$qc$groups[[1L]]$tabs[[2L]]$routes <- "DIA"
    expectGuiConfigRefused(bad, "qc$groups[[1]]$tabs[[2]]$routes")
    bad <- prot
    bad$qc$groups[[2L]]$id <- "peptide_qc"
    expectGuiConfigRefused(bad, "qc$groups$id")
    # An unnamed group is the only group.
    bad <- prot
    bad$qc$groups[[1L]]["id"] <- list(NULL)
    expectGuiConfigRefused(bad, "qc$groups")
    # Per-route values name the omic's routes.
    bad <- prot
    bad$qc$tabset_id <- bad$qc$tabset_id[c("DIA", "TMT")]
    expectGuiConfigRefused(bad, "qc$tabset_id")
    bad <- prot
    bad$qc$complete_state <- c(bad$qc$complete_state, GENOMICS = "protein_replicate_filtered")
    expectGuiConfigRefused(bad, "qc$complete_state")
    bad <- prot
    bad$qc$auto_init <- NA
    expectGuiConfigRefused(bad, "qc$auto_init")
    # The function fields are checked by type.
    for (field in guiConfigQcFinalizeFunctions) {
        bad <- lipid
        bad$qc[[field]] <- "a function name"
        expectGuiConfigRefused(bad, paste0("qc$", field))
    }
})

# ------------------------------------------------------------------------------------------------
# KIT-105 W1: the norm group (kit-101 design 4.7) for lipid and metab; NULL for proteomics
# (KIT-111). KIT-110 W3: trx's group (its own blocks below).
# ------------------------------------------------------------------------------------------------

guiConfigNormRequired <- c("module_id", "panels", "plot_tabs", "stage_labels", "assay_slots", "ruv_heading",
    "ruv_help", "run", "correlation", "pre_qc_on_entry", "export", "itsd_table", "testids",
    "skip_correlation_flags", "final_pca", "managed_correlation_notify")
# The option sections and each one's required options (the shared panel library's order).
guiConfigNormPanels <- list(aesthetics = character(), itsd = "aggregation", log2 = "offset",
    between_sample = c("methods", "default"), ruv_iii = c("modes", "default", "percentage", "k", "noun"),
    actions = character(), export = "button")
# Wired to the core by KIT-105 W4 for lipid and metab; required since KIT-105 W6.
guiConfigNormCoreFns <- list(run = "in_process", correlation = "in_process", export = c("collect", "save"))
guiConfigNormWired <- c(run.in_process = "runWorkflowOmicsInProcessNormalisation",
    correlation.in_process = "runWorkflowOmicsInProcessCorrelation", export.collect = "collectWorkflowOmicsNormSession",
    export.save = "saveWorkflowOmicsNormSessionArtifacts")

test_that("the norm group holds today's values for lipid and metab, trx has one (KIT-110) and proteomics' is NULL", {
    configs <- lapply(stats::setNames(guiConfigOmics, guiConfigOmics), omicGuiConfig)
    expect_true("norm" %in% names(configs$proteomics))
    expect_null(configs$proteomics$norm)
    # KIT-110 W3: trx's group is no longer NULL, and it validates.
    expect_false(is.null(configs$transcriptomics$norm))
    expect_identical(validateOmicGuiConfig(configs$transcriptomics), configs$transcriptomics)
    lipid <- configs$lipidomics$norm
    metab <- configs$metabolomics$norm
    for (omic in c("lipidomics", "metabolomics")) {
        norm <- configs[[omic]]$norm
        expect_true(all(guiConfigNormRequired %in% names(norm)), info = omic)
        # The tab's module id, as the shell's.
        tab <- Filter(function(tab) identical(tab$status_key, "normalization"), configs[[omic]]$shell$tabs)
        expect_identical(norm$module_id, tab[[1L]]$module_id, info = omic)
        expect_identical(norm$module_id, "norm")
        # The option sections in order, and their options (mod_lipid_norm_ui_helpers.R:1-264;
        # mod_metab_norm_ui.R:25-289).
        expect_identical(names(norm$panels), names(guiConfigNormPanels), info = omic)
        expect_identical(norm$panels$itsd$aggregation, c("Median (recommended)" = "median", "Mean" = "mean", "Sum" = "sum"))
        expect_identical(norm$panels$log2$offset, 1)
        expect_identical(norm$panels$between_sample$methods, c("Cyclic Loess" = "cyclicloess", "Quantile" = "quantile",
            "Scale (Median Absolute Values)" = "scale", "None" = "none"))
        expect_identical(norm$panels$between_sample$default, "cyclicloess")
        expect_identical(norm$panels$ruv_iii[c("modes", "default", "percentage", "k")], list(
            modes = c("Automatic (recommended)" = "automatic", "Manual (advanced users)" = "manual", "Skip RUV" = "skip"),
            default = "automatic", percentage = 5, k = 2))
        expect_identical(norm$panels$export$button, "export_session")
        # The QC tabs, the stage grid and the RUV tab's texts (buildLipidNormQcTabsetPanel, :285-307).
        expect_identical(norm$plot_tabs, c("itsd", "ruv", "pca", "density", "rle", "correlation", "correlation_filter", "log"))
        expect_identical(norm$stage_labels, c("Post-Filtering", "Post-Normalization", "RUV-Corrected"))
        expect_identical(norm$assay_slots, 2L)
        expect_identical(norm$ruv_heading, "RUV Parameter Optimization Results")
        expect_identical(norm$ruv_help,
            "Canonical correlation plots and optimization summary per assay. Each assay is optimized independently.")
        # The managed run and correlation are the core's; the tab opens with the pre-normalisation QC.
        expect_identical(norm$run$managed, runWorkflowOmicsManagedNormalisation)
        expect_identical(norm$correlation$managed, runWorkflowOmicsManagedCorrelation)
        expect_true(norm$pre_qc_on_entry)
    }
    # Per omic: the RUV noun, the test ids and the export (P-18) and D7 values.
    expect_identical(lipid$panels$ruv_iii$noun, "Lipids")
    expect_identical(metab$panels$ruv_iii$noun, "Metabolites")
    expect_identical(lipid$testids, c(method = "lipid-norm-method", ruv_mode = "lipid-norm-ruv-mode",
        run = "lipid-norm-run", export = "lipid-norm-export-session"))
    expect_identical(metab$testids, c(method = "metab-norm-method", ruv_mode = "metab-norm-ruv-mode",
        run = "metab-norm-run", export = "metab-norm-export-session"))
    expect_identical(lipid$export[c("prefix", "counts_field", "samples")],
        list(prefix = "lipid_", counts_field = "lipid_counts", samples = "design"))
    expect_identical(metab$export[c("prefix", "counts_field", "samples")],
        list(prefix = "metab_", counts_field = "metabolite_counts", samples = "columns"))
    expect_identical(lipid$skip_correlation_flags, "before_commit")
    expect_identical(metab$skip_correlation_flags, "after_commit")
    expect_identical(lipid$final_pca, list(preparing_size = 5, error_size = 4))
    expect_identical(metab$final_pca, list(preparing_size = 5, error_size = 5))
    # KIT-105 W5 (D7): the export's progress display and summary time, and the managed correlation's failure notice.
    expect_identical(lipid$export[c("progress", "summary_time")], list(progress = list(gather = "Gathering data and publishing outputs...",
        complete = "Export complete", steps = FALSE), summary_time = "session"))
    expect_identical(metab$export[c("progress", "summary_time")], list(progress = list(gather = "Gathering data...",
        complete = character(0), steps = TRUE), summary_time = "written"))
    expect_false(lipid$managed_correlation_notify)
    expect_true(metab$managed_correlation_notify)
    # The ITSD selection table calls the omic's core builder with its id-column argument (D8).
    testthat::local_mocked_bindings(
        buildLipidItsdSelectionTable = function(assay_data, lipid_id_col, annotation_cols = NULL, is_patterns = NULL) {
            list("lipid", assay_data, lipid_id_col, annotation_cols)
        },
        buildItsdSelectionTable = function(assay_data, metabolite_id_col, annotation_cols = NULL, is_patterns = NULL) {
            list("metab", assay_data, metabolite_id_col, annotation_cols)
        })
    expect_identical(lipid$itsd_table("assay", "lipid_id", "lipid"), list("lipid", "assay", "lipid_id", "lipid"))
    expect_identical(metab$itsd_table("assay", "database_identifier", "Name"),
        list("metab", "assay", "database_identifier", "Name"))
})

test_that("each norm field, required or conditional, is refused by class when missing or mistyped", {
    for (omic in c("lipidomics", "metabolomics")) {
        config <- omicGuiConfig(omic)
        for (field in guiConfigNormRequired) {
            removed <- config
            removed$norm[[field]] <- NULL
            expectGuiConfigRefused(removed, paste0("norm$", field))
            mistyped <- config
            mistyped$norm[field] <- list(NA)
            expectGuiConfigRefused(mistyped, paste0("norm$", field))
        }
    }
    lipid <- omicGuiConfig("lipidomics")
    # A group that is neither NULL nor a named list.
    bad <- lipid
    bad$norm <- "norm"
    expectGuiConfigRefused(bad, "norm")
    # Each option section's options.
    for (panel in names(guiConfigNormPanels)) {
        for (option in guiConfigNormPanels[[panel]]) {
            removed <- lipid
            removed$norm$panels[[panel]][[option]] <- NULL
            expectGuiConfigRefused(removed, paste0("norm$panels$", panel, "$", option))
        }
    }
    # The sections come from the panel library, in its order; the defaults are among the choices.
    bad <- lipid
    bad$norm$panels <- c(bad$norm$panels, list(tmm = list()))
    expectGuiConfigRefused(bad, "norm$panels")
    bad <- lipid
    bad$norm$panels <- bad$norm$panels[c("itsd", "aesthetics", "log2", "between_sample", "ruv_iii", "actions", "export")]
    expectGuiConfigRefused(bad, "norm$panels")
    bad <- lipid
    bad$norm$panels$between_sample$default <- "tmm"
    expectGuiConfigRefused(bad, "norm$panels$between_sample$default")
    bad <- lipid
    bad$norm$panels$ruv_iii$default <- "always"
    expectGuiConfigRefused(bad, "norm$panels$ruv_iii$default")
    # The QC tabs from the frame's set, once each and in order; three stage captions; whole assay rows.
    bad <- lipid
    bad$norm$plot_tabs <- c(bad$norm$plot_tabs, "tmm")
    expectGuiConfigRefused(bad, "norm$plot_tabs")
    bad <- lipid
    bad$norm$plot_tabs <- rev(bad$norm$plot_tabs)
    expectGuiConfigRefused(bad, "norm$plot_tabs")
    bad <- lipid
    bad$norm$stage_labels <- bad$norm$stage_labels[1:2]
    expectGuiConfigRefused(bad, "norm$stage_labels")
    bad <- lipid
    bad$norm$assay_slots <- 1.5
    expectGuiConfigRefused(bad, "norm$assay_slots")
    # The managed functions are required, and so (since KIT-105 W6) are the in-process twins and the
    # export's collector and saver, each wired to the core.
    for (field in c("run", "correlation")) {
        removed <- lipid
        removed$norm[[field]]$managed <- NULL
        expectGuiConfigRefused(removed, paste0("norm$", field, "$managed"))
    }
    for (field in names(guiConfigNormCoreFns)) {
        for (entry in guiConfigNormCoreFns[[field]]) {
            for (omic in c("lipidomics", "metabolomics")) {
                expect_identical(omicGuiConfig(omic)$norm[[field]][[entry]],
                    get(guiConfigNormWired[[paste(field, entry, sep = ".")]], envir = asNamespace("MultiScholaR")),
                    info = paste(omic, field, entry))
            }
            unwired <- lipid
            unwired$norm[[field]][[entry]] <- NULL
            expectGuiConfigRefused(unwired, paste0("norm$", field, "$", entry))
            bad <- lipid
            bad$norm[[field]][[entry]] <- "a function name"
            expectGuiConfigRefused(bad, paste0("norm$", field, "$", entry))
            wired <- lipid
            wired$norm[[field]][[entry]] <- function(...) NULL
            expect_identical(validateOmicGuiConfig(wired), wired, info = paste(field, entry))
        }
    }
    # The export's values and choices, the test ids and the D7 values.
    for (field in c("prefix", "counts_field", "samples")) {
        removed <- lipid
        removed$norm$export[[field]] <- NULL
        expectGuiConfigRefused(removed, paste0("norm$export$", field))
    }
    bad <- lipid
    bad$norm$export$samples <- "assay"
    expectGuiConfigRefused(bad, "norm$export$samples")
    # KIT-105 W5: the export's progress display (each of its fields) and summary time.
    for (omic in c("lipidomics", "metabolomics")) {
        config <- omicGuiConfig(omic)
        for (field in c("progress", "summary_time")) {
            removed <- config
            removed$norm$export[[field]] <- NULL
            expectGuiConfigRefused(removed, paste0("norm$export$", field))
            mistyped <- config
            mistyped$norm$export[field] <- list(NA)
            expectGuiConfigRefused(mistyped, paste0("norm$export$", field))
        }
        for (field in c("gather", "complete", "steps")) {
            removed <- config
            removed$norm$export$progress[[field]] <- NULL
            expectGuiConfigRefused(removed, paste0("norm$export$progress$", field))
            mistyped <- config
            mistyped$norm$export$progress[field] <- list(NA)
            expectGuiConfigRefused(mistyped, paste0("norm$export$progress$", field))
        }
    }
    bad <- lipid
    bad$norm$export$summary_time <- "export"
    expectGuiConfigRefused(bad, "norm$export$summary_time")
    for (field in c("method", "ruv_mode", "run", "export")) {
        removed <- lipid
        removed$norm$testids <- removed$norm$testids[names(removed$norm$testids) != field]
        expectGuiConfigRefused(removed, paste0("norm$testids$", field))
    }
    bad <- lipid
    bad$norm$skip_correlation_flags <- "never"
    expectGuiConfigRefused(bad, "norm$skip_correlation_flags")
    for (field in c("preparing_size", "error_size")) {
        removed <- lipid
        removed$norm$final_pca[[field]] <- NULL
        expectGuiConfigRefused(removed, paste0("norm$final_pca$", field))
    }
    bad <- lipid
    bad$norm$itsd_table <- "buildLipidItsdSelectionTable"
    expectGuiConfigRefused(bad, "norm$itsd_table")
})

# ------------------------------------------------------------------------------------------------
# KIT-110 W3: trx's norm group (D2, D4, D7) and the validator's per-tab and per-panel fields (D8).
# ------------------------------------------------------------------------------------------------

# D8: the fields required only with the tab that uses them, and the test ids by panel.
guiConfigNormByTab <- list(itsd = "itsd_table",
    correlation_filter = c("correlation", "skip_correlation_flags", "final_pca", "managed_correlation_notify"))
guiConfigNormTestidPanels <- list(method = "between_sample", ruv_mode = c("ruv_iii", "ruv_grid"), export = "export")

test_that("KIT-110 W3: trx's norm group holds today's trx tab values (mod_trx_norm.R, mod_transcriptomics.R at 38d9d25e)", {
    norm <- omicGuiConfig("transcriptomics")$norm
    expect_identical(norm$module_id, "norm")
    # D2: TMM/logCPM, the RUV grid and the actions; the grid's defaults (mod_trx_norm.R:26-47).
    expect_identical(names(norm$panels), c("tmm_logcpm", "ruv_grid", "actions"))
    expect_identical(norm$panels$tmm_logcpm, list())
    expect_identical(norm$panels$ruv_grid, list(modes = c("Automatic (recommended)" = "automatic",
        "Manual (advanced users)" = "manual"), default = "automatic", percentage_min = 5, percentage_max = 30,
        percentage_step = 5, k_penalty_weight = 0.05, percentage = 10, k = 2))
    expect_identical(norm$panels$actions, list())
    # The profile's tabs, captions and RUV texts (mod_transcriptomics.R:48-52), one assay row.
    expect_identical(norm$plot_tabs, c("ruv", "pca", "rle", "log"))
    expect_identical(norm$stage_labels, c("Raw logCPM", "TMM logCPM", "RUVs-corrected logCPM (display only)"))
    expect_identical(norm$assay_slots, 1L)
    expect_identical(norm$ruv_heading, "RUVs Control Selection (canonical correlation)")
    expect_identical(norm$ruv_help, paste("The CCA panels and the grid of percentages: the controls are inferred negative",
        "controls (the least significant genes in a preliminary group GLM)."))
    # D4: trx's own handlers, the core's.
    expect_identical(norm$run, list(managed = runTrxNormManaged, in_process = runTrxNormInProcess,
        settings = trxNormGridParameters, handler = runTrxNormalisation))
    expect_identical(norm$reset, resetTrxNormalisation)
    # D3, D7: the fixed assay, the image stem, the stage prefixes, the registered rows (P-21), no pre-QC on entry.
    expect_false(norm$pre_qc_on_entry)
    expect_identical(norm$assay_names, "counts")
    expect_identical(norm$image_stem, "trx")
    expect_identical(norm$stage_prefixes, c(pre_norm = "raw", post_norm = "tmm", ruv_corrected = "ruvs"))
    expect_identical(norm$output_assay_slots, 2L)
    expect_identical(norm$testids, c(ruv_mode = "trx-norm-ruv-mode", run = "trx-norm-run", reset = "trx-norm-reset"))
    # No tab or panel needs the rest.
    expect_identical(intersect(c("correlation", "export", "itsd_table", "final_pca", "skip_correlation_flags",
        "managed_correlation_notify"), names(norm)), character())
})

test_that("KIT-110 W3 (D8): the validator requires a norm field with the tab or panel that uses it", {
    trx <- omicGuiConfig("transcriptomics")
    lipid <- omicGuiConfig("lipidomics")
    # Always required, for trx too.
    for (field in c("module_id", "panels", "plot_tabs", "stage_labels", "assay_slots", "ruv_heading", "ruv_help", "run",
        "pre_qc_on_entry", "testids")) {
        removed <- trx
        removed$norm[[field]] <- NULL
        expectGuiConfigRefused(removed, paste0("norm$", field))
    }
    for (field in c("managed", "in_process")) {
        removed <- trx
        removed$norm$run[[field]] <- NULL
        expectGuiConfigRefused(removed, paste0("norm$run$", field))
    }
    removed <- trx
    removed$norm$testids <- removed$norm$testids[names(removed$norm$testids) != "run"]
    expectGuiConfigRefused(removed, "norm$testids$run")
    # Each conditional field is refused when its tab is present (lipid; trx with the tab added), and trx validates
    # without it.
    for (tab in names(guiConfigNormByTab)) {
        for (field in guiConfigNormByTab[[tab]]) {
            expect_false(field %in% names(trx$norm), info = field)
            removed <- lipid
            removed$norm[[field]] <- NULL
            expectGuiConfigRefused(removed, paste0("norm$", field))
            # Trx with the tab added and the tab's other fields (lipid's values): this one is refused.
            with_tab <- trx
            with_tab$norm$plot_tabs <- intersect(c("itsd", "ruv", "pca", "rle", "correlation_filter", "log"),
                c(trx$norm$plot_tabs, tab))
            others <- setdiff(guiConfigNormByTab[[tab]], field)
            with_tab$norm[others] <- lipid$norm[others]
            expectGuiConfigRefused(with_tab, paste0("norm$", field))
            with_tab$norm[field] <- lipid$norm[field]
            expect_identical(validateOmicGuiConfig(with_tab), with_tab, info = field)
        }
    }
    with_export <- trx
    with_export$norm$panels$export <- list(button = "export_session")
    with_export$norm$testids[["export"]] <- "trx-norm-export-session"
    expectGuiConfigRefused(with_export, "norm$export")
    # The test ids by panel: method with between_sample, ruv_mode with ruv_iii or ruv_grid, export with export.
    for (id in names(guiConfigNormTestidPanels)) {
        removed <- lipid
        removed$norm$testids <- removed$norm$testids[names(removed$norm$testids) != id]
        expectGuiConfigRefused(removed, paste0("norm$testids$", id))
    }
    removed <- trx
    removed$norm$testids <- removed$norm$testids[names(removed$norm$testids) != "ruv_mode"]
    expectGuiConfigRefused(removed, "norm$testids$ruv_mode")
    # The optional fields, checked by type when present.
    optional <- list(c("run", "settings"), c("run", "handler"), "reset", "assay_names", "image_stem", "stage_prefixes",
        "output_assay_slots")
    for (path in optional) {
        removed <- trx
        mistyped <- trx
        if (length(path) == 2L) {
            removed$norm[[path[[1L]]]][[path[[2L]]]] <- NULL
            mistyped$norm[[path[[1L]]]][path[[2L]]] <- list(NA)
        } else {
            removed$norm[[path]] <- NULL
            mistyped$norm[path] <- list(NA)
        }
        expect_identical(validateOmicGuiConfig(removed), removed, info = paste(path, collapse = "$"))
        expectGuiConfigRefused(mistyped, paste0("norm$", paste(path, collapse = "$")))
    }
    bad <- trx
    bad$norm$image_stem <- c("trx", "counts")
    expectGuiConfigRefused(bad, "norm$image_stem")
    bad <- trx
    bad$norm$output_assay_slots <- 1.5
    expectGuiConfigRefused(bad, "norm$output_assay_slots")
    bad <- trx
    bad$norm$stage_prefixes <- c(raw = "raw", tmm = "tmm", ruvs = "ruvs")
    expectGuiConfigRefused(bad, "norm$stage_prefixes")
    bad <- trx
    bad$norm$testids[["reset"]] <- ""
    expectGuiConfigRefused(bad, "norm$testids$reset")
    removed <- trx
    removed$norm$testids <- removed$norm$testids[names(removed$norm$testids) != "reset"]
    expect_identical(validateOmicGuiConfig(removed), removed)
    # The two new option sections: each of the grid's options; the default among its modes; the library's order.
    for (option in c("modes", "default", "percentage_min", "percentage_max", "percentage_step", "k_penalty_weight",
        "percentage", "k")) {
        removed <- trx
        removed$norm$panels$ruv_grid[[option]] <- NULL
        expectGuiConfigRefused(removed, paste0("norm$panels$ruv_grid$", option))
    }
    bad <- trx
    bad$norm$panels$ruv_grid$default <- "skip"
    expectGuiConfigRefused(bad, "norm$panels$ruv_grid$default")
    bad <- trx
    bad$norm$panels <- bad$norm$panels[c("ruv_grid", "tmm_logcpm", "actions")]
    expectGuiConfigRefused(bad, "norm$panels")
})

# ------------------------------------------------------------------------------------------------
# KIT-106 W1: the da group (D1) and its validator (D8), written from the ticket's table at the start
# commit b75561de.
# ------------------------------------------------------------------------------------------------

test_that("KIT-106 W1 (D1): the da group holds each omic's DA tab values at the start commit", {
    configs <- lapply(stats::setNames(guiConfigOmics, guiConfigOmics), function(omic) omicGuiConfig(omic)$da)
    lipid <- configs$lipidomics
    metab <- configs$metabolomics
    trx <- configs$transcriptomics
    prot <- configs$proteomics
    # The tab and status key: the shell's DA tab's, for every omic.
    for (omic in guiConfigOmics) {
        tab <- Filter(function(tab) tab$component %in% c("lipid_da", "metab_da", "trx_da", "prot_da"),
            omicGuiConfig(omic)$shell$tabs)[[1L]]
        expect_identical(c(configs[[omic]]$module_id, configs[[omic]]$status_key), c(tab$module_id, tab$status_key), info = omic)
    }
    expect_identical(c(lipid$module_id, lipid$status_key), c("de", "differential_analysis"))
    expect_identical(c(prot$module_id, prot$status_key), c("differential_expression", "differential_expression"))
    expect_identical(prot$status_keys, c("differential_abundance", "enrichment_analysis"))
    # The test ids (trx has no load section).
    for (prefix in c("lipid", "metab", "prot")) {
        da <- list(lipid = lipid, metab = metab, prot = prot)[[prefix]]
        expect_identical(da$testids, c(load_session = paste0(prefix, "-da-load-session"), run = paste0(prefix, "-da-run"),
            download = paste0(prefix, "-da-download"), results_tabs = paste0(prefix, "-da-results-tabs")), info = prefix)
    }
    expect_identical(trx$testids, c(run = "trx-da-run", download = "trx-da-download", results_tabs = "trx-da-results-tabs"))
    # The model: the formula input (lipid, metab, proteomics), or trx's model shown with its realised formula.
    for (da in list(lipid, metab)) {
        expect_identical(da$model, list(editable = TRUE, label = "Model Formula:", default = "~ 0 + group",
            help = "Formula fetched from S4 object @args if available"))
    }
    expect_identical(prot$model$help, "Formula fetched from S4 object @args")
    expect_false(trx$model$editable)
    expect_identical(trx$model[c("heading", "output")], list(heading = "Model", output = "formula_display"))
    expect_identical(trx$model$formula, trxRealisedFormula)
    # The threshold and fold-change inputs: today's values; trx's labels, and its fold change is not tested.
    numbers <- list(value = 0.05, min = 0.001, max = 0.2, step = 0.005)
    for (da in configs) expect_identical(da$threshold[names(numbers)], numbers)
    for (da in configs) expect_identical(da$lfc[c("value", "min", "max", "step")], list(value = 0, min = 0, max = 2, step = 0.1))
    expect_identical(c(lipid$threshold$label, metab$threshold$label, prot$threshold$label), rep("Q-value Threshold:", 3))
    expect_identical(c(lipid$lfc$label, metab$lfc$label, prot$lfc$label), rep("Log Fold-Change Cutoff:", 3))
    expect_identical(trx$threshold$label, "FDR (BH) threshold:")
    expect_identical(trx$lfc$label, "Fold-change line (display only; no test):")
    expect_identical(vapply(configs, function(da) da$lfc$tested, TRUE),
        c(proteomics = TRUE, lipidomics = TRUE, metabolomics = TRUE, transcriptomics = FALSE))
    # The table's columns, the summary and the nouns.
    expect_identical(lipid$columns, list(id = "lipid_id", label = "lipid_name", fold_change = "logFC", pvalue = "raw_pvalue",
        significance = "fdr_qvalue"))
    expect_identical(metab$columns[c("id", "label")], list(id = "metabolite_id", label = "metabolite_name"))
    expect_identical(trx$columns, list(id = "GeneID", label = "GeneName", fold_change = "logFC", pvalue = "PValue",
        significance = "FDR", significance_label = "FDR (BH)"))
    expect_identical(prot$columns, list(id_slot = "protein_id_column", fold_change = "log2FC", pvalue = "raw_pvalue",
        significance = "fdr_qvalue"))
    expect_identical(lapply(configs, `[[`, "summary"), list(proteomics = list(noun = "genes", significance_label = "q"),
        lipidomics = list(noun = "lipids", significance_label = "Q"), metabolomics = list(noun = "metabolites", significance_label = "Q"),
        transcriptomics = list(noun = "genes", significance_label = "FDR (BH)")))
    expect_identical(lipid$nouns, c(top_n = "Top N Lipids:", scaling = "Row (lipid) scaling", labels = "Show Lipid Labels",
        cluster = "lipids"))
    expect_identical(metab$nouns, c(top_n = "Top N Metabolites:", scaling = "Row (metabolite) scaling",
        labels = "Show Metabolite Labels", cluster = "metabolites"))
    expect_identical(prot$nouns, c(top_n = "Top N Genes:", scaling = "Row (gene) scaling", labels = "Show Gene Labels",
        label_top = "Label Top Proteins", cluster = "proteins"))
    expect_null(trx$nouns)
    # The tabs, the assay and table-filter rules.
    expect_identical(lipid$tabs, c(volcano = "Volcano Plot", heatmap = "Heatmap", table = "DA Results Table"))
    expect_identical(trx$tabs, c(volcano = "Volcano Plot", table = "DA Results Table"))
    expect_identical(vapply(configs, `[[`, "", "assays"),
        c(proteomics = "none", lipidomics = "multi", metabolomics = "multi", transcriptomics = "single"))
    expect_identical(vapply(configs, `[[`, "", "table_filter"),
        c(proteomics = "store_page", lipidomics = "all_only", metabolomics = "all_or_combined", transcriptomics = "all_only"))
    # The views: the static volcano, the heatmap with its tree-cut defaults, the table's defaults.
    expect_identical(lipid$views$static_volcano, list(output = "volcano_static", plot = generateLipidDAVolcanoStatic,
        show_labels = TRUE, n_labels = 15))
    expect_identical(metab$views$static_volcano$plot, generateMetabDAVolcanoStatic)
    expect_identical(trx$views$static_volcano, list(output = "volcano_static", plot = trxDaStaticVolcano,
        show_labels = TRUE, n_labels = 15))
    expect_identical(prot$views$static_volcano, list(output = "volcano_plot_static", plot = generateProtDAVolcanoStatic,
        show_labels = TRUE, n_labels = 10))
    expect_identical(lipid$views$heatmap, list(kernel = generateLipidDAHeatmap, labels_argument = "show_lipid_names",
        tree_cut = list(default = "none", n_clusters = 4, cut_height = 10, min_cluster_size = 5), save_prefix = "lipid_",
        save_note = TRUE))
    expect_identical(metab$views$heatmap[c("kernel", "labels_argument", "save_prefix", "save_note")],
        list(kernel = generateMetabDAHeatmap, labels_argument = "show_metabolite_names", save_prefix = "metab_", save_note = TRUE))
    expect_identical(metab$views$heatmap$tree_cut, lipid$views$heatmap$tree_cut)
    expect_identical(prot$views$heatmap, list(kernel = generateProtDAHeatmap, labels_argument = "show_gene_names",
        tree_cut = list(default = "k_clusters", n_clusters = 4, cut_height = 0.5, min_cluster_size = 3), save_prefix = "prot_",
        save_note = FALSE))
    expect_null(trx$views$heatmap)
    for (da in configs) expect_identical(da$views$table, list(significance = "significant", max_rows = 1000))
    # The run (proteomics: its plugin), the source, the session load and the preflight.
    for (da in list(lipid, metab, trx)) expect_identical(da$run$managed, runWorkflowOmicsManagedDa)
    expect_identical(c(lipid$run$in_process, metab$run$in_process, trx$run$in_process),
        c(runLipidReferenceDaAnalysis, runMetabReferenceDaAnalysis, runTrxReferenceDaAnalysis))
    expect_identical(lipid$source, list(read_metadata = readLipidWorkflowMetadata, class = "LipidomicsAssayData",
        no_data = "No lipidomics data loaded. Please load a filtered session first."))
    expect_identical(metab$source$read_metadata, readMetabWorkflowMetadata)
    # KIT-106 W2 merge: trx's source read, model text, static volcano, contrast entries and export report are wired to
    # the core (utils_workflow_trx_da.R), and lipid's preflight (func_omics_design_shared.R).
    expect_identical(trx$source, list(class = "TranscriptomicsAssayData", read_metadata = trxDaReadMetadata,
        no_data = "Run the normalisation first: the DA runs on the normalised (RUVs) state."))
    expect_identical(trx$model$text, trxDaModelText)
    expect_identical(trx$views$static_volcano$plot, trxDaStaticVolcano)
    expect_identical(c(trx$contrast_entries, trx$contrast_entry), c(trxDaContrastEntries, trxDaContrastEntry))
    expect_identical(trx$run$export_report, trxDaExport)
    expect_identical(lipid$preflight, validateLipidDaAnalysisPreflight)
    expect_identical(lipid$session_load, list(file = "lipid_filtered_session_data_latest.rds", assay_rule = "setequal",
        assay_omic_type = "lipidomics",
        assay_message = "Loaded lipid session assay_names do not match S4 lipid_data assays: session [%s], S4 [%s]"))
    expect_identical(metab$session_load, list(file = "metab_filtered_session_data_latest.rds", assay_rule = "identical",
        assay_message = "Session assay_names do not match current_s4_object assays: session=%s; object=%s",
        assay_check_first = TRUE))
    expect_null(trx$session_load)
    expect_identical(c(lipid$preflight_first, metab$preflight_first, trx$preflight_first), c(FALSE, TRUE, FALSE))
    # KIT-106 W4: metab's preflight adapts the frame's call (the state read for the run, the contrasts, the formula) to
    # validateMetabDesignDaPreflight's arguments, as the start commit's tab called it; without a formula it is not called.
    expect_identical(names(formals(metab$preflight)), c("currentS4", "contrastsTbl", "formulaString"))
    preflight_args <- NULL
    testthat::local_mocked_bindings(validateMetabDesignDaPreflight = function(...) {
        preflight_args <<- list(...)
        list(valid = FALSE, errors = "refused")
    })
    object <- module_ci_metab_da_object()
    contrasts <- module_ci_metab_da_contrasts("multi_group")
    expect_identical(metab$preflight(currentS4 = object, contrastsTbl = contrasts, formulaString = NULL), list(valid = TRUE))
    expect_null(preflight_args)
    expect_identical(metab$preflight(currentS4 = object, contrastsTbl = contrasts, formulaString = "~ 0 + group"),
        list(valid = FALSE, errors = "refused"))
    expect_identical(preflight_args, list(designMatrix = object@design_matrix, assayList = object@metabolite_data,
        contrastsTbl = contrasts, formulaString = "~ 0 + group", sampleIdCol = object@sample_id, groupCol = object@group_id,
        requireContrasts = TRUE))
    # The texts, the download and proteomics' plugins.
    for (da in list(lipid, metab)) {
        expect_identical(da$texts, c(status_ready = "Data loaded. Ready to run analysis.",
            status_waiting = "Waiting for data.\nClick 'Load Filtered Session' to begin.",
            contrasts_none = "No contrasts defined.\nLoad a filtered session or define contrasts in the design tab."))
    }
    expect_identical(trx$texts[c("status_ready", "status_waiting")], rep("Ready. The DA runs on the normalised (RUVs) state.", 2) |>
        stats::setNames(c("status_ready", "status_waiting")))
    expect_identical(prot$texts, c(status_waiting = "[WAITING] Waiting for analysis...\nClick 'Run DA Analysis' to start",
        contrasts_none = "No contrasts available.\nComplete normalization and\ncorrelation filtering first."))
    expect_identical(prot$status_text, diaDaStatusText)
    expect_identical(lapply(configs, `[[`, "download"), list(proteomics = list(stem = "DA_results_", extension = "zip"),
        lipidomics = list(stem = "lipidomics_da_results_", extension = "csv"),
        metabolomics = list(stem = "metabolomics_da_results_", extension = "csv"),
        transcriptomics = list(stem = "transcriptomics_da_results_", extension = "csv")))
    expect_identical(prot$plugins, c(init = "da_server_init_handlers", load = "da_server_load_session_handler",
        run = "da_server_run_analysis_handler", archive = "prot_da_session_archive"))
    expect_null(lipid$plugins)
    # Each named plugin handler is today's server's (the archive chooser is W5's to name).
    for (handler in prot$plugins[c("init", "load", "run")]) expect_true(is.function(get(handler, envir = asNamespace("MultiScholaR"))))
})

test_that("KIT-106 W1 (D8): each da field, required or conditional, is refused by class when missing or mistyped; the optional ones by type", {
    lipid <- omicGuiConfig("lipidomics")
    trx <- omicGuiConfig("transcriptomics")
    prot <- omicGuiConfig("proteomics")
    refused <- function(config, path, value = NULL, mistyped = list(NA)) {
        parts <- strsplit(path, "$", fixed = TRUE)[[1L]]
        removed <- config
        removed[["da"]][[parts]] <- NULL
        expectGuiConfigRefused(removed, paste0("da$", path))
        bad <- config
        bad[["da"]][[parts]] <- mistyped[[1L]]
        expectGuiConfigRefused(bad, paste0("da$", path))
    }
    # The always-required fields' members.
    for (path in c("model$editable", "model$help", "model$label", "model$default", "threshold$label", "threshold$value",
        "threshold$min", "threshold$max", "threshold$step", "lfc$label", "lfc$value", "lfc$tested", "columns$fold_change",
        "columns$pvalue", "columns$significance", "summary$noun", "summary$significance_label", "views$static_volcano",
        "views$table", "views$static_volcano$output", "views$static_volcano$show_labels", "views$static_volcano$n_labels",
        "views$table$significance", "views$table$max_rows", "download$stem", "download$extension",
        # KIT-106 W6: required since the W2 merge wired them for every omic that uses them.
        "views$static_volcano$plot", "source$read_metadata")) {
        refused(lipid, path)
    }
    # KIT-106 W6: a shown model's text and realised formula are required (trx's).
    for (path in c("model$heading", "model$output", "model$text", "model$formula")) refused(trx, path)
    for (key in c("status_waiting", "contrasts_none")) {
        bad <- lipid
        bad$da$texts <- bad$da$texts[names(bad$da$texts) != key]
        expectGuiConfigRefused(bad, paste0("da$texts$", key))
    }
    # The choices and the tabs from the frame's set in its order.
    for (field in c("assays", "table_filter")) {
        bad <- lipid
        bad$da[[field]] <- "every"
        expectGuiConfigRefused(bad, paste0("da$", field))
    }
    bad <- lipid
    bad$da$tabs <- rev(bad$da$tabs)
    expectGuiConfigRefused(bad, "da$tabs")
    bad <- lipid
    bad$da$tabs <- bad$da$tabs[c("heatmap", "table")]
    expectGuiConfigRefused(bad, "da$tabs")
    bad <- lipid
    bad$da$tabs <- c(bad$da$tabs, summary = "Summary")
    expectGuiConfigRefused(bad, "da$tabs")
    # The identifier column, by its name or by the object's slot (one of the two).
    bad <- lipid
    bad$da$columns$id <- NULL
    expectGuiConfigRefused(bad, "da$columns")
    bad <- prot
    bad$da$columns$id <- "uniprot_acc"
    expectGuiConfigRefused(bad, "da$columns")
    # The heatmap view with the heatmap tab: required for lipid and proteomics; trx has neither, and validates.
    for (config in list(lipid, prot)) {
        for (path in c("views$heatmap", "views$heatmap$kernel", "views$heatmap$labels_argument", "views$heatmap$tree_cut",
            "views$heatmap$save_prefix", "views$heatmap$save_note", "views$heatmap$tree_cut$default",
            "views$heatmap$tree_cut$n_clusters", "views$heatmap$tree_cut$cut_height", "views$heatmap$tree_cut$min_cluster_size")) {
            refused(config, path)
        }
    }
    bad <- lipid
    bad$da$views$heatmap$tree_cut$default <- "always"
    expectGuiConfigRefused(bad, "da$views$heatmap$tree_cut$default")
    with_tab <- trx
    with_tab$da$tabs <- c(volcano = "Volcano Plot", heatmap = "Heatmap", table = "DA Results Table")
    expectGuiConfigRefused(with_tab, "da$views$heatmap")
    with_tab$da$views$heatmap <- lipid$da$views$heatmap
    expect_identical(validateOmicGuiConfig(with_tab), with_tab)
    # The session load with the load section (its test id), unless the load plugin renders it.
    for (path in c("session_load", "session_load$file", "session_load$assay_rule", "session_load$assay_message")) refused(lipid, path)
    bad <- lipid
    bad$da$session_load$assay_rule <- "subset"
    expectGuiConfigRefused(bad, "da$session_load$assay_rule")
    with_load <- trx
    with_load$da$testids[["load_session"]] <- "trx-da-load-session"
    expectGuiConfigRefused(with_load, "da$session_load")
    expect_null(prot$da$session_load)
    # The run and the source, unless a plugin replaces them.
    for (path in c("run", "run$managed", "run$in_process", "source", "source$class", "source$no_data")) {
        refused(lipid, path)
        refused(trx, path)
    }
    without_plugins <- prot
    without_plugins$da$plugins <- without_plugins$da$plugins[c("load", "archive")]
    expectGuiConfigRefused(without_plugins, "da$run")
    without_plugins$da$run <- lipid$da$run
    expectGuiConfigRefused(without_plugins, "da$source")
    without_plugins$da$source <- lipid$da$source
    expect_identical(validateOmicGuiConfig(without_plugins), without_plugins)
    bad <- prot
    bad$da$plugins[["export"]] <- "prot_da_export"
    expectGuiConfigRefused(bad, "da$plugins")
    # The optional fields, checked by type when present (the function fields W2 wires among them).
    optional <- list(lipid = c("status_keys", "nouns", "status_text", "preflight", "preflight_first", "contrast_entries",
        "contrast_entry", "plugins", "model$text", "model$formula", "run$export_report",
        "session_load$assay_omic_type", "columns$label", "columns$significance_label",
        "texts$status_ready"))
    for (path in optional$lipid) {
        parts <- strsplit(path, "$", fixed = TRUE)[[1L]]
        removed <- lipid
        bad <- lipid
        if (identical(path, "texts$status_ready")) {
            removed$da$texts <- removed$da$texts[names(removed$da$texts) != "status_ready"]
            bad$da$texts[["status_ready"]] <- ""
        } else {
            removed[["da"]][[parts]] <- NULL
            bad[["da"]][[parts]] <- if (identical(path, "plugins")) list(NA) else NA
        }
        expect_identical(validateOmicGuiConfig(removed), removed, info = path)
        expectGuiConfigRefused(bad, paste0("da$", path))
    }
})

# ------------------------------------------------------------------------------------------------
# KIT-107 W1: the enrichment group (D1) and its validator (D9), written from the ticket's table at the start
# commit 8e7294b6.
# ------------------------------------------------------------------------------------------------

test_that("KIT-107 W1 (D1): the enrichment group holds proteomics' enrichment tab values at the start commit; NULL for the other omics", {
    for (omic in c("lipidomics", "metabolomics", "transcriptomics")) expect_null(omicGuiConfig(omic)$enrichment, info = omic)
    prot <- omicGuiConfig("proteomics")
    enrichment <- prot$enrichment
    # The tab: the shell's enrichment tab's module id, status key and value.
    tab <- Filter(function(tab) identical(tab$component, "prot_enrich"), prot$shell$tabs)[[1L]]
    expect_identical(c(enrichment$module_id, enrichment$status_key, enrichment$tab_value),
        c(tab$module_id, tab$status_key, tab$value))
    expect_identical(c(enrichment$module_id, enrichment$status_key, enrichment$tab_value),
        c("enrichment_analysis", "enrichment_analysis", "enrichment"))
    # The declaration and its run, the DA source.
    expect_identical(enrichment$operation, diaEnrichmentOperation)
    expect_identical(enrichment$run, runDIAReferenceEnrichment)
    expect_identical(enrichment$da_source, list(manifest = isDIAResultManifest, components = c("design_matrix", "protein_id_table"),
        read = createDAResultsForEnrichment))
    expect_identical(enrichment$operation()$name, "enrichment")
    expect_identical(enrichment$operation()$product_type, "dia-enrichment-results")
    # The identifier column, the cutoffs, the correction, the excluded IEA flag, the taxon.
    expect_identical(enrichment$id_column, list(slot = "protein_id_column", default = "uniprot_acc", gprofiler_override = "gene_name"))
    expect_identical(enrichment$cutoffs, list(up = list(value = 0, min = 0, max = 5, step = 0.1),
        down = list(value = 0, min = 0, max = 5, step = 0.1), q = list(value = 0.05, min = 0.001, max = 0.2, step = 0.005)))
    expect_identical(enrichment$correction, list(choices = c(`g:SCS (Conservative)` = "gSCS", `FDR (Less Conservative)` = "fdr"),
        selected = "gSCS"))
    expect_false(enrichment$exclude_iea)
    expect_identical(enrichment$taxon, list(default = "9606", prefill_from = "taxon_id"))
    # The result tabs (STRING-DB the placeholder), the organism filter's toggle, the completed state.
    expect_identical(enrichment$services, c("gprofiler2", "clusterprofiler", "stringdb"))
    expect_identical(enrichment$organism_filter$mixed_species_from, "mixed_species_analysis")
    expect_identical(enrichment$state_name, "enrichment_completed")
    # The declaration's parameters are the declaration's, not repeated here (D1); no labels, no da_columns.
    expect_false(any(c("parameters", "background", "labels", "da_columns", "tab") %in% names(enrichment)))
    # The values the UI renders today are these (mod_prot_enrich.R at the start commit).
    defaults <- omicsUiInputDefaults(mod_prot_enrich_ui("x"))
    expect_identical(c(defaults[["x-up_cutoff"]], defaults[["x-down_cutoff"]], defaults[["x-q_cutoff"]]),
        c(enrichment$cutoffs$up$value, enrichment$cutoffs$down$value, enrichment$cutoffs$q$value))
    expect_identical(defaults[["x-correction_method"]], enrichment$correction$selected)
    expect_identical(as.character(defaults[["x-organism_taxid"]]), enrichment$taxon$default)
})

test_that("KIT-107 W1 (D9): each enrichment field is refused by class when missing or mistyped; services and the organism filter are checked; the W2 function fields are optional, whole when present", {
    prot <- omicGuiConfig("proteomics")
    refused <- function(config, path, mistyped = list(NA)) {
        parts <- strsplit(path, "$", fixed = TRUE)[[1L]]
        removed <- config
        removed[["enrichment"]][[parts]] <- NULL
        expectGuiConfigRefused(removed, paste0("enrichment$", path))
        bad <- config
        bad[["enrichment"]][[parts]] <- mistyped[[1L]]
        expectGuiConfigRefused(bad, paste0("enrichment$", path))
    }
    # The required fields and their members (organism_filter may be NULL: its members only).
    for (path in c(setdiff(guiConfigEnrichmentRequired, "organism_filter"), "da_source$manifest", "da_source$components",
        "da_source$read", "id_column$slot", "id_column$default", "id_column$gprofiler_override", "cutoffs$up$value",
        "cutoffs$up$min", "cutoffs$up$max", "cutoffs$up$step", "cutoffs$down$value", "cutoffs$q$step", "correction$choices",
        "correction$selected", "taxon$default", "taxon$prefill_from", "organism_filter$mixed_species_from")) {
        refused(prot, path)
    }
    bad <- prot
    bad$enrichment$organism_filter <- "mixed_species_analysis"
    expectGuiConfigRefused(bad, "enrichment$organism_filter")
    # A wrong type for a list member, or a reference given by name.
    bad <- prot
    bad$enrichment$operation <- "diaEnrichmentOperation"
    expectGuiConfigRefused(bad, "enrichment$operation")
    bad <- prot
    bad$enrichment$cutoffs$q$value <- "0.05"
    expectGuiConfigRefused(bad, "enrichment$cutoffs$q$value")
    # The cutoffs name exactly up, down and q.
    bad <- prot
    bad$enrichment$cutoffs$q <- NULL
    expectGuiConfigRefused(bad, "enrichment$cutoffs")
    bad <- prot
    bad$enrichment$cutoffs$p <- bad$enrichment$cutoffs$q
    expectGuiConfigRefused(bad, "enrichment$cutoffs")
    # The correction's selection is one of its choices.
    bad <- prot
    bad$enrichment$correction$selected <- "bonferroni"
    expectGuiConfigRefused(bad, "enrichment$correction$selected")
    # The services are the tab's result tabs, once each.
    bad <- prot
    bad$enrichment$services <- c("gprofiler2", "gsea")
    expectGuiConfigRefused(bad, "enrichment$services")
    bad <- prot
    bad$enrichment$services <- c("gprofiler2", "gprofiler2")
    expectGuiConfigRefused(bad, "enrichment$services")
    # The organism filter: NULL is valid (no filter section); its W2 function fields are optional, by type when present.
    without_filter <- prot
    without_filter$enrichment["organism_filter"] <- list(NULL)
    expect_identical(validateOmicGuiConfig(without_filter), without_filter)
    for (field in c("mapping", "apply")) {
        with_ref <- prot
        with_ref$enrichment$organism_filter[[field]] <- identity
        expect_identical(validateOmicGuiConfig(with_ref), with_ref, info = field)
        bad <- prot
        bad$enrichment$organism_filter[[field]] <- "resolveProtEnrichOrganismMapping"
        expectGuiConfigRefused(bad, paste0("enrichment$organism_filter$", field))
    }
    # The W2 function groups: absent at W1; when present, whole (both members references).
    members <- list(organism = c("supported", "method"), annotation = c("resolve", "match"), kernel = c("select", "test_double"))
    for (group in names(members)) {
        expect_null(prot$enrichment[[group]], info = group)
        with_group <- prot
        with_group$enrichment[[group]] <- stats::setNames(list(identity, identity), members[[group]])
        expect_identical(validateOmicGuiConfig(with_group), with_group, info = group)
        for (member in members[[group]]) {
            bad <- with_group
            bad$enrichment[[group]][[member]] <- NULL
            expectGuiConfigRefused(bad, paste0("enrichment$", group, "$", member))
            bad <- with_group
            bad$enrichment[[group]][[member]] <- "identity"
            expectGuiConfigRefused(bad, paste0("enrichment$", group, "$", member))
        }
        bad <- prot
        bad$enrichment[[group]] <- "identity"
        expectGuiConfigRefused(bad, paste0("enrichment$", group))
    }
    # A non-proteomics omic may carry the group; NULL is valid, a malformed group is refused.
    lipid <- omicGuiConfig("lipidomics")
    bad <- lipid
    bad$enrichment <- list(module_id = "enrichment")
    expectGuiConfigRefused(bad, "enrichment$status_key")
    bad <- lipid
    bad$enrichment <- "enrichment"
    expectGuiConfigRefused(bad, "enrichment")
})
