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
        "success_text", "preview", "process")
)

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

test_that("each omic's configuration validates, with the shell, binding and import groups", {
    for (omic in guiConfigOmics) {
        config <- omicGuiConfig(omic)
        expect_identical(names(config), c("shell", "binding", "import"), info = omic)
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
