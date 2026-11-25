# KIT-263: a proteomics study_parameters.txt prints only proteomics' configuration sections. The writer
# `createWorkflowArgsFromConfig` takes `omic_type`; for "proteomics" it leaves the four sections of config.ini's
# Metabolomics and Transcriptomics blocks out of both printed blocks (VAL-279). Other omics, and a call without
# `omic_type`, print every section as before. The writer is called through `do.call` with `omic_type` added only when
# the function has it, so the same file runs on the base and on the head and fails on base by its asserted value.

kit263Config <- function() {
    readConfigFile(system.file("config", "config.ini", package = "MultiScholaR"), create_workers = FALSE)
}

#' Three proteins in four runs and two groups, with the shipped configuration as `@args`, as the design import
#' builds a proteomics object (R/utils_workflow_dia_state.R:1507-1516).
kit263Protein <- function(config) {
    ProteinQuantitativeData(
        protein_quant_table = data.frame(Protein.Ids = c("P1", "P2", "P3"), S1 = c(10, 11, 12), S2 = c(11, 12, 13),
            S3 = c(20, 21, 22), S4 = c(21, 22, 23), stringsAsFactors = FALSE),
        design_matrix = data.frame(Run = c("S1", "S2", "S3", "S4"), group = c("G1", "G1", "G2", "G2"),
            stringsAsFactors = FALSE),
        sample_id = "Run", group_id = "group", protein_id_column = "Protein.Ids", args = config
    )
}

#' The study parameters file's lines for `omic_type` (`NULL`: the argument is not passed).
kit263Lines <- function(omic_type = NULL, config = kit263Config()) {
    dir <- withr::local_tempdir(.local_envir = parent.frame())
    args <- list(workflow_name = "kit263", description = "d", source_dir_path = dir,
        final_s4_object = kit263Protein(config), workflow_data = list(config_list = config))
    if (!is.null(omic_type) && "omic_type" %in% names(formals(createWorkflowArgsFromConfig))) args$omic_type <- omic_type
    utils::capture.output(file <- do.call(createWorkflowArgsFromConfig, args))
    readLines(file.path(dir, "study_parameters.txt"), warn = FALSE)
}

kit263Foreign <- function(lines) grep("metabolit|transcript", lines, ignore.case = TRUE, perl = TRUE, value = TRUE)

test_that("(KIT-263 T1) a proteomics study_parameters.txt lists no metabolomics or transcriptomics setting", {
    lines <- kit263Lines("proteomics")
    foreign <- kit263Foreign(lines)
    cat(sprintf("\nKIT-263 T1: %d foreign lines\n%s\n", length(foreign), paste(foreign, collapse = "\n")))
    expect_identical(length(foreign), 0L)
})

test_that("(KIT-263 T2) other omics and a call without omic_type print every section, as before", {
    config <- kit263Config()
    held <- intersect(names(config), c("metaboliteIntensityFiltering", "removeMetabolitesWithOnlyOneReplicate",
        "transcriptIntensityFiltering", "removeTranscriptsWithOnlyOneReplicate"))
    # Two printed blocks, each with one header and one line per key, for each non-empty section of the four.
    expected <- sum(vapply(config[held], function(section) 2L * (1L + length(section)), 0L))
    expect_gt(expected, 0L)
    expect_identical(length(kit263Foreign(kit263Lines("metabolomics"))), expected)
    expect_identical(length(kit263Foreign(kit263Lines(NULL))), expected)
})

test_that("(KIT-263 T3) proteomics' own sections still print in both blocks", {
    lines <- kit263Lines("proteomics")
    start <- grep("Parameters from Final S4 Object:", lines, fixed = TRUE)[1L]
    additional <- grep("Additional Configuration Parameters:", lines, fixed = TRUE)[1L]
    expect_false(is.na(start))
    expect_false(is.na(additional))
    expect_true(all(c("[removeRowsWithMissingValuesPercent]", "[deAnalysisParameters]") %in% trimws(lines[start:additional])))
    expect_true("removeRowsWithMissingValuesPercent:" %in% trimws(lines[additional:length(lines)]))
})

test_that("(KIT-263 T4) the deny-list names the sections under config.ini's Metabolomics and Transcriptomics banners", {
    ini <- readLines(system.file("config", "config.ini", package = "MultiScholaR"), warn = FALSE)
    first <- grep("^; --- Metabolomics ---", ini)
    expect_length(first, 1L)
    banners <- grep("^; ====", ini)
    last <- min(banners[banners > first]) - 1L
    headers <- sub("^\\[(.*)\\]\\s*$", "\\1", grep("^\\[.*\\]\\s*$", ini[first:last], value = TRUE))
    expect_identical(headers, studyParameterForeignSections$proteomics,
        info = "update studyParameterForeignSections in R/func_general_filemgmt_study_parameter_writer.R")
})

# The Summary save passes the omic to the writer (the transcriptomics case of test-omics-summary-core.R:593-613; the
# file's helpers are copied, never sourced).
kit263Dirs <- function(omic, env = parent.frame()) {
    base <- withr::local_tempdir(.local_envir = env)
    paths <- list(base_dir = base, source_dir = file.path(base, "scripts", omic), integration_dir = file.path(base, "integration"),
        results_summary_dir = file.path(base, "results_summary"))
    for (dir in paths) dir.create(dir, recursive = TRUE, showWarnings = FALSE)
    stats::setNames(list(paths), omic)
}

test_that("(KIT-263 T5) saveWorkflowSummaryParameters passes the omic type to the writer", {
    config <- omicGuiConfig("transcriptomics")$summary
    dirs <- kit263Dirs("transcriptomics")
    paths <- dirs$transcriptomics
    config$final_state$context <- function(workflowData, catFn) list(finalS4Object = NULL, contrastsTbl = NULL,
        dataStateUsed = "trx_ruvs")
    config$integration$write <- function(workflowData, path) writeLines("integration", path)
    seen <- NULL
    args <- function(...) {
        seen <<- list(...)
        file.path(paths$source_dir, "study_parameters.txt")
    }
    manager <- new.env()
    manager$.context <- list(snapshot = list(exports = list(trx_ruviii_export = list(reference = "r"))))
    saveWorkflowSummaryParameters(config, list(state_manager = manager), dirs, "transcriptomics", "kit263", "d",
        notify = function(message, type, duration) NULL, done = function(value) NULL, failed = function(error) stop(error),
        createWorkflowArgsFn = args, writeAnalysisRecordFn = function(manager, path) invisible(path))
    expect_identical(seen$omic_type, "transcriptomics")
})

# ------------------------------------------------------------------------------------------------
# KIT-263 declared diagnostic (coordinator's ruling after R1): printing only, no expectation. (A) the T1 fixture's whole
# written text; (B) the text the application's own path writes for a proteomics project: `saveWorkflowSummaryParameters`
# (the Summary tab's core, default writer) on the module-ci proteomics bus (`helper-module-ci-prot-summary.R`'s bus, as `omicSummaryW1ProtSetup` builds it), B1 as that fixture is, B2 with the shipped
# configuration as the bus's config_list and the state object's @args, as the design import makes them.
# ------------------------------------------------------------------------------------------------

kit263PrintFile <- function(label, lines) {
    cat(sprintf("\n==== KIT-263 DIAG %s: %d lines ====\n", label, length(lines)))
    cat(sprintf("%3d| %s", seq_along(lines), lines), sep = "\n")
    cat(sprintf("==== end %s ====\n", label))
}

kit263PrintSections <- function(label, args, config_list) {
    cat(sprintf("\n---- KIT-263 DIAG %s: objects handed to the writer ----\n", label))
    cat("class(@args): ", paste(class(args), collapse = ","), "; length ", length(args), "\n")
    cat("names(@args): ", paste(names(args), collapse = " "), "\n")
    cat("class(config_list): ", paste(class(config_list), collapse = ","), "; length ", length(config_list), "\n")
    cat("names(config_list): ", paste(names(config_list), collapse = " "), "\n")
    cat("S4 block iterates over (names(@args)): ", paste(names(args), collapse = " "), "\n")
    cat("Additional block iterates over (names(config_list) then the @args names): ",
        paste(union(names(config_list), names(args)), collapse = " "), "\n")
    direct <- tryCatch(formatConfigList(config_list), error = function(e) paste("formatConfigList error:", conditionMessage(e)))
    cat("formatConfigList(config_list) alone returns ", length(direct), " elements; first lines:\n")
    cat(utils::head(unlist(strsplit(direct, "\n")), 12), sep = "\n")
}

test_that("(KIT-263 DIAG A) the T1 fixture's whole written file", {
    config <- kit263Config()
    kit263PrintSections("A", kit263Protein(config)@args, config)
    kit263PrintFile("A", kit263Lines("proteomics", config))
    succeed()
})

test_that("(KIT-263 DIAG B) the Summary save's own file for a proteomics project", {
    config <- omicGuiConfig("proteomics")$summary
    shipped <- kit263Config()
    for (variant in c("B1", "B2")) {
        # The construction of `omicSummaryW1ProtSetup` (test-omic-summary-module.R:812-858), without its tab server.
        bus <- module_ci_prot_summary_workflow(module_ci_prot_summary_branch_config("DIA"))
        project_dirs <- module_ci_prot_summary_paths()
        withr::defer(unlink(project_dirs$proteomics$base_dir, recursive = TRUE))
        if (identical(variant, "B2")) {
            object <- module_ci_prot_summary_object(module_ci_prot_summary_branch_config("DIA"))
            object@args <- shipped
            bus$state_manager <- module_ci_prot_summary_state_manager(object)
            bus$config_list <- shipped
        }
        final <- tryCatch(workflowSummaryFinalState(config, bus, "proteomics", function(...) NULL)$finalS4Object,
            error = function(e) NULL)
        kit263PrintSections(variant, if (isS4(final)) final@args else NULL, bus$config_list)
        err <- NULL
        saveWorkflowSummaryParameters(config, bus, project_dirs, "proteomics", "kit263", "d",
            notify = function(message, type, duration) cat("notify:", type, message, "\n"), done = function(value) NULL,
            failed = function(error) err <<- error, writeAnalysisRecordFn = function(manager, path) invisible(path),
            catFn = function(...) NULL)
        if (!is.null(err)) cat("save failed: ", conditionMessage(err), "\n")
        file <- file.path(project_dirs$proteomics$source_dir, "study_parameters.txt")
        kit263PrintFile(variant, if (file.exists(file)) readLines(file, warn = FALSE) else "<no study_parameters.txt>")
    }
    succeed()
})

# ------------------------------------------------------------------------------------------------
# KIT-263 declared diagnostic D2 (the coordinator's ruling after D1): printing only, no expectation. The application's own
# proteomics save with an OWNER (the `omicSummaryW1ProtSetup(owner = TRUE)` construction, test-omic-summary-module.R:
# 826-841, without its tab server) whose protein object carries the shipped configuration as @args and the bus the same
# as config_list, as the design import makes them; and a lipidomics owner's save (`module_ci_lipid_summary_workflow`),
# section headers only.
# ------------------------------------------------------------------------------------------------

test_that("(KIT-263 DIAG D2-prot) the Summary save's file for a proteomics owner", {
    config <- omicGuiConfig("proteomics")$summary
    shipped <- kit263Config()
    cfg <- module_ci_prot_summary_branch_config("DIA")
    bus <- module_ci_prot_summary_workflow(cfg)
    project_dirs <- module_ci_prot_summary_paths()
    paths <- project_dirs$proteomics
    withr::defer(unlink(paths$base_dir, recursive = TRUE))
    design <- module_ci_prot_da_design()
    rownames(design) <- NULL
    quant <- module_ci_prot_da_table()
    rownames(quant) <- NULL
    object <- module_ci_prot_da_object(workflow_type = cfg$workflow_type, use_limpa = cfg$use_limpa,
        report_template = cfg$report_template, design = design, protein_quant_table = quant, args = shipped)
    store_root <- file.path(paths$base_dir, "store")
    dir.create(store_root, showWarnings = FALSE)
    manager <- newDIAStateAdapter(newWorkflowEntityStore(store_root, "kit263-prot", cfg$workflow_type),
        WorkflowState$new(audit_enabled = FALSE))
    withr::defer(closeDIAStateAdapter(manager))
    manager$.context$reference_states <- TRUE
    manager$.context$experiment_paths <- paths
    manager$saveState("correlation_filtered", object, list(), "kit263 fixture state")
    bus$state_manager <- manager
    bus$config_list <- shipped
    final <- tryCatch(workflowSummaryFinalState(config, bus, "proteomics", function(...) NULL)$finalS4Object,
        error = function(e) e)
    cat("\n---- KIT-263 DIAG D2-prot: the final object handed to the writer ----\n")
    cat("class(final): ", paste(class(final), collapse = ","), "; isSavedWorkflowReference: ",
        tryCatch(isSavedWorkflowReference(final), error = function(e) NA), "\n")
    args <- tryCatch(readDIAStateMetadata(manager, final)$data$slots$args, error = function(e) NULL)
    if (is.null(args) && isS4(final)) args <- final@args
    kit263PrintSections("D2-prot", args, bus$config_list)
    err <- NULL
    saveWorkflowSummaryParameters(config, bus, project_dirs, "proteomics", "kit263", "d",
        notify = function(message, type, duration) cat("notify:", type, message, "\n"), done = function(value) NULL,
        failed = function(error) err <<- error, writeAnalysisRecordFn = function(manager, path) invisible(path),
        catFn = function(...) NULL)
    if (!is.null(err)) cat("save failed: ", conditionMessage(err), "\n")
    file <- file.path(paths$source_dir, "study_parameters.txt")
    kit263PrintFile("D2-prot", if (file.exists(file)) readLines(file, warn = FALSE) else "<no study_parameters.txt>")
    succeed()
})

test_that("(KIT-263 DIAG D2-lipid) the Summary save's section headers for a lipidomics owner", {
    config <- omicGuiConfig("lipidomics")$summary
    bus <- module_ci_lipid_summary_workflow(module_ci_lipid_summary_config("pos_neg"))
    project_dirs <- module_ci_lipid_summary_paths()
    paths <- project_dirs$lipidomics
    withr::defer(unlink(paths$base_dir, recursive = TRUE))
    err <- NULL
    saveWorkflowSummaryParameters(config, bus, project_dirs, "lipidomics", "kit263", "d",
        notify = function(message, type, duration) cat("notify:", type, message, "\n"), done = function(value) NULL,
        failed = function(error) err <<- error, writeAnalysisRecordFn = function(manager, path) invisible(path),
        catFn = function(...) NULL)
    if (!is.null(err)) cat("save failed: ", conditionMessage(err), "\n")
    file <- file.path(paths$source_dir, "study_parameters.txt")
    lines <- if (file.exists(file)) readLines(file, warn = FALSE) else "<no study_parameters.txt>"
    cat("\n==== KIT-263 DIAG D2-lipid: ", length(lines), " lines; headers (bracketed or ending in a colon at column 1) ====\n")
    cat(grep("^\\[|^[A-Za-z][^ ]*:$|^[A-Z][A-Za-z ]*:$", lines, value = TRUE), sep = "\n")
    cat("foreign lines: ", length(kit263Foreign(lines)), "\n==== end D2-lipid ====\n")
    succeed()
})
