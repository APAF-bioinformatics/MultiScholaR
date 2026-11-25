# KIT-263: a proteomics study_parameters.txt prints only proteomics' configuration sections. The writer
# `createWorkflowArgsFromConfig` takes `omic_type`; for "proteomics" it leaves the four sections of config.ini's
# Metabolomics and Transcriptomics blocks out of both printed blocks (VAL-279). Other omics, and a call without
# `omic_type`, print every section as before. The writer is called through `do.call` with `omic_type` added only when
# the function has it, so the same file runs on the base and on the head and fails on base by its asserted value.

kit263Config <- function() {
    readConfigFile(system.file("config", "config.ini", package = "MultiScholaR"), create_workers = FALSE)
}

#' The application's own construction of a proteomics owner (`omicSummaryW1ProtSetup(owner = TRUE)`,
#' test-omic-summary-module.R:826-841, without its tab server; KIT-263 D2, run 45681727116): the design state saved in a
#' DIA state owner, whose protein object carries the shipped configuration as `@args`, as the design import makes them
#' (R/utils_workflow_dia_state.R:1507-1516), and the bus (`workflow_data`) the same as `config_list`. The writer gets the
#' saved reference of the owner's head, as the Summary save hands it.
kit263Owner <- function(config, env = parent.frame()) {
    cfg <- module_ci_prot_summary_branch_config("DIA")
    bus <- module_ci_prot_summary_workflow(cfg)
    project_dirs <- module_ci_prot_summary_paths()
    paths <- project_dirs$proteomics
    withr::defer(unlink(paths$base_dir, recursive = TRUE), envir = env)
    design <- module_ci_prot_da_design()
    rownames(design) <- NULL
    quant <- module_ci_prot_da_table()
    rownames(quant) <- NULL
    object <- module_ci_prot_da_object(workflow_type = cfg$workflow_type, use_limpa = cfg$use_limpa,
        report_template = cfg$report_template, design = design, protein_quant_table = quant, args = config)
    store_root <- file.path(paths$base_dir, "store")
    dir.create(store_root, showWarnings = FALSE)
    manager <- newDIAStateAdapter(newWorkflowEntityStore(store_root, "kit263-prot", cfg$workflow_type),
        WorkflowState$new(audit_enabled = FALSE))
    withr::defer(closeDIAStateAdapter(manager), envir = env)
    manager$.context$reference_states <- TRUE
    manager$.context$experiment_paths <- paths
    manager$saveState("correlation_filtered", object, list(), "kit263 fixture state")
    bus$state_manager <- manager
    bus$config_list <- config
    final <- workflowSummaryFinalState(omicGuiConfig("proteomics")$summary, bus, "proteomics", function(...) NULL)$finalS4Object
    list(final = final, bus = bus, dir = paths$source_dir)
}

#' The study parameters file's lines for `omic_type` (`NULL`: the argument is not passed).
kit263Lines <- function(omic_type = NULL, config = kit263Config()) {
    owner <- kit263Owner(config)
    args <- list(workflow_name = "kit263", description = "d", source_dir_path = owner$dir,
        final_s4_object = owner$final, workflow_data = owner$bus)
    if (!is.null(omic_type) && "omic_type" %in% names(formals(createWorkflowArgsFromConfig))) args$omic_type <- omic_type
    utils::capture.output(do.call(createWorkflowArgsFromConfig, args))
    readLines(file.path(owner$dir, "study_parameters.txt"), warn = FALSE)
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

