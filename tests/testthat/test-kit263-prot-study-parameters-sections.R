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
