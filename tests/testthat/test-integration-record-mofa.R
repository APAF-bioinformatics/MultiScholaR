# KIT-142 end to end with MOFA2: import, sample matching with an uploaded match table, training, results saved, the
# record and the report saved through the Record tab, the project closed, reopened and the analysis repeated from its
# record. The tier-1 corpus is copied into a project with its transcriptomics samples renamed (`S-<name>`), as in
# test-integration-results-mofa.R, so the uploaded table is part of what is replayed. Skipped where MOFA2 or a Python
# with mofapy2 is absent (the package's own CI); run on the MOFA2 lane with the bundle's Python. Declared before the run:
#   - the repeated model equals the recorded one within 1e-5 relative (INT-6; on one computer it is expected to be
#     identical, and the largest difference is reported);
#   - the record names MOFA2 1.22.0 and mofapy2 0.7.3 (INT-5), double precision, the uploaded match and the six data
#     sets, and the report (Word, rendered by the package's template) holds the methods paragraph with those values;
#   - MULTISCHOLAR_INTEGRATION_EVIDENCE_OUT set: the methods paragraph, the study parameters, the report, the record and
#     the comparison are written there for the hand-off.

intRecMofaData <- function() normalizePath(testthat::test_path("..", "testdata", "integration"), mustWork = FALSE)
intRecMofaOmics <- c("lipidomics", "metabolomics", "proteomics", "transcriptomics")

intRecMofaSkip <- function() {
    testthat::skip_if_not(dir.exists(file.path(intRecMofaData(), "corpus")), "the tier-1 corpus is not present")
    testthat::skip_if_not(workflowOwnerCapable(), "the integration store needs qs2, DuckDB and Arrow")
    capability <- integrationMofaCapability()
    testthat::skip_if_not(isTRUE(capability$available), capability$reason)
}

intRecMofaProject <- function(.local_envir = parent.frame()) {
    root <- withr::local_tempdir(.local_envir = .local_envir)
    results_dir <- file.path(root, "results", "integration")
    paths <- list(base_dir = root, integration_dir = file.path(root, "integration"),
        source_dir = file.path(root, "scripts", "integration"), results_dir = results_dir,
        multiomic_inputs_dir = file.path(results_dir, "multiomic_inputs"), mofa_tables_dir = file.path(results_dir, "mofa", "tables"),
        mofa_plots_dir = file.path(results_dir, "mofa", "plots"),
        study_report_dir = file.path(root, "results_summary", "integration", "Study_report"))
    dir.create(paths$source_dir, recursive = TRUE)
    for (omic in intRecMofaOmics) {
        dir <- file.path(paths$integration_dir, omic)
        dir.create(dir, recursive = TRUE)
        stopifnot(all(file.copy(list.files(file.path(intRecMofaData(), "corpus", omic), full.names = TRUE), dir)))
    }
    dir <- file.path(paths$integration_dir, "transcriptomics")
    input <- readIntegrationInputWithMatrices(dir)
    read <- input$manifest
    manifest <- list(omic_type = read$omic_type, workflow = read$workflow, matrix_role = read$matrix_role, scale = read$scale,
        correction = NULL, source_state = list(name = "kit142-renamed"), sample_id_column = read$sample_id_column,
        group_column = read$group_column, technical_replicate_column = read$technical_replicate_column,
        biological_sample_column = read$biological_sample_column)
    samples <- input$samples
    samples[[manifest$biological_sample_column]] <- paste0("S-", samples[[manifest$biological_sample_column]])
    writeIntegrationInput(dir, manifest, samples, lapply(input$views, `[[`, "matrix"), lapply(input$views, `[[`, "features"))
    paths
}

test_that("KIT-142: with MOFA2, a saved integration is reopened, repeated from its record to the same model, and reported", {
    intRecMofaSkip()
    paths <- intRecMofaProject()
    evidence <- Sys.getenv("MULTISCHOLAR_INTEGRATION_EVIDENCE_OUT", unset = "")
    bus <- do.call(shiny::reactiveValues, omicWorkflowInitialBus(integrationGuiConfig()$shell))
    withr::defer(if (!is.null(shiny::isolate(bus$state_manager))) closeWorkflowStateOwner(shiny::isolate(bus$state_manager)))
    shiny::testServer(mod_integration_import_server, args = list(workflow_data = bus, paths = paths), {
        session$flushReact()
        session$setInputs(views = bus$integration_scan$views$view, import = 1)
    })
    manager <- shiny::isolate(bus$state_manager)
    edited <- integrationSampleMatchTemplate(integrationImportedObject(manager))
    rows <- edited$omic == "transcriptomics"
    edited$biological_sample[rows] <- sub("^S-", "", edited$saved_biological_sample[rows])
    upload <- file.path(paths$base_dir, "edited.tsv")
    utils::write.table(edited, upload, sep = "\t", quote = FALSE, row.names = FALSE)
    shiny::testServer(mod_integration_alignment_server, args = list(workflow_data = bus, paths = paths), {
        session$setInputs(match_upload = data.frame(name = "edited.tsv", size = file.size(upload), type = "text/tab-separated-values",
            datapath = upload, stringsAsFactors = FALSE))
        session$setInputs(confirm = 1)
    })
    shiny::testServer(mod_integration_model_server, args = list(workflow_data = bus, paths = paths,
        withProgressFn = function(message, value, expr) force(expr)), {
        session$setInputs(run = 1)
        expect_identical(bus$integration_model$status, "computed", info = output$run_status$html)
    })
    shiny::testServer(mod_integration_results_server, args = list(workflow_data = bus, paths = paths), {
        session$flushReact()
        session$setInputs(save = 1)
        expect_false(is.null(bus$integration_results), info = output$save_status$html)
    })
    shiny::testServer(mod_integration_record_server, args = list(workflow_data = bus, paths = paths,
        withProgressFn = function(message, value, expr) force(expr)), {
        session$setInputs(save = 1)
        expect_false(is.null(bus$integration_record), info = output$record_status$html)
    })
    saved <- shiny::isolate(bus$integration_record)
    expect_identical(basename(saved$files), c("analysis_record.json", "analysis_record.xlsx", "study_parameters.txt",
        "integration_report.docx", "model.hdf5", "sample_match.tsv"))
    expect_true(all(file.size(saved$files) > 0))
    record <- jsonlite::read_json(saved$files[[1L]], simplifyVector = FALSE)
    block <- record$analyses$integration_mofa_model$integration
    # KIT-148: the model file MOFA2 wrote, its checksum the record's, without the training data (no `data` group); MOFA2
    # loads it; the uploaded match as the runs it moved.
    expect_identical(digest::digest(file = saved$files[[5L]], algo = "sha256"), block$model$model_file_sha256[[1L]])
    groups <- rhdf5::h5ls(saved$files[[5L]], recursive = FALSE)$name
    message(sprintf("KIT-148 model.hdf5: %d bytes; top-level groups %s", as.integer(file.size(saved$files[[5L]])),
        paste(groups, collapse = ", ")))
    expect_false("data" %in% groups)
    expect_true(all(c("expectations", "variance_explained") %in% groups))
    loaded <- suppressWarnings(MOFA2::load_model(saved$files[[5L]]))
    expect_identical(ncol(MOFA2::get_factors(loaded)[[1L]]), 5L)
    # KIT-239 (L8, a guard): the Limitations' fixed sentence holds for the model MOFA2 trained: every likelihood Gaussian
    # (the record's MOFA2 options, decoded) and one group of samples (`get_factors` gives one matrix per group).
    opts <- decodeWorkflowRecordValue(block$model$mofa2_options)
    message(sprintf("KIT-239 L8: likelihoods %s; groups %d", paste(unlist(opts$model$likelihoods), collapse = ", "),
        length(MOFA2::get_factors(loaded))))
    expect_true(length(opts$model$likelihoods) > 0L && all(unlist(opts$model$likelihoods) == "gaussian"))
    expect_identical(length(MOFA2::get_factors(loaded)), 1L)
    expect_identical(nrow(utils::read.delim(saved$files[[6L]], colClasses = "character")), sum(rows))
    expect_identical(block$model$versions$MOFA2, "1.22.0")
    expect_identical(block$model$versions$mofapy2, "0.7.3")
    expect_identical(block$model$precision, "double")
    expect_identical(block$matching$rule, "user_table")
    methods <- saved$methods
    message("KIT-142 methods paragraph (draft, from the record):\n", methods)
    for (fragment in c("MOFA2 1.22.0 (mofapy2 0.7.3, Python 3.12.16", "in double precision", "6 data sets from 4 omics",
        sprintf("assigned %d runs to another biological sample", sum(rows)))) {
        expect_true(grepl(fragment, methods, fixed = TRUE), info = fragment)
    }
    # The Word report holds the methods paragraph with the recorded values.
    unpacked <- withr::local_tempdir()
    utils::unzip(saved$files[[4L]], files = "word/document.xml", exdir = unpacked)
    document <- paste(readLines(file.path(unpacked, "word", "document.xml"), warn = FALSE), collapse = "")
    words <- gsub("<[^>]+>", "", document)
    for (fragment in c(integrationRecordTexts()[["report_title"]], "MOFA2 1.22.0", "mofapy2 0.7.3", "double precision",
        sprintf("Training settled after %d rounds", as.integer(block$model$iterations[[1L]])))) {
        expect_true(grepl(fragment, words, fixed = TRUE), info = fragment)
    }
    # Closed, reopened, repeated from the record.
    root <- dirname(manager$.context$store$directory)
    recorded_sha <- block$model$model_file_sha256
    closeWorkflowStateOwner(manager)
    started <- Sys.time()
    replay <- replayIntegrationProject(root, file.path(paths$base_dir, "replay"))
    seconds <- as.numeric(difftime(Sys.time(), started, units = "secs"))
    table <- replay$comparison$tables
    message(sprintf("KIT-142 replay from the record: %.1f s; largest relative difference %.3g (identical: %s); per table: %s",
        seconds, replay$comparison$largest_relative, replay$comparison$identical,
        paste(sprintf("%s %.3g", table$table, table$largest_relative), collapse = "; ")))
    message(sprintf("KIT-142 model file sha256: recorded %s, repeated %s", recorded_sha, replay$replayed$model_file$sha256))
    expect_true(replay$comparison$equal)
    expect_lte(replay$comparison$largest_relative, 1e-5)
    expect_true(all(replay$versions$same))
    expect_identical(replay$parameters$sample_match, block$matching$sample_match_sha256[[1L]])
    if (nzchar(evidence)) {
        dir.create(evidence, recursive = TRUE, showWarnings = FALSE)
        writeLines(methods, file.path(evidence, "methods-paragraph.txt"))
        file.copy(saved$files, evidence, overwrite = TRUE)
        utils::write.table(table, file.path(evidence, "replay-comparison.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
        utils::write.table(replay$versions, file.path(evidence, "replay-versions.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
        writeLines(sprintf("replay seconds %.1f; largest relative %.17g; identical %s; model sha256 recorded %s repeated %s",
            seconds, replay$comparison$largest_relative, replay$comparison$identical, recorded_sha,
            replay$replayed$model_file$sha256), file.path(evidence, "replay-summary.txt"))
    }
})
