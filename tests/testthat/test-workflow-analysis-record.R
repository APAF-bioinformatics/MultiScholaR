# KIT-027: the analysis record derived from head metadata.

analysisRecordFixture <- function(root) {
    manager <- newDIAStateAdapter(newWorkflowEntityStore(root, "record-test", "DIA"), WorkflowState$new(audit_enabled = FALSE))
    manager$setWorkflowType("DIA")
    manager$.context$reference_states <- TRUE
    peptide <- makeWorkflowCheckpointPeptide(readWorkflowCheckpointRecipe("small"))
    peptide@args$filterMinNumPeptidesPerSample <- list(peptides_per_sample_cutoff = 500,
        inclusion_list = paste0("sample_", seq_len(40)))
    manager$saveState("raw_data_s4", peptide, list(), "imported")
    filtered <- peptide
    filtered@args$srlQvalueProteotypicPeptideClean <- list(qvalue_threshold = 0.05, global_qvalue_threshold = 0.01,
        global_pg_qvalue_threshold = 0.01, confidence_provenance_mode = "diann_main_report")
    entry <- workflowOperationEntry("peptide_qvalue_filter")
    setDIAStepProvenance(manager, "peptide_qvalue_filter", "whole", entry$parameters)
    manager$saveState("qvalue_filtered", filtered, list(proteotypic_only = TRUE), "q-value filter")
    fasta <- file.path(root, "proteins.fasta")
    writeLines(c(">sp|P1|ONE", "MPEPTIDE"), fasta)
    recordDIAWorkflowInputs(manager, c(fasta = fasta))
    protein <- makeWorkflowCheckpointProtein(readWorkflowCheckpointRecipe("small"))
    protein@args$chooseBestProteinAccession <- list(delim = ";", aggregation_method = "mean")
    entry <- workflowOperationEntry("protein_accession_cleanup")
    setDIAStepProvenance(manager, "protein_accession_cleanup", "whole", entry$parameters)
    manager$saveState("protein_accession_cleaned", protein, list(delimiter = ";", aggregation_method = "mean",
        fasta_sha256 = digest::digest(file = fasta, algo = "sha256")), "accession cleanup")
    recordWorkflowExternalResource(manager, "string-db", "12.0", "api_version", query_summary = "fixture")
    list(manager = manager, peptide = peptide, fasta = fasta)
}

test_that("the record has one step per state with declared parameters from args or config", {
    fixture <- analysisRecordFixture(withr::local_tempdir())
    manager <- fixture$manager
    withr::defer(closeDIAStateAdapter(manager))
    reads <- manager$.context$full_s4_reads
    record <- buildWorkflowAnalysisRecord(manager)
    expect_identical(manager$.context$full_s4_reads, reads)
    check <- validateWorkflowAnalysisRecord(record)
    expect_true(check$valid)
    expect_false(any(grepl("declared parameter not recorded", check$incomplete)))
    history <- unlist(manager$.context$snapshot$state_history)
    expect_length(record$steps, length(history))
    expect_identical(vapply(record$steps, `[[`, "", "state"), history)
    for (step in record$steps) {
        expect_identical(step$operation, manager$.context$snapshot$step_provenance[[step$state]]$operation)
        expect_true(nzchar(step$step_label))
    }
    by_state <- stats::setNames(record$steps, history)
    parameters <- function(state) stats::setNames(by_state[[state]]$parameters,
        vapply(by_state[[state]]$parameters, `[[`, "", "name"))
    qvalue <- parameters("qvalue_filtered")
    expect_setequal(names(qvalue), workflowOperationEntry("peptide_qvalue_filter")$parameters)
    expect_identical(qvalue$qvalue_threshold$provenance, "user_set")
    expect_identical(qvalue$qvalue_threshold$source, "args")
    expect_identical(qvalue$global_qvalue_threshold$provenance, "default")
    expect_identical(qvalue$proteotypic_only$source, "config")
    expect_identical(qvalue$proteotypic_only$provenance, "default")
    accession <- parameters("protein_accession_cleaned")
    expect_identical(accession$fasta_sha256$source, "config")
    expect_identical(accession$fasta_sha256$provenance, "input")
    fasta_input <- Filter(function(input) identical(input$role, "fasta"), record$inputs)[[1]]
    expect_identical(decodeWorkflowRecordValue(accession$fasta_sha256$value), fasta_input$sha256)
    expect_identical(accession$delimiter$provenance, "default")
    expect_identical(record$external_resources[[1]]$release, "12.0")
    expect_true(grepl("worker", record$limits$external_resources))
    expect_true(any(grepl("controller-process", check$incomplete)))
})

test_that("every atomic argument survives the JSON round trip at full length and the record is stable", {
    fixture <- analysisRecordFixture(withr::local_tempdir())
    manager <- fixture$manager
    withr::defer(closeDIAStateAdapter(manager))
    path <- withr::local_tempfile(fileext = ".json")
    writeWorkflowAnalysisRecord(manager, path, now = as.POSIXct("2026-10-17 10:00:00", tz = "UTC"))
    parsed <- jsonlite::read_json(path, simplifyVector = FALSE)
    expect_true(validateWorkflowAnalysisRecord(parsed)$valid)
    legacy <- Filter(function(step) identical(step$state, "raw_data_s4"), parsed$steps)[[1]]
    args <- fixture$peptide@args
    compare <- function(original, encoded, path) {
        if (is.data.frame(original) || is.null(original)) return(invisible())
        if (is.list(original)) {
            for (name in names(original)) compare(original[[name]], encoded$value[[match(name, unlist(encoded$names))]], c(path, name))
            return(invisible())
        }
        expect_identical(decodeWorkflowRecordValue(encoded), original, label = paste(path, collapse = "/"))
    }
    for (parameter in legacy$parameters) compare(args[[parameter$name]], parameter$value, parameter$name)
    inclusion <- Filter(function(p) identical(p$name, "filterMinNumPeptidesPerSample"), legacy$parameters)[[1]]
    expect_length(decodeWorkflowRecordValue(inclusion$value)$inclusion_list, 40L)
    again <- withr::local_tempfile(fileext = ".json")
    writeWorkflowAnalysisRecord(manager, again, now = as.POSIXct("2026-10-17 10:00:00", tz = "UTC"))
    expect_identical(readLines(path), readLines(again))
})

test_that("a store without provenance or inputs gives legacy steps and reports them as incomplete", {
    manager <- newDIAStateAdapter(newWorkflowEntityStore(withr::local_tempdir(), "record-legacy", "DIA"), WorkflowState$new(audit_enabled = FALSE))
    withr::defer(closeDIAStateAdapter(manager))
    manager$.context$reference_states <- TRUE
    manager$saveState("raw_data_s4", makeWorkflowCheckpointPeptide(readWorkflowCheckpointRecipe("small")), list(), "imported")
    record <- buildWorkflowAnalysisRecord(manager)
    check <- validateWorkflowAnalysisRecord(record)
    expect_true(check$valid)
    expect_true(any(grepl("legacy step: raw_data_s4", check$incomplete)))
    expect_true("no inputs recorded" %in% check$incomplete)
})

test_that("KIT-101 F-15: a completed enrichment state on the head is recorded and the earlier steps and design are unchanged", {
    fixture <- analysisRecordFixture(withr::local_tempdir())
    manager <- fixture$manager
    withr::defer(closeDIAStateAdapter(manager))
    now <- as.POSIXct("2026-11-07 10:00:00", tz = "UTC")
    before <- buildWorkflowAnalysisRecord(manager, now)
    # As the proteomics enrichment module commits it: the declared result class, checkpointed on the adapter.
    results <- methods::new("EnrichmentResults", contrasts = tibble::tibble(contrasts = "groupB-groupA"),
        enrichment_data = list(`groupB-groupA` = list(up = data.frame(term_id = "GO:0008150", p_value = 0.01))),
        enrichment_plots = list(), enrichment_plotly = list(), enrichment_summaries = list())
    manager$saveState("enrichment_completed", results, list(q_cutoff = 0.05), "fixture enrichment")
    expect_identical(manager$current_state, "enrichment_completed")
    reads <- manager$.context$full_s4_reads
    after <- buildWorkflowAnalysisRecord(manager, now)
    expect_identical(manager$.context$full_s4_reads, reads)
    expect_true(validateWorkflowAnalysisRecord(after)$valid)
    expect_identical(vapply(after$steps, `[[`, "", "state"), unlist(manager$.context$snapshot$state_history))
    expect_identical(after$steps[seq_along(before$steps)], before$steps)
    expect_identical(after$steps[[length(after$steps)]]$description, "fixture enrichment")
    expect_false(is.null(before$design))
    expect_identical(after$design, before$design)
    expect_identical(after[setdiff(names(after), "steps")], before[setdiff(names(before), "steps")])
})

test_that("the summary page writes analysis_record.json after study_parameters.txt for an on-disk owner", {
    calls <- character()
    record_args <- NULL
    manager <- structure(list(), class = "DIAStateAdapter")
    workflow <- list(state_manager = manager, contrasts_tbl = NULL, config_list = NULL)
    source_dir <- withr::local_tempdir()
    # KIT-108 W4: the core's save under proteomics' group, its final state replaced as before.
    config <- omicGuiConfig("proteomics")$summary
    config$final_state$context <- function(...) list(finalS4Object = NULL)
    saveWorkflowSummaryParameters(config, workflow, projectDirs = list(proteomics = list(source_dir = source_dir,
            integration_dir = withr::local_tempdir(), base_dir = source_dir)), omicType = "proteomics",
        label = "x", description = "d", notify = function(...) NULL, done = function(...) NULL, failed = function(...) NULL,
        createWorkflowArgsFn = function(...) { calls <<- c(calls, "study_parameters"); "study_parameters.txt" },
        writeAnalysisRecordFn = function(manager, path) { calls <<- c(calls, "analysis_record"); record_args <<- path },
        saveRDSFn = function(...) NULL, catFn = function(...) NULL)
    expect_identical(calls[1:2], c("study_parameters", "analysis_record"))
    expect_identical(record_args, file.path(source_dir, "analysis_record.json"))
})

test_that("KIT-027 audit PH-01/02/03: declared values come from provenance, both correlation operations and the DA settings are recorded", {
    fixture <- analysisRecordFixture(withr::local_tempdir())
    manager <- fixture$manager
    withr::defer(closeDIAStateAdapter(manager))
    protein <- makeWorkflowCheckpointProtein(readWorkflowCheckpointRecipe("small"))
    norm <- workflowOperationEntry("normalization-ruv")
    requested <- list(norm_method = "cyclicloess", ruv_mode = "automatic", auto_percentage_min = 1, auto_percentage_max = 20,
        separation_metric = "max_difference", k_penalty_weight = 0.5, adaptive_k_penalty = TRUE, ruv_percentage = 5,
        ruv_k = 5, grouping = "group")
    setDIAStepProvenance(manager, "normalization-ruv", "whole", norm$parameters, values = requested)
    # The commit's config records the chosen k as `ruv_k` (the optimisation picked 3).
    manager$saveState("ruv_corrected", protein, list(norm_method = "cyclicloess", ruv_mode = "automatic", ruv_k = 3,
        percentage_as_neg_ctrl = 5), "ruv")
    setDIAStepProvenance(manager, "pearson-pairs", "whole", workflowOperationEntry("pearson-pairs")$parameters,
        values = list(group = "group", tech_rep_remove_regex = "pool"))
    setDIAStepProvenance(manager, "pearson-selection", "whole", workflowOperationEntry("pearson-selection")$parameters,
        values = list(pairs = strrep("a", 64), threshold = 0.5, correlations = data.frame(x = 1)))
    manager$saveState("correlation_filtered", protein, list(min_pearson_correlation_threshold = 0.5), "correlation")
    manifest <- list(schema = "dia.da_results.v1", source = getDIAStateReference(manager, "correlation_filtered"),
        contrasts = data.frame(contrasts = "B-A", friendly_names = "B_vs_A", full_format = "B_vs_A=B-A"),
        results = list(), parameters = list(formula = "~ 0 + group", qvalue = 0.05, lfc = 0, weights = NA),
        model = list(group_id = "group", technical_replicate_id = NA_character_, block_reason = "no_duplicates"))
    snapshot <- manager$.context$snapshot
    snapshot$da_results <- saveDIAStateValue(manager$.context, manifest, "da-manifest", codec = "rds")
    publishDIAStateHead(manager$.context, snapshot)
    record <- buildWorkflowAnalysisRecord(manager)
    check <- validateWorkflowAnalysisRecord(record)
    expect_true(check$valid)
    expect_false(any(grepl("declared parameter not recorded", check$incomplete)))
    by_state <- stats::setNames(record$steps, vapply(record$steps, `[[`, "", "state"))
    named <- function(parameters) stats::setNames(parameters, vapply(parameters, `[[`, "", "name"))
    ruv <- named(by_state$ruv_corrected$parameters)
    expect_setequal(names(ruv), norm$parameters)
    expect_identical(ruv$ruv_k$source, "provenance")
    expect_identical(decodeWorkflowRecordValue(ruv$ruv_k$value), 5)
    expect_identical(ruv$ruv_k$provenance, "default")
    selection <- named(by_state$correlation_filtered$parameters)
    expect_identical(decodeWorkflowRecordValue(selection$threshold$value), 0.5)
    expect_identical(by_state$correlation_filtered$preceded_by[[1]]$operation, "pearson-pairs")
    pairs <- named(by_state$correlation_filtered$preceded_by[[1]]$parameters)
    expect_identical(pairs$tech_rep_remove_regex$provenance, "default")
    # A table among the operation's inputs is never carried into the head.
    expect_null(manager$.context$snapshot$step_provenance$correlation_filtered$values$correlations)
    da <- record$analyses$differential_abundance
    expect_identical(da$source_state, "correlation_filtered")
    da_parameters <- named(da$parameters)
    expect_identical(da_parameters$formula$provenance, "default")
    expect_identical(da_parameters$qvalue$provenance, "default")
    expect_identical(da_parameters$lfc$provenance, "default")
    expect_identical(da$model$block_reason, "no_duplicates")
})

test_that("KIT-027 audit PH-04/07/08: backfilled steps, exact default comparison and a visible record failure", {
    fixture <- analysisRecordFixture(withr::local_tempdir())
    manager <- fixture$manager
    withr::defer(closeDIAStateAdapter(manager))
    snapshot <- manager$.context$snapshot
    snapshot$step_provenance$raw_data_s4$backfilled <- TRUE
    snapshot$step_provenance$raw_data_s4$committed_at <- "2030-01-01T00:00:00.000+0000"
    publishDIAStateHead(manager$.context, snapshot)
    step <- Filter(function(x) identical(x$state, "raw_data_s4"), buildWorkflowAnalysisRecord(manager)$steps)[[1]]
    expect_true(is.na(step$software$sha))
    expect_true(is.na(step$defaults_sha256))
    expect_false(identical(step$committed_at, "2030-01-01T00:00:00.000+0000"))
    expect_true(workflowParameterMatchesDefault("(\\\\d+),x", "(\\\\d+),x"))
    expect_true(workflowParameterMatchesDefault(c("pdf", "png"), "pdf,png"))
    expect_false(workflowParameterMatchesDefault(0.0100000001, "0.01"))
    expect_true(workflowParameterMatchesDefault(0.01, "0.01"))
    notified <- NULL
    # KIT-108 W4: the core's save under proteomics' group, its final state replaced as before.
    config <- omicGuiConfig("proteomics")$summary
    config$final_state$context <- function(...) list(finalS4Object = NULL)
    saveWorkflowSummaryParameters(config, list(state_manager = structure(list(), class = "DIAStateAdapter")),
        projectDirs = list(proteomics = list(source_dir = withr::local_tempdir(), integration_dir = withr::local_tempdir(),
            base_dir = withr::local_tempdir())), omicType = "proteomics", label = "x", description = "d",
        notify = function(message, ...) notified <<- c(notified, message), done = function(...) NULL,
        failed = function(...) NULL, createWorkflowArgsFn = function(...) "study_parameters.txt",
        writeAnalysisRecordFn = function(...) stop("disk full"), saveRDSFn = function(...) NULL, catFn = function(...) NULL)
    expect_true(any(grepl("analysis record .* not written: disk full", notified)))
})

test_that("KIT-028: replay refuses mismatched inputs, a store beyond its preconditions and a legacy step", {
    root <- withr::local_tempdir()
    manager <- newDIAStateAdapter(newWorkflowEntityStore(root, "replay-test", "DIA"), WorkflowState$new(audit_enabled = FALSE))
    manager$.context$reference_states <- TRUE
    fasta <- file.path(root, "proteins.fasta")
    writeLines(c(">sp|P1|ONE", "MPEPTIDE"), fasta)
    recordDIAWorkflowInputs(manager, c(fasta = fasta))
    record <- jsonlite::fromJSON(jsonlite::toJSON(buildWorkflowAnalysisRecord(manager), auto_unbox = TRUE, digits = I(17),
        null = "null", na = "null"), simplifyVector = FALSE)
    closeDIAStateAdapter(manager)
    expect_identical(vapply(record$steps, `[[`, "", "operation"), "import")
    other <- file.path(withr::local_tempdir(), "other.fasta")
    writeLines(">sp|X|OTHER", other)
    scratch <- withr::local_tempdir()
    # Inputs are checked before the store is opened.
    expect_error(replayWorkflowAnalysisRecord(record, file.path(root, "absent"), c(fasta = other), list(), scratch),
        class = "workflow_replay_inputs")
    expect_error(replayWorkflowAnalysisRecord(record, root, character(), list(), scratch), class = "workflow_replay_inputs")
    # A legacy step after the preconditions stops replay, naming the step.
    legacy <- record
    legacy$steps[[2]] <- list(order = 2, state = "qvalue_filtered", operation = "legacy", step_label = "qvalue_filtered",
        parameters = list(), products = list())
    expect_error(replayWorkflowAnalysisRecord(legacy, root, c(fasta = fasta), list(), scratch),
        regexp = "qvalue_filtered", class = "workflow_replay_unsupported")
    # A record whose import state is not the store's is refused before any step.
    moved <- record
    moved$steps[[1]]$result_entity_id <- strrep("0", 64)
    expect_error(replayWorkflowAnalysisRecord(moved, root, c(fasta = fasta), list(), scratch), class = "workflow_replay_precondition")
})

docxText <- function(path) {
    directory <- withr::local_tempdir(.local_envir = parent.frame())
    utils::unzip(path, files = "word/document.xml", exdir = directory)
    xml <- paste(readLines(file.path(directory, "word", "document.xml"), warn = FALSE), collapse = "")
    gsub("<[^>]+>", "", xml)
}

test_that("KIT-029: the workbook has one parameters row per record parameter and matches the record", {
    fixture <- analysisRecordFixture(withr::local_tempdir())
    manager <- fixture$manager
    withr::defer(closeDIAStateAdapter(manager))
    json <- withr::local_tempfile(fileext = ".json")
    writeWorkflowAnalysisRecord(manager, json)
    record <- jsonlite::read_json(json, simplifyVector = FALSE)
    long <- record
    long$steps[[2]]$parameters[[length(long$steps[[2]]$parameters) + 1L]] <- list(name = "long_vector", source = "args",
        value = encodeWorkflowRecordValue(seq_len(300) + 0.5), provenance = "input")
    xlsx <- withr::local_tempfile(fileext = ".xlsx")
    writeWorkflowAnalysisRecordWorkbook(long, xlsx)
    parameters <- openxlsx::read.xlsx(xlsx, "parameters")
    expected <- sum(vapply(long$steps, function(step) length(step$parameters) +
        sum(vapply(step$preceded_by %||% list(), function(p) length(p$parameters), integer(1))), integer(1)))
    expect_identical(nrow(parameters), as.integer(expected))
    expect_identical(sort(openxlsx::read.xlsx(xlsx, "inputs")$sha256), sort(vapply(long$inputs, function(i) i$sha256 %||% NA_character_, "")[
        !is.na(vapply(long$inputs, function(i) i$sha256 %||% NA_character_, ""))]))
    expect_identical(openxlsx::read.xlsx(xlsx, "external_resources")$release, "12.0")
    expect_identical(openxlsx::read.xlsx(xlsx, "software")$value[[1]], "MultiScholaR")
    long_value <- parameters$value[parameters$name == "long_vector"]
    expect_identical(length(strsplit(long_value, ", ", fixed = TRUE)[[1]]), 300L)
})

test_that("KIT-029: the record renders to DOCX and HTML with every parameter name, input hash and the full long vector", {
    skip_if_not(rmarkdown::pandoc_available())
    skip_if(!nzchar(system.file("reports", "shared", "analysis_record_SI.rmd", package = "MultiScholaR")))
    fixture <- analysisRecordFixture(withr::local_tempdir())
    manager <- fixture$manager
    withr::defer(closeDIAStateAdapter(manager))
    json <- withr::local_tempfile(fileext = ".json")
    writeWorkflowAnalysisRecord(manager, json)
    record <- jsonlite::read_json(json, simplifyVector = FALSE)
    record$steps[[2]]$parameters[[length(record$steps[[2]]$parameters) + 1L]] <- list(name = "long_vector", source = "args",
        value = encodeWorkflowRecordValue(seq_len(300) + 0.5), provenance = "input")
    jsonlite::write_json(record, json, auto_unbox = TRUE, digits = I(17), null = "null", na = "null")
    out <- withr::local_tempdir()
    # Relative paths, as a caller in the project directory passes them (render changes directory first).
    withr::local_dir(dirname(json))
    docx <- renderWorkflowAnalysisRecord(basename(json), out, "docx")
    html <- renderWorkflowAnalysisRecord(json, out, "html")
    texts <- list(docx = docxText(docx), html = paste(readLines(html, warn = FALSE), collapse = "\n"))
    names_all <- unique(unlist(lapply(record$steps, function(step) vapply(step$parameters, `[[`, "", "name"))))
    hashes <- Filter(Negate(is.null), lapply(record$inputs, `[[`, "sha256"))
    for (kind in names(texts)) {
        text <- texts[[kind]]
        for (name in names_all) expect_true(grepl(name, text, fixed = TRUE), label = paste(kind, name))
        for (hash in hashes) expect_true(grepl(hash, text, fixed = TRUE), label = paste(kind, hash))
        expect_true(grepl("300.5", text, fixed = TRUE) && grepl("1.5, 2.5, 3.5", text, fixed = TRUE), label = paste(kind, "long vector"))
    }
})

test_that("KIT-029: report generation exports the record bundle for an on-disk owner and shows a failure", {
    exported <- NULL
    notified <- NULL
    workflow <- list(state_manager = structure(list(), class = "DIAStateAdapter"))
    # KIT-108 W4: the core's report entry under proteomics' group (runProtSummaryReportProgress deleted): the route's
    # template replaced and present in the project, the render replaced (a report written, or none); TRUE when the
    # report is generated, as before.
    dirs <- list(proteomics = list(base_dir = withr::local_tempdir(), study_report_dir = withr::local_tempdir()))
    dir.create(file.path(dirs$proteomics$base_dir, "scripts", "proteomics"), recursive = TRUE)
    writeLines("template", file.path(dirs$proteomics$base_dir, "scripts", "proteomics", "t.rmd"))
    rendered <- file.path(dirs$proteomics$study_report_dir, "t_proteomics_x.docx")
    writeLines("report", rendered)
    config <- omicGuiConfig("proteomics")$summary
    config$report$template <- function(...) list(templateFilename = "t.rmd")
    run <- function(generated, export) !is.null(renderWorkflowSummaryReport(config, workflow, dirs, "proteomics", "x",
        notify = function(message, ...) notified <<- c(notified, message), progress = function(...) NULL,
        renderFn = function(...) if (generated) rendered, exportRecordFn = export, catFn = function(...) NULL)$path)
    expect_true(run(TRUE, function(manager, directory) exported <<- directory))
    expect_identical(exported, dirs$proteomics$study_report_dir)
    exported <- NULL
    expect_false(run(FALSE, function(manager, directory) exported <<- directory))
    expect_null(exported)
    expect_true(run(TRUE, function(...) stop("pandoc missing")))
    expect_true(any(grepl("bundle was not exported: pandoc missing", notified)))
})

test_that("KIT-033: study_parameters.txt on the record path carries every value in full and the record's step parameters", {
    fixture <- analysisRecordFixture(withr::local_tempdir())
    manager <- fixture$manager
    withr::defer(closeDIAStateAdapter(manager))
    source_dir <- withr::local_tempdir()
    path <- suppressMessages(capture.output(file <- createWorkflowArgsFromConfig("record-test", source_dir_path = source_dir,
        final_s4_object = getDIAStateReference(manager, "qvalue_filtered"),
        workflow_data = list(state_manager = manager, config_list = list())), type = "output"))
    text <- readLines(file.path(source_dir, "study_parameters.txt"))
    expect_false(any(grepl("omitted|\\[list object\\]|\\.\\.\\.", text)))
    # The 40-sample inclusion list is written in full.
    inclusion <- grep("^  inclusion_list = ", text, value = TRUE)
    expect_identical(length(strsplit(sub("^  inclusion_list = ", "", inclusion[[1]]), ", ")[[1]]), 40L)
    expect_true("Analysis Record:" %in% text)
    record <- buildWorkflowAnalysisRecord(manager)
    for (step in Filter(function(s) !identical(s$operation, "legacy") && length(s$parameters), record$steps)) {
        for (parameter in step$parameters) expect_true(sprintf("  %s [%s] = %s", parameter$name, parameter$provenance,
            workflowRecordValueText(parameter$value)) %in% text, label = paste(step$state, parameter$name))
    }
})

test_that("KIT-033: an owner without a record keeps the current writer", {
    source_dir <- withr::local_tempdir()
    object <- makeWorkflowCheckpointPeptide(readWorkflowCheckpointRecipe("small"))
    object@args$filterMinNumPeptidesPerSample <- list(inclusion_list = paste0("s", 1:40))
    suppressMessages(capture.output(createWorkflowArgsFromConfig("plain", source_dir_path = source_dir, final_s4_object = object,
        workflow_data = list(state_manager = WorkflowState$new(audit_enabled = FALSE), config_list = list())), type = "output"))
    text <- readLines(file.path(source_dir, "study_parameters.txt"))
    expect_false("Analysis Record:" %in% text)
    expect_true(any(grepl("^  inclusion_list = character vector \\[40 values", text)))
})

test_that("KIT-033: the report takes the record's values and names the defaults it shows", {
    record <- list(steps = list(list(operation = "minimum_peptides_per_sample", parameters = list(
        list(name = "peptides_per_sample_cutoff", value = encodeWorkflowRecordValue(300))))),
        analyses = list(differential_abundance = list(parameters = list(list(name = "qvalue", value = encodeWorkflowRecordValue(0.1))))))
    path <- withr::local_tempfile(fileext = ".json")
    jsonlite::write_json(record, path, auto_unbox = TRUE, digits = I(17), null = "null", na = "null")
    params <- list(min_peptides_per_sample = 500, s4_de_q_val_thresh = 0.05, proportion_samples_below_cutoff = 0.5,
        missing_values_cutoff = 33.333, qvalue_threshold = 0.01)
    study <- list(`Workflow Parameters` = list(`Parameters from Final S4 Object` = list(
        removeRowsWithMissingValuesPercent = list(groupwise_percentage_cutoff = "60"),
        srlQvalueProteotypicPeptideClean = list(qvalue_threshold = "0.01"))))
    overlaid <- workflowRecordReportParameters(path, params, study)
    expect_identical(overlaid$min_peptides_per_sample, 300)
    expect_identical(overlaid$s4_de_q_val_thresh, 0.1)
    expect_true(grepl("Parameter not recorded; default shown: proportion_samples_below_cutoff = 0.5", overlaid$record_defaults_note, fixed = TRUE))
    expect_false(grepl("missing_values_cutoff|qvalue_threshold", overlaid$record_defaults_note))
    expect_identical(workflowRecordReportParameters(file.path(dirname(path), "absent.json"), params, study), params)
})

test_that("KIT-043: the LFQ/TMT templates' variables take the protein intensity filter and the DA settings from the record", {
    record <- list(steps = list(list(operation = "protein_intensity_filter", parameters = list(
        list(name = "proteins_intensity_cutoff_percentile", value = encodeWorkflowRecordValue(2)),
        list(name = "min_reps_per_group", value = encodeWorkflowRecordValue(2L)),
        list(name = "min_groups", value = encodeWorkflowRecordValue(2L))))),
        analyses = list(differential_abundance = list(parameters = list(
            list(name = "qvalue", value = encodeWorkflowRecordValue(0.1)), list(name = "lfc", value = encodeWorkflowRecordValue(0.5))))),
        software = list(version = "0.5.0", sha = "abc123", build = "def456", r_version = "R 4.6.1", blas = "openblas"))
    path <- withr::local_tempfile(fileext = ".json")
    jsonlite::write_json(record, path, auto_unbox = TRUE, digits = I(17), null = "null", na = "null")
    # The LFQ template's extracted parameters: template defaults for the filter and the DA settings.
    params <- list(proteins_intensity_cutoff_percentile = 1, groupwise_percentage_cutoff = 60, max_groups_percentage_cutoff = 60,
        ui_de_q_value_threshold = 0.05, ui_de_log_fold_change_cutoff = 0, s4_de_q_val_thresh = 0.05)
    study <- list(`Workflow Parameters` = list(`Parameters from Final S4 Object` = list(
        removeRowsWithMissingValuesPercent = list(groupwise_percentage_cutoff = "60"))))
    overlaid <- workflowRecordReportParameters(path, params, study)
    expect_identical(overlaid$proteins_intensity_cutoff_percentile, 2)
    expect_null(overlaid$min_reps_per_group)
    expect_identical(overlaid$ui_de_q_value_threshold, 0.1)
    expect_identical(overlaid$ui_de_log_fold_change_cutoff, 0.5)
    expect_identical(overlaid$s4_de_q_val_thresh, 0.1)
    expect_identical(overlaid$software_version, "0.5.0")
    expect_identical(overlaid$source_commit, "abc123")
    expect_true(isTRUE(overlaid$record_present))
    # The cross-group cutoff is neither recorded nor in the file: the methods text names the default it shows, once
    # (audit F-01: the source map holds every variable once).
    expect_identical(overlaid$record_defaults_note, " Parameter not recorded; default shown: max_groups_percentage_cutoff = 60.")
    expect_identical(anyDuplicated(names(workflowRecordReportSources)), 0L)
    # The DIA-NN template's parameters (its parser's names): the same record gives its DA text the record's q-value
    # and fold-change cutoff (audit F-02, declared), and no note it did not have before.
    diann <- list(missing_values_cutoff = 60, max_groups_percentage_cutoff = 0, intensity_cutoff_percentile = 1,
        ui_de_q_value_threshold = 0.05, ui_de_log_fold_change_cutoff = 0, s4_de_q_val_thresh = 0.05)
    diann_study <- list(`Workflow Parameters` = list(`Parameters from Final S4 Object` = list(
        removeRowsWithMissingValuesPercent = list(groupwise_percentage_cutoff = "60", max_groups_percentage_cutoff = "0"),
        peptideIntensityFiltering = list(peptides_intensity_cutoff_percentile = "1"))))
    overlaid_diann <- workflowRecordReportParameters(path, diann, diann_study)
    expect_identical(overlaid_diann$ui_de_q_value_threshold, 0.1)
    expect_identical(overlaid_diann$ui_de_log_fold_change_cutoff, 0.5)
    expect_identical(overlaid_diann$record_defaults_note, "")
})

test_that("KIT-043: study_parameters.txt on a protein-only owner is rendered from the record (audit F-04)", {
    root <- withr::local_tempdir()
    manager <- newDIAStateAdapter(newWorkflowEntityStore(root, "protein-only", "LFQ"), WorkflowState$new(audit_enabled = FALSE))
    withr::defer(closeDIAStateAdapter(manager))
    manager$.context$reference_states <- TRUE
    manager$setWorkflowType("LFQ")
    protein <- makeWorkflowCheckpointProtein(readWorkflowCheckpointRecipe("small"))
    setDIAStepProvenance(manager, "design_builder", fixed = list(pivot = "protein x run", transform = "log2(x + 1)"))
    manager$saveState("protein_s4_initial", protein, list(), "design")
    entry <- workflowOperationEntry("protein_intensity_filter")
    values <- list(strict_mode = FALSE, min_reps_per_group = 2L, min_groups = 2L, proteins_intensity_cutoff_percentile = 2L)
    setDIAStepProvenance(manager, "protein_intensity_filter", "whole", entry$parameters, values = values)
    manager$saveState("protein_intensity_filtered", protein, values, "intensity filter")
    source_dir <- withr::local_tempdir()
    invisible(suppressMessages(capture.output(createWorkflowArgsFromConfig("protein-only", source_dir_path = source_dir,
        final_s4_object = getDIAStateReference(manager), workflow_data = list(state_manager = manager, config_list = list())), type = "output")))
    text <- readLines(file.path(source_dir, "study_parameters.txt"))
    expect_true("Analysis Record:" %in% text)
    expect_false(any(grepl("omitted|\\[list object\\]|\\.\\.\\.", text)))
    record <- buildWorkflowAnalysisRecord(manager)
    expect_identical(record$workflow, "LFQ")
    expect_false(any(vapply(record$steps, function(step) identical(step$operation, "legacy"), logical(1))))
    for (step in Filter(function(s) length(s$parameters), record$steps)) for (parameter in step$parameters)
        expect_true(sprintf("  %s [%s] = %s", parameter$name, parameter$provenance, workflowRecordValueText(parameter$value)) %in% text,
            label = paste(step$state, parameter$name))
    expect_true(any(grepl("proteins_intensity_cutoff_percentile [user_set] = 2", text, fixed = TRUE)))
})

test_that("KIT-038: fixed constants are listed for every registered step, a legacy step alone dumps its arguments", {
    fixture <- analysisRecordFixture(withr::local_tempdir())
    manager <- fixture$manager
    withr::defer(closeDIAStateAdapter(manager))
    rollup <- workflowOperationEntry("precursor_rollup")
    expect_identical(rollup$declaration$fixed, list(backend = "P14"))
    # A registered step that declares no parameters (the P14 rollup) lists its fixed constants and nothing else.
    setDIAStepProvenance(manager, "precursor_rollup", "p14", rollup$parameters, fixed = rollup$declaration$fixed)
    manager$saveState("rolled_up", fixture$peptide, list(), "rollup")
    record <- buildWorkflowAnalysisRecord(manager)
    by_state <- stats::setNames(record$steps, vapply(record$steps, `[[`, "", "state"))
    rows <- by_state$rolled_up$parameters
    expect_length(rows, 1L)
    expect_identical(rows[[1]]$name, "backend")
    expect_identical(rows[[1]]$provenance, "fixed")
    expect_identical(decodeWorkflowRecordValue(rows[[1]]$value), "P14")
    # The legacy raw_data_s4 step still lists its arguments whole; the provenance vocabulary is closed and validated.
    raw <- by_state$raw_data_s4$parameters
    expect_true(length(raw) > 1L)
    expect_true(all(vapply(raw, `[[`, "", "provenance") == "unknown"))
    expect_true(validateWorkflowAnalysisRecord(record)$valid)
    expect_true(all(unlist(lapply(record$steps, function(step) vapply(step$parameters, `[[`, "", "provenance"))) %in%
        workflowRecordProvenanceValues))
    expect_false(any(unlist(lapply(record$steps, function(step) vapply(step$parameters, `[[`, "", "provenance"))) == "no_default"))
    record$steps[[2]]$parameters[[1]]$provenance <- "no_default"
    expect_true(validateWorkflowAnalysisRecord(record)$valid)
    record$steps[[2]]$parameters[[1]]$provenance <- "bogus"
    check <- validateWorkflowAnalysisRecord(record)
    expect_false(check$valid)
    expect_true(any(grepl("provenance vocabulary", check$errors)))
    # Every step carries outcomes with its dimensions; the record and each step name the build.
    expect_true(all(vapply(record$steps, function(step) is.list(step$outcomes) && !is.null(step$outcomes$rows), logical(1))))
    expect_identical(record$software$build, manager$.context$snapshot$step_provenance$rolled_up$build)
    expect_identical(by_state$rolled_up$software$build, record$software$build)
    for (name in c("blas", "lapack", "os_release", "locale", "duckdb_threads")) expect_false(is.null(record$software[[name]]))
    expect_identical(record$software$duckdb_threads$rollup, 8L)
})

test_that("KIT-038: the normalisation step states the realised RUV outcome and a requested NA is not_set", {
    fixture <- analysisRecordFixture(withr::local_tempdir())
    manager <- fixture$manager
    withr::defer(closeDIAStateAdapter(manager))
    protein <- makeWorkflowCheckpointProtein(readWorkflowCheckpointRecipe("small"))
    norm <- workflowOperationEntry("normalization-ruv")
    requested <- list(norm_method = "cyclicloess", ruv_mode = "automatic", auto_percentage_min = 1, auto_percentage_max = 20,
        separation_metric = "max_difference", k_penalty_weight = 0.5, adaptive_k_penalty = TRUE, ruv_percentage = 5,
        ruv_k = NA, grouping = "group")
    # The optimisation product carries the control index; the commit's config carries the realised k and percentage.
    controls <- c(rep(TRUE, 7), rep(FALSE, 13))
    product <- saveDIAStateValue(manager$.context, list(best_k = 3, controls = controls, percentage = 5), "normalization-science", codec = "rds")
    checkpoint <- checkpointDIAExecutionResult(manager, protein)
    checkpoint$products <- list(normalization_science = product)
    setDIAStepProvenance(manager, "normalization-ruv", "whole", norm$parameters, values = requested)
    commitDIAReferenceState(manager, "ruv_corrected", checkpoint, list(norm_method = "cyclicloess", ruv_mode = "automatic",
        ruv_k = 3, percentage_as_neg_ctrl = 5), "ruv")
    record <- buildWorkflowAnalysisRecord(manager)
    expect_true(validateWorkflowAnalysisRecord(record)$valid)
    by_state <- stats::setNames(record$steps, vapply(record$steps, `[[`, "", "state"))
    named <- function(parameters) stats::setNames(parameters, vapply(parameters, `[[`, "", "name"))
    ruv <- named(by_state$ruv_corrected$parameters)
    expect_identical(ruv$ruv_k$provenance, "not_set")
    expect_identical(ruv$ruv_mode$provenance, "no_shipped_default")
    outcomes <- by_state$ruv_corrected$outcomes
    expect_true(isTRUE(outcomes$ruv_applied))
    expect_identical(outcomes$k, 3)
    expect_identical(outcomes$negative_control_percentage, 5)
    expect_identical(outcomes$negative_control_count, 7L)
    # A skipped RUV commits `normalized`: applied FALSE, no k, no count; a skip flag in config surfaces as an outcome.
    setDIAStepProvenance(manager, "normalization-ruv", "whole", norm$parameters, values = requested)
    manager$saveState("normalized", protein, list(norm_method = "cyclicloess", ruv_mode = "skip", ruv_k = NA,
        percentage_as_neg_ctrl = NA, ruv_applied = FALSE, skipped = TRUE), "skipped")
    record <- buildWorkflowAnalysisRecord(manager)
    by_state <- stats::setNames(record$steps, vapply(record$steps, `[[`, "", "state"))
    skipped <- by_state$normalized$outcomes
    expect_false(skipped$ruv_applied)
    expect_true(is.na(skipped$k))
    expect_true(is.na(skipped$negative_control_count))
    expect_true(isTRUE(skipped$skipped))
})

test_that("KIT-039: the report overlay takes the software and the realised RUV outcome from the record, tolerating their absence", {
    record <- list(software = list(version = "0.5.0", sha = "abc123", build = "deadbeef", r_version = "R 4.6.1", blas = "openblas"),
        steps = list(list(operation = "normalization-ruv", parameters = list(), outcomes = list(ruv_applied = TRUE, k = 2,
            negative_control_percentage = 10, negative_control_count = 242))),
        analyses = list())
    path <- withr::local_tempfile(fileext = ".json")
    jsonlite::write_json(record, path, auto_unbox = TRUE, digits = I(17), null = "null", na = "null")
    params <- list(ruv_applied = FALSE, best_k = 5, best_percentage = "1", num_neg_ctrl = "3", commit = "N/A")
    overlaid <- workflowRecordReportParameters(path, params, list())
    expect_identical(overlaid$software_version, "0.5.0")
    expect_identical(overlaid$source_commit, "abc123")
    expect_identical(overlaid$build, "deadbeef")
    expect_true(overlaid$ruv_applied)
    expect_equal(overlaid$best_k, 2)
    expect_equal(overlaid$best_percentage, 10)
    expect_equal(overlaid$num_neg_ctrl, 242)
    expect_true(isTRUE(overlaid$record_present))
    # A KIT-033-vintage record (no software fields beyond the parser's, no outcomes) does not fail.
    old <- list(steps = list(list(operation = "normalization-ruv", parameters = list())), analyses = list())
    jsonlite::write_json(old, path, auto_unbox = TRUE, digits = I(17), null = "null", na = "null")
    tolerant <- workflowRecordReportParameters(path, params, list())
    expect_identical(tolerant$software_version, "not recorded")
    expect_false(tolerant$ruv_applied)
    expect_equal(tolerant$best_k, 5)
})

test_that("KIT-039: the record tables carry the step sentence and a label for every parameter", {
    fixture <- analysisRecordFixture(withr::local_tempdir())
    manager <- fixture$manager
    withr::defer(closeDIAStateAdapter(manager))
    record <- buildWorkflowAnalysisRecord(manager)
    tables <- workflowAnalysisRecordTables(record)
    expect_true("description" %in% names(tables$steps))
    expect_identical(tables$steps$description[tables$steps$state == "qvalue_filtered"], "q-value filter")
    expect_true("label" %in% names(tables$parameters))
    declared <- tables$parameters[tables$parameters$provenance != "unknown", ]
    expect_true(all(nzchar(declared$label)))
    expect_identical(declared$label[declared$name == "qvalue_threshold"][[1]], "Run-level q-value threshold")
    expect_true("build" %in% tables$software$field)
    expect_true("BLAS" %in% tables$software$field)
})

test_that("a DIA_LIMPA record names the dpcDE model and replays its protein stage, covariates and DA from the record (KIT-053)", {
    skip_if_not_installed("limpa")
    skip_if_not_installed("mixOmics")
    root <- withr::local_tempdir()
    manager <- newDIAStateAdapter(newWorkflowEntityStore(root, "limpa-record", "DIA"), WorkflowState$new(audit_enabled = FALSE))
    withr::defer(closeDIAStateAdapter(manager))
    manager$.context$reference_states <- TRUE
    manager$setWorkflowType("DIA_LIMPA")
    rolled <- rollUpPrecursorToPeptide(makeWorkflowCheckpointPeptide(readWorkflowCheckpointRecipe("small")), core_utilisation = NA)
    quantified <- suppressMessages(proteinMissingValueImputationLimpa(log2TransformPeptideMatrix(calcPeptideMatrix(rolled)),
        dpc_slope = 0.8, verbose = FALSE, chunk = 1000, engine = "limpa"))
    # The record's preconditions are the import and design states; this fixture store starts at the DPC state, so the
    # DPC state stands in for the design precondition (the replay checks the states it holds, not their names).
    setDIAStepProvenance(manager, "design_builder")
    manager$saveState("protein_dpc_quant", quantified, list(), "fixture")
    paths <- list(source_dir = file.path(root, "src"), protein_qc_dir = file.path(root, "pqc"), da_output_dir = file.path(root, "da"),
        publication_graphs_dir = file.path(root, "pub"), time_dir = file.path(root, "time"))
    for (p in paths) dir.create(p, recursive = TRUE, showWarnings = FALSE)
    manager$.context$experiment_paths <- paths
    workflow <- list2env(list(state_manager = manager, qc_params = list(protein_qc = list()), config_list = list()))
    contrasts <- data.frame(contrasts = "treated.vs.control", full_format = "treated.vs.control=grouptreated-groupcontrol",
        friendly_names = "treated.vs.control", stringsAsFactors = FALSE)
    withDIAStateScope(manager, function() {
        v1 <- computeDIAReferenceProteinNormalisation(manager, getDIAStateReference(manager), list(norm_method = "none"))
        setDIAStepProvenance(manager, "protein_between_sample_normalisation", "whole", "norm_method", NULL, values = list(norm_method = "none"))
        finishDIAReferenceProteinStep(workflow, v1, publish = FALSE)
        v2 <- computeDIAReferenceProteinPeptideCountFilter(manager, getDIAStateReference(manager),
            list(num_peptides_per_protein_thresh = 1, num_peptidoforms_per_protein_thresh = 2))
        setDIAStepProvenance(manager, "protein_peptide_count_filter", "whole", c("num_peptides_per_protein_thresh", "num_peptidoforms_per_protein_thresh"), NULL, values = v2$config[1:2])
        finishDIAReferenceProteinStep(workflow, v2, publish = FALSE)
        source <- getDIAStateReference(manager)
        pairs <- runSavedOperation(manager, diaPearsonPairsOperation(), source, list(group = "group", tech_rep_remove_regex = "pool"))
        selection <- runSavedOperation(manager, diaPearsonSelectionOperation(), source,
            list(pairs = pairs$reference$entity_id, threshold = 0.5, correlations = pairs$value), dependencies = list(pairs$reference))
        commitDIAReferenceState(manager, "correlation_filtered", readDIAStateMetadata(manager, selection$value)$data,
            list(min_pearson_correlation_threshold = 0.5), "correlation")
        parameters <- resolveDIARuvCovariatesParameters(list(), "group")
        run <- suppressMessages(suppressWarnings(runSavedOperation(manager, diaRuvCovariatesOperation(), getDIAStateReference(manager), parameters)))
        final <- readDIAStateMetadata(manager, run$value$covariates)$data
        final$products$ruv_covariates_science <- run$reference
        saved <- run$value
        # The configuration the route's runner (and the replay) commits, so the state summaries compare like for like.
        commitDIAReferenceState(manager, "ruv_covariates", final, list(ruv_mode = "covariates", ruv_k = saved$best_k, percentage_as_neg_ctrl = saved$percentage,
            control_count = saved$control_count, separation_score = saved$ruv_result$best_separation_score,
            composite_score = saved$ruv_result$best_composite_score, max_acceptable_k = saved$ruv_result$max_acceptable_k,
            replicate_structure = saved$replicate_structure, block_reason = saved$block_reason, replicate_blocks = saved$blocks,
            ruv_grouping_variable = parameters$ruv_grouping_variable, ruv_applied = TRUE,
            outcomes = list(k = saved$best_k, percentage_as_neg_ctrl = saved$percentage, control_count = saved$control_count,
                separation_score = saved$ruv_result$best_separation_score, composite_score = saved$ruv_result$best_composite_score,
                replicate_structure = saved$replicate_structure, block_reason = saved$block_reason, ruv_mode = "covariates")),
            sprintf("RUV-III factors appended to the design as covariates (W1..W%d); the quantified matrix unchanged", saved$best_k))
        da_data <- new.env(parent = emptyenv()); da_data$current_s4_object <- getDIAStateReference(manager); da_data$contrasts_available <- contrasts$contrasts
        workflow$contrasts_tbl <- contrasts
        manifest <- suppressMessages(suppressWarnings(runDIAReferenceDaAnalysis(list(formula_string = "~ 0 + group", da_q_val_thresh = 0.05, treat_lfc_cutoff = 0),
            da_data, workflow, paths, contrasts = contrasts, publish = FALSE)))
        expect_identical(manifest$model$model, "dpcDE")
        expect_identical(manifest$model$route, "limpa")
        expect_identical(manifest$model$formula_used, paste("~ 0 + group +", paste(paste0("W", seq_len(run$value$best_k)), collapse = " + ")))
        expect_identical(manifest$model$k, run$value$best_k)
        expect_false(manifest$model$treat_applied)
    })
    json <- file.path(root, "analysis_record.json")
    writeWorkflowAnalysisRecord(manager, json)
    record <- jsonlite::read_json(json, simplifyVector = FALSE)
    expect_true(validateWorkflowAnalysisRecord(record)$valid)
    expect_identical(record$workflow, "DIA_LIMPA")
    expect_identical(record$analyses$differential_abundance$model$model, "dpcDE")
    expect_identical(vapply(record$steps, function(s) s$operation, ""), c("import", "design_builder", "protein_between_sample_normalisation",
        "protein_peptide_count_filter", "pearson-selection", "protein_ruv_covariates"))
    overlay <- workflowRecordReportParameters(json, list(s4_formula_string = "~ 0 + group"), list())
    expect_identical(overlay$da_model, "dpcDE")
    expect_identical(overlay$limpa_norm_method, "none")
    ruv_step <- Filter(function(s) identical(s$state, "ruv_covariates"), record$steps)[[1L]]
    expect_identical(overlay$ruvcov_k, ruv_step$outcomes$k)
    expect_identical(overlay$limpa_min_peptides_per_protein, 1)
    # Replay from the DPC state on a copy of the store: every later state hash-identical, the DA refitted.
    origin <- file.path(root, ".workflow-entities")
    copy <- withr::local_tempdir()
    stopifnot(system2("cp", c("-a", file.path(root, "."), copy)) == 0L)
    h <- new.env(parent = globalenv())
    sys.source(file.path(testthat::test_path("..", ".."), "tools", "profiling", "workflow_kit.R"), h)
    api <- function(name) get(name, asNamespace("MultiScholaR"))
    prepared <- h$kitPrepareReplayProject(copy, c("initial", "protein_dpc_quant"), api)
    scratch <- withr::local_tempdir()
    replay <- h$kitReplayJourney(record, prepared, root, list(), paths, scratch, api)
    expect_true(all(c("protein_normalised", "peptide_count_filtered", "correlation_filtered", "ruv_covariates") %in% replay$states))
    expect_identical(as.character(replay$mismatches), character())
    expect_true(any(vapply(replay$per_step, function(x) identical(x$operation, "DA-contrast-output"), logical(1))))
    # KIT-056: the mode is a recorded parameter (covariates here) and reaches the overlay; a record written before the
    # parameter existed (KIT-053's) still replays as covariates: the same steps with the parameter removed.
    expect_identical(overlay$ruvcov_mode, "covariates")
    expect_identical(overlay$da_ruv_mode, "covariates")
    legacy <- record
    legacy$steps[[6L]]$parameters <- Filter(function(parameter) !identical(parameter$name, "ruv_mode"), legacy$steps[[6L]]$parameters)
    expect_false("ruv_mode" %in% vapply(legacy$steps[[6L]]$parameters, `[[`, "", "name"))
    copy_legacy <- withr::local_tempdir()
    stopifnot(system2("cp", c("-a", file.path(root, "."), copy_legacy)) == 0L)
    prepared_legacy <- h$kitPrepareReplayProject(copy_legacy, c("initial", "protein_dpc_quant"), api)
    replay_legacy <- h$kitReplayJourney(legacy, prepared_legacy, root, list(), paths, withr::local_tempdir(), api)
    expect_true("ruv_covariates" %in% replay_legacy$states)
    expect_identical(as.character(replay_legacy$mismatches), character())
})

test_that("a DIA_LIMPA record in the RUV correction mode commits ruv_corrected, names the mode in the model and replays exact (KIT-056)", {
    skip_if_not_installed("limpa")
    skip_if_not_installed("mixOmics")
    root <- withr::local_tempdir()
    manager <- newDIAStateAdapter(newWorkflowEntityStore(root, "limpa-record-corr", "DIA"), WorkflowState$new(audit_enabled = FALSE))
    withr::defer(closeDIAStateAdapter(manager))
    manager$.context$reference_states <- TRUE
    manager$setWorkflowType("DIA_LIMPA")
    rolled <- rollUpPrecursorToPeptide(makeWorkflowCheckpointPeptide(readWorkflowCheckpointRecipe("small")), core_utilisation = NA)
    quantified <- suppressMessages(proteinMissingValueImputationLimpa(log2TransformPeptideMatrix(calcPeptideMatrix(rolled)),
        dpc_slope = 0.8, verbose = FALSE, chunk = 1000, engine = "limpa"))
    setDIAStepProvenance(manager, "design_builder")
    manager$saveState("protein_dpc_quant", quantified, list(), "fixture")
    paths <- list(source_dir = file.path(root, "src"), protein_qc_dir = file.path(root, "pqc"), da_output_dir = file.path(root, "da"),
        publication_graphs_dir = file.path(root, "pub"), time_dir = file.path(root, "time"))
    for (p in paths) dir.create(p, recursive = TRUE, showWarnings = FALSE)
    manager$.context$experiment_paths <- paths
    workflow <- list2env(list(state_manager = manager, qc_params = list(protein_qc = list()), config_list = list()))
    contrasts <- data.frame(contrasts = "treated.vs.control", full_format = "treated.vs.control=grouptreated-groupcontrol",
        friendly_names = "treated.vs.control", stringsAsFactors = FALSE)
    norm_data <- new.env(parent = emptyenv())
    withDIAStateScope(manager, function() {
        v1 <- computeDIAReferenceProteinNormalisation(manager, getDIAStateReference(manager), list(norm_method = "none"))
        setDIAStepProvenance(manager, "protein_between_sample_normalisation", "whole", "norm_method", NULL, values = list(norm_method = "none"))
        finishDIAReferenceProteinStep(workflow, v1, publish = FALSE)
        v2 <- computeDIAReferenceProteinPeptideCountFilter(manager, getDIAStateReference(manager),
            list(num_peptides_per_protein_thresh = 1, num_peptidoforms_per_protein_thresh = 2))
        setDIAStepProvenance(manager, "protein_peptide_count_filter", "whole", c("num_peptides_per_protein_thresh", "num_peptidoforms_per_protein_thresh"), NULL, values = v2$config[1:2])
        finishDIAReferenceProteinStep(workflow, v2, publish = FALSE)
        source <- getDIAStateReference(manager)
        pairs <- runSavedOperation(manager, diaPearsonPairsOperation(), source, list(group = "group", tech_rep_remove_regex = "pool"))
        selection <- runSavedOperation(manager, diaPearsonSelectionOperation(), source,
            list(pairs = pairs$reference$entity_id, threshold = 0.5, correlations = pairs$value), dependencies = list(pairs$reference))
        commitDIAReferenceState(manager, "correlation_filtered", readDIAStateMetadata(manager, selection$value)$data,
            list(min_pearson_correlation_threshold = 0.5), "correlation")
        norm_data$correlation_filtered_obj <- getDIAStateReference(manager)
        # The route's runner in correction mode (as the tab calls it), publication to the fixture paths.
        accepted <- suppressMessages(suppressWarnings(runDIAReferenceRuvCovariates(list(ruv_mode = "correction", color_variable = NULL, shape_variable = NULL),
            workflow, norm_data, paths, "proteomics", "fixture", function(reference) NULL, function() "group")))
        expect_identical(manager$current_state, "ruv_corrected")
        corrected <- readDIAExecutionComponents(manager)
        expect_false(any(grepl("^W[0-9]+$", names(corrected@design_matrix))))
        expect_identical(corrected@args$ruvCovariates$ruv_mode, "correction")
        expect_identical(readDIAStateMetadata(manager)$config$ruv_mode, "correction")
        da_data <- new.env(parent = emptyenv()); da_data$current_s4_object <- getDIAStateReference(manager); da_data$contrasts_available <- contrasts$contrasts
        workflow$contrasts_tbl <- contrasts
        manifest <- suppressMessages(suppressWarnings(runDIAReferenceDaAnalysis(list(formula_string = "~ 0 + group", da_q_val_thresh = 0.05, treat_lfc_cutoff = 0),
            da_data, workflow, paths, contrasts = contrasts, publish = FALSE)))
        expect_identical(manifest$model$model, "dpcDE")
        expect_identical(manifest$model$formula_used, "~ 0 + group")
        expect_identical(manifest$model$k, 0L)
        expect_identical(manifest$model$ruv_mode, "correction")
    })
    json <- file.path(root, "analysis_record.json")
    writeWorkflowAnalysisRecord(manager, json)
    record <- jsonlite::read_json(json, simplifyVector = FALSE)
    expect_true(validateWorkflowAnalysisRecord(record)$valid)
    ruv_step <- Filter(function(s) identical(s$state, "ruv_corrected"), record$steps)[[1L]]
    expect_identical(ruv_step$operation, "protein_ruv_covariates")
    expect_identical(ruv_step$outcomes$ruv_mode, "correction")
    expect_true("ruv_mode" %in% vapply(ruv_step$parameters, `[[`, "", "name"))
    overlay <- workflowRecordReportParameters(json, list(s4_formula_string = "~ 0 + group"), list())
    expect_identical(overlay$ruvcov_mode, "correction")
    expect_identical(overlay$da_ruv_mode, "correction")
    expect_identical(overlay$da_formula_used, "~ 0 + group")
    copy <- withr::local_tempdir()
    stopifnot(system2("cp", c("-a", file.path(root, "."), copy)) == 0L)
    h <- new.env(parent = globalenv())
    sys.source(file.path(testthat::test_path("..", ".."), "tools", "profiling", "workflow_kit.R"), h)
    api <- function(name) get(name, asNamespace("MultiScholaR"))
    prepared <- h$kitPrepareReplayProject(copy, c("initial", "protein_dpc_quant"), api)
    replay <- h$kitReplayJourney(record, prepared, root, list(), paths, withr::local_tempdir(), api)
    expect_true("ruv_corrected" %in% replay$states)
    expect_identical(as.character(replay$mismatches), character())
})

# KIT-107 F-16 (D4, D5, D6 test 3): the enrichment step carries its provenance and replays after the DA.

#' The proteomics `enrichment` group with its function fields (KIT-107 D1): until the W2 merge wires them in
#' utils_omic_gui_config.R they are optional; this fills an absent one with the value the merge wires, and changes
#' nothing once they are wired.
recordEnrichmentConfig <- function(config) {
    config[["organism"]] <- config[["organism"]] %||% list(supported = buildProtEnrichSupportedOrganisms,
        method = resolveProtEnrichAnalysisMethod)
    config$annotation <- config$annotation %||% list(resolve = resolveProtEnrichUniprotAnnotations,
        match = resolveProtEnrichAnnotationMatching)
    config$kernel <- config$kernel %||% list(select = resolveProtEnrichProcessEnrichmentsFn,
        test_double = runProtEnrichDeterministicProcessEnrichments)
    config$organism_filter$mapping <- config$organism_filter$mapping %||% resolveProtEnrichOrganismMapping
    config$organism_filter$apply <- config$organism_filter$apply %||% applyProtEnrichOrganismFilter
    config
}

#' The enrichment product a store's head roots (schema `dia_enrichment_product_v1`), hydrated.
recordEnrichmentProduct <- function(manager) {
    for (reference in manager$.context$snapshot$products) {
        value <- tryCatch(hydrateSavedWorkflowEntity(manager$.context$store, reference)$value, error = function(e) NULL)
        if (is.list(value) && identical(value$schema, "dia_enrichment_product_v1")) return(value)
    }
    NULL
}

test_that("KIT-107 F-16: a DIA record's enrichment step replays after the DA to the same enrichment data, summaries and product (D6 test 3)", {
    withr::local_envvar(MULTISCHOLAR_TEST_MODE = "true")
    enrichment <- omicGuiProteomicsEnrichment
    testthat::local_mocked_bindings(omicGuiProteomicsEnrichment = function() recordEnrichmentConfig(enrichment()))
    root <- withr::local_tempdir()
    manager <- newDIAStateAdapter(newWorkflowEntityStore(root, "enrichment-record", "DIA"), WorkflowState$new(audit_enabled = FALSE))
    withr::defer(closeDIAStateAdapter(manager))
    manager$.context$reference_states <- TRUE
    manager$setWorkflowType("DIA")
    # W1 (b)'s fixture state (MCI-009's object, row names dropped), standing in for the design precondition.
    design <- module_ci_prot_da_design()
    rownames(design) <- NULL
    quant <- module_ci_prot_da_table()
    rownames(quant) <- NULL
    setDIAStepProvenance(manager, "design_builder")
    manager$saveState("correlation_filtered", module_ci_prot_da_object(design = design, protein_quant_table = quant), list(), "fixture")
    paths <- list(source_dir = file.path(root, "source"), results_dir = file.path(root, "results"),
        da_output_dir = file.path(root, "results", "da_proteins"), pathway_dir = file.path(root, "results", "pathway_enrichment"),
        publication_graphs_dir = file.path(root, "publication_graphs"))
    for (p in paths) dir.create(p, recursive = TRUE, showWarnings = FALSE)
    manager$.context$experiment_paths <- paths
    contrasts <- module_ci_prot_da_contrasts("two_group")
    workflow <- list2env(list(state_manager = manager, contrasts_tbl = contrasts, tab_status = list(),
        qc_params = list(protein_qc = list()), config_list = list()))
    origin <- NULL
    withDIAStateScope(manager, function() {
        da_data <- new.env(parent = emptyenv())
        da_data$current_s4_object <- getDIAStateReference(manager)
        da_data$contrasts_available <- contrasts$contrasts
        manifest <- suppressMessages(suppressWarnings(runDIAReferenceDaAnalysis(list(formula_string = "~ 0 + group",
            da_q_val_thresh = 0.05, treat_lfc_cutoff = 0), da_data, workflow, paths, contrasts = contrasts, publish = FALSE)))
        # The DA tables as the GUI's DA publication writes them (its one body), then the enrichment as the tab runs it.
        suppressMessages(suppressWarnings(produceDIADaPublication(manager, manifest, data.frame(), 0.05,
            paths[c("da_output_dir", "publication_graphs_dir")], save_assets = FALSE)))
        utils::capture.output(origin <<- suppressMessages(suppressWarnings(runWorkflowEnrichment(workflow,
            omicGuiConfig("proteomics")$enrichment, manifest, list(selected_contrast = "B_vs_A", organism_taxid = "9606",
                up_cutoff = 0, down_cutoff = 0, q_cutoff = 0.05, correction_method = "gSCS", enable_organism_filter = FALSE),
            paths))))
    })
    expect_identical(manager$current_state, "enrichment_completed")
    json <- file.path(root, "analysis_record.json")
    writeWorkflowAnalysisRecord(manager, json)
    record <- jsonlite::read_json(json, simplifyVector = FALSE)
    check <- validateWorkflowAnalysisRecord(record)
    expect_true(check$valid)
    expect_false(any(grepl("enrichment_completed", check$incomplete, fixed = TRUE)))
    expect_identical(vapply(record$steps, function(s) s$operation, ""), c("import", "design_builder", "enrichment"))
    step <- record$steps[[3L]]
    expect_identical(step$variant, "whole")
    expect_identical(vapply(step$parameters, `[[`, "", "name"), diaEnrichmentOperation()$parameters)
    expect_identical(vapply(step$parameters, `[[`, "", "source"), rep("provenance", 9L))
    expect_false(is.null(record$analyses$differential_abundance))
    state <- withDIAStateScope(manager, function() manager$getState("enrichment_completed"))
    product <- recordEnrichmentProduct(manager)
    expect_identical(product$parameters, origin$parameters[setdiff(names(origin$parameters), "results")])
    expect_identical(product$results, origin$parameters$results)
    # Replay on a copy cut to the import and design states: the DA refitted, then the enrichment on the double.
    copy <- withr::local_tempdir()
    stopifnot(system2("cp", c("-a", file.path(root, "."), copy)) == 0L)
    h <- new.env(parent = globalenv())
    sys.source(file.path(testthat::test_path("..", ".."), "tools", "profiling", "workflow_kit.R"), h)
    api <- function(name) get(name, asNamespace("MultiScholaR"))
    prepared <- h$kitPrepareReplayProject(copy, c("initial", "correlation_filtered"), api)
    seen <- character()
    replayed <- list()
    utils::capture.output(replay <- suppressMessages(suppressWarnings(replayWorkflowAnalysisRecord(record, prepared,
        character(), paths, withr::local_tempdir(), on_step = function(step, manager) {
            seen <<- c(seen, step$operation)
            if (identical(step$operation, "enrichment")) {
                replayed$state <<- manager$getState("enrichment_completed")
                replayed$product <<- recordEnrichmentProduct(manager)
                replayed$provenance <<- manager$.context$snapshot$step_provenance$enrichment_completed
                replayed$tables <<- lapply(replayed$product$results, function(id)
                    describeSavedWorkflowEntity(manager$.context$store, id)$identity)
            }
        }))))
    expect_identical(replay$states, c("initial", "correlation_filtered", "enrichment_completed"))
    expect_identical(seen, c("DA-contrast-output", "enrichment"))
    # The KIT-016 measure: the enrichment data and summaries, and the product's parameters and results; the plot slots
    # carry random widget ids and are not compared. `results` holds the DA tables' entity ids, which are store-local (an
    # entity id hashes its payload's file name): the replayed tables are compared by their content identity.
    expect_identical(digest::digest(replayed$state@enrichment_data), digest::digest(state@enrichment_data))
    expect_identical(digest::digest(replayed$state@enrichment_summaries), digest::digest(state@enrichment_summaries))
    expect_identical(replayed$product$parameters, product$parameters)
    expect_identical(names(replayed$product$results), names(product$results))
    expect_identical(replayed$tables, lapply(product$results, function(id) describeSavedWorkflowEntity(manager$.context$store, id)$identity))
    expect_false(is.null(replayed$tables[[1L]]$input))
    expect_identical(replayed$product$service, product$service)
    expect_identical(replayed$provenance$operation, "enrichment")
    recorded <- manager$.context$snapshot$step_provenance$enrichment_completed$values
    keep <- setdiff(names(recorded), "results")
    expect_identical(replayed$provenance$values[keep], recorded[keep])
    expect_identical(names(replayed$provenance$values$results), names(recorded$results))
    # A record whose enrichment step is legacy (a store written before F-16) is still refused by name.
    legacy <- record
    legacy$steps[[3L]]$operation <- "legacy"
    legacy$steps[[3L]]$parameters <- list()
    copy_legacy <- withr::local_tempdir()
    stopifnot(system2("cp", c("-a", file.path(root, "."), copy_legacy)) == 0L)
    prepared_legacy <- h$kitPrepareReplayProject(copy_legacy, c("initial", "correlation_filtered"), api)
    utils::capture.output(expect_error(suppressMessages(suppressWarnings(replayWorkflowAnalysisRecord(legacy, prepared_legacy,
        character(), paths, withr::local_tempdir()))), regexp = "enrichment_completed", class = "workflow_replay_unsupported"))
    # A recorded value the replay does not produce refuses the replay, naming the parameter.
    changed <- record
    expect_identical(changed$steps[[3L]]$parameters[[7L]]$name, "protein_id_column")
    changed$steps[[3L]]$parameters[[7L]]$value$value <- "not_the_id_column"
    copy_changed <- withr::local_tempdir()
    stopifnot(system2("cp", c("-a", file.path(root, "."), copy_changed)) == 0L)
    prepared_changed <- h$kitPrepareReplayProject(copy_changed, c("initial", "correlation_filtered"), api)
    utils::capture.output(expect_error(suppressMessages(suppressWarnings(replayWorkflowAnalysisRecord(changed, prepared_changed,
        character(), paths, withr::local_tempdir()))), regexp = "protein_id_column", class = "workflow_replay_invalid"))
    # A record with an enrichment step and no DA analysis is refused before any step runs.
    no_da <- record
    no_da$analyses$differential_abundance <- NULL
    expect_error(replayWorkflowAnalysisRecord(no_da, prepared, character(), paths, withr::local_tempdir()),
        regexp = "enrichment_completed", class = "workflow_replay_unsupported")
})
