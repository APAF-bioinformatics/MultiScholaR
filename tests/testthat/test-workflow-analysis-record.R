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

test_that("the summary page writes analysis_record.json after study_parameters.txt for an on-disk owner", {
    calls <- character()
    record_args <- NULL
    manager <- structure(list(), class = "DIAStateAdapter")
    workflow <- list(state_manager = manager, contrasts_tbl = NULL, config_list = NULL)
    source_dir <- withr::local_tempdir()
    completeProtSummaryWorkflowArgsSave(output = list(), values = list(), projectDirs = list(proteomics = list(source_dir = source_dir,
            integration_dir = withr::local_tempdir(), base_dir = source_dir)),
        workflowData = workflow, experimentLabel = "x", description = "d",
        resolveFinalS4StateFn = function(...) list(finalS4Object = NULL),
        createWorkflowArgsFn = function(...) { calls <<- c(calls, "study_parameters"); "study_parameters.txt" },
        writeAnalysisRecordFn = function(manager, path) { calls <<- c(calls, "analysis_record"); record_args <<- path },
        saveRDSFn = function(...) NULL, renderTextFn = function(...) NULL, showNotificationFn = function(...) NULL,
        catFn = function(...) NULL)
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
    completeProtSummaryWorkflowArgsSave(output = list(), values = list(), projectDirs = list(proteomics = list(
            source_dir = withr::local_tempdir(), integration_dir = withr::local_tempdir(), base_dir = withr::local_tempdir())),
        workflowData = list(state_manager = structure(list(), class = "DIAStateAdapter")), experimentLabel = "x", description = "d",
        resolveFinalS4StateFn = function(...) list(finalS4Object = NULL), createWorkflowArgsFn = function(...) "study_parameters.txt",
        writeAnalysisRecordFn = function(...) stop("disk full"), saveRDSFn = function(...) NULL, renderTextFn = function(...) NULL,
        showNotificationFn = function(message, ...) notified <<- c(notified, message), catFn = function(...) NULL)
    expect_true(any(grepl("analysis record .* not written: disk full", notified)))
})
