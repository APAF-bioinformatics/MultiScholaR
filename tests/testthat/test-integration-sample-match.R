# KIT-063 item A (KIT-139 follow-up): the sample match edited by the user. Samples are matched across omics by the
# biological-sample column each omic saved; a user whose omics name one biological sample differently uploads a match
# table (per omic, a run or a saved name, and its biological sample). Declared before the run:
#   - without a table, or with the downloaded table uploaded unchanged, the parameters and the model's identity are
#     those of before (the same saved model is reused);
#   - on the tier-1 corpus with the transcriptomics samples renamed, the uploaded table (by run, or by saved name)
#     gives views identical to those of the corpus as composed;
#   - a run mapped twice and a run without a biological sample are refused; unknown omics, runs and names are counted
#     and ignored; a biological sample given two groups is refused by the alignment, as without a table;
#   - the match is part of the model's identity, recorded with the model and replayed.
# MOFA2 is replaced by a stand-in trainer; the end-to-end training with an uploaded table is
# test-integration-results-mofa.R on the MOFA2 lane.

intMatchCorpus <- function() normalizePath(testthat::test_path("..", "testdata", "integration", "corpus"), mustWork = FALSE)
intMatchOmics <- c("lipidomics", "metabolomics", "proteomics", "transcriptomics")

intMatchTrainer <- function(views, options) {
    samples <- colnames(views[[1L]])
    k <- options$num_factors
    labels <- paste0("Factor", seq_len(k))
    list(factors = matrix(seq_len(length(samples) * k) / 10, length(samples), k, dimnames = list(samples, labels)),
        weights = lapply(views, function(v) matrix(0.5, nrow(v), k, dimnames = list(rownames(v), labels))),
        r2_per_factor = matrix(1, k, length(views), dimnames = list(labels, names(views))),
        r2_total = stats::setNames(rep(k, length(views)), names(views)), elbo = c(-20, -10, NaN), iterations = 2L,
        options = options, versions = list(MOFA2 = "stand-in"), messages = character(),
        model_file = list(bytes = as.raw(1:3), sha256 = strrep("b", 64L), size = 3))
}

# A project whose integration folder holds the corpus; with `rename`, the transcriptomics input is written again with
# its biological samples renamed (`S-<name>`), as an omic that names its samples differently would save them.
intMatchProject <- function(rename = FALSE, .local_envir = parent.frame()) {
    testthat::skip_if_not(dir.exists(intMatchCorpus()), "the tier-1 corpus is not present")
    testthat::skip_if_not(workflowOwnerCapable(), "the integration store needs qs2, DuckDB and Arrow")
    root <- withr::local_tempdir(.local_envir = .local_envir)
    paths <- list(base_dir = root, integration_dir = file.path(root, "integration"),
        source_dir = file.path(root, "scripts", "integration"), results_dir = file.path(root, "results", "integration"),
        multiomic_inputs_dir = file.path(root, "results", "integration", "multiomic_inputs"))
    dir.create(paths$source_dir, recursive = TRUE)
    for (omic in intMatchOmics) {
        dir <- file.path(paths$integration_dir, omic)
        dir.create(dir, recursive = TRUE)
        stopifnot(all(file.copy(list.files(file.path(intMatchCorpus(), omic), full.names = TRUE), dir)))
    }
    if (rename) {
        dir <- file.path(paths$integration_dir, "transcriptomics")
        input <- readIntegrationInput(dir)
        read <- input$manifest
        manifest <- list(omic_type = read$omic_type, workflow = read$workflow, matrix_role = read$matrix_role, scale = read$scale,
            correction = NULL, source_state = list(name = "kit063-renamed"), sample_id_column = read$sample_id_column,
            group_column = read$group_column, technical_replicate_column = read$technical_replicate_column,
            biological_sample_column = read$biological_sample_column)
        samples <- input$samples
        column <- manifest$biological_sample_column
        samples[[column]] <- paste0("S-", samples[[column]])
        writeIntegrationInput(dir, manifest, samples, lapply(input$views, `[[`, "matrix"), lapply(input$views, `[[`, "features"))
    }
    list(root = root, paths = paths)
}

intMatchStore <- function(paths, name = "store", .local_envir = parent.frame()) {
    manager <- newIntegrationStateOwner(file.path(paths$source_dir, name))
    withr::defer(closeWorkflowStateOwner(manager), envir = .local_envir)
    scan <- integrationScanProject(paths$integration_dir)
    integrationCommitImport(manager, integrationImportRequest(scan, scan$views$view))
    manager
}

intMatchBus <- function(manager, paths) {
    bus <- do.call(shiny::reactiveValues, omicWorkflowInitialBus(integrationGuiConfig()$shell))
    shiny::isolate({
        bus$state_manager <- manager
        scan <- integrationScanProject(paths$integration_dir)
        bus$integration_import <- list(views = scan$views$view, manifests = integrationStateManifests(integrationImportedObject(manager)))
    })
    bus
}

# The match that undoes the renaming: every transcriptomics run back to its corpus name.
intMatchUndo <- function(object, by = c("run", "name")) {
    by <- match.arg(by)
    template <- integrationSampleMatchTemplate(object)
    if (by == "run") {
        rows <- template$omic == "transcriptomics"
        template$biological_sample[rows] <- sub("^S-", "", template$saved_biological_sample[rows])
        return(template)
    }
    names <- unique(template$saved_biological_sample[template$omic == "transcriptomics"])
    data.frame(omic = "transcriptomics", run = "", saved_biological_sample = names, biological_sample = sub("^S-", "", names),
        stringsAsFactors = FALSE)
}

test_that("KIT-063 A: without a table, or with the table as downloaded, the parameters and the model are those of before", {
    project <- intMatchProject()
    manager <- intMatchStore(project$paths)
    object <- integrationImportedObject(manager)
    plain <- integrationResolveMofaParameters(object)
    expect_identical(names(plain), integrationMofaParameterNames())
    expect_identical(integrationResolveMofaParameters(object, list(sample_match_table = NULL)), plain)
    # The downloaded table, one row per run, uploaded unchanged: nothing is matched differently, nothing is added.
    template <- integrationSampleMatchTemplate(object)
    expect_identical(names(template), c("omic", "run", "saved_biological_sample", "biological_sample"))
    expect_identical(nrow(template), nrow(object@design_matrix))
    expect_identical(template$saved_biological_sample, template$biological_sample)
    unchanged <- integrationResolveSampleMatch(object, template)
    expect_null(unchanged$table)
    expect_identical(unchanged$report$runs_changed, 0L)
    expect_identical(unchanged$report$runs_mapped, nrow(template))
    expect_identical(integrationResolveMofaParameters(object, list(sample_match_table = template)), plain)
    # The declaration and its identity: the registry's, unchanged.
    expect_identical(integrationMofaDeclaration(plain), workflowOperationEntry(INTEGRATION_MODEL_OPERATION)$declaration)
    expect_identical(integrationMofaRecordedNames(plain), integrationMofaParameterNames())
    first <- runIntegrationModel(manager, trainer = intMatchTrainer)
    expect_identical(first$status, "computed")
    again <- runIntegrationModel(manager, list(sample_match_table = template), trainer = intMatchTrainer)
    expect_identical(again$status, "hit")
    expect_identical(again$reference$entity_id, first$reference$entity_id)
    expect_identical(names(manager$.context$snapshot$exports$integration_mofa_model$parameters), integrationMofaParameterNames())
    # The preview with the unchanged table equals the preview without one.
    views <- names(object@views)
    without <- integrationAlignmentPreview(object, views)
    with <- integrationAlignmentPreview(object, views, template)
    expect_identical(with$parameters, without$parameters)
    expect_identical(with$aligned, without$aligned)
    expect_identical(with$match$report$runs_changed, 0L)
})

test_that("KIT-063 A: an omic that names its samples differently is matched by the uploaded table, by run or by name", {
    corpus <- intMatchProject()
    reference <- integrationImportedObject(intMatchStore(corpus$paths))
    expected <- integrationAlignViews(reference, integrationResolveMofaParameters(reference))
    project <- intMatchProject(rename = TRUE)
    manager <- intMatchStore(project$paths)
    object <- integrationImportedObject(manager)
    # As saved, the renamed transcriptomics samples do not match the others: none is shared by every data set.
    expect_error(integrationAlignViews(object, integrationResolveMofaParameters(object)), class = "integration_no_shared_samples")
    by_run <- integrationResolveSampleMatch(object, intMatchUndo(object, "run"))
    by_name <- integrationResolveSampleMatch(object, intMatchUndo(object, "name"))
    runs <- sum(object@design_matrix$omic == "transcriptomics")
    expect_identical(by_run$report$runs_changed, runs)
    expect_identical(length(by_run$table$run), runs)
    expect_identical(unique(by_run$table$omic), "transcriptomics")
    expect_false(any(startsWith(by_run$table$biological_sample, "S-")))
    # One canonical match whichever way it was written.
    expect_identical(by_name$table, by_run$table)
    expect_identical(by_name$sha256, by_run$sha256)
    expect_identical(by_run$sha256, digest::digest(by_run$table, algo = "sha256"))
    parameters <- integrationResolveMofaParameters(object, list(sample_match_table = intMatchUndo(object, "name")))
    expect_identical(names(parameters), c(integrationMofaParameterNames(), "sample_match", "sample_match_table"))
    expect_identical(parameters$sample_match, by_run$sha256)
    aligned <- integrationAlignViews(object, parameters)
    expect_equal(aligned$views, expected$views, tolerance = 1e-12)
    expect_identical(aligned$samples, expected$samples)
    expect_identical(aligned$view_summary, expected$view_summary)
    # The canonical form is accepted as it is (the record's and the worker's), and its checksum is checked.
    expect_identical(integrationResolveMofaParameters(object, list(sample_match_table = by_run$table,
        sample_match = by_run$sha256)), parameters)
    expect_error(integrationResolveMofaParameters(object, list(sample_match_table = by_run$table,
        sample_match = strrep("0", 64L))), class = "integration_parameters_invalid")
    expect_error(integrationResolveMofaParameters(object, list(sample_match = by_run$sha256)),
        class = "integration_parameters_invalid")
    # The table offered for download with the match applied, and the match saved beside the design.
    current <- integrationSampleMatchTemplate(object, by_run$table)
    expect_identical(current$biological_sample[current$omic == "transcriptomics"],
        sub("^S-", "", current$saved_biological_sample[current$omic == "transcriptomics"]))
    saved <- utils::read.delim(integrationSaveSampleMatch(by_run$table, project$paths$multiomic_inputs_dir),
        stringsAsFactors = FALSE)
    expect_identical(as.list(saved), by_run$table)
})

test_that("KIT-063 A: a run mapped twice or left without a sample is refused; unknown rows are counted and ignored", {
    project <- intMatchProject()
    object <- integrationImportedObject(intMatchStore(project$paths))
    texts <- integrationGuiTexts()
    template <- integrationSampleMatchTemplate(object)
    twice <- rbind(template, template[1L, ])
    expect_error(integrationResolveSampleMatch(object, twice), class = "integration_sample_match_invalid",
        regexp = "mapped more than once")
    # A run row and a name row that both reach one run.
    name_row <- data.frame(omic = template$omic[[1L]], run = "", saved_biological_sample = template$saved_biological_sample[[1L]],
        biological_sample = "X", stringsAsFactors = FALSE)
    expect_error(integrationResolveSampleMatch(object, rbind(template, name_row)), class = "integration_sample_match_invalid")
    empty <- template
    empty$biological_sample[[3L]] <- ""
    expect_error(integrationResolveSampleMatch(object, empty), class = "integration_sample_match_invalid",
        regexp = "give no biological sample")
    expect_error(integrationResolveSampleMatch(object, template[c("omic", "run")]), class = "integration_sample_match_invalid",
        regexp = "needs the columns")
    unknown <- rbind(template, data.frame(omic = c("glycomics", "proteomics", "proteomics", "proteomics"),
        run = c("g1", "no-such-run", "", ""), saved_biological_sample = c("", "", "no-such-name", ""),
        biological_sample = c("A", "B", "C", "D"), stringsAsFactors = FALSE))
    report <- integrationResolveSampleMatch(object, unknown)$report
    expect_identical(report[c("unknown_omics", "unknown_runs", "unknown_names", "rows_empty", "runs_changed")],
        list(unknown_omics = 1L, unknown_runs = 1L, unknown_names = 1L, rows_empty = 1L, runs_changed = 0L))
    expect_identical(integrationSampleMatchFacts(report, texts), c(texts[["match_unchanged"]],
        sprintf(texts[["match_unknown_omics"]], 1L), sprintf(texts[["match_unknown_runs"]], 1L),
        sprintf(texts[["match_unknown_names"]], 1L), sprintf(texts[["match_rows_empty"]], 1L)))
    # A biological sample given two groups: refused by the alignment, as without a table.
    conflict <- template
    wt <- which(conflict$omic == "proteomics" & startsWith(conflict$saved_biological_sample, "WT_"))[[1L]]
    conflict$biological_sample[[wt]] <- "KO_1"
    expect_error(integrationAlignmentPreview(object, names(object@views), conflict), class = "integration_group_conflict")
    # The refusal in plain words.
    refusal <- tryCatch(integrationResolveSampleMatch(object, twice), error = function(error) error)
    expect_identical(integrationRefusalText(refusal, texts)$sentence, texts[["refusal_sample_match_invalid"]])
})

test_that("KIT-063 A: the match is part of the model's identity, recorded with the model and replayed", {
    project <- intMatchProject(rename = TRUE)
    manager <- intMatchStore(project$paths)
    object <- integrationImportedObject(manager)
    match <- intMatchUndo(object, "run")
    run <- runIntegrationModel(manager, list(sample_match_table = match), trainer = intMatchTrainer)
    expect_identical(run$status, "computed")
    checked <- integrationResolveSampleMatch(object, match)
    export <- manager$.context$snapshot$exports$integration_mofa_model
    expect_identical(names(export$parameters), c(integrationMofaParameterNames(), "sample_match", "sample_match_table"))
    expect_identical(export$parameters$sample_match, checked$sha256)
    expect_identical(export$parameters$sample_match_table, checked$table)
    expect_identical(nrow(run$value$samples), 20L)
    # Another match (one run moved to another biological sample of its group) is another model.
    other <- match
    design <- object@design_matrix
    group <- design$group[match(paste(other$omic, other$run), paste(design$omic, design$run))]
    rows <- which(other$omic == "transcriptomics")
    first <- rows[[1L]]
    target <- rows[group[rows] == group[[first]] & other$biological_sample[rows] != other$biological_sample[[first]]][[1L]]
    other$biological_sample[[first]] <- other$biological_sample[[target]]
    moved <- runIntegrationModel(manager, list(sample_match_table = other), trainer = intMatchTrainer)
    expect_false(identical(moved$reference$entity_id, run$reference$entity_id))
    rerun <- runIntegrationModel(manager, list(sample_match_table = match), trainer = intMatchTrainer)
    expect_identical(rerun$status, "hit")
    # The record lists the match; replayed from the record on another store of the same import, the same model.
    json <- file.path(project$root, "analysis_record.json")
    writeWorkflowAnalysisRecord(manager, json)
    record <- jsonlite::read_json(json, simplifyVector = FALSE)
    analysis <- record$analyses$integration_mofa_model
    expect_true("sample_match" %in% vapply(analysis$parameters, function(p) p$name, ""))
    replayed_store <- intMatchStore(project$paths, "replay")
    testthat::local_mocked_bindings(integrationTrainMofa = intMatchTrainer,
        integrationMofaCapability = function() list(available = TRUE, reason = NA_character_), .package = "MultiScholaR")
    replayed <- replayIntegrationModel(replayed_store, analysis)
    replayed_export <- replayed_store$.context$snapshot$exports$integration_mofa_model
    expect_identical(replayed_export$parameters$sample_match, checked$sha256)
    expect_identical(replayed_export$parameters$sample_match_table, checked$table)
    expect_identical(replayed$value$factors, run$value$factors)
    expect_identical(replayed$value$samples, run$value$samples)
})

test_that("KIT-063 A: the Sample Matching tab offers the table, takes the edited one back and confirms it for the model", {
    project <- intMatchProject(rename = TRUE)
    manager <- intMatchStore(project$paths)
    texts <- integrationGuiTexts()
    bus <- intMatchBus(manager, project$paths)
    object <- integrationImportedObject(manager)
    upload <- file.path(project$root, "edited.tsv")
    utils::write.table(intMatchUndo(object, "run"), upload, sep = "\t", quote = FALSE, row.names = FALSE)
    ui <- as.character(htmltools::renderTags(mod_integration_alignment_ui("a"))$html)
    for (id in c("integration-alignment-match-download", "integration-alignment-match-upload", "integration-alignment-match-clear")) {
        expect_match(ui, sprintf('data-testid="%s"', id), fixed = TRUE)
    }
    expect_match(ui, texts[["match_intro"]], fixed = TRUE)
    shiny::testServer(mod_integration_alignment_server, args = list(workflow_data = bus, paths = project$paths), {
        # As saved: the renamed samples share nothing with the other omics.
        expect_s3_class(session$returned(), "integration_alignment_refused")
        session$setInputs(match_upload = data.frame(name = "edited.tsv", size = file.size(upload), type = "text/tab-separated-values",
            datapath = upload, stringsAsFactors = FALSE))
        preview <- session$returned()
        expect_false(inherits(preview, "integration_alignment_refused"))
        runs <- sum(object@design_matrix$omic == "transcriptomics")
        expect_identical(preview$match$report$runs_changed, runs)
        status <- output$match_status$html
        expect_match(status, texts[["match_active"]], fixed = TRUE)
        expect_match(status, sprintf(texts[["match_changed"]], runs), fixed = TRUE)
        expect_identical(preview$aligned$shared_samples, 19L)
        session$setInputs(confirm = 1)
        confirmed <- bus$integration_alignment
        expect_identical(confirmed$parameters$sample_match, preview$parameters$sample_match)
        expect_true(file.exists(confirmed$sample_match_file))
        expect_match(paste(integrationModelMatchingSettings(confirmed$parameters), collapse = "\n"),
            sprintf(texts[["setting_sample_match"]], runs), fixed = TRUE)
        # A broken table is refused in words; the match in use stays.
        broken <- file.path(project$root, "broken.tsv")
        writeLines(c("omic\trun", "transcriptomics\tx"), broken)
        session$setInputs(match_upload = data.frame(name = "broken.tsv", size = file.size(broken), type = "text/plain",
            datapath = broken, stringsAsFactors = FALSE))
        expect_match(output$match_status$html, sprintf(texts[["match_refused"]], texts[["refusal_sample_match_invalid"]]),
            fixed = TRUE)
        expect_identical(session$returned()$parameters$sample_match, preview$parameters$sample_match)
        # Back to the names as saved: the confirmation is cleared, the samples no longer match.
        session$setInputs(match_clear = 1)
        expect_null(bus$integration_alignment)
        expect_s3_class(session$returned(), "integration_alignment_refused")
    })
    # The Model tab trains on the confirmed match.
    bus2 <- intMatchBus(manager, project$paths)
    shiny::testServer(mod_integration_alignment_server, args = list(workflow_data = bus2, paths = project$paths), {
        session$setInputs(match_upload = data.frame(name = "edited.tsv", size = file.size(upload), type = "text/tab-separated-values",
            datapath = upload, stringsAsFactors = FALSE))
        session$setInputs(confirm = 1)
    })
    runs <- sum(object@design_matrix$omic == "transcriptomics")
    shiny::testServer(mod_integration_model_server, args = list(workflow_data = bus2, paths = project$paths,
        probe = function() list(available = TRUE, reason = NA_character_), trainer = intMatchTrainer,
        withProgressFn = function(message, value, expr) force(expr)), {
        session$setInputs(run = 1)
        model <- bus2$integration_model
        expect_identical(model$status, "computed", info = output$run_status$html)
        expect_identical(model$parameters$sample_match, shiny::isolate(bus2$integration_alignment$parameters$sample_match))
        expect_match(output$matching_settings$html, sprintf(texts[["setting_sample_match"]], runs), fixed = TRUE)
    })
})
