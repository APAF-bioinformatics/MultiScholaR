# KIT-142: the integration analysis record, its replay and its report, without MOFA2. The project is the KIT-137 tier-1
# corpus copied into `<project>/integration/<omic>/`, with two changes written by the contract's writer: the
# transcriptomics samples renamed (`S-<name>`, so an uploaded match table is needed, KIT-063 item A) and the proteomics
# input given a corrected matrix with the groups as replicate sets (KIT-153: the model takes the uncorrected matrices by
# default, so the claim limit of KIT-060 D10 does not apply and the uncorrected matrices' caution does). The model is
# the stored double-precision reference (tests/testdata/integration/mofa_reference) returned by a stand-in trainer that
# reports its own software versions, so that the record and the methods paragraph can be seen to read them. The
# journey runs through the tabs' server logic: import, an uploaded match confirmed, training, results saved, record
# saved. Declared before the run:
#   - the record holds every input (omic, data sets, manifest checksum equal to the file's, matrix role, correction),
#     the confirmed match (the uploaded table, the saved sample table's checksum), every setting with its status
#     (decided, provisional default, set by the user, from the data; the uploaded match is "user_set", not "unknown"),
#     the versions, the model's checksum, the saved results' settings and files (checksums equal to the files'), MOFA2's
#     warnings in words and the claim limit;
#   - the methods paragraph and the study parameters read every value from the record (a record with other versions or
#     settings gives other words; none is hard-coded), the same from the record as built and as read back from JSON;
#   - the record, its workbook, the study parameters and the report are written through a short folder and checked: a
#     report the renderer did not write stops the save by name;
#   - reopened and replayed from its record, the analysis gives the same model (largest relative difference 0 on one
#     computer) and the same results; a changed or missing input, MOFA2 absent and other software versions are said in
#     plain words (refusals; the versions a warning).
# The same journey with MOFA2 training in the job, the report rendered, is test-integration-record-mofa.R (MOFA2 lane).

intRecData <- function() normalizePath(testthat::test_path("..", "testdata", "integration"), mustWork = FALSE)
intRecOmics <- c("lipidomics", "metabolomics", "proteomics", "transcriptomics")
intRecVersions <- list(R = "4.6.1", MOFA2 = "9.8.7", reticulate = "1.2.3", python = "3.12.16 (main, stand-in)",
    mofapy2 = "6.5.4", numpy = "1", scipy = "1", h5py = "1")
intRecOverallLevel <- paste("WARNING: Factor(s) 2 are strongly correlated with the average expression of features for",
    "at least one of your omics. Such factors appear when there are differences in the total 'levels' between your",
    "samples, *sometimes* because of poor normalisation in the preprocessing steps.")

intRecStored <- function(name) utils::read.delim(file.path(intRecData(), "mofa_reference", paste0(name, ".tsv")),
    check.names = FALSE, stringsAsFactors = FALSE)

# The stored reference as `integrationTrainMofa` returns a model, with the given versions and warnings.
intRecTrainer <- function(versions = intRecVersions, messages = intRecOverallLevel) {
    function(views, options) {
        factors <- intRecStored("factors")
        f <- as.matrix(factors[-1L])
        rownames(f) <- factors$biological_sample
        f <- f[colnames(views[[1L]]), , drop = FALSE]
        w <- intRecStored("weights")
        weights <- lapply(stats::setNames(names(views), names(views)), function(view) {
            rows <- w[w$view == view, , drop = FALSE]
            m <- as.matrix(rows[-(1:2)])
            rownames(m) <- paste(view, rows$feature_id, sep = ":")
            m[rownames(views[[view]]), , drop = FALSE]
        })
        v <- intRecStored("variance_explained")
        r2 <- matrix(v$r2_percent, nrow = length(unique(v$factor)), dimnames = list(unique(v$factor), unique(v$view)))
        r2 <- r2[, names(views), drop = FALSE]
        list(factors = f, weights = weights, r2_per_factor = r2, r2_total = colSums(r2), elbo = c(-2, -1), iterations = 186L,
            options = options, versions = versions, messages = messages,
            model_file = list(bytes = as.raw(1:3), sha256 = strrep("d", 64L), size = 3))
    }
}

# Rewrite one omic's input with the contract's writer: its samples renamed, or its matrix marked as corrected.
intRecRewrite <- function(dir, rename = FALSE, correction = NULL) {
    input <- readIntegrationInputWithMatrices(dir)
    read <- input$manifest
    manifest <- list(omic_type = read$omic_type, workflow = read$workflow,
        matrix_role = read$matrix_role, scale = read$scale, correction = correction,
        source_state = list(name = "kit142"), sample_id_column = read$sample_id_column, group_column = read$group_column,
        technical_replicate_column = read$technical_replicate_column, biological_sample_column = read$biological_sample_column)
    samples <- input$samples
    if (rename) samples[[manifest$biological_sample_column]] <- paste0("S-", samples[[manifest$biological_sample_column]])
    # KIT-153: a correction comes with its corrected matrix (here the same values), beside the normalised one.
    matrices <- lapply(input$views, `[[`, "matrix")
    # KIT-184: the corrected matrix equals the normalised one here, so no removed factor can be determined; the writer's
    # classed warning is expected and silenced (test-integration-removed-factors.R, K3, asserts it).
    suppressWarnings(writeIntegrationInput(dir, manifest, samples, matrices, lapply(input$views, `[[`, "features"),
        corrected = if (is.null(correction)) list() else matrices), classes = "integration_input_removed_factors_undetermined")
}

intRecProject <- function(.local_envir = parent.frame()) {
    testthat::skip_if_not(dir.exists(file.path(intRecData(), "corpus")), "the tier-1 corpus is not present")
    testthat::skip_if_not(workflowOwnerCapable(), "the integration store needs qs2, DuckDB and Arrow")
    root <- withr::local_tempdir(.local_envir = .local_envir)
    results_dir <- file.path(root, "results", "integration")
    paths <- list(base_dir = root, integration_dir = file.path(root, "integration"),
        source_dir = file.path(root, "scripts", "integration"), results_dir = results_dir,
        multiomic_inputs_dir = file.path(results_dir, "multiomic_inputs"), mofa_tables_dir = file.path(results_dir, "mofa", "tables"),
        mofa_plots_dir = file.path(results_dir, "mofa", "plots"),
        study_report_dir = file.path(root, "results_summary", "integration", "Study_report"))
    dir.create(paths$source_dir, recursive = TRUE)
    for (omic in intRecOmics) {
        dir <- file.path(paths$integration_dir, omic)
        dir.create(dir, recursive = TRUE)
        stopifnot(all(file.copy(list.files(file.path(intRecData(), "corpus", omic), full.names = TRUE), dir)))
    }
    intRecRewrite(file.path(paths$integration_dir, "transcriptomics"), rename = TRUE)
    intRecRewrite(file.path(paths$integration_dir, "proteomics"), correction = list(state = "ruv_corrected",
        method = "RUVIII_C_Varying", replicate_structure = "group", grouping_variable = "group", k = 2L))
    list(root = root, paths = paths)
}

# A stand-in report: a file holding the methods paragraph the record gives.
intRecRender <- function(record_path, output_dir, format) {
    path <- file.path(output_dir, paste0("integration_report.", format))
    writeLines(c("stand-in report", integrationMethodsText(jsonlite::read_json(record_path, simplifyVector = FALSE))), path)
    path
}

# The journey through the tabs: import, the uploaded match confirmed, training, results saved, record saved.
intRecJourney <- function(project, renderFn = intRecRender, .local_envir = parent.frame()) {
    paths <- project$paths
    bus <- do.call(shiny::reactiveValues, omicWorkflowInitialBus(integrationGuiConfig()$shell))
    withr::defer(if (!is.null(shiny::isolate(bus$state_manager))) closeWorkflowStateOwner(shiny::isolate(bus$state_manager)),
        envir = .local_envir)
    shiny::testServer(mod_integration_import_server, args = list(workflow_data = bus, paths = paths), {
        session$flushReact()
        session$setInputs(views = bus$integration_scan$views$view, import = 1)
    })
    manager <- shiny::isolate(bus$state_manager)
    edited <- integrationSampleMatchTemplate(integrationImportedObject(manager))
    rows <- edited$omic == "transcriptomics"
    edited$biological_sample[rows] <- sub("^S-", "", edited$saved_biological_sample[rows])
    upload <- file.path(project$root, "edited.tsv")
    utils::write.table(edited, upload, sep = "\t", quote = FALSE, row.names = FALSE)
    shiny::testServer(mod_integration_alignment_server, args = list(workflow_data = bus, paths = paths), {
        session$setInputs(match_upload = data.frame(name = "edited.tsv", size = file.size(upload),
            type = "text/tab-separated-values", datapath = upload, stringsAsFactors = FALSE))
        session$setInputs(confirm = 1)
    })
    shiny::testServer(mod_integration_model_server, args = list(workflow_data = bus, paths = paths, trainer = intRecTrainer(),
        probe = function() list(available = TRUE, reason = NA_character_),
        withProgressFn = function(message, value, expr) force(expr)), {
        session$setInputs(run = 1)
    })
    shiny::testServer(mod_integration_results_server, args = list(workflow_data = bus, paths = paths), {
        session$flushReact()
        session$setInputs(save = 1)
    })
    shiny::testServer(mod_integration_record_server, args = list(workflow_data = bus, paths = paths, renderFn = renderFn,
        withProgressFn = function(message, value, expr) force(expr)), {
        session$setInputs(save = 1)
    })
    list(bus = bus, manager = manager, paths = paths, runs_moved = sum(rows))
}

intRecSha <- function(path) digest::digest(file = path, algo = "sha256")

test_that("KIT-142: the record holds every input, the confirmed match, each setting's status, the model and the saved results", {
    project <- intRecProject()
    journey <- intRecJourney(project)
    manager <- journey$manager
    bus <- journey$bus
    expect_false(is.null(shiny::isolate(bus$integration_results)))
    record <- buildWorkflowAnalysisRecord(manager)
    expect_true(validateWorkflowAnalysisRecord(record)$valid)
    analysis <- record$analyses$integration_mofa_model
    block <- analysis$integration
    expect_identical(block$schema, "multischolar.integration_record.v1")
    expect_true(block$complete)
    # Every input: omic, data sets, manifest checksum (the file's), matrix role, correction.
    inputs <- block$inputs
    expect_identical(vapply(inputs, `[[`, "", "omic"), intRecOmics)
    for (input in inputs) {
        manifest <- file.path(journey$paths$integration_dir, input$omic, "integration_input.json")
        expect_identical(input$manifest_sha256, intRecSha(manifest), info = input$omic)
        expect_identical(input$views_in_model, input$views, info = input$omic)
    }
    proteomics <- inputs[[3L]]
    # KIT-153: the proteomics input holds a corrected matrix; the model took the uncorrected one (the default).
    expect_identical(proteomics$matrix_role, "normalised")
    expect_identical(proteomics$matrix_used, "uncorrected")
    expect_true(proteomics$corrected_available)
    expect_identical(proteomics$correction$method, "RUVIII_C_Varying")
    expect_false(proteomics$group_corrected)
    expect_false(inputs[[4L]]$group_corrected)
    expect_false(inputs[[4L]]$corrected_available)
    expect_identical(inputs[[4L]]$matrix_role, "normalised")
    # The confirmed match: the uploaded table, the saved sample table and its checksum.
    matching <- block$matching
    expect_identical(matching$rule, "user_table")
    expect_identical(matching$runs_changed, journey$runs_moved)
    export <- manager$.context$snapshot$exports$integration_mofa_model
    expect_identical(matching$table, export$parameters$sample_match_table)
    expect_identical(matching$sample_match_sha256, export$parameters$sample_match)
    confirmed <- shiny::isolate(bus$integration_alignment)
    expect_identical(matching$confirmation$design_sha256, intRecSha(confirmed$file))
    expect_identical(matching$confirmation$sample_match_sha256, intRecSha(confirmed$sample_match_file))
    expect_identical(nrow(matching$samples), 20L)
    # Each setting's status; the uploaded match is the user's (its provenance was "unknown").
    status <- stats::setNames(vapply(analysis$parameters, `[[`, "", "status"), vapply(analysis$parameters, `[[`, "", "name"))
    provenance <- stats::setNames(vapply(analysis$parameters, `[[`, "", "provenance"), names(status))
    expect_identical(unname(provenance[c("sample_match", "sample_match_table")]), c("user_set", "user_set"))
    expect_false("unknown" %in% provenance)
    expect_identical(unname(status[c("input_manifests", "views")]), c("input", "input"))
    expect_identical(status[["precision"]], "decided")
    expect_identical(status[["input_matrix"]], "decided")
    expect_identical(unname(status[c("num_factors", "seed", "technical_replicates")]), rep("provisional_default", 3L))
    # KIT-151: decided by the user (INT-12, INT-15, INT-24).
    expect_identical(unname(status[c("feature_selection", "convergence_mode", "scaling", "scale_views")]), rep("decided", 4L))
    # The model: versions, checksum, training.
    model <- block$model
    expect_identical(model$versions, intRecVersions)
    expect_identical(model$model_file_sha256, strrep("d", 64L))
    expect_identical(model$entity_id, export$reference$entity_id)
    expect_true(model$converged)
    expect_identical(model$stopped_by, "convergence_rule")
    expect_identical(c(model$factors_kept, model$factors_requested), c(5L, 5L))
    expect_identical(model$precision, "double")
    # The saved results: their settings with a status, every file with the file's checksum.
    results <- block$results
    expect_identical(results$label, "integration-results")
    expect_identical(vapply(results$parameters, `[[`, "", "name"), names(integrationResultsParameterSpec()))
    expect_identical(unique(vapply(results$parameters, `[[`, "", "status")), "provisional_default")
    saved <- shiny::isolate(bus$integration_results)
    expect_identical(results$tables$file, saved$tables)
    expect_identical(results$tables$sha256, unname(vapply(file.path(saved$tables_dir, saved$tables), intRecSha, "")))
    expect_identical(results$figures$sha256, unname(vapply(file.path(saved$plots_dir, saved$figures), intRecSha, "")))
    expect_true(all(results$association$status == "tested"))
    # MOFA2's warning in words; KIT-153: no claim limit on the uncorrected matrices, and their caution.
    expect_true("overall_level" %in% vapply(block$advisories, `[[`, "", "kind"))
    expect_false(block$claim_limit$applies)
    expect_length(block$claim_limit$group_corrected, 0L)
    expect_true(block$claim_limit$uncorrected)
    expect_true(integrationGuiTexts()[["claim_uncorrected"]] %in% block$claim_limit$text)
    # The omics' record fields are unchanged: the analysis keeps its registry fields.
    expect_identical(analysis$operation, "integration_mofa_model")
    expect_identical(analysis$outcomes$model_sha256, strrep("d", 64L))
})

test_that("KIT-142: the methods paragraph and the study parameters read the record's values, as built and from JSON", {
    project <- intRecProject()
    journey <- intRecJourney(project)
    manager <- journey$manager
    record <- buildWorkflowAnalysisRecord(manager)
    json <- file.path(project$root, "record.json")
    writeWorkflowAnalysisRecord(manager, json)
    read <- jsonlite::read_json(json, simplifyVector = FALSE)
    t <- integrationRecordTexts()
    methods <- integrationMethodsText(record)
    expect_identical(integrationMethodsText(read), methods)
    expect_identical(integrationStudyParameterLines(read), integrationStudyParameterLines(record))
    match_sha <- record$analyses$integration_mofa_model$integration$matching$sample_match_sha256
    for (fragment in c(
        "6 data sets from 4 omics were combined",
        "proteomics (proteomics; normalised, before correction for unwanted variation; the omic's corrected values (applied: the approach of earlier versions (RUV-III-C; replicate sets: ",
        "transcriptomics (transcriptomics; normalised, not corrected for unwanted variation)",
        "MOFA2 9.8.7 (mofapy2 6.5.4, Python 3.12.16, R 4.6.1), in double precision",
        sprintf(t[["methods_matching_user"]], journey$runs_moved, match_sha),
        t[["methods_replicates_average"]], t[["methods_missing_keep_missing"]],
        # KIT-151 (INT-12): the default selection, its observed share and its rule.
        sprintf(t[["methods_sparse_selected"]], 10L, 20L, "0.5"),
        paste0(sprintf(t[["methods_selection_variance_above_trend"]], 1000L, 1000L), ", ", t[["methods_scaling_feature_zscore"]]),
        # KIT-151 (INT-15): slow, a limit of 10,000 rounds, and how training stopped.
        "asked for 5 factors (convergence mode slow, at most 10000 training rounds, seed 20261119",
        "Training settled after 186 rounds, stopped by its convergence rule, and kept 5 factors; 20 biological samples",
        # KIT-151 (INT-14): five models, one per seed; the stand-in trains the same factors from every seed.
        sprintf(t[["methods_seeds"]], 5L, 20261119L, 20261123L, 20261119L,
            paste(sprintf("Factor%d 1.000", 1:5), collapse = ", ")),
        sprintf(t[["methods_association"]], "Benjamini-Hochberg"),
        t[["methods_claim_none"]], t[["methods_claim_uncorrected"]],
        "Settings that were provisional defaults at the time of the analysis: ", "Number of factors (5)")) {
        expect_true(grepl(fragment, methods, fixed = TRUE), info = fragment)
    }
    # Nothing hard-coded: the versions are the record's, and precision is decided, so not listed as provisional.
    expect_false(grepl("1.22.0", methods, fixed = TRUE))
    expect_false(grepl("Numeric precision (", methods, fixed = TRUE))
    expect_false(grepl(t[["methods_user_set"]] |> sub(pattern = "%s.", replacement = "", fixed = TRUE), methods, fixed = TRUE))
    lines <- integrationStudyParameterLines(record)
    expect_identical(lines[[1L]], t[["report_title"]])
    expect_true(all(unname(t[c("section_inputs", "section_matching", "section_settings", "section_model", "section_notes",
        "section_methods")]) %in% lines))
    expect_identical(lines[[length(lines) - 1L]], methods)
    # Other settings give other words: a model with 3 factors and seed 7 (its results not saved yet).
    confirmed <- shiny::isolate(journey$bus$integration_alignment)
    # KIT-151 (INT-17): a model trained with the maths library pinned says so.
    pin <- list(applied = TRUE, variables = as.list(INTEGRATION_MATHS_PIN), reason = NA_character_)
    runIntegrationModel(manager, utils::modifyList(confirmed$parameters, list(num_factors = 3L, seed = 7L)), trainer = intRecTrainer(
        versions = utils::modifyList(intRecVersions, list(MOFA2 = "1.0.0", maths_pin = pin))))
    other <- integrationMethodsText(buildWorkflowAnalysisRecord(manager))
    expect_true(grepl(sprintf(t[["methods_pin_applied"]],
        "OPENBLAS_NUM_THREADS=1, OMP_NUM_THREADS=1, MKL_NUM_THREADS=1, OPENBLAS_CORETYPE=Haswell"), other, fixed = TRUE))
    expect_false(grepl(sprintf(t[["methods_pin_applied"]], ""), methods, fixed = TRUE))
    expect_true(grepl("asked for 3 factors", other, fixed = TRUE))
    expect_true(grepl("seed 7;", other, fixed = TRUE))
    expect_true(grepl("MOFA2 1.0.0 (", other, fixed = TRUE))
    expect_true(grepl(sprintf(t[["methods_user_set"]], "Number of factors (3); Random seed (7)"), other, fixed = TRUE))
    expect_false(grepl(sprintf(t[["methods_association"]], "Benjamini-Hochberg"), other, fixed = TRUE))
})

test_that("KIT-142: the record, workbook, study parameters and report are saved through a short folder and checked", {
    project <- intRecProject()
    journey <- intRecJourney(project)
    bus <- journey$bus
    saved <- shiny::isolate(bus$integration_record)
    expect_identical(shiny::isolate(bus$tab_status$integration_record), "complete")
    expect_identical(basename(saved$files), c("analysis_record.json", "analysis_record.xlsx", "study_parameters.txt",
        "integration_report.docx", "model.hdf5", "sample_match.tsv"))
    expect_identical(normalizePath(dirname(saved$files[[1L]]), winslash = "/"),
        normalizePath(journey$paths$study_report_dir, winslash = "/"))
    expect_true(all(file.size(saved$files) > 0))
    read <- jsonlite::read_json(saved$files[[1L]], simplifyVector = FALSE)
    # KIT-148: the model file (the stand-in trainer's three bytes) with the record's checksum, and the uploaded match.
    expect_identical(readBin(saved$files[[5L]], "raw", 16L), as.raw(1:3))
    export <- journey$manager$.context$snapshot$exports$integration_mofa_model
    match <- utils::read.delim(saved$files[[6L]], colClasses = "character", stringsAsFactors = FALSE)
    expect_identical(as.list(match), export$parameters$sample_match_table)
    expect_identical(nrow(match), journey$runs_moved)
    expect_identical(digest::digest(as.list(match), algo = "sha256"), export$parameters$sample_match)
    expect_identical(saved$methods, integrationMethodsText(read))
    expect_identical(readLines(saved$files[[3L]]), integrationStudyParameterLines(read))
    expect_identical(readLines(saved$files[[4L]])[[2L]], saved$methods)
    expect_true(all(INTEGRATION_RECORD_SHEETS %in% openxlsx::getSheetNames(saved$files[[2L]])))
    settings <- openxlsx::read.xlsx(saved$files[[2L]], sheet = "integration_settings", check.names = FALSE)
    expect_true(all(c("provisional default", "decided", "set by the user", "from the data") %in%
        settings[[integrationRecordTexts()[["col_status"]]]]))
    # A report the renderer did not write stops the save by name; nothing is lost silently.
    nothing <- function(record_path, output_dir, format) invisible(NULL)
    directory <- file.path(project$root, "again")
    expect_error(integrationWriteRecordBundle(journey$manager, directory, renderFn = nothing), class = "workflow_figure_unwritten")
    # Without a trained model there is nothing to record.
    empty <- newIntegrationStateOwner(file.path(project$root, "empty"))
    withr::defer(closeWorkflowStateOwner(empty))
    importIntegrationInputs(empty, file.path(journey$paths$integration_dir, intRecOmics))
    expect_error(integrationWriteRecordBundle(empty, directory, renderFn = intRecRender), class = "integration_record_no_model")
    # The omics' record tables carry no integration sheet.
    expect_length(integrationRecordTables(list(omic_type = "proteomics", analyses = list())), 0L)
})

test_that("KIT-142: the report renders from the record with the package's template", {
    testthat::skip_if_not(rmarkdown::pandoc_available(), "pandoc is not available")
    project <- intRecProject()
    journey <- intRecJourney(project)
    bundle <- integrationWriteRecordBundle(journey$manager, file.path(project$root, "html"), format = "html")
    # pandoc wraps long lines in the HTML it writes, so a sentence can span lines: compare with whitespace collapsed.
    page <- gsub("\\s+", " ", paste(readLines(bundle$files[[4L]], warn = FALSE), collapse = " "))
    expect_match(page, integrationRecordTexts()[["report_title"]], fixed = TRUE)
    expect_match(page, "MOFA2 9.8.7", fixed = TRUE)
    expect_match(page, "Training settled after 186 rounds", fixed = TRUE)
})

test_that("KIT-142: reopened and replayed from its record, the same model; refusals and warnings in plain words", {
    project <- intRecProject()
    journey <- intRecJourney(project)
    manager <- journey$manager
    root <- dirname(manager$.context$store$directory)
    # KIT-148: the recorded model as held outside its store (the weights read from it before it is closed).
    recorded <- integrationModelComparable(manager, withWorkflowStateScope(manager, function()
        hydrateSavedWorkflowEntity(manager$.context$store, manager$.context$snapshot$exports$integration_mofa_model$reference)$value))
    closeWorkflowStateOwner(manager)
    t <- integrationRecordTexts()
    scratch <- file.path(project$root, "replay")
    replay <- suppressWarnings(replayIntegrationProject(root, scratch, trainer = intRecTrainer(), current_versions = intRecVersions))
    expect_true(replay$comparison$equal)
    expect_true(replay$comparison$identical)
    expect_identical(replay$comparison$largest_relative, 0)
    expect_identical(replay$replayed$factors, recorded$factors)
    expect_identical(replay$parameters$sample_match, replay$record$analyses$integration_mofa_model$integration$matching$sample_match_sha256)
    expect_true(all(replay$versions$same))
    recorded_results <- replay$record$analyses$integration_mofa_model$integration$results
    expect_identical(replay$results$association$r2_percent, recorded_results$association$r2_percent)
    expect_identical(replay$results$parameters, integrationResolveResultsParameters(list(), replay$results$factors))
    # From the record read back from JSON, the same.
    json <- file.path(project$root, "record.json")
    jsonlite::write_json(replay$record, json, auto_unbox = TRUE, digits = I(17), null = "null", na = "null")
    from_json <- suppressWarnings(replayIntegrationRecord(jsonlite::read_json(json, simplifyVector = FALSE), scratch,
        recorded = recorded, trainer = intRecTrainer(), current_versions = intRecVersions))
    expect_identical(from_json$comparison$largest_relative, 0)
    # Other software versions: a warning in words, the comparison kept.
    expect_warning(other <- replayIntegrationRecord(replay$record, scratch, recorded = recorded, trainer = intRecTrainer(),
        current_versions = utils::modifyList(intRecVersions, list(MOFA2 = "1.0.0"))), class = "integration_replay_versions_differ")
    expect_identical(other$versions$same, c(FALSE, TRUE, TRUE, TRUE))
    expect_true(other$comparison$equal)
    # MOFA2 absent: refused before anything is trained.
    absent <- tryCatch(replayIntegrationRecord(replay$record, scratch, probe = function() list(available = FALSE, reason = "none")),
        error = function(error) error)
    expect_s3_class(absent, "integration_mofa_unavailable")
    expect_identical(conditionMessage(absent), t[["replay_mofa_unavailable"]])
    # An input not where it was: refused, and accepted at its new place.
    moved <- file.path(project$root, "moved", "proteomics")
    dir.create(moved, recursive = TRUE)
    file.copy(list.files(file.path(journey$paths$integration_dir, "proteomics"), full.names = TRUE), moved)
    unlink(file.path(journey$paths$integration_dir, "proteomics"), recursive = TRUE)
    missing <- tryCatch(replayIntegrationRecord(replay$record, scratch, trainer = intRecTrainer()), error = function(error) error)
    expect_s3_class(missing, "integration_replay_input_missing")
    expect_match(conditionMessage(missing), "The saved data of proteomics is no longer at", fixed = TRUE)
    found <- suppressWarnings(replayIntegrationRecord(replay$record, scratch, recorded = recorded, trainer = intRecTrainer(),
        manifests = c(proteomics = file.path(moved, "integration_input.json"))))
    expect_identical(found$comparison$largest_relative, 0)
    # An input changed since the record (its samples renamed once more): refused.
    intRecRewrite(file.path(journey$paths$integration_dir, "transcriptomics"), rename = TRUE)
    changed <- tryCatch(replayIntegrationRecord(replay$record, scratch, trainer = intRecTrainer(),
        manifests = c(proteomics = file.path(moved, "integration_input.json"))), error = function(error) error)
    expect_s3_class(changed, "integration_replay_input_changed")
    expect_identical(conditionMessage(changed), sprintf(t[["replay_input_changed"]], "transcriptomics"))
    expect_identical(integrationRefusalText(changed)$sentence, conditionMessage(changed))
    # A record without an integration model.
    expect_error(replayIntegrationRecord(list(analyses = list()), scratch), class = "integration_replay_no_model")
})

test_that("KIT-142: the comparison tells a different model apart", {
    # KIT-148: models as the comparison holds them (`integrationModelComparable`: the weights as one table).
    a <- list(factors = data.frame(biological_sample = "x", Factor1 = 1, Factor2 = 2, stringsAsFactors = FALSE),
        weights = data.frame(view = "v", feature_id = "f", Factor1 = 0.5, Factor2 = -0.5, stringsAsFactors = FALSE),
        variance_explained = data.frame(factor = c("Factor1", "Factor2"), view = "v", r2_percent = 1, stringsAsFactors = FALSE),
        variance_total = data.frame(view = "v", r2_percent = 2, stringsAsFactors = FALSE),
        samples = data.frame(biological_sample = "x", stringsAsFactors = FALSE), view_summary = data.frame(view = "v"))
    b <- a
    b$factors$Factor1 <- 1 + 1e-9
    near <- integrationCompareModels(a, b)
    expect_true(near$equal)
    expect_false(near$identical)
    expect_equal(near$largest_relative, 1e-9 / 2, tolerance = 1e-6)
    b$factors$Factor1 <- 1.5
    expect_false(integrationCompareModels(a, b)$equal)
    b <- a
    b$factors$biological_sample <- "y"
    expect_identical(integrationCompareModels(a, b)$largest_relative, Inf)
})

test_that("KIT-246: the Methods paragraph's scaling sentence states the simulations' numbers", {
    # INT-31 corrects INT-24's evidence sentence: the z-score's sentence said "weak factors" (undefined) and "otherwise
    # made no practical difference" (KIT-151 d: inconclusive at 20 samples without missing values). One journey with
    # the default scaling; the paragraph as built and as read back from JSON. E2 and E5 (the owner's words and the
    # agreement with the Model tab's help) are added with the fix.
    old <- paste("each scaled to mean 0 and standard deviation 1 after the selection (which used the unscaled values);",
        "in simulation, with few samples and missing values this scaling improved the recovery of weak factors, and",
        "otherwise made no practical difference.")
    project <- intRecProject()
    journey <- intRecJourney(project)
    record <- buildWorkflowAnalysisRecord(journey$manager)
    json <- file.path(project$root, "record.json")
    writeWorkflowAnalysisRecord(journey$manager, json)
    read <- jsonlite::read_json(json, simplifyVector = FALSE)
    t <- integrationRecordTexts()
    # E1: every record prints the table's sentence, and not the old one.
    for (methods in list(built = integrationMethodsText(record), json = integrationMethodsText(read))) {
        expect_true(grepl(t[["methods_scaling_feature_zscore"]], methods, fixed = TRUE))
        expect_false(grepl("weak factors", methods, fixed = TRUE))
        expect_false(grepl("no practical difference", methods, fixed = TRUE))
        expect_false(grepl(old, methods, fixed = TRUE))
    }
    # G1: a record whose scaling is "none" keeps its own sentence (the value replaced in the record's encoded form).
    parameters <- read$analyses$integration_mofa_model$parameters
    hit <- vapply(parameters, \(parameter) identical(integrationRecordScalar(parameter$name), "scaling"), logical(1))
    expect_identical(sum(hit), 1L)
    parameters[[which(hit)]]$value <- encodeWorkflowRecordValue("none")
    read$analyses$integration_mofa_model$parameters <- parameters
    unscaled <- integrationMethodsText(read)
    expect_true(grepl(t[["methods_scaling_none"]], unscaled, fixed = TRUE))
    expect_true(grepl("without per-feature scaling.", unscaled, fixed = TRUE))
    expect_false(grepl(t[["methods_scaling_feature_zscore"]], unscaled, fixed = TRUE))
})
