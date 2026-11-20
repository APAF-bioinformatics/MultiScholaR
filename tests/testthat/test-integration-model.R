# KIT-062: the integration owner and the MOFA2 declaration, without Python. Two small integration inputs written by
# the contract's writer (KIT-061): proteomics (one run per biological sample, A_1..A_3, B_1, B_2; P4 observed once)
# and lipidomics (two technical replicates per biological sample, A_3 absent; ITSD1 an internal standard). MOFA2 is
# replaced by a stand-in trainer where a model is needed; the end-to-end training is test-integration-mofa-model.R.

intModelInput <- function(root, omic, views, samples, standards = list()) {
    matrices <- lapply(views, function(m) data.frame(feature_id = rownames(m), m, check.names = FALSE,
        stringsAsFactors = FALSE))
    features <- lapply(stats::setNames(names(views), names(views)), function(view) {
        table <- data.frame(feature_id = rownames(views[[view]]), display_name = tolower(rownames(views[[view]])),
            stringsAsFactors = FALSE)
        if (!is.null(standards[[view]])) table$is_internal_standard <- standards[[view]]
        table
    })
    manifest <- list(omic_type = omic, workflow = "TEST", matrix_role = "normalised", scale = "log2", correction = NULL,
        source_state = list(name = "kit062-test"), sample_id_column = "Run", group_column = "group",
        technical_replicate_column = "tech_rep_group", biological_sample_column = "tech_rep_group")
    dir <- file.path(root, "integration", omic)
    writeIntegrationInput(dir, manifest, samples, matrices, features)
    dir
}

intModelSamples <- function(runs, biological) {
    data.frame(Run = runs, group = sub("_.*$", "", biological), replicates = as.integer(sub("^.*_", "", biological)),
        tech_rep_group = biological, stringsAsFactors = FALSE)
}

intModelProteomics <- function(root, shift = 0, biological = c("A_1", "A_2", "A_3", "B_1", "B_2")) {
    runs <- paste0("p", seq_along(biological))
    values <- rbind(P1 = c(10, 11, 12, 13, 14), P2 = c(5, 7, 6, 9, 8), P3 = c(1, 1.5, 2, 2.5, 3.5),
        P4 = c(NA, NA, 4, NA, NA)) + shift
    colnames(values) <- runs
    intModelInput(root, "proteomics", list(proteomics = values), intModelSamples(runs, biological))
}

intModelLipidomics <- function(root) {
    biological <- c("A_1", "A_1", "A_2", "A_2", "B_1", "B_1", "B_2", "B_2")
    runs <- paste0("l", seq_along(biological))
    values <- rbind(L1 = c(20, 22, 21, 23, 25, 27, 26, 28), L2 = c(3, NA, 4, 5, 6, 7, NA, NA),
        ITSD1 = c(15, 15, 15, 15, 15, 15, 15, 15))
    colnames(values) <- runs
    intModelInput(root, "lipidomics", list(lipidomics_LCMS_Pos = values), intModelSamples(runs, biological),
        standards = list(lipidomics_LCMS_Pos = c(FALSE, FALSE, TRUE)))
}

intModelFakeTrainer <- function(views, options) {
    samples <- colnames(views[[1L]])
    k <- options$num_factors
    labels <- paste0("Factor", seq_len(k))
    list(factors = matrix(seq_len(length(samples) * k) / 10, length(samples), k, dimnames = list(samples, labels)),
        weights = lapply(views, function(v) matrix(0.5, nrow(v), k, dimnames = list(rownames(v), labels))),
        r2_per_factor = matrix(1, k, length(views), dimnames = list(labels, names(views))),
        r2_total = stats::setNames(rep(k, length(views)), names(views)), elbo = c(-20, -10, NaN), iterations = 2L,
        options = options, versions = list(MOFA2 = "stand-in"), messages = character(),
        model_file = list(bytes = raw(0), sha256 = NA_character_, size = 0))
}

intModelInputs <- function() {
    root <- withr::local_tempdir(.local_envir = parent.frame())
    list(root = root, proteomics = intModelProteomics(root), lipidomics = intModelLipidomics(root))
}

test_that("KIT-062: the inputs state aligns the views on the biological sample, averages runs and drops standards", {
    inputs <- intModelInputs()
    read <- integrationViewsFromInputs(c(inputs$proteomics, inputs$lipidomics))
    object <- read$object
    # The omics in name order; the manifests by omic.
    expect_identical(names(object@views), c("lipidomics_LCMS_Pos", "proteomics"))
    expect_identical(names(read$manifests), c("lipidomics", "proteomics"))
    expect_identical(unname(integrationStateManifests(object)),
        unname(vapply(read$manifests, function(path) digest::digest(file = path, algo = "sha256"), character(1))))
    expect_identical(names(object@design_matrix), c("omic", "run", "group", "biological_sample"))
    expect_identical(object@inputs$lipidomics$biological_sample_column, "tech_rep_group")
    parameters <- integrationResolveMofaParameters(object)
    expect_identical(names(parameters), integrationMofaParameterNames())
    expect_identical(parameters$views, c("lipidomics_LCMS_Pos", "proteomics"))
    expect_warning(aligned <- integrationAlignViews(object, parameters), class = "integration_few_samples")
    lipid <- aligned$views$lipidomics_LCMS_Pos
    proteomics <- aligned$views$proteomics
    # Samples of every view kept (D4): A_3 is missing in lipidomics only.
    expect_identical(colnames(proteomics), c("A_1", "A_2", "A_3", "B_1", "B_2"))
    expect_identical(colnames(lipid), colnames(proteomics))
    expect_true(all(is.na(lipid[, "A_3"])))
    # Technical replicates averaged (INT-2), a missing run ignored; the internal standard dropped (INT-3).
    expect_identical(rownames(lipid), c("L1", "L2"))
    expect_equal(unname(lipid["L1", c("A_1", "A_2", "B_1", "B_2")]), c(21, 22, 26, 27))
    expect_equal(unname(lipid["L2", c("A_1", "B_2")]), c(3, NA))
    # P4 is observed in one sample only: dropped.
    expect_identical(rownames(proteomics), c("P1", "P2", "P3"))
    summary <- aligned$view_summary
    expect_identical(summary$internal_standards_dropped, c(1L, 0L))
    expect_identical(summary$sparse_dropped, c(0L, 1L))
    expect_identical(summary$features_model, c(2L, 3L))
    expect_identical(aligned$shared_samples, 4L)
    expect_length(aligned$notes, 1L)
    samples <- aligned$samples
    expect_identical(samples$group, c("A", "A", "A", "B", "B"))
    expect_identical(samples$views_present, c(2L, 2L, 1L, 2L, 2L))
    expect_identical(samples$lipidomics_LCMS_Pos_runs, c("l1;l2", "l3;l4", "", "l5;l6", "l7;l8"))
    expect_identical(samples$proteomics_runs, c("p1", "p2", "p3", "p4", "p5"))
    # The choices: standards kept; only the samples every view has; the most variable feature; per-feature Z-score.
    kept <- suppressWarnings(integrationAlignViews(object, integrationResolveMofaParameters(object, list(internal_standards = "keep"))))
    expect_identical(rownames(kept$views$lipidomics_LCMS_Pos), c("L1", "L2", "ITSD1"))
    complete <- suppressWarnings(integrationAlignViews(object, integrationResolveMofaParameters(object,
        list(missing_samples = "drop_incomplete"))))
    expect_identical(colnames(complete$views$proteomics), c("A_1", "A_2", "B_1", "B_2"))
    top <- suppressWarnings(integrationAlignViews(object, integrationResolveMofaParameters(object,
        list(feature_selection = "top_variance", features_per_view = 1L))))
    expect_identical(rownames(top$views$proteomics), "P1")
    scaled <- suppressWarnings(integrationAlignViews(object, integrationResolveMofaParameters(object,
        list(scaling = "feature_zscore"))))
    z <- scaled$views$proteomics
    expect_equal(unname(rowMeans(z, na.rm = TRUE)), c(0, 0, 0))
    expect_equal(unname(apply(z, 1L, stats::sd, na.rm = TRUE)), c(1, 1, 1))
})

test_that("KIT-062: the inputs and parameters are refused in plain words with classed conditions", {
    inputs <- intModelInputs()
    expect_error(integrationViewsFromInputs(inputs$proteomics), class = "integration_too_few_omics")
    other <- withr::local_tempdir()
    expect_error(integrationViewsFromInputs(c(inputs$proteomics, intModelProteomics(other))),
        class = "integration_input_duplicate")
    object <- integrationViewsFromInputs(c(inputs$proteomics, inputs$lipidomics))$object
    expect_error(integrationResolveMofaParameters(object, list(views = "proteomics")), class = "integration_too_few_views")
    expect_error(integrationResolveMofaParameters(object, list(views = c("proteomics", "nothing"))),
        class = "integration_parameters_invalid")
    expect_error(integrationResolveMofaParameters(object, list(scaling = "vst")), class = "integration_parameters_invalid")
    expect_error(integrationResolveMofaParameters(object, list(num_factors = 2.5)), class = "integration_parameters_invalid")
    expect_error(integrationResolveMofaParameters(object, list(features_per_view = 10L)), class = "integration_parameters_invalid")
    expect_error(integrationResolveMofaParameters(object, list(feature_selection = "top_variance")),
        class = "integration_parameters_invalid")
    expect_error(integrationResolveMofaParameters(object, list(input_manifests = c(lipidomics = "a", proteomics = "b"))),
        class = "integration_input_stale")
    # Every refusal is also integration_refused.
    expect_error(integrationResolveMofaParameters(object, list(views = "proteomics")), class = "integration_refused")
    # Disjoint biological samples: nothing to align on.
    disjoint <- withr::local_tempdir()
    apart <- integrationViewsFromInputs(c(intModelProteomics(disjoint, biological = c("C_1", "C_2", "C_3", "D_1", "D_2")),
        intModelLipidomics(disjoint)))$object
    expect_error(integrationAlignViews(apart, integrationResolveMofaParameters(apart)), class = "integration_no_shared_samples")
    # The provisional parameters are one table, each row naming what settles it.
    table <- integrationMofaParameterTable()
    expect_identical(table$parameter, setdiff(integrationMofaParameterNames(), c("input_manifests", "views")))
    expect_true(all(nzchar(table$settled_by)))
    # Every method choice is the workbook's to settle, except the precision, decided by the user (INT-6).
    expect_true(all(grepl("workbook", table$settled_by[table$parameter != "precision"])))
    expect_identical(table$default[table$parameter == "precision"], "double")
    expect_match(table$settled_by[table$parameter == "precision"], "INT-6")
    expect_identical(integrationResolveMofaParameters(object)$precision, "double")
    expect_identical(integrationMofaOptions(integrationResolveMofaParameters(object, list(precision = "single")))$precision, "single")
    expect_error(integrationResolveMofaParameters(object, list(precision = "half")), class = "integration_parameters_invalid")
    # The package loads without MOFA2; the probe then says so.
    if (!requireNamespace("MOFA2", quietly = TRUE)) {
        capability <- integrationMofaCapability()
        expect_false(capability$available)
        expect_match(capability$reason, "MOFA2")
    }
})

test_that("KIT-062: the owner commits the inputs state, the identity follows the inputs and every parameter, and the model is accepted and recorded", {
    inputs <- intModelInputs()
    root <- file.path(inputs$root, "project")
    manager <- newIntegrationStateOwner(root)
    withr::defer(closeWorkflowStateOwner(manager))
    importIntegrationInputs(manager, c(inputs$proteomics, inputs$lipidomics))
    expect_identical(manager$getHistory(), c("initial", "integration_inputs"))
    roles <- manager$.context$snapshot$file_dependencies
    expect_setequal(names(roles), c("integration_input_lipidomics", "integration_input_proteomics"))
    expect_identical(roles$integration_input_proteomics$sha256,
        digest::digest(file = file.path(inputs$proteomics, "integration_input.json"), algo = "sha256"))
    source <- getWorkflowStateReference(manager, "integration_inputs")
    identity <- function(parameters) withWorkflowStateScope(manager, function() {
        workflowComponentScientificIdentity(manager, source, "integration_mofa_model", parameters[integrationMofaParameterNames()])
    })
    object <- withWorkflowStateScope(manager, function() readWorkflowExecutionComponents(manager, source, fields = "inputs"))
    defaults <- integrationResolveMofaParameters(object)
    # One form whatever the caller wrote: 5 and 5L give one identity.
    expect_identical(identity(defaults), identity(integrationResolveMofaParameters(object, list(num_factors = 5))))
    changed <- defaults
    changed$num_factors <- 3L
    expect_false(identical(identity(defaults), identity(changed)))
    changed <- defaults
    changed$input_manifests[["proteomics"]] <- strrep("0", 64L)
    expect_false(identical(identity(defaults), identity(changed)))
    changed <- defaults
    changed$precision <- "single"
    expect_false(identical(identity(defaults), identity(changed)))
    changed <- defaults
    changed$views <- rev(defaults$views)
    expect_false(identical(identity(defaults), identity(changed)))
    # Unavailable MOFA2 or Python: refused before training, nothing left pending.
    expect_error(runIntegrationModel(manager, probe = function() list(available = FALSE, reason = "no Python here")),
        class = "integration_mofa_unavailable", regexp = "no Python here")
    expect_null(manager$.context$pending_provenance)
    # A stand-in trainer: computed, accepted under the exports, then reused for the same identity.
    expect_warning(run <- runIntegrationModel(manager, trainer = intModelFakeTrainer), class = "integration_few_samples")
    expect_identical(run$status, "computed")
    expect_identical(manager$getHistory(), c("initial", "integration_inputs"))
    expect_null(manager$.context$pending_provenance)
    export <- manager$.context$snapshot$exports$integration_mofa_model
    expect_identical(export$reference$entity_id, run$reference$entity_id)
    expect_identical(names(export$parameters), integrationMofaParameterNames())
    value <- run$value
    expect_identical(value$schema, "multischolar.integration_mofa_model.v1")
    expect_identical(names(value$factors), c("biological_sample", paste0("Factor", 1:5)))
    expect_setequal(value$weights$feature_id, c("L1", "L2", "P1", "P2", "P3"))
    expect_identical(value$elbo[c("iterations", "evaluations", "last")], list(iterations = 2L, evaluations = 2L, last = -10))
    again <- runIntegrationModel(manager, trainer = intModelFakeTrainer)
    expect_identical(again$status, "hit")
    # The record: the import step, the inputs by role, the model as an analysis with its parameters and outcomes.
    record <- buildWorkflowAnalysisRecord(manager)
    expect_true(validateWorkflowAnalysisRecord(record)$valid)
    expect_identical(record$omic_type, "integration")
    expect_identical(vapply(record$steps, function(step) step$operation, ""), c("import", "integration_import"))
    expect_setequal(vapply(record$inputs, function(input) input$role, ""),
        c("integration_input_lipidomics", "integration_input_proteomics"))
    analysis <- record$analyses$integration_mofa_model
    expect_identical(analysis$operation, "integration_mofa_model")
    expect_identical(analysis$source_state, "integration_inputs")
    expect_identical(vapply(analysis$parameters, function(p) p$name, ""), integrationMofaParameterNames())
    expect_identical(unique(vapply(analysis$parameters[-(1:2)], function(p) p$provenance, "")), "no_shipped_default")
    expect_identical(analysis$outcomes$views, c("lipidomics_LCMS_Pos", "proteomics"))
    expect_identical(analysis$outcomes$factors, 5L)
    # An input changed on disk since the import is refused.
    intModelProteomics(inputs$root, shift = 1)
    expect_error(runIntegrationModel(manager, trainer = intModelFakeTrainer), class = "integration_input_stale")
})

test_that("KIT-062: the replay checks the integration inputs by sha256 and takes the import as its precondition", {
    inputs <- intModelInputs()
    root <- file.path(inputs$root, "project")
    manager <- newIntegrationStateOwner(root)
    importIntegrationInputs(manager, c(inputs$proteomics, inputs$lipidomics))
    closeWorkflowStateOwner(manager)
    prepared <- file.path(inputs$root, "prepared")
    dir.create(prepared)
    file.copy(list.files(root, all.files = TRUE, full.names = TRUE, no.. = TRUE), prepared, recursive = TRUE)
    identity <- readSavedStoreIdentity(file.path(root, ".workflow-entities", "store.json"))
    manager <- newWorkflowStateOwner(newWorkflowEntityStore(root, identity$project_id, identity$workflow_id), reopen = TRUE)
    withr::defer(closeWorkflowStateOwner(manager))
    # A reopened integration store is the integration's.
    expect_identical(manager$.context$omic_type, "integration")
    suppressWarnings(runIntegrationModel(manager, trainer = intModelFakeTrainer))
    json <- file.path(inputs$root, "analysis_record.json")
    writeWorkflowAnalysisRecord(manager, json)
    record <- jsonlite::read_json(json, simplifyVector = FALSE)
    manifests <- c(integration_input_lipidomics = file.path(inputs$lipidomics, "integration_input.json"),
        integration_input_proteomics = file.path(inputs$proteomics, "integration_input.json"))
    # A manifest that is not the recorded one: refused before any step.
    wrong <- manifests
    wrong[["integration_input_proteomics"]] <- manifests[["integration_input_lipidomics"]]
    expect_error(replayWorkflowAnalysisRecord(record, prepared, wrong, list(), withr::local_tempdir()),
        class = "workflow_replay_inputs")
    # A store without the import: the precondition fails.
    empty <- file.path(inputs$root, "empty")
    closeWorkflowStateOwner(newIntegrationStateOwner(empty))
    expect_error(replayWorkflowAnalysisRecord(record, empty, manifests, list(), withr::local_tempdir()),
        class = "workflow_replay_precondition")
    # The right inputs on the prepared store pass both checks and reach the model, which needs MOFA2 (the end-to-end
    # file replays it where MOFA2 and mofapy2 are present).
    skip_if(isTRUE(integrationMofaCapability()$available), "MOFA2 is available: the replay is the end-to-end file's")
    expect_error(suppressWarnings(replayWorkflowAnalysisRecord(record, prepared, manifests, list(), withr::local_tempdir())),
        class = "integration_mofa_unavailable")
})

# KIT-147: the streamed import (contract v2). Each data set is streamed row group by row group from its Parquet file
# into the store; the state is checkpointed from zero-row prototypes and its `views` components swapped for the saved
# tables before the commit.

test_that("KIT-147: the streamed import gives the inputs state the in-memory reader builds, value for value", {
    inputs <- intModelInputs()
    dirs <- c(inputs$proteomics, inputs$lipidomics)
    memory <- integrationViewsFromInputs(dirs)$object
    manager <- newIntegrationStateOwner(file.path(inputs$root, "project"))
    withr::defer(closeWorkflowStateOwner(manager))
    importIntegrationInputs(manager, dirs)
    source <- getWorkflowStateReference(manager, "integration_inputs")
    stored <- withWorkflowStateScope(manager, function() readWorkflowExecutionComponents(manager, source, fields = TRUE))
    for (slot in c("views", "features")) {
        expect_identical(lapply(methods::slot(stored, slot), as.list), lapply(methods::slot(memory, slot), as.list), info = slot)
    }
    expect_identical(as.list(stored@design_matrix), as.list(memory@design_matrix))
    expect_identical(stored@inputs, memory@inputs)
    # The checkpoint took the swapped references: each view is a saved table of its file's rows and columns.
    spec <- withWorkflowStateScope(manager, function() readWorkflowStateMetadata(manager, source)$data)
    expect_true(all(vapply(spec$components$views, isWorkflowTableReference, logical(1))))
    expect_identical(lapply(spec$components$views, function(table) as.integer(table$metadata$shape)),
        list(lipidomics_LCMS_Pos = c(3L, 9L), proteomics = c(4L, 6L)))
    expect_identical(spec$dimensions$rows, 7L)
    # The same files imported again reuse the saved tables.
    importIntegrationInputs(manager, dirs)
    again <- withWorkflowStateScope(manager, function() readWorkflowStateMetadata(manager,
        getWorkflowStateReference(manager, "integration_inputs"))$data)
    expect_identical(lapply(again$components$views, `[[`, "entity_id"), lapply(spec$components$views, `[[`, "entity_id"))
    # The alignment on the stored state is the alignment on the in-memory one.
    parameters <- integrationResolveMofaParameters(memory)
    expect_identical(suppressWarnings(integrationAlignViews(stored, parameters)),
        suppressWarnings(integrationAlignViews(memory, parameters)))
})

test_that("KIT-147: a data set of more than 64 runs and 10,000 features imports with no text reader and no whole matrix", {
    root <- withr::local_tempdir()
    runs <- sprintf("big_%03d", seq_len(70L))
    biological <- sprintf("%s_%d", rep(c("A", "B"), each = 35L), rep(seq_len(35L), 2L))
    n <- 10050L
    values <- matrix(withr::with_seed(1147L, stats::rnorm(n * length(runs), 20, 2)), n, length(runs),
        dimnames = list(sprintf("BIG%05d", seq_len(n)), runs))
    values[cbind(1:50, 1:50)] <- NA
    big <- intModelInput(root, "proteomics", list(proteomics = values), intModelSamples(runs, biological))
    lipid <- intModelLipidomics(root)
    # A whole-file text reader or a whole-matrix read fails the import; each row group read is traced.
    pages <- integer()
    original <- readIntegrationInputRowGroup
    testthat::local_mocked_bindings(read_tsv = function(...) stop("a whole-file text reader was called"), .package = "readr")
    testthat::local_mocked_bindings(vroom = function(...) stop("a whole-file text reader was called"), .package = "vroom")
    testthat::local_mocked_bindings(
        readIntegrationInputRowGroup = function(...) {
            page <- original(...)
            pages <<- c(pages, nrow(page))
            page
        },
        readIntegrationInputMatrix = function(...) stop("a whole matrix was read"))
    project <- file.path(root, "project")
    manager <- newIntegrationStateOwner(project)
    importIntegrationInputs(manager, c(big, lipid))
    # The largest object read is one row group of at most 10,000 rows; every row was read once.
    expect_identical(pages, c(3L, 10000L, 50L))
    source <- getWorkflowStateReference(manager, "integration_inputs")
    spec <- withWorkflowStateScope(manager, function() readWorkflowStateMetadata(manager, source)$data)
    expect_identical(as.integer(spec$components$views$proteomics$metadata$shape), c(n, 71L))
    expected <- as.list(data.frame(feature_id = rownames(values), values, check.names = FALSE, stringsAsFactors = FALSE))
    stored <- withWorkflowStateScope(manager, function() readWorkflowExecutionComponents(manager, source,
        fields = "views/proteomics"))
    expect_identical(as.list(stored@views$proteomics), expected)
    # The model's path runs on the streamed state (a stand-in trainer).
    run <- suppressWarnings(runIntegrationModel(manager, trainer = intModelFakeTrainer))
    expect_identical(run$status, "computed")
    expect_identical(run$value$view_summary$features_input, c(3L, n))
    closeWorkflowStateOwner(manager)
    # The swapped references survive the store's reopening.
    identity <- readSavedStoreIdentity(file.path(project, ".workflow-entities", "store.json"))
    reopened <- newWorkflowStateOwner(newWorkflowEntityStore(project, identity$project_id, identity$workflow_id), reopen = TRUE)
    withr::defer(closeWorkflowStateOwner(reopened))
    again <- withWorkflowStateScope(reopened, function() readWorkflowExecutionComponents(reopened,
        getWorkflowStateReference(reopened, "integration_inputs"), fields = "views/proteomics"))
    expect_identical(as.list(again@views$proteomics), expected)
})
