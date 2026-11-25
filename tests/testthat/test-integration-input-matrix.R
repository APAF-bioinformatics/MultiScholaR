# KIT-153 (methods ledger INT-23): each omic's integration input holds, per data set, the normalised matrix before
# correction for unwanted variation (RUV) and, where the omic applied RUV, the RUV-corrected matrix beside it; the model
# takes the uncorrected matrices by default, the corrected ones selectable (`input_matrix`), refused for a data set that
# has none. The choice enters the model's identity, the claim limit, the record, the Methods paragraph and the replay.
# An input saved by the earlier build (one matrix per data set, no stage) is refused in plain words. The omics' writers
# are tested in their own files (test-workflow-integration-input.R, test-integration-input-lipid-metab.R,
# test-integration-input-transcriptomics.R, test-kit129-run-path.R). MOFA2 is replaced by a stand-in trainer.

intMatTable <- function(m) data.frame(feature_id = rownames(m), m, check.names = FALSE, stringsAsFactors = FALSE, row.names = NULL)

# One omic's input: one data set `view`, one run per biological sample, its corrected matrix and correction if given.
intMatInput <- function(root, omic, view, values, corrected = NULL, correction = NULL) {
    biological <- colnames(values)
    runs <- paste0(substr(omic, 1L, 1L), "_", biological)
    colnames(values) <- runs
    if (!is.null(corrected)) colnames(corrected) <- runs
    samples <- data.frame(Run = runs, group = sub("_.*$", "", biological), tech_rep_group = biological, stringsAsFactors = FALSE)
    features <- data.frame(feature_id = rownames(values), display_name = tolower(rownames(values)), stringsAsFactors = FALSE)
    manifest <- list(omic_type = omic, workflow = "TEST", matrix_role = "normalised", scale = "log2",
        normalisation = list(method = "test", source = list(name = "kit153-test")), correction = correction,
        source_state = list(name = "kit153-test"), sample_id_column = "Run", group_column = "group",
        technical_replicate_column = "tech_rep_group", biological_sample_column = "tech_rep_group")
    dir <- file.path(root, "integration", omic)
    writeIntegrationInput(dir, manifest, samples, stats::setNames(list(intMatTable(values)), view),
        stats::setNames(list(features), view),
        corrected = if (is.null(corrected)) list() else stats::setNames(list(intMatTable(corrected)), view))
    dir
}

intMatValues <- function(features, offset) {
    biological <- c("A_1", "A_2", "A_3", "B_1", "B_2", "B_3")
    m <- outer(seq_along(features), seq_along(biological), function(i, j) offset + i + j / 7 + (j > 3) * i / 3 + ((i * j) %% 5) / 11)
    dimnames(m) <- list(features, biological)
    m
}

# The corrected matrix: the normalised one less a per-sample shift (an unwanted factor removed).
intMatCorrect <- function(m) sweep(m, 2L, c(0.4, -0.3, 0.2, -0.1, 0.35, -0.25), "-")

intMatCorrection <- list(state = "ruv_corrected", id = "C1", method = "RUVIII_C_Varying", k = 1L, replicate_structure = "group",
    grouping_variable = "group")

intMatInputs <- function(lipid_corrected = TRUE, .local_envir = parent.frame()) {
    root <- withr::local_tempdir(.local_envir = .local_envir)
    proteomics <- intMatValues(c("P1", "P2", "P3", "P4"), 10)
    lipidomics <- intMatValues(c("L1", "L2", "L3"), 20)
    list(root = root, proteomics_values = proteomics, lipidomics_values = lipidomics,
        proteomics = intMatInput(root, "proteomics", "proteomics", proteomics, intMatCorrect(proteomics), intMatCorrection),
        lipidomics = intMatInput(root, "lipidomics", "lipidomics_LCMS_Pos", lipidomics,
            if (isTRUE(lipid_corrected)) intMatCorrect(lipidomics), if (isTRUE(lipid_corrected)) intMatCorrection))
}

intMatTrainer <- function(views, options) {
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

intMatOwner <- function(inputs, name = "project", .local_envir = parent.frame()) {
    manager <- newIntegrationStateOwner(file.path(inputs$root, name))
    withr::defer(closeWorkflowStateOwner(manager), envir = .local_envir)
    importIntegrationInputs(manager, c(inputs$proteomics, inputs$lipidomics))
    manager
}

# The values as saved (no feature selection, no per-feature scaling), whatever the defaults, so that they can be compared.
intMatRun <- function(manager, parameters = list()) {
    suppressWarnings(runIntegrationModel(manager, utils::modifyList(list(feature_selection = "none", scaling = "none"), parameters),
        trainer = intMatTrainer))
}

# The values handed to the model for one data set (the product's aligned values), as a features x samples matrix.
intMatAligned <- function(manager, run, view) {
    table <- integrationModelViewTable(manager, run$value, "aligned_values", view)
    m <- as.matrix(table[-1L])
    rownames(m) <- sub("^[^:]*:", "", table$feature_id)
    m
}

test_that("KIT-153: a data set's normalised matrix and its RUV-corrected matrix are written, named and read back", {
    root <- withr::local_tempdir()
    values <- intMatValues(c("P1", "P2", "P3"), 10)
    values[2L, 3L] <- NA_real_
    corrected <- intMatCorrect(values)
    dir <- intMatInput(root, "proteomics", "proteomics", values, corrected, intMatCorrection)
    # KIT-184: with a corrected matrix and a correction of k = 1, the removed factor is written beside it.
    expect_setequal(list.files(dir), c("integration_input.json", "matrix_proteomics.parquet", "features_proteomics.parquet",
        "corrected_proteomics.parquet", "removed_proteomics.parquet", "samples.parquet"))
    input <- readIntegrationInput(dir)
    manifest <- input$manifest
    expect_identical(manifest$schema, "multischolar.integration_input.v2")
    # The data set's matrix file is the normalised matrix (its role and stage); the corrected one is named beside it.
    expect_identical(manifest$matrix_role, "normalised")
    expect_identical(manifest$views[[1]][c("matrix_file", "matrix_stage", "corrected_matrix_file", "corrected_row_groups")],
        list(matrix_file = "matrix_proteomics.parquet", matrix_stage = "normalised",
            corrected_matrix_file = "corrected_proteomics.parquet", corrected_row_groups = 1L))
    expect_identical(manifest$normalisation$method, "test")
    expect_identical(manifest$correction$method, "RUVIII_C_Varying")
    expect_identical(manifest$files$corrected_proteomics.parquet,
        digest::digest(file = file.path(dir, "corrected_proteomics.parquet"), algo = "sha256"))
    expect_true(input$views$proteomics$corrected)
    expect_identical(input$views$proteomics$corrected_row_groups, 1L)
    # Each matrix comes back exactly, missing values included, in the samples' order.
    runs <- input$samples$Run
    expected <- intMatTable(values)
    names(expected)[-1L] <- runs
    expect_identical(integrationInputMatrix(dir, "proteomics", input), expected)
    expected_corrected <- intMatTable(corrected)
    names(expected_corrected)[-1L] <- runs
    expect_identical(integrationInputMatrix(dir, "proteomics", input, "ruv_corrected"), expected_corrected)
    view <- openIntegrationInputView(dir, "proteomics", input, "ruv_corrected")
    expect_identical(view$sha256, manifest$files$corrected_proteomics.parquet)
    expect_identical(view$columns, c("feature_id", runs))
    # Without RUV there is no corrected matrix, and none can be opened.
    plain <- intMatInput(withr::local_tempdir(), "proteomics", "proteomics", values)
    plain_input <- readIntegrationInput(plain)
    expect_false(plain_input$views$proteomics$corrected)
    expect_null(plain_input$manifest$views[[1]]$corrected_matrix_file)
    expect_false(file.exists(file.path(plain, "corrected_proteomics.parquet")))
    expect_error(openIntegrationInputView(plain, "proteomics", plain_input, "ruv_corrected"), class = "integration_input_invalid",
        regexp = "no RUV-corrected matrix")
    # The writer refuses a matrix role other than normalised, a corrected matrix without its correction and one whose
    # features are not the data set's in its order.
    write <- function(manifest, corrected) {
        writeIntegrationInput(file.path(withr::local_tempdir(), "proteomics"), manifest, input$samples,
            list(proteomics = expected), list(proteomics = input$views$proteomics$features), corrected = list(proteomics = corrected))
    }
    wrong <- manifest[c("omic_type", "workflow", "matrix_role", "scale", "correction", "source_state", "sample_id_column",
        "group_column", "technical_replicate_column")]
    wrong$matrix_role <- "corrected"
    expect_error(write(wrong, expected_corrected), class = "integration_input_invalid", regexp = "must be normalised")
    wrong$matrix_role <- "normalised"
    wrong$correction <- NULL
    expect_error(write(wrong, expected_corrected), class = "integration_input_invalid", regexp = "correction")
    wrong$correction <- intMatCorrection
    expect_error(write(wrong, expected_corrected[c(2L, 1L, 3L), ]), class = "integration_input_invalid", regexp = "feature ids in order")
    # The reader refuses a corrected file under another name and a corrected matrix whose correction was removed.
    manifest_file <- file.path(dir, "integration_input.json")
    rewrite <- function(change) {
        intMatInput(root, "proteomics", "proteomics", values, corrected, intMatCorrection)
        json <- change(jsonlite::read_json(manifest_file, simplifyVector = FALSE))
        jsonlite::write_json(json, manifest_file, auto_unbox = TRUE, null = "null", na = "null", digits = NA, pretty = TRUE)
    }
    rewrite(function(json) {
        file.copy(file.path(dir, "corrected_proteomics.parquet"), file.path(dir, "other.parquet"))
        json$files$other.parquet <- json$files$corrected_proteomics.parquet
        json$views[[1]]$corrected_matrix_file <- "other.parquet"
        json
    })
    expect_error(readIntegrationInput(dir), class = "integration_input_invalid", regexp = "corrected_proteomics.parquet")
    rewrite(function(json) {
        json$correction <- NULL
        json
    })
    expect_error(readIntegrationInput(dir), class = "integration_input_invalid", regexp = "correction")
})

test_that("KIT-153: an input saved by the earlier build (one matrix per data set, no stage) is refused in plain words", {
    root <- withr::local_tempdir()
    values <- intMatValues(c("P1", "P2", "P3"), 10)
    dir <- intMatInput(file.path(root, "project"), "proteomics", "proteomics", values)
    manifest_file <- file.path(dir, "integration_input.json")
    earlier <- function(role, correction = NULL) {
        json <- jsonlite::read_json(manifest_file, simplifyVector = FALSE)
        json$views <- lapply(json$views, function(view) {
            view$matrix_stage <- NULL
            view
        })
        json$normalisation <- NULL
        json$matrix_role <- role
        json$correction <- correction
        jsonlite::write_json(json, manifest_file, auto_unbox = TRUE, null = "null", na = "null", digits = NA, pretty = TRUE)
    }
    # As the earlier build wrote an omic that applied RUV (its one matrix the corrected one), and one that did not.
    for (role in c("corrected", "normalised")) {
        intMatInput(file.path(root, "project"), "proteomics", "proteomics", values)
        earlier(role, if (identical(role, "corrected")) intMatCorrection)
        expect_error(readIntegrationInput(dir), class = "integration_input_outdated", regexp = "earlier development version",
            info = role)
        expect_error(readIntegrationInput(dir), class = "integration_input_outdated", regexp = "proteomics Summary again", info = role)
        expect_error(readIntegrationInput(dir), class = "integration_input_invalid", info = role)
        expect_error(openIntegrationInputView(dir, "proteomics"), class = "integration_input_outdated", info = role)
    }
    # The Import tab's scan lists it as refused, with that sentence.
    scan <- integrationScanProject(file.path(root, "project", "integration"))
    expect_identical(scan$omics$status, "refused")
    expect_match(scan$omics$reason, "earlier development version", fixed = TRUE)
})

test_that("KIT-153: the model takes the uncorrected matrices by default; the corrected ones are selectable, change the identity and are refused where a data set has none", {
    testthat::skip_if_not(workflowOwnerCapable(), "the integration store needs qs2, DuckDB and Arrow")
    inputs <- intMatInputs()
    manager <- intMatOwner(inputs)
    source <- getWorkflowStateReference(manager, "integration_inputs")
    # The import keeps both matrices of each data set as tables of the store.
    spec <- withWorkflowStateScope(manager, function() readWorkflowStateMetadata(manager, source)$data)
    expect_setequal(names(spec$components$corrected_views), c("lipidomics_LCMS_Pos", "proteomics"))
    expect_true(all(vapply(spec$components$corrected_views, isWorkflowTableReference, logical(1))))
    object <- withWorkflowStateScope(manager, function() readWorkflowExecutionComponents(manager, source, fields = "inputs"))
    expect_setequal(integrationCorrectedViewNames(object), c("lipidomics_LCMS_Pos", "proteomics"))
    # The default is the uncorrected matrix, decided by INT-23.
    expect_identical(integrationMofaParameterSpec()$input_matrix$default, "uncorrected")
    expect_identical(integrationMofaParameterSpec()$input_matrix$choices, c("uncorrected", "ruv_corrected"))
    defaults <- integrationResolveMofaParameters(object)
    expect_identical(defaults$input_matrix, "uncorrected")
    expect_true("input_matrix" %in% integrationMofaParameterNames())
    expect_error(integrationResolveMofaParameters(object, list(input_matrix = "raw")), class = "integration_parameters_invalid")
    # The choice enters the model's identity.
    identity <- function(parameters) withWorkflowStateScope(manager, function() {
        workflowComponentScientificIdentity(manager, source, "integration_mofa_model", parameters[integrationMofaParameterNames()])
    })
    chosen <- integrationResolveMofaParameters(object, list(input_matrix = "ruv_corrected"))
    expect_false(identical(identity(defaults), identity(chosen)))
    # The default model is trained on the normalised values, the chosen one on the corrected values.
    uncorrected <- intMatRun(manager)
    expect_identical(uncorrected$status, "computed")
    expect_identical(manager$.context$snapshot$exports$integration_mofa_model$parameters$input_matrix, "uncorrected")
    expect_equal(unname(intMatAligned(manager, uncorrected, "proteomics")), unname(inputs$proteomics_values), tolerance = 0)
    corrected <- intMatRun(manager, list(input_matrix = "ruv_corrected"))
    expect_identical(corrected$status, "computed")
    expect_false(identical(corrected$reference$entity_id, uncorrected$reference$entity_id))
    expect_identical(manager$.context$snapshot$exports$integration_mofa_model$parameters$input_matrix, "ruv_corrected")
    expect_equal(unname(intMatAligned(manager, corrected, "proteomics")), unname(intMatCorrect(inputs$proteomics_values)),
        tolerance = 1e-12)
    expect_equal(unname(intMatAligned(manager, corrected, "lipidomics_LCMS_Pos")), unname(intMatCorrect(inputs$lipidomics_values)),
        tolerance = 1e-12)
    # The whole-matrix alignment (the partitioned one is the default) reads the same corrected tables (a second owner,
    # so that the model is computed again).
    other <- intMatOwner(inputs, "whole")
    whole <- intMatRun(other, list(input_matrix = "ruv_corrected", variant = "whole"))
    expect_identical(whole$status, "computed")
    for (view in c("lipidomics_LCMS_Pos", "proteomics")) {
        expect_identical(intMatAligned(other, whole, view), intMatAligned(manager, corrected, view), info = view)
    }
    # A data set without a corrected matrix: the corrected choice is refused in plain words, the default trains.
    lacking <- intMatInputs(lipid_corrected = FALSE)
    lacking_manager <- intMatOwner(lacking)
    refusal <- tryCatch(intMatRun(lacking_manager, list(input_matrix = "ruv_corrected")), error = function(error) error)
    expect_s3_class(refusal, "integration_input_matrix_unavailable")
    expect_s3_class(refusal, "integration_parameters_invalid")
    expect_match(conditionMessage(refusal), "lipidomics_LCMS_Pos", fixed = TRUE)
    expect_identical(integrationRefusalText(refusal)$sentence, integrationGuiTexts()[["refusal_input_matrix_unavailable"]])
    expect_identical(intMatRun(lacking_manager)$status, "computed")
    # The Import tab says which inputs hold both matrices.
    scan <- integrationScanProject(file.path(lacking$root, "integration"))
    expect_identical(scan$omics$values[match(c("lipidomics", "proteomics"), scan$omics$omic)],
        unname(integrationGuiTexts()[c("values_normalised", "values_both")]))
})

test_that("KIT-153: the claim limit, the record, the Methods paragraph and the replay's parameters follow the matrix the model used", {
    testthat::skip_if_not(workflowOwnerCapable(), "the integration store needs qs2, DuckDB and Arrow")
    inputs <- intMatInputs()
    manager <- intMatOwner(inputs)
    gui <- integrationGuiTexts()
    texts <- integrationRecordTexts()
    object <- integrationImportedObject(manager, "inputs")
    views <- c("lipidomics_LCMS_Pos", "proteomics")
    # The claim limit: no data set group-corrected under the uncorrected matrices; each under the corrected ones.
    expect_identical(unname(integrationGroupCorrectedViews(object, views, "uncorrected")), c(FALSE, FALSE))
    expect_identical(unname(integrationGroupCorrectedViews(object, views, "ruv_corrected")), c(TRUE, TRUE))
    expect_match(paste(integrationClaimLimitText(list(claim = list(group_corrected = character(), corrections = c(a = "x"),
        uncorrected = TRUE)), gui), collapse = " "), gui[["claim_uncorrected"]], fixed = TRUE)
    record_of <- function() {
        record <- buildWorkflowAnalysisRecord(manager)
        expect_true(validateWorkflowAnalysisRecord(record)$valid)
        record
    }
    parameter <- function(analysis, name) Filter(function(p) identical(p$name, name), analysis$parameters)[[1L]]
    # The default: the record names the uncorrected matrices, the caution, and no claim limit.
    intMatRun(manager)
    record <- record_of()
    analysis <- record$analyses$integration_mofa_model
    block <- analysis$integration
    expect_identical(vapply(block$inputs, function(input) input$matrix_used, ""), c("uncorrected", "uncorrected"))
    expect_true(all(vapply(block$inputs, function(input) input$corrected_available, logical(1))))
    expect_false(block$claim_limit$applies)
    expect_true(block$claim_limit$uncorrected)
    expect_true(gui[["claim_uncorrected"]] %in% block$claim_limit$text)
    expect_identical(parameter(analysis, "input_matrix")$status, "decided")
    methods <- integrationMethodsText(record)
    expect_true(grepl(texts[["methods_claim_uncorrected"]], methods, fixed = TRUE))
    expect_true(grepl(sprintf("proteomics (proteomics; %s)", sprintf(texts[["methods_input_uncorrected_available"]],
        integrationCorrectionText(object@inputs$proteomics$manifest, gui))), methods, fixed = TRUE))
    tables <- integrationRecordTables(record)
    expect_identical(unique(tables$integration_inputs[[texts[["col_values"]]]]), texts[["values_uncorrected_available"]])
    expect_true("input_matrix" %in% tables$integration_settings[[texts[["col_name"]]]])
    # The corrected matrices chosen: the claim limit for both data sets, no caution, the choice set by the user.
    intMatRun(manager, list(input_matrix = "ruv_corrected"))
    record <- record_of()
    analysis <- record$analyses$integration_mofa_model
    block <- analysis$integration
    expect_identical(vapply(block$inputs, function(input) input$matrix_used, ""), c("ruv_corrected", "ruv_corrected"))
    expect_true(block$claim_limit$applies)
    expect_identical(block$claim_limit$group_corrected, views)
    expect_false(block$claim_limit$uncorrected)
    expect_identical(parameter(analysis, "input_matrix")$status, "user_set")
    methods <- integrationMethodsText(record)
    expect_true(grepl(sprintf(texts[["methods_claim_limit"]], paste(views, collapse = ", ")), methods, fixed = TRUE))
    expect_false(grepl(texts[["methods_claim_uncorrected"]], methods, fixed = TRUE))
    expect_true(grepl(sprintf("proteomics (proteomics; %s)", sprintf(texts[["methods_input_corrected"]],
        integrationCorrectionText(object@inputs$proteomics$manifest, gui))), methods, fixed = TRUE))
    expect_identical(unique(integrationRecordTables(record)$integration_inputs[[texts[["col_values"]]]]), texts[["values_corrected"]])
    # The replay takes the recorded choice.
    replayed <- workflowReplayParameters(analysis$parameters, workflowOperationEntry(INTEGRATION_MODEL_OPERATION)$parameters)
    expect_identical(replayed$input_matrix, "ruv_corrected")
    # Read back from JSON, the record says the same.
    json <- file.path(inputs$root, "record.json")
    writeWorkflowAnalysisRecord(manager, json)
    read <- jsonlite::read_json(json, simplifyVector = FALSE)
    expect_identical(integrationMethodsText(read), methods)
})
