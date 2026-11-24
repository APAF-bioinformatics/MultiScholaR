# KIT-061: the per-omic integration input (the contract, utils_workflow_integration_input.R) and the proteomics writer
# at the Summary hook (utils_workflow_dia_summary.R). KIT-147: contract v2, Parquet tables; the matrix is read through
# `openIntegrationInputView()` (helper-integration-input.R).

integrationInputFixture <- function() {
    samples <- data.frame(Run = c("s1", "s2", "s3", "s4"), group = c("A", "A", "B", "B"),
        tech_rep_group = c("A_1", "A_2", "B_1", "B_2"), stringsAsFactors = FALSE)
    matrix <- data.frame(feature_id = c("P1", "P2", "P3"), s1 = c(20.125, 21.5, NA), s2 = c(20.25, 21.375, 18.0625),
        s3 = c(22.1, 21.4, 18.3), s4 = c(22.3, 21.6, 1 / 3), stringsAsFactors = FALSE)
    features <- data.frame(feature_id = c("P1", "P2", "P3"), display_name = c("GENE1", "P2", "GENE3"),
        stringsAsFactors = FALSE)
    manifest <- list(omic_type = "proteomics", workflow = "DIA", matrix_role = "normalised", scale = "log2",
        correction = NULL, source_state = list(name = "correlation_filtered", entity_id = strrep("a", 64),
            sha256 = strrep("b", 64)), sample_id_column = "Run", group_column = "group",
        technical_replicate_column = "tech_rep_group")
    list(manifest = manifest, samples = samples, matrix = matrix, features = features)
}

writeIntegrationInputFixture <- function(dir, fixture = integrationInputFixture(), corrected = list()) {
    writeIntegrationInput(dir, fixture$manifest, fixture$samples, list(proteomics = fixture$matrix),
        list(proteomics = fixture$features), now = as.POSIXct("2026-11-17 10:00:00", tz = "UTC"), corrected = corrected)
}

# KIT-153: one state's protein table (`feature_id` first), and the written matrix against it on the matrix's rows and
# columns, value for value.
integrationLimpaTable <- function(manager, reference) {
    table <- withDIAStateScope(manager, function() {
        readDIAExecutionComponents(manager, reference, fields = c("protein_quant_table", "design_matrix"))@protein_quant_table
    })
    names(table)[[1]] <- "feature_id"
    as.data.frame(table, stringsAsFactors = FALSE)
}

integrationLimpaExpectSameTable <- function(written, table) {
    at <- match(written$feature_id, as.character(table$feature_id))
    expect_false(anyNA(at))
    expect_identical(unname(as.matrix(written[-1L])), unname(as.matrix(table[at, names(written)[-1L], drop = FALSE])))
}

test_that("KIT-061: the contract round trip writes the fixed files and reads them back validated", {
    root <- withr::local_tempdir()
    dir <- file.path(root, "integration", "proteomics")
    fixture <- integrationInputFixture()
    path <- writeIntegrationInputFixture(dir, fixture)
    expect_identical(path, file.path(dir, "integration_input.json"))
    expect_setequal(list.files(dir), c("integration_input.json", "matrix_proteomics.parquet", "features_proteomics.parquet",
        "samples.parquet"))
    # Nothing is left beside the folder (the staging folder is renamed into place).
    expect_identical(list.files(file.path(root, "integration"), all.files = TRUE, no.. = TRUE), "proteomics")
    input <- readIntegrationInput(dir)
    manifest <- input$manifest
    expect_identical(manifest$schema, "multischolar.integration_input.v2")
    expect_identical(manifest$omic_type, "proteomics")
    expect_identical(manifest$workflow, "DIA")
    expect_identical(manifest$matrix_role, "normalised")
    expect_null(manifest$correction)
    expect_identical(manifest$source_state$entity_id, strrep("a", 64))
    expect_identical(manifest$software$package, "MultiScholaR")
    expect_identical(manifest$software$version, as.character(utils::packageVersion("MultiScholaR")))
    expect_identical(manifest$written_at, "2026-11-17T10:00:00Z")
    expect_identical(manifest$views[[1]][c("name", "matrix_file", "features_file", "rows", "columns", "row_groups")],
        list(name = "proteomics", matrix_file = "matrix_proteomics.parquet", features_file = "features_proteomics.parquet",
            rows = 3L, columns = 4L, row_groups = 1L))
    expect_identical(manifest$parquet, list(format_version = "2.6", compression = "snappy", row_group_rows = 10000L,
        column_types = list(feature_id = "string", runs = "double", samples = "string", features = "string")))
    for (name in names(manifest$files)) {
        expect_identical(manifest$files[[name]], digest::digest(file = file.path(dir, name), algo = "sha256"))
    }
    expect_setequal(names(manifest$files), c("matrix_proteomics.parquet", "features_proteomics.parquet", "samples.parquet"))
    # The values, missing values included, come back exactly.
    expect_identical(integrationInputMatrix(dir, "proteomics", input), fixture$matrix)
    expect_identical(input$views$proteomics[c("rows", "columns", "row_groups")], list(rows = 3L, columns = 4L, row_groups = 1L))
    expect_identical(input$views$proteomics$features, fixture$features)
    expect_identical(input$samples, fixture$samples)
    # A second write replaces the folder whole (KIT-153: with a corrected matrix beside the normalised one).
    fixture$manifest$correction <- list(id = "C1", k = 2L)
    writeIntegrationInputFixture(dir, fixture, corrected = list(proteomics = fixture$matrix))
    expect_identical(readIntegrationInput(dir)$manifest$correction$k, 2L)
    expect_identical(list.files(file.path(root, "integration"), all.files = TRUE, no.. = TRUE), "proteomics")
})

test_that("KIT-061 (K061-2): a failed publish or move-aside leaves the previous input intact", {
    root <- withr::local_tempdir()
    dir <- file.path(root, "integration", "proteomics")
    writeIntegrationInputFixture(dir)
    before <- readIntegrationInput(dir)
    fixture <- integrationInputFixture()
    fixture$manifest$correction <- list(id = "C1", k = 2L)
    write <- function(rename) writeIntegrationInput(dir, fixture$manifest, fixture$samples, list(proteomics = fixture$matrix),
        list(proteomics = fixture$features), renameFn = rename, corrected = list(proteomics = fixture$matrix))
    # The staging folder cannot be renamed into place: the previous input is restored.
    calls <- 0L
    failPublish <- function(from, to) {
        calls <<- calls + 1L
        if (identical(to, dir) && grepl("^\\.proteomics-", basename(from)) && !grepl("previous", basename(from))) FALSE
        else file.rename(from, to)
    }
    expect_error(write(failPublish), class = "integration_input_invalid", regexp = "previous input is kept")
    expect_identical(calls, 3L)
    expect_identical(readIntegrationInput(dir), before)
    expect_identical(list.files(dirname(dir), all.files = TRUE, no.. = TRUE), "proteomics")
    # The previous input cannot be moved aside: nothing is replaced.
    expect_error(write(function(from, to) FALSE), class = "integration_input_invalid", regexp = "moved aside")
    expect_identical(readIntegrationInput(dir), before)
    expect_identical(list.files(dirname(dir), all.files = TRUE, no.. = TRUE), "proteomics")
})

test_that("KIT-061 (K061-6): the reader ties the files to the fixed names, the omic to its folder and refuses duplicate views and dot names", {
    root <- withr::local_tempdir()
    dir <- file.path(root, "proteomics")
    manifest_file <- file.path(dir, "integration_input.json")
    rewrite <- function(change) {
        writeIntegrationInputFixture(dir)
        manifest <- jsonlite::read_json(manifest_file, simplifyVector = FALSE)
        manifest <- change(manifest)
        jsonlite::write_json(manifest, manifest_file, auto_unbox = TRUE, null = "null", na = "null", digits = NA,
            pretty = TRUE)
    }
    # The matrix file renamed in the manifest (with its own sha256): refused, not read.
    rewrite(function(manifest) {
        file.copy(file.path(dir, "matrix_proteomics.parquet"), file.path(dir, "other.parquet"))
        manifest$files$other.parquet <- manifest$files$matrix_proteomics.parquet
        manifest$views[[1]]$matrix_file <- "other.parquet"
        manifest
    })
    expect_error(readIntegrationInput(dir), class = "integration_input_invalid", regexp = "matrix_proteomics.parquet")
    rewrite(function(manifest) {
        manifest$views[[1]]$features_file <- "samples.parquet"
        manifest
    })
    expect_error(readIntegrationInput(dir), class = "integration_input_invalid", regexp = "features_proteomics.parquet")
    # The omic is the folder's.
    rewrite(function(manifest) {
        manifest$omic_type <- "lipidomics"
        manifest
    })
    expect_error(readIntegrationInput(dir), class = "integration_input_invalid", regexp = "folder")
    # A view listed twice.
    rewrite(function(manifest) {
        manifest$views <- c(manifest$views, manifest$views)
        manifest
    })
    expect_error(readIntegrationInput(dir), class = "integration_input_invalid", regexp = "unique")
    # "." and ".." as file names or a view name.
    for (name in c(".", "..")) {
        rewrite(function(manifest) {
            manifest$files[[name]] <- strrep("0", 64)
            manifest
        })
        expect_error(readIntegrationInput(dir), class = "integration_input_invalid")
        rewrite(function(manifest) {
            manifest$views[[1]]$name <- name
            manifest
        })
        expect_error(readIntegrationInput(dir), class = "integration_input_invalid")
    }
})

test_that("KIT-061 (K061-4): the accession is the first member's, without its isoform suffix; the feature id stays the key", {
    ids <- c("P1-2;Q2", "P3;P4-1", "P5-10", "P6", "sp|X-Y")
    annotation <- data.frame(Entry = c("P1", "P3", "P5"), gene_names = c("ONE", "THREE X", "FIVE"),
        stringsAsFactors = FALSE)
    features <- protIntegrationInputFeatures(ids, annotation)
    expect_identical(features$feature_id, ids)
    expect_identical(features$uniprot_accession, c("P1", "P3", "P5", "P6", "sp|X-Y"))
    expect_identical(features$gene_names, c("ONE", "THREE X", "FIVE", NA, NA))
    expect_identical(features$display_name, c("ONE", "THREE", "FIVE", "P6", "sp|X-Y"))
})

test_that("KIT-061: a tampered file, a wrong shape, a missing manifest and inconsistent ids are refused", {
    dir <- file.path(withr::local_tempdir(), "proteomics")
    writeIntegrationInputFixture(dir)
    matrix_file <- file.path(dir, "matrix_proteomics.parquet")
    original <- integrationInputMatrix(dir, "proteomics")
    # One byte changed: the file no longer matches its sha256.
    bytes <- readBin(matrix_file, "raw", file.size(matrix_file))
    bytes[[64L]] <- xor(bytes[[64L]], as.raw(1L))
    writeBin(bytes, matrix_file)
    expect_error(readIntegrationInput(dir), class = "integration_input_invalid", regexp = "sha256")
    # One sample column dropped and the hash updated: the shape no longer matches the manifest.
    writeIntegrationInputParquet(original[-5L], matrix_file)
    manifest_file <- file.path(dir, "integration_input.json")
    manifest <- jsonlite::read_json(manifest_file, simplifyVector = FALSE)
    manifest$files$matrix_proteomics.parquet <- digest::digest(file = matrix_file, algo = "sha256")
    jsonlite::write_json(manifest, manifest_file, auto_unbox = TRUE, null = "null", na = "null", digits = NA)
    expect_error(readIntegrationInput(dir), class = "integration_input_invalid", regexp = "is not 3 x 4")
    # A corrected matrix without its correction (KIT-153: the corrected file beside the normalised one).
    fixture <- integrationInputFixture()
    fixture$manifest$correction <- list(id = "C1", k = 2L)
    writeIntegrationInputFixture(dir, fixture, corrected = list(proteomics = fixture$matrix))
    manifest <- jsonlite::read_json(manifest_file, simplifyVector = FALSE)
    manifest$correction <- NULL
    jsonlite::write_json(manifest, manifest_file, auto_unbox = TRUE, null = "null", na = "null", digits = NA)
    expect_error(readIntegrationInput(dir), class = "integration_input_invalid", regexp = "correction")
    unlink(manifest_file)
    expect_error(readIntegrationInput(dir), class = "integration_input_invalid", regexp = "no integration input manifest")
    # The writer refuses a matrix whose samples are not the design's, and features out of the matrix's order.
    fixture <- integrationInputFixture()
    fixture$samples <- fixture$samples[1:3, ]
    expect_error(writeIntegrationInputFixture(dir, fixture), class = "integration_input_invalid", regexp = "samples")
    fixture <- integrationInputFixture()
    fixture$features <- fixture$features[c(2, 1, 3), ]
    expect_error(writeIntegrationInputFixture(dir, fixture), class = "integration_input_invalid", regexp = "features")
    fixture <- integrationInputFixture()
    fixture$matrix$feature_id[[2]] <- "P1"
    expect_error(writeIntegrationInputFixture(dir, fixture), class = "integration_input_invalid", regexp = "unique")
})

test_that("KIT-147: the v2 round trip is bit-identical: doubles exact, missing values kept, columns in the samples' order, row groups of 10,000", {
    dir <- file.path(withr::local_tempdir(), "proteomics")
    n <- 25001L
    runs <- c("r3", "r1", "r2")
    samples <- data.frame(Run = runs, group = c("A", "B", "A"), tech_rep_group = c("A_2", "B_1", "A_1"),
        replicates = c(2L, 1L, 1L), stringsAsFactors = FALSE)
    values <- withr::with_seed(147L, stats::rnorm(3L * n, 20, 3))
    # Values the v1 text round trip could not keep: below 1e-148, subnormal, -0, NaN beside NA, the extremes.
    special <- c(NA, NaN, -0, 1 / 3, 1e-300, .Machine$double.xmin / 3, .Machine$double.xmax, -Inf, Inf, 2^-1074)
    values[seq_along(special)] <- special
    # The matrix's columns given in another order than the samples table's.
    matrix <- data.frame(feature_id = sprintf("F%06d", seq_len(n)), r1 = values[seq_len(n)], r2 = values[n + seq_len(n)],
        r3 = values[2L * n + seq_len(n)], stringsAsFactors = FALSE)
    features <- data.frame(feature_id = matrix$feature_id, display_name = c(NA, "", matrix$feature_id[-(1:2)]),
        stringsAsFactors = FALSE)
    manifest <- integrationInputFixture()$manifest
    writeIntegrationInput(dir, manifest, samples, list(proteomics = matrix), list(proteomics = features))
    input <- readIntegrationInput(dir)
    expected <- matrix[c("feature_id", runs)]
    view <- openIntegrationInputView(dir, "proteomics", input)
    expect_identical(view$columns, c("feature_id", runs))
    expect_identical(c(view$rows, view$row_groups), c(n, 3L))
    expect_identical(input$manifest$views[[1]]$row_groups, 3L)
    whole <- view$read()
    expect_identical(names(whole), c("feature_id", runs))
    # Bit for bit (num.eq = FALSE tells -0 from 0; NA and NaN stay apart).
    expect_true(identical(whole, expected, num.eq = FALSE))
    expect_identical(whole$r1[seq_along(special)], special)
    expect_true(identical(whole$r1[[3L]], -0, num.eq = FALSE))
    expect_true(is.na(whole$r1[[1L]]) && !is.nan(whole$r1[[1L]]))
    expect_true(is.nan(whole$r1[[2L]]))
    # Row group by row group, the same rows.
    pages <- lapply(seq_len(view$row_groups), view$read_row_group)
    expect_identical(vapply(pages, nrow, integer(1)), c(10000L, 10000L, 5001L))
    expect_true(identical(do.call(rbind, pages)$r3, expected$r3, num.eq = FALSE))
    # The samples and features as written, as text; a missing and an empty text kept apart.
    expect_identical(input$samples, data.frame(Run = runs, group = c("A", "B", "A"), tech_rep_group = c("A_2", "B_1", "A_1"),
        replicates = c("2", "1", "1"), stringsAsFactors = FALSE))
    expect_identical(input$views$proteomics$features, features)
    # Written twice, the same bytes.
    again <- file.path(withr::local_tempdir(), "proteomics")
    writeIntegrationInput(again, manifest, samples, list(proteomics = matrix), list(proteomics = features),
        now = as.POSIXct(input$manifest$written_at, format = "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"))
    for (file in names(input$manifest$files)) {
        expect_identical(digest::digest(file = file.path(again, file), algo = "sha256"), input$manifest$files[[file]], info = file)
    }
})

test_that("KIT-147: a v1 folder (TSV tables) is refused in plain words: it was saved by an earlier development build", {
    root <- withr::local_tempdir()
    dir <- file.path(root, "integration", "proteomics")
    dir.create(dir, recursive = TRUE)
    fixture <- integrationInputFixture()
    utils::write.table(fixture$matrix, file.path(dir, "matrix_proteomics.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
    manifest <- c(fixture$manifest, list(schema = "multischolar.integration_input.v1", views = list(list(name = "proteomics",
        matrix_file = "matrix_proteomics.tsv", features_file = "features_proteomics.tsv", rows = 3L, columns = 4L)),
        samples_file = "samples.tsv", files = list(matrix_proteomics.tsv = digest::digest(file = file.path(dir,
        "matrix_proteomics.tsv"), algo = "sha256"))))
    jsonlite::write_json(manifest, file.path(dir, "integration_input.json"), auto_unbox = TRUE, null = "null", digits = NA)
    expect_error(readIntegrationInput(dir), class = "integration_input_outdated", regexp = "earlier development version")
    expect_error(readIntegrationInput(dir), class = "integration_input_outdated", regexp = "proteomics Summary again")
    expect_error(readIntegrationInput(dir), class = "integration_input_invalid")
    expect_error(openIntegrationInputView(dir, "proteomics"), class = "integration_input_outdated")
    # The Import tab's scan lists it as refused, with that sentence.
    scan <- integrationScanProject(file.path(root, "integration"))
    expect_identical(scan$omics$status, "refused")
    expect_match(scan$omics$reason, "earlier development version", fixed = TRUE)
})

# The proteomics writer, through the Summary hook (`saveWorkflowSummaryParameters` under proteomics' group).

integrationSummarySave <- function(manager, root, workflow = list2env(list(state_manager = manager,
    contrasts_tbl = NULL, config_list = NULL))) {
    paths <- list(source_dir = file.path(root, "scripts", "proteomics"), integration_dir = file.path(root, "integration"),
        base_dir = root)
    dir.create(paths$source_dir, recursive = TRUE, showWarnings = FALSE)
    shown <- character()
    saved <- NULL
    failure <- NULL
    saveWorkflowSummaryParameters(omicGuiConfig("proteomics")$summary, workflow,
        projectDirs = list(proteomics = paths), omicType = "proteomics", label = "kit061", description = "d",
        notify = function(message, ...) shown <<- c(shown, message), done = function(value) saved <<- value,
        failed = function(error) failure <<- error, createWorkflowArgsFn = function(...) "study_parameters.txt",
        writeAnalysisRecordFn = function(...) NULL, catFn = function(...) NULL)
    if (!is.null(failure)) stop(failure)
    list(shown = shown, saved = saved, dir = file.path(paths$integration_dir, "proteomics"),
        integration_dir = paths$integration_dir)
}

integrationProteinOwner <- function(root, workflow_type) {
    manager <- newDIAStateAdapter(newWorkflowEntityStore(root, "kit061", "DIA"), WorkflowState$new(audit_enabled = FALSE))
    manager$setWorkflowType(workflow_type)
    manager$.context$reference_states <- TRUE
    manager
}

test_that("KIT-061: a classic route with RUV writes the corrected matrix after the correlation filter, and no archive or S4 object", {
    root <- withr::local_tempdir()
    manager <- integrationProteinOwner(root, "DIA")
    withr::defer(closeDIAStateAdapter(manager))
    protein <- makeWorkflowCheckpointProtein(readWorkflowCheckpointRecipe("small"))
    controls <- seq_len(nrow(protein@protein_quant_table)) <= 7L
    # KIT-153: the normalisation's product names the normalised matrix its RUV read (here other values than the
    # corrected state's, so that the two files can be told apart).
    normalised <- protein
    samples_all <- setdiff(names(normalised@protein_quant_table), normalised@protein_id_column)
    normalised@protein_quant_table[samples_all] <- normalised@protein_quant_table[samples_all] - 0.5
    manager$saveState("protein_replicate_filtered", normalised, list(), "normalised stand-in")
    product <- saveDIAStateValue(manager$.context, list(best_k = 3, controls = controls, percentage = 5,
        parameters = list(grouping = "group"), normalized = getDIAStateReference(manager, "protein_replicate_filtered")),
        "normalization-science", codec = "rds")
    checkpoint <- checkpointDIAExecutionResult(manager, protein)
    checkpoint$products <- list(normalization_science = product)
    commitDIAReferenceState(manager, "ruv_corrected", checkpoint, list(norm_method = "cyclicloess", ruv_mode = "automatic",
        ruv_k = 3, percentage_as_neg_ctrl = 5), "ruv")
    # The correlation filter drops one sample.
    filtered <- protein
    filtered@protein_quant_table <- filtered@protein_quant_table[, names(filtered@protein_quant_table) != "run_0005"]
    filtered@design_matrix <- filtered@design_matrix[filtered@design_matrix$Run != "run_0005", ]
    manager$saveState("correlation_filtered", filtered, list(min_pearson_correlation_threshold = 0.5), "correlation")
    # An earlier save's archive is left as it was; the writer adds nothing beside its folder.
    # KIT-184 K4: these stand-in states differ by a constant per protein (0.5), so no removed direction can be determined:
    # the writer says so with its classed warning and writes no factors file.
    expect_warning(run <- integrationSummarySave(manager, root), class = "integration_input_removed_factors_undetermined")
    expect_false(file.exists(file.path(run$dir, "removed_proteomics.parquet")))
    expect_identical(readIntegrationInput(run$dir)$manifest$views[[1L]]$removed_factors$status, "not_determined")
    expect_identical(run$saved$integrationFile, file.path("proteomics", "integration_input.json"))
    expect_true("Saved Integration input" %in% run$shown)
    expect_identical(list.files(run$integration_dir, all.files = TRUE, no.. = TRUE), "proteomics")
    input <- readIntegrationInput(run$dir)
    manifest <- input$manifest
    expect_identical(manifest$workflow, "DIA")
    # KIT-153: the role of the data set's matrix file, the normalised one; the corrected matrix beside it.
    expect_identical(manifest$matrix_role, "normalised")
    expect_true(input$views$proteomics$corrected)
    expect_identical(manifest$normalisation$source$entity_id, getDIAStateReference(manager, "protein_replicate_filtered")$entity_id)
    expect_identical(manifest$source_state$name, "correlation_filtered")
    expect_identical(manifest$source_state$entity_id, getDIAStateReference(manager, "correlation_filtered")$entity_id)
    expect_null(manifest$source_state$product)
    # KIT-062: the biological sample's column, by the other omics' rule (the declared technical-replicate column).
    expect_identical(manifest$biological_sample_column, "tech_rep_group")
    expect_true(manifest$biological_sample_column %in% names(input$samples))
    expect_identical(manifest$correction[c("state", "id", "method", "ruv_mode", "k", "percentage_as_neg_ctrl",
        "control_count", "grouping_variable")], list(state = "ruv_corrected", id = "C1", method = "RUVIII_C_Varying",
        ruv_mode = "automatic", k = 3L, percentage_as_neg_ctrl = 5L, control_count = 7L, grouping_variable = "group"))
    expect_identical(unlist(manifest$correction$grouping_levels), c("control", "treated"))
    expected <- filtered@protein_quant_table
    names(expected)[[1]] <- "feature_id"
    expect_identical(integrationInputMatrix(run$dir, "proteomics", input, "ruv_corrected"), expected)
    # KIT-153: the data set's own matrix is the normalised one, on the corrected matrix's proteins and samples.
    uncorrected <- normalised@protein_quant_table[names(filtered@protein_quant_table)]
    names(uncorrected)[[1]] <- "feature_id"
    expect_identical(integrationInputMatrix(run$dir, "proteomics", input), uncorrected)
    expect_identical(input$samples$Run, setdiff(protein@design_matrix$Run, "run_0005"))
    expect_identical(input$samples$group, as.character(filtered@design_matrix$group))
    # No UniProt table: each display name is the feature id.
    expect_identical(input$views$proteomics$features$display_name, expected$feature_id)
})

test_that("KIT-153: classic proteomics: the data set's matrix is the normalised matrix RUV-III-C corrected; the corrected one is beside it", {
    root <- withr::local_tempdir()
    manager <- integrationProteinOwner(root, "DIA")
    withr::defer(closeDIAStateAdapter(manager))
    manager$saveState("protein_replicate_filtered", makeWorkflowCheckpointProtein(readWorkflowCheckpointRecipe("small")), list(),
        "source")
    parameters <- list(norm_method = "quantile", ruv_mode = "manual", auto_percentage_min = NULL, auto_percentage_max = NULL,
        separation_metric = NULL, k_penalty_weight = NULL, adaptive_k_penalty = NULL, ruv_percentage = 50, ruv_k = 1,
        grouping = "group")
    run <- NULL
    utils::capture.output(run <- suppressWarnings(suppressMessages(withDIAStateScope(manager, function() {
        run <- computeDIAReferenceNormalization(manager, getDIAStateReference(manager), parameters)
        commit <- diaNormalizationCommit(manager, parameters, run)
        commitDIAReferenceState(manager, commit$state_name, commit$checkpoint, commit$config, commit$description)
        run
    }))))
    expect_identical(manager$current_state, "ruv_corrected")
    saved_dir <- integrationSummarySave(manager, root)$dir
    input <- readIntegrationInputWithMatrices(saved_dir)
    manifest <- input$manifest
    expect_identical(manifest$matrix_role, "normalised")
    expect_true(input$views$proteomics$corrected)
    expect_identical(manifest$correction[c("id", "method")], list(id = "C1", method = "RUVIII_C_Varying"))
    expect_equal(manifest$correction$k, run$value$best_k)
    # KIT-184 K4 (gate G1): the written factors are what RUV-III-C removed.
    expectIntegrationRemovedFactors(saved_dir, "proteomics", run$value$best_k, input = input)
    # The normalisation named: the method and the normalised checkpoint the RUV step read.
    expect_identical(manifest$normalisation$method, "quantile")
    expect_identical(manifest$normalisation$source[c("name", "entity_id")],
        list(name = "normalized", entity_id = run$value$normalized$entity_id))
    read <- function(reference) {
        withDIAStateScope(manager, function() readDIAExecutionComponents(manager, reference))
    }
    table <- function(object) {
        out <- as.data.frame(object@protein_quant_table, stringsAsFactors = FALSE)
        names(out)[[1]] <- "feature_id"
        out
    }
    normalised <- read(run$value$normalized)
    # The data set's matrix: the normalised one, on the corrected matrix's proteins and runs.
    integrationLimpaExpectSameTable(input$views$proteomics$matrix, table(normalised))
    expect_identical(input$views$proteomics$matrix$feature_id, input$views$proteomics$corrected_matrix$feature_id)
    # The corrected one: the state's table.
    integrationLimpaExpectSameTable(input$views$proteomics$corrected_matrix, table(read(run$value$final)))
    # And it is the matrix that entered RUV: RUV-III-C with the step's k and controls, then its missing-value cleanup, run
    # again on it gives the corrected matrix.
    again <- suppressWarnings(suppressMessages(removeRowsWithMissingValuesPercent(theObject = ruvIII_C_Varying(normalised,
        ruv_grouping_variable = "group", ruv_number_k = run$value$best_k, ctrl = run$value$controls))))
    again <- table(again)
    at <- match(input$views$proteomics$corrected_matrix$feature_id, again$feature_id)
    expect_false(anyNA(at))
    expect_equal(unname(as.matrix(input$views$proteomics$corrected_matrix[-1L])),
        unname(as.matrix(again[at, names(input$views$proteomics$corrected_matrix)[-1L]])), tolerance = 1e-10)
})

# KIT-184 K4 (gate G1): the classic RUV-III-C step run with a given number of factors on a protein fixture.
integrationRuvStep <- function(manager, protein, k) {
    manager$saveState("protein_replicate_filtered", protein, list(), "source")
    parameters <- list(norm_method = "quantile", ruv_mode = "manual", auto_percentage_min = NULL, auto_percentage_max = NULL,
        separation_metric = NULL, k_penalty_weight = NULL, adaptive_k_penalty = NULL, ruv_percentage = 50, ruv_k = k,
        grouping = "group")
    run <- NULL
    utils::capture.output(run <- suppressWarnings(suppressMessages(withDIAStateScope(manager, function() {
        run <- computeDIAReferenceNormalization(manager, getDIAStateReference(manager), parameters)
        commit <- diaNormalizationCommit(manager, parameters, run)
        commitDIAReferenceState(manager, commit$state_name, commit$checkpoint, commit$config, commit$description)
        run
    }))))
    run
}

test_that("KIT-184 K4: classic proteomics with two factors and missing values: the written factors are what was removed", {
    root <- withr::local_tempdir()
    manager <- integrationProteinOwner(root, "DIA")
    withr::defer(closeDIAStateAdapter(manager))
    run <- integrationRuvStep(manager, makeWorkflowCheckpointProtein(readWorkflowCheckpointRecipe("small")), 2)
    dir <- integrationSummarySave(manager, root)$dir
    input <- readIntegrationInput(dir)
    expect_equal(run$value$best_k, 2)
    expect_equal(input$manifest$correction$k, 2)
    expectIntegrationRemovedFactors(dir, "proteomics", 2L, input = input)
    # Proteins with a missing run were corrected one by one and are left out of the factors (gate G2).
    expect_lt(input$manifest$views[[1L]]$removed_factors$features, input$views$proteomics$rows)
})

test_that("KIT-184 K4: after the correlation filter removed a run, the removed factors are approximate, with their share", {
    root <- withr::local_tempdir()
    manager <- integrationProteinOwner(root, "DIA")
    withr::defer(closeDIAStateAdapter(manager))
    protein <- makeWorkflowCheckpointProtein(readWorkflowCheckpointRecipe("small"))
    samples <- setdiff(names(protein@protein_quant_table), protein@protein_id_column)
    removed_run <- samples[[5L]]
    # Otherwise complete proteins made missing in the run the filter will remove.
    complete <- which(stats::complete.cases(protein@protein_quant_table[samples]))
    expect_gte(length(complete), 6L)
    protein@protein_quant_table[complete[1:3], removed_run] <- NA_real_
    run <- integrationRuvStep(manager, protein, 2)
    final <- withDIAStateScope(manager, function() readDIAExecutionComponents(manager, run$value$final))
    filtered <- final
    filtered@protein_quant_table <- filtered@protein_quant_table[, names(filtered@protein_quant_table) != removed_run]
    filtered@design_matrix <- filtered@design_matrix[filtered@design_matrix$Run != removed_run, ]
    manager$saveState("correlation_filtered", filtered, list(min_pearson_correlation_threshold = 0.5), "correlation")
    dir <- integrationSummarySave(manager, root)$dir
    input <- readIntegrationInput(dir)
    expect_false(removed_run %in% input$samples$Run)
    expectIntegrationRemovedFactors(dir, "proteomics", 2L, status = "approximate", input = input)
    record <- input$manifest$views[[1L]]$removed_factors
    expect_gt(record$rank, 2L)
    expect_lt(record$share, 1)
    expect_true(file.exists(file.path(dir, "removed_proteomics.parquet")))
})

test_that("KIT-061: a classic route without RUV writes the normalised matrix and names the UniProt gene as display name", {
    root <- withr::local_tempdir()
    manager <- integrationProteinOwner(root, "LFQ")
    withr::defer(closeDIAStateAdapter(manager))
    protein <- makeWorkflowCheckpointProtein(readWorkflowCheckpointRecipe("small"))
    manager$saveState("normalized", protein, list(norm_method = "cyclicloess", ruv_mode = "skip", ruv_applied = FALSE),
        "skipped")
    manager$saveState("correlation_filtered", protein, list(min_pearson_correlation_threshold = 0.5), "correlation")
    ids <- protein@protein_quant_table$Protein.Group
    uniprot <- data.frame(Entry = ids[1:2], gene_names = c("ALPHA BETA", "GAMMA:x"), stringsAsFactors = FALSE)
    workflow <- list2env(list(state_manager = manager, contrasts_tbl = NULL, config_list = NULL, uniprot_dat_cln = uniprot))
    input <- readIntegrationInput(integrationSummarySave(manager, root, workflow)$dir)
    expect_identical(input$manifest$workflow, "LFQ")
    expect_identical(input$manifest$matrix_role, "normalised")
    expect_null(input$manifest$correction)
    # KIT-184 K4: without RUV there are no removed factors: no file, no field.
    expect_null(input$manifest$views[[1L]]$removed_factors)
    expect_false("removed_proteomics.parquet" %in% names(input$manifest$files))
    features <- input$views$proteomics$features
    expect_identical(features$display_name, c("ALPHA", "GAMMA", ids[-(1:2)]))
    expect_identical(features$gene_names[1:2], c("ALPHA BETA", "GAMMA:x"))
    expect_true(all(is.na(features$gene_names[-(1:2)])))
})

test_that("KIT-061 (K061-1): after a classic re-run with RUV skipped the head's normalised matrix is written, not the old correction", {
    root <- withr::local_tempdir()
    manager <- integrationProteinOwner(root, "DIA")
    withr::defer(closeDIAStateAdapter(manager))
    protein <- makeWorkflowCheckpointProtein(readWorkflowCheckpointRecipe("small"))
    product <- saveDIAStateValue(manager$.context, list(best_k = 3, controls = rep(TRUE, 7), percentage = 5,
        parameters = list(grouping = "group")), "normalization-science", codec = "rds")
    checkpoint <- checkpointDIAExecutionResult(manager, protein)
    checkpoint$products <- list(normalization_science = product)
    commitDIAReferenceState(manager, "ruv_corrected", checkpoint, list(norm_method = "cyclicloess", ruv_mode = "automatic",
        ruv_k = 3, percentage_as_neg_ctrl = 5), "ruv")
    manager$saveState("correlation_filtered", protein, list(min_pearson_correlation_threshold = 0.5), "correlation")
    # The normalisation run again with RUV skipped (its name is new, so it is appended after the correlation filter).
    skipped <- protein
    skipped@protein_quant_table[[2]] <- skipped@protein_quant_table[[2]] + 1
    manager$saveState("normalized", skipped, list(norm_method = "cyclicloess", ruv_mode = "skip", ruv_applied = FALSE),
        "skipped")
    manifest <- readIntegrationInput(integrationSummarySave(manager, root)$dir)$manifest
    expect_identical(manifest$source_state$name, "normalized")
    expect_identical(manifest$source_state$entity_id, getDIAStateReference(manager, "normalized")$entity_id)
    expect_identical(manifest$matrix_role, "normalised")
    expect_null(manifest$correction)
    expect_null(manifest$views[[1L]]$removed_factors)
    expect_false("removed_proteomics.parquet" %in% names(manifest$files))
    # Then the correlation filter on it: its own lineage has no RUV state.
    manager$saveState("correlation_filtered", skipped, list(min_pearson_correlation_threshold = 0.5), "correlation")
    dir <- integrationSummarySave(manager, root)$dir
    input <- readIntegrationInput(dir)
    expect_identical(input$manifest$source_state$name, "correlation_filtered")
    expect_identical(input$manifest$source_state$entity_id, getDIAStateReference(manager, "correlation_filtered")$entity_id)
    expect_identical(input$manifest$matrix_role, "normalised")
    expect_null(input$manifest$correction)
    expect_identical(integrationInputMatrix(dir, "proteomics", input)[[2]], skipped@protein_quant_table[[2]])
    # A completed enrichment on the head (a result-class state, as the enrichment module commits it): the walk reads
    # through it to the correlation filter.
    results <- methods::new("EnrichmentResults", contrasts = tibble::tibble(contrasts = "treated-control"),
        enrichment_data = list(`treated-control` = list(up = data.frame(term_id = "GO:0008150", p_value = 0.01))),
        enrichment_plots = list(), enrichment_plotly = list(), enrichment_summaries = list())
    manager$saveState("enrichment_completed", results, list(q_cutoff = 0.05), "fixture enrichment")
    expect_identical(manager$current_state, "enrichment_completed")
    after <- readIntegrationInput(integrationSummarySave(manager, root)$dir)$manifest
    expect_identical(after$source_state$name, "correlation_filtered")
    expect_identical(after$source_state$entity_id, input$manifest$source_state$entity_id)
    expect_identical(after$matrix_role, "normalised")
})

test_that("KIT-061: a session with no final matrix writes nothing and says so", {
    root <- withr::local_tempdir()
    manager <- integrationProteinOwner(root, "DIA")
    withr::defer(closeDIAStateAdapter(manager))
    manager$saveState("imputed", makeWorkflowCheckpointProtein(readWorkflowCheckpointRecipe("small")), list(), "imputed")
    # K061-3: an earlier save's input is removed, so a stale input cannot validate.
    writeIntegrationInputFixture(file.path(root, "integration", "proteomics"))
    run <- integrationSummarySave(manager, root)
    expect_false(any(grepl("Saved Integration", run$shown)))
    expect_identical(run$saved$integrationFile, "not written (no final proteomics matrix yet)")
    expect_false(dir.exists(run$dir))
})

# DIA_LIMPA (D14): the route's runner on the KIT-053/056 fixture (test-workflow-analysis-record.R), in both modes.
integrationLimpaRuv <- function(root, ruv_mode) {
    manager <- integrationProteinOwner(root, "DIA_LIMPA")
    rolled <- rollUpPrecursorToPeptide(makeWorkflowCheckpointPeptide(readWorkflowCheckpointRecipe("small")),
        core_utilisation = NA)
    quantified <- suppressMessages(proteinMissingValueImputationLimpa(log2TransformPeptideMatrix(calcPeptideMatrix(rolled)),
        dpc_slope = 0.8, verbose = FALSE, chunk = 1000, engine = "limpa"))
    setDIAStepProvenance(manager, "design_builder")
    manager$saveState("protein_dpc_quant", quantified, list(), "fixture")
    paths <- list(source_dir = file.path(root, "src"), protein_qc_dir = file.path(root, "pqc"),
        da_output_dir = file.path(root, "da"), publication_graphs_dir = file.path(root, "pub"),
        time_dir = file.path(root, "time"))
    for (p in paths) dir.create(p, recursive = TRUE, showWarnings = FALSE)
    manager$.context$experiment_paths <- paths
    workflow <- list2env(list(state_manager = manager, qc_params = list(protein_qc = list()), config_list = list()))
    norm_data <- new.env(parent = emptyenv())
    correlation <- NULL
    withDIAStateScope(manager, function() {
        v1 <- computeDIAReferenceProteinNormalisation(manager, getDIAStateReference(manager), list(norm_method = "none"))
        setDIAStepProvenance(manager, "protein_between_sample_normalisation", "whole", "norm_method", NULL,
            values = list(norm_method = "none"))
        finishDIAReferenceProteinStep(workflow, v1, publish = FALSE)
        v2 <- computeDIAReferenceProteinPeptideCountFilter(manager, getDIAStateReference(manager),
            list(num_peptides_per_protein_thresh = 1, num_peptidoforms_per_protein_thresh = 2))
        setDIAStepProvenance(manager, "protein_peptide_count_filter", "whole", c("num_peptides_per_protein_thresh",
            "num_peptidoforms_per_protein_thresh"), NULL, values = v2$config[1:2])
        finishDIAReferenceProteinStep(workflow, v2, publish = FALSE)
        source <- getDIAStateReference(manager)
        pairs <- runSavedOperation(manager, diaPearsonPairsOperation(), source, list(group = "group",
            tech_rep_remove_regex = "pool"))
        selection <- runSavedOperation(manager, diaPearsonSelectionOperation(), source, list(pairs = pairs$reference$entity_id,
            threshold = 0.5, correlations = pairs$value), dependencies = list(pairs$reference))
        commitDIAReferenceState(manager, "correlation_filtered", readDIAStateMetadata(manager, selection$value)$data,
            list(min_pearson_correlation_threshold = 0.5), "correlation")
        norm_data$correlation_filtered_obj <- getDIAStateReference(manager)
        correlation <<- norm_data$correlation_filtered_obj
        # K061-1: several modes in turn are the tab's mode switches.
        for (mode in ruv_mode) suppressMessages(suppressWarnings(runDIAReferenceRuvCovariates(list(ruv_mode = mode,
            color_variable = NULL, shape_variable = NULL), workflow, norm_data, paths, "proteomics", "fixture",
            function(reference) NULL, function() "group")))
    })
    list(manager = manager, adjusted = file.path(paths$protein_qc_dir, "ruv_adjusted_protein_matrix.tsv"),
        correlation = correlation)
}

integrationLimpaExpected <- function(path) {
    id <- names(readr::read_tsv(path, n_max = 0L, show_col_types = FALSE))[[1]]
    types <- readr::cols(.default = readr::col_double())
    types$cols[[id]] <- readr::col_character()
    expected <- as.data.frame(readr::read_tsv(path, col_types = types, progress = FALSE))
    names(expected)[[1]] <- "feature_id"
    expected
}

integrationStateTable <- function(manager) {
    withDIAStateScope(manager, function() {
        readDIAExecutionComponents(manager, fields = c("protein_quant_table", "design_matrix"))@protein_quant_table
    })
}

test_that("KIT-061 (D14): DIA_LIMPA in covariates mode writes the RUV-adjusted product, not the state's uncorrected table", {
    skip_if_not_installed("limpa")
    skip_if_not_installed("mixOmics")
    root <- withr::local_tempdir()
    fixture <- integrationLimpaRuv(root, "covariates")
    manager <- fixture$manager
    withr::defer(closeDIAStateAdapter(manager))
    expect_identical(manager$current_state, "ruv_covariates")
    saved_dir <- integrationSummarySave(manager, root)$dir
    input <- readIntegrationInputWithMatrices(saved_dir)
    manifest <- input$manifest
    expect_identical(manifest$workflow, "DIA_LIMPA")
    expect_identical(manifest$matrix_role, "normalised")
    expect_identical(manifest$source_state$name, "ruv_covariates")
    expect_identical(manifest$source_state$product$type, "ruv-adjusted-science")
    expect_identical(manifest$correction$id, "C2")
    expect_identical(manifest$correction$ruv_mode, "covariates")
    expect_equal(manifest$correction$k, readDIAStateMetadata(manager)$config$ruv_k)
    expect_identical(manifest$correction$grouping_variable, "group")
    # The corrected matrix is the published adjusted matrix (the step's own file), and differs from the state's table.
    expect_equal(input$views$proteomics$corrected_matrix, integrationLimpaExpected(fixture$adjusted), tolerance = 1e-12)
    state <- integrationStateTable(manager)
    expect_false(isTRUE(all.equal(unname(as.matrix(input$views$proteomics$corrected_matrix[-1])), unname(as.matrix(state[-1])))))
    # KIT-153: the data set's own matrix is the table the RUV step read: the state's own (uncorrected) table, which is the
    # correlation filter's (the step's source).
    expect_identical(manifest$normalisation$source$name, "ruv_step_input")
    names(state)[[1]] <- "feature_id"
    integrationLimpaExpectSameTable(input$views$proteomics$matrix, state)
    integrationLimpaExpectSameTable(input$views$proteomics$matrix, integrationLimpaTable(manager, fixture$correlation))
    # The design is the state's, with its W columns.
    expect_true(all(paste0("W", seq_len(manifest$correction$k)) %in% names(input$samples)))
    # KIT-184 K4 (gate G1): the written factors are what was removed, and they span the design's W1..Wk.
    factors <- expectIntegrationRemovedFactors(saved_dir, "proteomics", manifest$correction$k, input = input)
    w <- vapply(paste0("W", seq_len(manifest$correction$k)), function(name) as.numeric(input$samples[[name]]),
        numeric(nrow(input$samples)))
    w <- scale(matrix(w, nrow = nrow(input$samples)), scale = FALSE)
    expect_lt(max(abs(w - factors %*% crossprod(factors, w))) / max(abs(w)), 1e-8)
})

test_that("KIT-061: DIA_LIMPA in correction mode writes the corrected state's table", {
    skip_if_not_installed("limpa")
    skip_if_not_installed("mixOmics")
    root <- withr::local_tempdir()
    fixture <- integrationLimpaRuv(root, "correction")
    manager <- fixture$manager
    withr::defer(closeDIAStateAdapter(manager))
    expect_identical(manager$current_state, "ruv_corrected")
    saved_dir <- integrationSummarySave(manager, root)$dir
    input <- readIntegrationInputWithMatrices(saved_dir)
    manifest <- input$manifest
    expect_identical(manifest$matrix_role, "normalised")
    expect_identical(manifest$source_state$name, "ruv_corrected")
    expect_null(manifest$source_state$product)
    expect_identical(manifest$correction$id, "C2")
    expect_identical(manifest$correction$ruv_mode, "correction")
    state <- integrationStateTable(manager)
    names(state)[[1]] <- "feature_id"
    expect_identical(input$views$proteomics$corrected_matrix, state)
    # The same fit as the covariates mode's adjusted matrix.
    expect_equal(input$views$proteomics$corrected_matrix, integrationLimpaExpected(fixture$adjusted), tolerance = 1e-12)
    expect_false(any(grepl("^W[0-9]+$", names(input$samples))))
    # KIT-184 K4 (gate G1): the written factors are what the correction removed.
    expectIntegrationRemovedFactors(saved_dir, "proteomics", manifest$correction$k, input = input)
    # KIT-153: the data set's own matrix is the table the RUV step read, the correlation filter's.
    integrationLimpaExpectSameTable(input$views$proteomics$matrix, integrationLimpaTable(manager, fixture$correlation))
})

test_that("KIT-061 (K061-1): DIA_LIMPA switched covariates -> correction -> covariates writes the head's adjusted product", {
    skip_if_not_installed("limpa")
    skip_if_not_installed("mixOmics")
    root <- withr::local_tempdir()
    fixture <- integrationLimpaRuv(root, c("covariates", "correction", "covariates"))
    manager <- fixture$manager
    withr::defer(closeDIAStateAdapter(manager))
    expect_identical(manager$current_state, "ruv_covariates")
    expect_true("ruv_corrected" %in% diaStateHistory(manager))
    input <- readIntegrationInputWithMatrices(integrationSummarySave(manager, root)$dir)
    manifest <- input$manifest
    expect_identical(manifest$source_state$name, "ruv_covariates")
    expect_identical(manifest$source_state$entity_id, getDIAStateReference(manager)$entity_id)
    expect_identical(manifest$source_state$product$type, "ruv-adjusted-science")
    expect_identical(manifest$correction$state, "ruv_covariates")
    expect_identical(manifest$correction$ruv_mode, "covariates")
    expect_equal(input$views$proteomics$corrected_matrix, integrationLimpaExpected(fixture$adjusted), tolerance = 1e-12)
    expect_true(all(paste0("W", seq_len(manifest$correction$k)) %in% names(input$samples)))
})

test_that("KIT-061 (K061-7): DIA_LIMPA correlation -> covariates writes the corrected product; a correlation re-apply then writes the unadjusted table, normalised", {
    skip_if_not_installed("limpa")
    skip_if_not_installed("mixOmics")
    root <- withr::local_tempdir()
    fixture <- integrationLimpaRuv(root, "covariates")
    manager <- fixture$manager
    withr::defer(closeDIAStateAdapter(manager))
    corrected <- readIntegrationInputWithMatrices(integrationSummarySave(manager, root)$dir)
    expect_identical(corrected$manifest$source_state$name, "ruv_covariates")
    expect_identical(corrected$manifest$matrix_role, "normalised")
    expect_identical(corrected$manifest$correction$id, "C2")
    expect_equal(corrected$views$proteomics$corrected_matrix, integrationLimpaExpected(fixture$adjusted), tolerance = 1e-12)
    # The correlation re-applied from its pre-RUV source (`diaCorrelationSource` on DIA_LIMPA), committed on the
    # covariates head as the tab commits it: its record's previous state is ruv_covariates.
    withDIAStateScope(manager, function() {
        commit <- diaCorrelationCommit(manager, fixture$correlation, 0.5, FALSE)
        commitDIAReferenceState(manager, commit$state_name, commit$checkpoint, commit$config, commit$description)
    })
    expect_identical(manager$current_state, "correlation_filtered")
    expect_identical(readDIAStateMetadata(manager)$previous_state, "ruv_covariates")
    input <- readIntegrationInputWithMatrices(integrationSummarySave(manager, root)$dir)
    manifest <- input$manifest
    expect_identical(manifest$source_state$name, "correlation_filtered")
    expect_identical(manifest$matrix_role, "normalised")
    expect_null(manifest$correction)
    expect_null(manifest$source_state$product)
    state <- integrationStateTable(manager)
    names(state)[[1]] <- "feature_id"
    expect_identical(input$views$proteomics$matrix, state)
    expect_false(input$views$proteomics$corrected)
    expect_false(isTRUE(all.equal(unname(as.matrix(input$views$proteomics$matrix[-1])),
        unname(as.matrix(corrected$views$proteomics$corrected_matrix[-1])))))
})

test_that("KIT-061 (D2): the normalisation tab's S4 compatibility download no longer copies its RDS to integration/", {
    root <- withr::local_tempdir()
    paths <- list(protein_qc_dir = file.path(root, "protein_qc"), integration_dir = file.path(root, "integration"))
    for (path in paths) dir.create(path, recursive = TRUE)
    testthat::local_mocked_bindings(downloadHandler = function(filename, content, ...) list(filename = filename,
        content = content), .package = "shiny")
    testthat::local_mocked_bindings(
        resolveProtNormCorrelationResultFilenames = function(...) list(rdsFilename = "normalised_results_cln_with_replicates.RDS"),
        exportDIAS4Compatibility = function(manager, path) saveRDS("compatibility object", path))
    handler <- protNormS4DownloadHandler(list(state_manager = NULL, ruv_optimization_result = NULL),
        list(ruv_optimization_result = NULL), paths, "proteomics", "kit061")
    # The download to the user is unchanged.
    expect_identical(handler$filename(), "proteomics_kit061_final_s4.RDS")
    delivered <- file.path(root, "download.RDS")
    handler$content(delivered)
    expect_identical(readRDS(delivered), "compatibility object")
    expect_true(file.exists(file.path(paths$protein_qc_dir, "normalised_results_cln_with_replicates.RDS")))
    # integration/ holds only the integration input.
    expect_identical(list.files(paths$integration_dir, all.files = TRUE, no.. = TRUE), character())
})
