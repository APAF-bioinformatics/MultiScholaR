# KIT-061: the per-omic integration input (the contract, utils_workflow_integration_input.R) and the proteomics writer
# at the Summary hook (utils_workflow_dia_summary.R).

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

writeIntegrationInputFixture <- function(dir, fixture = integrationInputFixture()) {
    writeIntegrationInput(dir, fixture$manifest, fixture$samples, list(proteomics = fixture$matrix),
        list(proteomics = fixture$features), now = as.POSIXct("2026-11-17 10:00:00", tz = "UTC"))
}

test_that("KIT-061: the contract round trip writes the fixed files and reads them back validated", {
    root <- withr::local_tempdir()
    dir <- file.path(root, "integration", "proteomics")
    fixture <- integrationInputFixture()
    path <- writeIntegrationInputFixture(dir, fixture)
    expect_identical(path, file.path(dir, "integration_input.json"))
    expect_setequal(list.files(dir), c("integration_input.json", "matrix_proteomics.tsv", "features_proteomics.tsv",
        "samples.tsv"))
    # Nothing is left beside the folder (the staging folder is renamed into place).
    expect_identical(list.files(file.path(root, "integration"), all.files = TRUE, no.. = TRUE), "proteomics")
    input <- readIntegrationInput(dir)
    manifest <- input$manifest
    expect_identical(manifest$schema, "multischolar.integration_input.v1")
    expect_identical(manifest$omic_type, "proteomics")
    expect_identical(manifest$workflow, "DIA")
    expect_identical(manifest$matrix_role, "normalised")
    expect_null(manifest$correction)
    expect_identical(manifest$source_state$entity_id, strrep("a", 64))
    expect_identical(manifest$software$package, "MultiScholaR")
    expect_identical(manifest$software$version, as.character(utils::packageVersion("MultiScholaR")))
    expect_identical(manifest$written_at, "2026-11-17T10:00:00Z")
    expect_identical(manifest$views[[1]][c("name", "matrix_file", "features_file", "rows", "columns")],
        list(name = "proteomics", matrix_file = "matrix_proteomics.tsv", features_file = "features_proteomics.tsv",
            rows = 3L, columns = 4L))
    for (name in names(manifest$files)) {
        expect_identical(manifest$files[[name]], digest::digest(file = file.path(dir, name), algo = "sha256"))
    }
    expect_setequal(names(manifest$files), c("matrix_proteomics.tsv", "features_proteomics.tsv", "samples.tsv"))
    # The values, missing values included, come back exactly.
    expect_identical(input$views$proteomics$matrix, fixture$matrix)
    expect_identical(input$views$proteomics$features, fixture$features)
    expect_identical(input$samples, fixture$samples)
    # A second write replaces the folder whole.
    fixture$manifest$matrix_role <- "corrected"
    fixture$manifest$correction <- list(id = "C1", k = 2L)
    writeIntegrationInputFixture(dir, fixture)
    expect_identical(readIntegrationInput(dir)$manifest$correction$k, 2L)
    expect_identical(list.files(file.path(root, "integration"), all.files = TRUE, no.. = TRUE), "proteomics")
})

test_that("KIT-061 (K061-2): a failed publish or move-aside leaves the previous input intact", {
    root <- withr::local_tempdir()
    dir <- file.path(root, "integration", "proteomics")
    writeIntegrationInputFixture(dir)
    before <- readIntegrationInput(dir)
    fixture <- integrationInputFixture()
    fixture$manifest$matrix_role <- "corrected"
    fixture$manifest$correction <- list(id = "C1", k = 2L)
    write <- function(rename) writeIntegrationInput(dir, fixture$manifest, fixture$samples, list(proteomics = fixture$matrix),
        list(proteomics = fixture$features), renameFn = rename)
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
        file.copy(file.path(dir, "matrix_proteomics.tsv"), file.path(dir, "other.tsv"))
        manifest$files$other.tsv <- manifest$files$matrix_proteomics.tsv
        manifest$views[[1]]$matrix_file <- "other.tsv"
        manifest
    })
    expect_error(readIntegrationInput(dir), class = "integration_input_invalid", regexp = "matrix_proteomics.tsv")
    rewrite(function(manifest) {
        manifest$views[[1]]$features_file <- "samples.tsv"
        manifest
    })
    expect_error(readIntegrationInput(dir), class = "integration_input_invalid", regexp = "features_proteomics.tsv")
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
    matrix_file <- file.path(dir, "matrix_proteomics.tsv")
    original <- readLines(matrix_file)
    # One value changed: the file no longer matches its sha256.
    writeLines(sub("21.5", "21.6", original, fixed = TRUE), matrix_file)
    expect_error(readIntegrationInput(dir), class = "integration_input_invalid", regexp = "sha256")
    # One sample column dropped and the hash updated: the shape no longer matches the manifest.
    dropped <- vapply(strsplit(original, "\t", fixed = TRUE), function(cells) paste(cells[-5L], collapse = "\t"), "")
    writeLines(dropped, matrix_file)
    manifest_file <- file.path(dir, "integration_input.json")
    manifest <- jsonlite::read_json(manifest_file, simplifyVector = FALSE)
    manifest$files$matrix_proteomics.tsv <- digest::digest(file = matrix_file, algo = "sha256")
    jsonlite::write_json(manifest, manifest_file, auto_unbox = TRUE, null = "null", na = "null", digits = NA)
    expect_error(readIntegrationInput(dir), class = "integration_input_invalid", regexp = "is not 3 x 4")
    # A corrected matrix without its correction.
    writeIntegrationInputFixture(dir)
    manifest <- jsonlite::read_json(manifest_file, simplifyVector = FALSE)
    manifest$matrix_role <- "corrected"
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
    product <- saveDIAStateValue(manager$.context, list(best_k = 3, controls = controls, percentage = 5,
        parameters = list(grouping = "group")), "normalization-science", codec = "rds")
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
    run <- integrationSummarySave(manager, root)
    expect_identical(run$saved$integrationFile, file.path("proteomics", "integration_input.json"))
    expect_true("Saved Integration input" %in% run$shown)
    expect_identical(list.files(run$integration_dir, all.files = TRUE, no.. = TRUE), "proteomics")
    input <- readIntegrationInput(run$dir)
    manifest <- input$manifest
    expect_identical(manifest$workflow, "DIA")
    expect_identical(manifest$matrix_role, "corrected")
    expect_identical(manifest$source_state$name, "correlation_filtered")
    expect_identical(manifest$source_state$entity_id, getDIAStateReference(manager, "correlation_filtered")$entity_id)
    expect_null(manifest$source_state$product)
    expect_identical(manifest$correction[c("state", "id", "method", "ruv_mode", "k", "percentage_as_neg_ctrl",
        "control_count", "grouping_variable")], list(state = "ruv_corrected", id = "C1", method = "RUVIII_C_Varying",
        ruv_mode = "automatic", k = 3L, percentage_as_neg_ctrl = 5L, control_count = 7L, grouping_variable = "group"))
    expect_identical(unlist(manifest$correction$grouping_levels), c("control", "treated"))
    expected <- filtered@protein_quant_table
    names(expected)[[1]] <- "feature_id"
    expect_identical(input$views$proteomics$matrix, expected)
    expect_identical(input$samples$Run, setdiff(protein@design_matrix$Run, "run_0005"))
    expect_identical(input$samples$group, as.character(filtered@design_matrix$group))
    # No UniProt table: each display name is the feature id.
    expect_identical(input$views$proteomics$features$display_name, expected$feature_id)
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
    # Then the correlation filter on it: its own lineage has no RUV state.
    manager$saveState("correlation_filtered", skipped, list(min_pearson_correlation_threshold = 0.5), "correlation")
    input <- readIntegrationInput(integrationSummarySave(manager, root)$dir)
    expect_identical(input$manifest$source_state$name, "correlation_filtered")
    expect_identical(input$manifest$source_state$entity_id, getDIAStateReference(manager, "correlation_filtered")$entity_id)
    expect_identical(input$manifest$matrix_role, "normalised")
    expect_null(input$manifest$correction)
    expect_identical(input$views$proteomics$matrix[[2]], skipped@protein_quant_table[[2]])
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
        # K061-1: several modes in turn are the tab's mode switches.
        for (mode in ruv_mode) suppressMessages(suppressWarnings(runDIAReferenceRuvCovariates(list(ruv_mode = mode,
            color_variable = NULL, shape_variable = NULL), workflow, norm_data, paths, "proteomics", "fixture",
            function(reference) NULL, function() "group")))
    })
    list(manager = manager, adjusted = file.path(paths$protein_qc_dir, "ruv_adjusted_protein_matrix.tsv"))
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
    input <- readIntegrationInput(integrationSummarySave(manager, root)$dir)
    manifest <- input$manifest
    expect_identical(manifest$workflow, "DIA_LIMPA")
    expect_identical(manifest$matrix_role, "corrected")
    expect_identical(manifest$source_state$name, "ruv_covariates")
    expect_identical(manifest$source_state$product$type, "ruv-adjusted-science")
    expect_identical(manifest$correction$id, "C2")
    expect_identical(manifest$correction$ruv_mode, "covariates")
    expect_equal(manifest$correction$k, readDIAStateMetadata(manager)$config$ruv_k)
    expect_identical(manifest$correction$grouping_variable, "group")
    # The matrix is the published adjusted matrix (the step's own file), and differs from the state's table.
    expect_equal(input$views$proteomics$matrix, integrationLimpaExpected(fixture$adjusted), tolerance = 1e-12)
    state <- integrationStateTable(manager)
    expect_false(isTRUE(all.equal(unname(as.matrix(input$views$proteomics$matrix[-1])), unname(as.matrix(state[-1])))))
    # The design is the state's, with its W columns.
    expect_true(all(paste0("W", seq_len(manifest$correction$k)) %in% names(input$samples)))
})

test_that("KIT-061: DIA_LIMPA in correction mode writes the corrected state's table", {
    skip_if_not_installed("limpa")
    skip_if_not_installed("mixOmics")
    root <- withr::local_tempdir()
    fixture <- integrationLimpaRuv(root, "correction")
    manager <- fixture$manager
    withr::defer(closeDIAStateAdapter(manager))
    expect_identical(manager$current_state, "ruv_corrected")
    input <- readIntegrationInput(integrationSummarySave(manager, root)$dir)
    manifest <- input$manifest
    expect_identical(manifest$matrix_role, "corrected")
    expect_identical(manifest$source_state$name, "ruv_corrected")
    expect_null(manifest$source_state$product)
    expect_identical(manifest$correction$id, "C2")
    expect_identical(manifest$correction$ruv_mode, "correction")
    state <- integrationStateTable(manager)
    names(state)[[1]] <- "feature_id"
    expect_identical(input$views$proteomics$matrix, state)
    # The same fit as the covariates mode's adjusted matrix.
    expect_equal(input$views$proteomics$matrix, integrationLimpaExpected(fixture$adjusted), tolerance = 1e-12)
    expect_false(any(grepl("^W[0-9]+$", names(input$samples))))
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
    input <- readIntegrationInput(integrationSummarySave(manager, root)$dir)
    manifest <- input$manifest
    expect_identical(manifest$source_state$name, "ruv_covariates")
    expect_identical(manifest$source_state$entity_id, getDIAStateReference(manager)$entity_id)
    expect_identical(manifest$source_state$product$type, "ruv-adjusted-science")
    expect_identical(manifest$correction$state, "ruv_covariates")
    expect_identical(manifest$correction$ruv_mode, "covariates")
    expect_equal(input$views$proteomics$matrix, integrationLimpaExpected(fixture$adjusted), tolerance = 1e-12)
    expect_true(all(paste0("W", seq_len(manifest$correction$k)) %in% names(input$samples)))
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
