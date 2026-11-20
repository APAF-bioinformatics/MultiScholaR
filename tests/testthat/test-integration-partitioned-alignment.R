# KIT-149: the integration model's alignment in row pages (`integrationAlignViewsPartitioned`, the declaration's
# `partitioned` variant) against the whole-matrix alignment (`integrationAlignViews`). Declared before the first run:
#   - on a fixture built here (audit ENG-06) of more than 10,000 features and more than 64 runs, the partitioned
#     alignment read from the store is `identical()` to the whole-matrix one under both scaling settings and both
#     missing-sample rules, with the store's pages of 10,000 rows and blocks of 64 columns exercised;
#   - the fixture plants sixty features of equal variance across the page boundary (rows 9,971 to 10,030) and a cap
#     (545) that cuts inside them: 500 features are more variable, so exactly the first 45 tied features in row order
#     are kept, 30 from the first page and 15 from the second;
#   - internal standards on both pages are dropped by their row position; sparse and constant features on both pages
#     are dropped and counted as the whole path counts them;
#   - the model trained through either variant gets the same views and gives the same product; the variant is not part
#     of the identity (a second variant on the same inputs is a hit).
# The default variant was switched to `partitioned` once this file had passed on a runner (plan item C; run-tests
# 44668558288).

intPartWriteInput <- function(root, omic, view, values, samples, standards = NULL) {
    matrix <- data.frame(feature_id = rownames(values), values, check.names = FALSE, stringsAsFactors = FALSE)
    features <- data.frame(feature_id = rownames(values), display_name = tolower(rownames(values)), stringsAsFactors = FALSE)
    if (!is.null(standards)) features$is_internal_standard <- standards
    manifest <- list(omic_type = omic, workflow = "TEST", matrix_role = "normalised", scale = "log2", correction = NULL,
        source_state = list(name = "kit149-test"), sample_id_column = "Run", group_column = "group",
        technical_replicate_column = "tech_rep_group", biological_sample_column = "tech_rep_group")
    dir <- file.path(root, "integration", omic)
    writeIntegrationInput(dir, manifest, samples, stats::setNames(list(matrix), view), stats::setNames(list(features), view))
    dir
}

intPartSamples <- function(runs, biological) {
    data.frame(Run = runs, group = sub("_.*$", "", biological), replicates = as.integer(sub("^.*_", "", biological)),
        tech_rep_group = biological, stringsAsFactors = FALSE)
}

# The fixture (ENG-06): lipidomics, 12,000 features x 70 runs (35 biological samples, two runs each), and proteomics,
# 40 features x 34 runs (one per sample; B_17 absent, so the two missing-sample rules differ).
intPartFixture <- function(.local_envir = parent.frame()) {
    testthat::skip_if_not(workflowOwnerCapable(), "the integration store needs qs2, DuckDB and Arrow")
    root <- withr::local_tempdir(.local_envir = .local_envir)
    biological <- sprintf("%s_%d", rep(c("A", "B"), c(18L, 17L)), c(seq_len(18L), seq_len(17L)))
    lipid_samples <- rep(biological, each = 2L)
    runs <- sprintf("l%03d", seq_along(lipid_samples))
    n <- 12000L
    tie <- 9971:10030
    sparse <- c(5L, 9000L, 9965L, 11500L)
    constant <- c(7L, 9500L, 10500L, 11900L)
    standards <- c(100L, 9960L, 10040L, 11000L)
    special <- c(tie, sparse, constant, standards)
    values <- withr::with_seed(1149L, {
        noise <- matrix(stats::rnorm(n * length(runs), 20, 1), n, length(runs))
        noise[sample.int(length(noise), floor(0.02 * length(noise)))] <- NA
        # 500 features far more variable than the rest (a per-sample effect), on both pages.
        top <- sort(c(sample(setdiff(1:9900, special), 300L), sample(setdiff(10101:11900, special), 200L)))
        effect <- matrix(stats::rnorm(length(top) * length(biological), 0, 10), length(top), length(biological))
        noise[top, ] <- 20 + effect[, match(lipid_samples, biological)] + stats::rnorm(length(top) * length(runs), 0, 0.1)
        list(values = noise, top = top)
    })
    top <- values$top
    values <- values$values
    # Sixty features of one variance across the page boundary: the same value per biological sample (17 or 23).
    level <- ifelse(seq_along(biological) %% 2L == 0L, 23, 17)[match(lipid_samples, biological)]
    values[tie, ] <- matrix(level, length(tie), length(runs), byrow = TRUE)
    # Sparse (one biological sample observed) and constant features on both pages.
    values[sparse, ] <- NA
    values[sparse, 1:2] <- 21
    values[constant, ] <- 20
    # Internal standards on both pages, more variable than the top features: the cap would take them if they were kept.
    values[standards, ] <- 20 + matrix(withr::with_seed(2149L, stats::rnorm(length(standards) * length(runs), 0, 30)),
        length(standards))
    dimnames(values) <- list(sprintf("LIP%05d", seq_len(n)), runs)
    flags <- seq_len(n) %in% standards
    lipid <- intPartWriteInput(root, "lipidomics", "lipidomics_LCMS_Pos", values, intPartSamples(runs, lipid_samples), flags)
    protein_samples <- setdiff(biological, "B_17")
    protein <- matrix(withr::with_seed(3149L, stats::rnorm(40L * length(protein_samples), 25, 2)), 40L,
        dimnames = list(sprintf("PRO%03d", 1:40), sprintf("p%02d", seq_along(protein_samples))))
    proteomics <- intPartWriteInput(root, "proteomics", "proteomics", protein, intPartSamples(colnames(protein), protein_samples))
    list(root = root, dirs = c(lipid, proteomics), n = n, runs = length(runs), tie = sprintf("LIP%05d", tie),
        top = sprintf("LIP%05d", top), standards = sprintf("LIP%05d", standards), cap = 545L)
}

intPartStore <- function(fixture, name = "project", .local_envir = parent.frame()) {
    manager <- newIntegrationStateOwner(file.path(fixture$root, name))
    withr::defer(closeWorkflowStateOwner(manager), envir = .local_envir)
    importIntegrationInputs(manager, fixture$dirs)
    manager
}

intPartTrainer <- function(seen = NULL) {
    function(views, options) {
        if (is.environment(seen)) seen$views <- views
        samples <- colnames(views[[1L]])
        k <- options$num_factors
        labels <- paste0("Factor", seq_len(k))
        list(factors = matrix(seq_len(length(samples) * k) / 10, length(samples), k, dimnames = list(samples, labels)),
            weights = lapply(views, function(v) matrix(rowMeans(v, na.rm = TRUE), nrow(v), k, dimnames = list(rownames(v), labels))),
            r2_per_factor = matrix(1, k, length(views), dimnames = list(labels, names(views))),
            r2_total = stats::setNames(rep(k, length(views)), names(views)), elbo = c(-20, -10, NaN), iterations = 2L,
            options = options, versions = list(MOFA2 = "stand-in"), messages = character(),
            model_file = list(bytes = as.raw(1:3), sha256 = strrep("c", 64L), size = 3))
    }
}

test_that("KIT-149 (ENG-06): over 10,000 features and 64 runs, ties across the page boundary and a cap inside them, the partitioned alignment is identical to the whole", {
    fixture <- intPartFixture()
    manager <- intPartStore(fixture)
    source <- getWorkflowStateReference(manager, "integration_inputs")
    spec <- withWorkflowStateScope(manager, function() readWorkflowStateMetadata(manager, source)$data)
    expect_identical(as.integer(spec$components$views$lipidomics_LCMS_Pos$metadata$shape), c(fixture$n, fixture$runs + 1L))
    whole_object <- withWorkflowStateScope(manager, function() readWorkflowExecutionComponents(manager, source, fields = TRUE))
    paged_object <- withWorkflowStateScope(manager, function() readWorkflowExecutionComponents(manager, source,
        fields = c("features", "design_matrix", "inputs")))
    # The partitioned path is given no values: its data sets are zero-row prototypes.
    expect_identical(vapply(paged_object@views, nrow, integer(1)), c(lipidomics_LCMS_Pos = 0L, proteomics = 0L))
    store_pages <- integrationStateViewPages(manager, source)
    pages <- list()
    traced <- function(view, columns, kernel) store_pages(view, columns, function(page, offset) {
        pages[[length(pages) + 1L]] <<- c(rows = nrow(page), columns = ncol(page), offset = offset)
        kernel(page, offset)
    })
    for (scaling in c("none", "feature_zscore")) {
        for (missing_samples in c("keep_missing", "drop_incomplete")) {
            info <- paste(scaling, missing_samples)
            parameters <- integrationResolveMofaParameters(whole_object, list(feature_selection = "top_variance",
                features_per_view = fixture$cap, scaling = scaling, missing_samples = missing_samples))
            whole <- integrationAlignViews(whole_object, parameters, prefix_rows = TRUE)
            pages <- list()
            partitioned <- integrationAlignViewsPartitioned(paged_object, parameters, traced, prefix_rows = TRUE)
            expect_identical(partitioned, whole, info = info)
            # Two passes over the large data set, each in two pages (10,000 and 2,000 rows) of every column (more than
            # 64, so read in 64-column blocks); one page per pass for the small one.
            big <- Filter(function(page) page[["columns"]] == fixture$runs + 1L, pages)
            expect_identical(unname(vapply(big, function(page) page[["rows"]], numeric(1))), c(10000, 2000, 10000, 2000), info = info)
            expect_identical(unname(vapply(big, function(page) page[["offset"]], numeric(1))), c(0, 10000, 0, 10000), info = info)
            # The cap cuts inside the tie: the 500 more variable features, then the first 45 tied ones in row order
            # (30 before the page boundary, 15 after it); the internal standards are not among them.
            kept <- sub("^lipidomics_LCMS_Pos:", "", rownames(partitioned$views$lipidomics_LCMS_Pos))
            expect_identical(length(kept), fixture$cap, info = info)
            expect_identical(kept[kept %in% fixture$tie], fixture$tie[1:45], info = info)
            expect_identical(sort(kept[!kept %in% fixture$tie]), sort(fixture$top), info = info)
            summary <- partitioned$view_summary
            expect_identical(summary$features_input, c(fixture$n, 40L), info = info)
            expect_identical(summary$internal_standards_dropped, c(4L, 0L), info = info)
            expect_identical(summary$sparse_dropped[[1L]], 4L, info = info)
            expect_identical(summary$constant_dropped[[1L]], if (scaling == "feature_zscore") 4L else 0L, info = info)
            expect_identical(summary$samples_present, if (missing_samples == "keep_missing") c(35L, 34L) else c(34L, 34L),
                info = info)
        }
    }
    # Without a cap every feature left is aligned, in row order, from both pages.
    parameters <- integrationResolveMofaParameters(whole_object, list(views = c("proteomics", "lipidomics_LCMS_Pos")))
    expect_identical(integrationAlignViewsPartitioned(paged_object, parameters, store_pages),
        integrationAlignViews(whole_object, parameters))
    # The counts alone (values = FALSE): the same sample table and summary, no values.
    counts <- integrationAlignViewsPartitioned(paged_object, parameters, store_pages, values = FALSE)
    expect_null(counts$views)
    expect_identical(counts[c("samples", "view_summary", "shared_samples", "notes")],
        integrationAlignViews(whole_object, parameters)[c("samples", "view_summary", "shared_samples", "notes")])
})

test_that("KIT-149: the model trained through either variant gets the same data sets and gives the same product; the variant is not identity", {
    fixture <- intPartFixture()
    settings <- list(feature_selection = "top_variance", features_per_view = fixture$cap)
    runs <- list()
    seen <- list()
    for (variant in c("whole", "partitioned")) {
        manager <- intPartStore(fixture, paste0("project-", variant))
        holder <- new.env()
        run <- runIntegrationModel(manager, c(settings, list(variant = variant)), trainer = intPartTrainer(holder))
        expect_identical(run$status, "computed", info = variant)
        seen[[variant]] <- holder$views
        value <- run$value
        tables <- lapply(stats::setNames(c("weights", "aligned_values"), c("weights", "aligned_values")), function(kind) {
            lapply(stats::setNames(names(value[[kind]]), names(value[[kind]])), function(view)
                integrationModelViewTable(manager, value, kind, view))
        })
        runs[[variant]] <- c(value[c("factors", "variance_explained", "variance_total", "samples", "view_summary",
            "shared_samples", "notes", "elbo")], tables)
        # The other variant on the same inputs and parameters is the same model: a hit.
        other <- setdiff(c("whole", "partitioned"), variant)
        again <- runIntegrationModel(manager, c(settings, list(variant = other)), trainer = intPartTrainer())
        expect_identical(again$status, "hit", info = variant)
        expect_identical(again$reference$entity_id, run$reference$entity_id, info = variant)
    }
    expect_identical(seen$partitioned, seen$whole)
    expect_identical(runs$partitioned, runs$whole)
    expect_identical(nrow(runs$whole$aligned_values$lipidomics_LCMS_Pos), fixture$cap)
    # A variant the registry does not declare is refused in words.
    expect_error(runIntegrationModel(manager, list(variant = "rows"), trainer = intPartTrainer()),
        class = "integration_parameters_invalid")
    entry <- workflowOperationEntry(INTEGRATION_MODEL_OPERATION)
    expect_identical(entry$variants, c("whole", "partitioned"))
    # The default, once the identity was shown on a runner (run-tests 44668558288).
    expect_identical(entry$default_variant, "partitioned")
})

# KIT-149, audit ENG-05: no data set is read whole where nothing needs its values. Every whole-table read of the store
# (`materializeWorkflowTable`, and its alias `materializeDIATable` the checkpoint reader calls) is traced; the data sets'
# tables are named by the inputs state. The tier-1 corpus (four omics, six data sets) is copied into a project.
intPartCorpusProject <- function(.local_envir = parent.frame()) {
    corpus <- normalizePath(testthat::test_path("..", "testdata", "integration", "corpus"), mustWork = FALSE)
    testthat::skip_if_not(dir.exists(corpus), "the tier-1 corpus is not present")
    testthat::skip_if_not(workflowOwnerCapable(), "the integration store needs qs2, DuckDB and Arrow")
    root <- withr::local_tempdir(.local_envir = .local_envir)
    paths <- list(base_dir = root, integration_dir = file.path(root, "integration"),
        source_dir = file.path(root, "scripts", "integration"), results_dir = file.path(root, "results", "integration"),
        multiomic_inputs_dir = file.path(root, "results", "integration", "multiomic_inputs"))
    dir.create(paths$source_dir, recursive = TRUE)
    for (omic in c("lipidomics", "metabolomics", "proteomics", "transcriptomics")) {
        dir <- file.path(paths$integration_dir, omic)
        dir.create(dir, recursive = TRUE)
        stopifnot(all(file.copy(list.files(file.path(corpus, omic), full.names = TRUE), dir)))
    }
    list(root = root, paths = paths)
}

intPartTraceReads <- function(.local_envir = parent.frame()) {
    reads <- new.env()
    reads$ids <- character()
    trace <- function(original) {
        # Forced here: a promise would be read after the mock is in place and call the mock itself.
        force(original)
        function(context, reference, ...) {
            reads$ids <- c(reads$ids, reference$entity_id)
            original(context, reference, ...)
        }
    }
    testthat::local_mocked_bindings(materializeWorkflowTable = trace(materializeWorkflowTable),
        materializeDIATable = trace(materializeDIATable), .env = .local_envir)
    reads
}

intPartViewTables <- function(manager) withWorkflowStateScope(manager, function() {
    source <- getWorkflowStateReference(manager, "integration_inputs")
    vapply(readWorkflowStateMetadata(manager, source)$data$components$views, function(table) table$entity_id, character(1))
})

test_that("KIT-149 (ENG-05): the GUI's import, the Sample Matching preview, the match file and the confirmation read no data set whole", {
    project <- intPartCorpusProject()
    bus <- do.call(shiny::reactiveValues, omicWorkflowInitialBus(integrationGuiConfig()$shell))
    withr::defer(if (!is.null(shiny::isolate(bus$state_manager))) closeWorkflowStateOwner(shiny::isolate(bus$state_manager)))
    reads <- intPartTraceReads()
    # Not named `scan`: the import module's own `scan` would shadow it inside testServer.
    all_views <- integrationScanProject(project$paths$integration_dir)$views$view
    shiny::testServer(mod_integration_import_server, args = list(workflow_data = bus, paths = project$paths), {
        session$setInputs(views = all_views, import = 1)
    })
    manager <- shiny::isolate(bus$state_manager)
    expect_identical(manager$getHistory(), c("initial", "integration_inputs"))
    tables <- intPartViewTables(manager)
    expect_length(tables, 6L)
    # The import action reads the names, the design and the manifests, not a data set.
    expect_identical(intersect(reads$ids, tables), character())
    texts <- integrationGuiTexts()
    shiny::testServer(mod_integration_alignment_server, args = list(workflow_data = bus, paths = project$paths), {
        preview <- session$returned()
        expect_false(inherits(preview, "integration_alignment_refused"))
        expect_null(preview$aligned$views)
        # The match file as offered, taken back unchanged, then the confirmation.
        downloaded <- output$match_download
        expect_true(file.exists(downloaded))
        session$setInputs(match_upload = data.frame(name = "match.tsv", size = file.size(downloaded),
            type = "text/tab-separated-values", datapath = downloaded, stringsAsFactors = FALSE))
        expect_match(output$match_status$html, texts[["match_unchanged"]], fixed = TRUE)
        session$setInputs(confirm = 1)
        expect_identical(bus$tab_status$integration_alignment, "complete")
        expect_identical(intersect(reads$ids, tables), character())
        # The numbers the tab shows are the whole-matrix alignment's.
        shown <- session$returned()
        object <- integrationImportedObject(manager)
        whole <- suppressWarnings(integrationAlignViews(object, shown$parameters))
        fields <- c("samples", "view_summary", "shared_samples", "notes")
        expect_identical(shown$aligned[fields], whole[fields])
        expect_identical(bus$integration_alignment$view_summary, whole$view_summary)
        expect_identical(integrationAlignmentFacts(shown), integrationAlignmentFacts(c(list(aligned = whole),
            shown[setdiff(names(shown), "aligned")])))
    })
})

test_that("KIT-149 (ENG-05): the model reads only the selected data sets whole (whole variant), and none whole (partitioned)", {
    project <- intPartCorpusProject()
    manager <- newIntegrationStateOwner(file.path(project$paths$source_dir, "store"))
    withr::defer(closeWorkflowStateOwner(manager))
    importIntegrationInputs(manager, file.path(project$paths$integration_dir,
        c("lipidomics", "metabolomics", "proteomics", "transcriptomics")))
    tables <- intPartViewTables(manager)
    selected <- c("proteomics", "transcriptomics")
    reads <- intPartTraceReads()
    whole <- runIntegrationModel(manager, list(views = selected, variant = "whole"), trainer = intPartTrainer())
    expect_identical(whole$status, "computed")
    expect_setequal(intersect(reads$ids, tables), unname(tables[selected]))
    reads$ids <- character()
    other <- newIntegrationStateOwner(file.path(project$paths$source_dir, "store-partitioned"))
    withr::defer(closeWorkflowStateOwner(other))
    importIntegrationInputs(other, file.path(project$paths$integration_dir,
        c("lipidomics", "metabolomics", "proteomics", "transcriptomics")))
    reads$ids <- character()
    partitioned <- runIntegrationModel(other, list(views = selected, variant = "partitioned"), trainer = intPartTrainer())
    expect_identical(partitioned$status, "computed")
    expect_identical(intersect(reads$ids, intPartViewTables(other)), character())
    expect_identical(partitioned$value$view_summary, whole$value$view_summary)
})
