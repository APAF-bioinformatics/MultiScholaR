# KIT-184 (VAL-53 (1), VAL-58, VAL-62): the unwanted factors a RUV correction removed, written with the integration
# input, carried in the inputs state and put into the design of the integration enrichment's tests that have a design
# (camera, camera_floor), on any model, by the setting `adjust_ruv_factors`, which is off unless chosen (the owner's
# decision VAL-160 after the counted simulation; before it, on by default for a model trained on RUV-corrected values).
# Declared before the run (ticket section 10, revision 1.4; each test's failing value on the base
# is stated beside it there):
#   K1  the removed factors of a hand case: exact, rank 2, orthonormal, centred, sign fixed, the samples table's order;
#   K2  approximate (a third direction), not determined (too few features, no k, equal matrices, a constant per feature);
#   K3  the contract: written, listed, read back bit for bit; malformed refused; an input without a correction as before;
#   K5  the covariates at test time: averaged as the values, the sample match applied, matched by name, a basis;
#   K6  exactness: the kernel equals limma and lm; equals the same model on the corrected values; differs from option 1;
#   K7  (i) (a) the declarations are base's; (ii) without the setting no model gets covariates: KIT-183's product;
#   K8  refusals in plain words; no test run in another's place;
#   K9  (b) the replay repeats the saved form;
#   K10 the record and Methods per entry; K11 the tab and the summary's sentences; K12 the workbook;
#   K14 counts, run_complete, the share of the factor the removed factors explain; K15 the test of record by n_eq;
#   K16 the four combinations of the setting, camera_floor, the permutation, the corpus check of 6.6.
# K3's base parts, K7 (i) (b), (c) and (iii) and K9 (a) read a store made by the base commit's code: they run only inside
# docs/qa/workflow-execution/kit/kit-184/base-store-job.R's runner job (VAL-64 (2)), which names it in K184_BASE_STORE.

# --- Fixtures ---

intRfData <- function() normalizePath(testthat::test_path("..", "testdata", "integration"), mustWork = FALSE)
intRfStored <- function(name) utils::read.delim(file.path(intRfData(), "mofa_reference", paste0(name, ".tsv")),
    check.names = FALSE, stringsAsFactors = FALSE)
intRfSettings <- function() integrationGseaSettings(p_adjust_method = "BH", min_size = 5L, max_size = 500L, p_cutoff = 1,
    tie_jitter = 0)
intRfTable <- function(m) data.frame(feature_id = rownames(m), m, check.names = FALSE, stringsAsFactors = FALSE,
    row.names = NULL)
intRfUndetermined <- "integration_input_removed_factors_undetermined"

# The hand case of K1: 30 features by 8 runs; corrected = normalised - W a - c with a known 8 by 2 W (not centred); 6
# features with a missing value.
intRfHand <- function() {
    runs <- sprintf("r%d", 1:8)
    w <- cbind(c(1.5, -0.5, 1.5, -0.5, 1.5, -0.5, 1.5, -0.5), c(-3, -2, -1, 0, 0, 1, 2, 3.5))
    a <- cbind(0.3 + 0.05 * seq_len(30), -0.2 + 0.03 * (seq_len(30) %% 7))
    level <- rep(c(0, 0.7), 15L)
    normalised <- 15 + withr::with_seed(20261841L, matrix(stats::rnorm(240), 30L, 8L))
    dimnames(normalised) <- list(sprintf("F%02d", 1:30), runs)
    corrected <- normalised - a %*% t(w) - level
    dimnames(corrected) <- dimnames(normalised)
    for (i in 1:6) {
        normalised[i, i] <- NA_real_
        corrected[i, i] <- NA_real_
    }
    list(runs = runs, w = w, normalised = normalised, corrected = corrected)
}

# One omic's input folder from two matrices (one run per biological sample), as the shared writer writes it.
intRfWrite <- function(root, normalised, corrected = NULL, correction = NULL, omic = "proteomics", view = omic) {
    runs <- colnames(normalised)
    samples <- data.frame(Run = runs, group = rep(c("A", "B"), length.out = length(runs)), tech_rep_group = runs,
        stringsAsFactors = FALSE)
    features <- data.frame(feature_id = rownames(normalised), display_name = tolower(rownames(normalised)),
        stringsAsFactors = FALSE)
    manifest <- list(omic_type = omic, workflow = "TEST", matrix_role = "normalised", scale = "log2",
        normalisation = list(method = "test", source = list(name = "kit184-test")), correction = correction,
        source_state = list(name = "kit184-test"), sample_id_column = "Run", group_column = "group",
        technical_replicate_column = "tech_rep_group", biological_sample_column = "tech_rep_group")
    dir <- file.path(root, "integration", omic)
    writeIntegrationInput(dir, manifest, samples, stats::setNames(list(intRfTable(normalised)), view),
        stats::setNames(list(features), view),
        corrected = if (is.null(corrected)) list() else stats::setNames(list(intRfTable(corrected)), view))
    dir
}

# Rewrite a written input's manifest (and, when given, one of its files), keeping every listed sha256 true: the reader
# then reaches its own checks of the removed factors.
intRfTamper <- function(dir, edit = identity, file = NULL, table = NULL) {
    path <- file.path(dir, INTEGRATION_INPUT_MANIFEST)
    manifest <- jsonlite::read_json(path, simplifyVector = FALSE)
    if (!is.null(file)) {
        writeIntegrationInputParquet(table, file.path(dir, file))
        manifest$files[[file]] <- digest::digest(file = file.path(dir, file), algo = "sha256")
    }
    manifest <- edit(manifest)
    jsonlite::write_json(manifest, path, auto_unbox = TRUE, pretty = TRUE, digits = NA, null = "null", na = "null")
    invisible(dir)
}

# An input as a version before KIT-184 saved it: the removed factors' file and fields taken out.
intRfStrip <- function(dir, views = NULL) {
    intRfTamper(dir, function(manifest) {
        for (i in seq_along(manifest$views)) {
            name <- manifest$views[[i]]$name
            if (!is.null(views) && !name %in% views) next
            file <- manifest$views[[i]]$removed_factors_file
            if (!is.null(file)) {
                manifest$files[[file]] <- NULL
                unlink(file.path(dir, file))
            }
            manifest$views[[i]]$removed_factors_file <- NULL
            manifest$views[[i]]$removed_factors <- NULL
        }
        manifest
    })
}

intRfTrainer <- function() {
    function(views, options) {
        factors <- intRfStored("factors")
        f <- as.matrix(factors[-1L])
        rownames(f) <- factors$biological_sample
        f <- f[colnames(views[[1L]]), , drop = FALSE]
        w <- intRfStored("weights")
        weights <- lapply(stats::setNames(names(views), names(views)), function(view) {
            rows <- w[w$view == view, , drop = FALSE]
            m <- as.matrix(rows[-(1:2)])
            rownames(m) <- paste(view, rows$feature_id, sep = ":")
            m[rownames(views[[view]]), , drop = FALSE]
        })
        v <- intRfStored("variance_explained")
        r2 <- matrix(v$r2_percent, nrow = length(unique(v$factor)), dimnames = list(unique(v$factor), unique(v$view)))
        r2 <- r2[, names(views), drop = FALSE]
        list(factors = f, weights = weights, r2_per_factor = r2, r2_total = colSums(r2), elbo = c(-2, -1), iterations = 186L,
            options = options, versions = list(MOFA2 = "stored reference"), messages = character(),
            model_file = list(bytes = as.raw(1:3), sha256 = strrep("d", 64L), size = 3))
    }
}

intRfKeggIds <- function(n) {
    ids <- sprintf("C%05d", 10000L + seq_len(n))
    pairs <- seq(2L, 20L, by = 2L)
    ids[pairs] <- ids[pairs - 1L]
    ids
}

intRfMapping <- function(sets) {
    pathway <- paste0("path:hsa", sprintf("%05d", seq_along(sets)))
    list(database = "KEGG", release = "KEGG Release fixture", sha256 = strrep("e", 64L), url = "fixture",
        downloaded_at = "2026-11-24T00:00:00+0000",
        table = data.frame(pathway = rep(pathway, lengths(sets)), compound = paste0("cpd:", unlist(sets, use.names = FALSE)),
            name = rep(names(sets), lengths(sets)), stringsAsFactors = FALSE))
}

# Two known directions over an omic's runs, the same for the runs of one biological sample and balanced over the two
# groups: the parity of the replicate number, and the replicate number itself.
intRfDirections <- function(samples) {
    replicate <- as.integer(sub("^.*_", "", as.character(samples$tech_rep_group)))
    w <- cbind(ifelse(replicate %% 2L == 0L, 1, -1), replicate - 5.5)
    rownames(w) <- as.character(samples$Run)
    w
}

# A corrected matrix made from a normalised one (ticket 6.6): less the two known directions times known per-feature
# coefficients, the level of every second feature shifted by 0.5. `constant = TRUE`: less 0.5 only (a difference that is
# a constant per feature, so no direction can be determined).
intRfCorrect <- function(table, w, constant = FALSE) {
    runs <- names(table)[-1L]
    m <- as.matrix(table[runs])
    i <- seq_len(nrow(m))
    if (constant) {
        m <- m - 0.5
    } else {
        a <- cbind(0.30 + 0.01 * (i %% 13L), 0.02 + 0.05 * ((i %% 5L) - 2))
        m <- m - a %*% t(w[runs, , drop = FALSE]) - ifelse(i %% 2L == 0L, 0.5, 0)
    }
    out <- data.frame(feature_id = table$feature_id, m, check.names = FALSE, stringsAsFactors = FALSE)
    names(out) <- names(table)
    out
}

# A project holding the corpus: metabolomics (both data sets given a `kegg_id` column) and proteomics rewritten with a
# corrected matrix and a correction of k = 2 (lipidomics and transcriptomics keep no correction).
#   metabolomics = "exact": both of its data sets corrected by the two known directions;
#   metabolomics = "mixed": metabolomics_GCMS as above with its removed factors stripped afterwards (an input saved by an
#     earlier version), metabolomics_LCMS_Pos corrected by a constant per feature (not determined).
intRfProject <- function(metabolomics = "exact", .local_envir = parent.frame()) {
    testthat::skip_if_not(dir.exists(file.path(intRfData(), "corpus")), "the tier-1 corpus is not present")
    testthat::skip_if_not(workflowOwnerCapable(), "the integration store needs qs2, DuckDB and Arrow")
    root <- withr::local_tempdir(.local_envir = .local_envir)
    results_dir <- file.path(root, "results", "integration")
    paths <- list(base_dir = root, integration_dir = file.path(root, "integration"),
        source_dir = file.path(root, "scripts", "integration"), results_dir = results_dir,
        mofa_tables_dir = file.path(results_dir, "mofa", "tables"), mofa_plots_dir = file.path(results_dir, "mofa", "plots"),
        study_report_dir = file.path(root, "results_summary", "integration", "Study_report"))
    dir.create(paths$source_dir, recursive = TRUE)
    for (omic in c("lipidomics", "metabolomics", "proteomics", "transcriptomics")) {
        dir <- file.path(paths$integration_dir, omic)
        dir.create(dir, recursive = TRUE)
        stopifnot(all(file.copy(list.files(file.path(intRfData(), "corpus", omic), full.names = TRUE), dir)))
    }
    for (omic in c("metabolomics", "proteomics")) {
        dir <- file.path(paths$integration_dir, omic)
        input <- readIntegrationInputWithMatrices(dir)
        read <- input$manifest
        manifest <- list(omic_type = read$omic_type, workflow = read$workflow, matrix_role = read$matrix_role,
            scale = read$scale, correction = list(state = "ruv_corrected", method = "kit184-test", k = 2L,
                replicate_structure = "group", grouping_variable = "group"),
            source_state = list(name = "kit184"), sample_id_column = read$sample_id_column, group_column = read$group_column,
            technical_replicate_column = read$technical_replicate_column,
            biological_sample_column = read$biological_sample_column)
        features <- lapply(input$views, `[[`, "features")
        matrices <- lapply(input$views, `[[`, "matrix")
        if (omic == "metabolomics") for (view in names(features)) features[[view]]$kegg_id <- intRfKeggIds(nrow(features[[view]]))
        w <- intRfDirections(input$samples)
        corrected <- lapply(stats::setNames(names(matrices), names(matrices)), function(view) intRfCorrect(matrices[[view]], w,
            constant = omic == "metabolomics" && metabolomics == "mixed" && view == "metabolomics_LCMS_Pos"))
        suppressWarnings(writeIntegrationInput(dir, manifest, input$samples, matrices, features, corrected = corrected),
            classes = intRfUndetermined)
        if (omic == "metabolomics" && metabolomics == "mixed") intRfStrip(dir, "metabolomics_GCMS")
    }
    paths
}

intRfViewsCorrected <- c("proteomics", "metabolomics_LCMS_Pos", "metabolomics_GCMS")

# A store on the project with the stand-in model: trained on the corrected values of the three corrected data sets
# (`corrected = TRUE`) or, by default, on the values before RUV of every data set.
intRfModel <- function(paths, corrected = FALSE, .local_envir = parent.frame()) {
    manager <- newIntegrationStateOwner(file.path(paths$source_dir, paste0("store-", basename(tempfile()))))
    withr::defer(closeWorkflowStateOwner(manager), envir = .local_envir)
    scan <- integrationScanProject(paths$integration_dir)
    integrationCommitImport(manager, integrationImportRequest(scan, scan$views$view))
    parameters <- if (corrected) list(input_matrix = "ruv_corrected", views = intRfViewsCorrected) else list()
    run <- runIntegrationModel(manager, parameters, trainer = intRfTrainer())
    list(manager = manager, run = run)
}

# A data set of the store averaged per biological sample over S_d, read whole from a slot of the inputs state.
intRfIndependent <- function(manager, value, view, slot = "views") {
    object <- integrationImportedObject(manager)
    table <- methods::slot(object, slot)[[view]]
    m <- as.matrix(table[-1L])
    rownames(m) <- as.character(table$feature_id)
    design <- object@design_matrix
    design <- design[design$omic == integrationViewOmics(object)[[view]], , drop = FALSE]
    key <- design$biological_sample[match(colnames(m), design$run)]
    samples <- as.character(value$factors$biological_sample)
    samples <- samples[samples %in% key]
    x <- vapply(samples, function(s) {
        out <- rowMeans(m[, key == s, drop = FALSE], na.rm = TRUE)
        out[is.nan(out)] <- NA_real_
        out
    }, numeric(nrow(m)))
    rownames(x) <- rownames(m)
    list(x = x, samples = samples, features = object@features[[view]], key = key, runs = colnames(m), raw = m)
}

# Three sets of ten compounds each from complete features on single-feature identifiers.
intRfSets <- function(manager, value, view = "metabolomics_GCMS") {
    independent <- intRfIndependent(manager, value, view)
    ids <- stats::setNames(independent$features$kegg_id, independent$features$feature_id)
    complete <- rownames(independent$x)[stats::complete.cases(independent$x)]
    pool <- intersect(names(ids)[!duplicated(ids) & !duplicated(ids, fromLast = TRUE)], complete)
    list(first = unname(ids[pool[1:10]]), second = unname(ids[pool[11:20]]), third = unname(ids[pool[21:30]]))
}

intRfRequest <- function(view = "metabolomics_GCMS", ...) c(list(factor = "Factor1", view = view, database = "kegg",
    id_column = "kegg_id"), list(...))
intRfEntries <- function(manager) manager$.context$snapshot$exports$integration_mofa_model$enrichment

# A small data set in memory (no store): 40 metabolites in 12 runs; `biological` names each run's biological sample
# (default: one run each), `samples` the model's order of them. m01 and m02 share an identifier; m05 lacks one run.
# `factors`: the removed factors by run (12 x k), kept in the inputs state in another run order; `record`: their record.
intRfUnit <- function(biological = sprintf("s%02d", 1:12), samples = unique(biological), factors = NULL,
    record = list(method = "difference_svd", k_declared = NCOL(factors), rank = NCOL(factors), k = NCOL(factors), share = 1,
        features = 30L, tolerance = 1e-8, status = "exact", reason = NULL), corrected = TRUE, match = NULL) {
    runs <- sprintf("r%02d", 1:12)
    n <- length(samples)
    f <- stats::setNames(as.numeric(scale(seq_len(n) + 1.7 * sin(2.1 * seq_len(n)))), samples)
    ids <- sprintf("m%02d", 1:40)
    kegg <- c("C00001", "C00001", sprintf("C%05d", 3:40))
    noise <- withr::with_seed(20261184L, matrix(stats::rnorm(40 * 12, sd = 0.3), 40L, 12L))
    mu <- withr::with_seed(8L, stats::runif(40, 14, 20))
    key <- biological
    if (!is.null(match)) key[match(match$run, runs)] <- match$biological_sample
    x <- mu + noise + outer(c(rep(0.6, 10L), rep(-0.5, 10L), rep(0, 20L)), unname(f[key]))
    x[5L, 2L] <- NA_real_
    dimnames(x) <- list(ids, runs)
    table <- data.frame(feature_id = ids, x, check.names = FALSE, stringsAsFactors = FALSE)
    features <- data.frame(feature_id = ids, display_name = ids, kegg = kegg, stringsAsFactors = FALSE)
    entry <- list(name = "metab", matrix_file = "matrix_metab.parquet", features_file = "features_metab.parquet", rows = 40L,
        columns = 12L, matrix_stage = "normalised", row_groups = 1L)
    if (corrected) entry$corrected_matrix_file <- "corrected_metab.parquet"
    input <- list(omic = "metabolomics", path = "unused", manifest_sha256 = strrep("b", 64L), views = "metab",
        manifest = list(omic_type = "metabolomics", correction = if (corrected) list(k = 2L), views = list(entry)))
    if (!is.null(record) && corrected) {
        input$manifest$views[[1L]]$removed_factors <- record
        if (!is.null(factors)) {
            input$manifest$views[[1L]]$removed_factors_file <- "removed_metab.parquet"
            dimnames(factors) <- list(runs, paste0("ruv_W", seq_len(ncol(factors))))
            input$removed_factors <- list(metab = c(record, list(sha256 = strrep("a", 64L),
                factors = factors[rev(runs), , drop = FALSE])))
        }
    }
    object <- methods::new("IntegrationViewsData", views = list(metab = table[0L, , drop = FALSE]),
        features = list(metab = features),
        design_matrix = data.frame(omic = "metabolomics", run = runs, group = rep(c("A", "B"), each = 6L),
            biological_sample = biological, stringsAsFactors = FALSE),
        inputs = list(metabolomics = input))
    value <- list(factors = data.frame(biological_sample = samples, Factor1 = unname(f), stringsAsFactors = FALSE),
        weights = list(metab = "unused"))
    averaged <- integrationAverageRuns(x, key)[, samples, drop = FALSE]
    correlation <- vapply(seq_len(40L), function(i) {
        ok <- !is.na(averaged[i, ])
        stats::cor(averaged[i, ok], f[ok])
    }, numeric(1))
    ranking <- list(factor = "Factor1", view = "metab", min_samples = 5L,
        direction = list(rule = "largest_absolute_weight_positive", sign = 1L, view = "metab", feature_id = "m01", weight = 1),
        counts = list(features = 40L, ranked = 40L, too_few_samples = 0L, no_variation = 0L, internal_standard = 0L),
        ranking = data.frame(feature_id = ids, n_samples = as.integer(rowSums(!is.na(averaged))), correlation = correlation,
            status = "ranked", stringsAsFactors = FALSE))
    mapping <- intRfMapping(list(up = sprintf("C%05d", 3:10), down = sprintf("C%05d", 11:20), null = sprintf("C%05d", 21:30),
        mixed = sprintf("C%05d", c(7:10, 31:36))))
    list(object = object, value = value, ranking = ranking, mapping = mapping, x = x, f = f, key = key, samples = samples,
        runs = runs, pages = function(view, columns, kernel) kernel(table[columns], 0L),
        model_parameters = list(input_matrix = "ruv_corrected", sample_match_table = match))
}

# Removed factors over the unit's 12 runs: two directions between runs (k = 2), or `k` random ones.
intRfUnitFactors <- function(k = 2L) {
    m <- withr::with_seed(20261842L, matrix(stats::rnorm(12L * k), 12L, k))
    qr.Q(qr(scale(m, scale = FALSE)))
}

# The unit's request resolved, and its product through the correlated declaration's own kernel.
intRfUnitResolve <- function(unit, ..., model_parameters = unit$model_parameters) {
    integrationResolveEnrichmentParameters(c(list(factor = "Factor1", view = "metab", database = "kegg", id_column = "kegg"),
        list(...)), unit$value, unit$object, length(unit$samples), model_parameters = model_parameters)
}
intRfUnitProduct <- function(unit, ..., model_parameters = unit$model_parameters, .local_envir = parent.frame()) {
    testthat::local_mocked_bindings(integrationStateViewPages = function(manager, source, ...) unit$pages,
        .env = .local_envir)
    p <- intRfUnitResolve(unit, ..., model_parameters = model_parameters)
    meta <- unit$mapping[c("database", "release", "sha256", "url", "downloaded_at")]
    integrationEnrichmentCorrelatedOperation()$kernel(unit$object, list(p = p, ranking_value = unit$ranking,
        mapping_table = unit$mapping$table, mapping_meta = meta, model_value = unit$value,
        model_parameters = model_parameters), list(manager = NULL, source = NULL))
}
intRfUnitValues <- function(unit, covariates = NULL) {
    p <- list(database = "kegg", omic = "metabolomics", id_type = "kegg", id_column = "kegg")
    mapped <- integrationEnrichmentMappedList(unit$ranking, p, unit$object@features$metab)
    integrationEnrichmentTestValues(unit$object, NULL, NULL, unit$value, unit$model_parameters, unit$ranking, mapped,
        "metab", pages = unit$pages, covariates = covariates)
}

# The largest residual of the columns of `a` projected onto the span of the columns of `b` (both centred first).
intRfResidual <- function(a, b) {
    a <- scale(as.matrix(a), scale = FALSE)
    b <- scale(as.matrix(b), scale = FALSE)
    max(abs(stats::lm.fit(b, a)$residuals))
}

# --- K1, K2: the removed factors of one data set ---

test_that("KIT-184 K1: the removed factors of a hand case are exact, orthonormal, centred and in the samples' order", {
    hand <- intRfHand()
    out <- integrationInputRemovedFactors(intRfTable(hand$normalised),
        intRfTable(hand$corrected)[c("feature_id", rev(hand$runs))], 2L, hand$runs)
    record <- out$record
    expect_identical(record$method, "difference_svd")
    expect_identical(record$status, "exact")
    expect_identical(record$rank, 2L)
    expect_identical(record$k, 2L)
    expect_identical(record$k_declared, 2L)
    expect_identical(record$features, 24L)
    expect_identical(record$tolerance, 1e-8)
    expect_lt(abs(record$share - 1), 1e-12)
    expect_identical(names(out$table), c("feature_id", hand$runs))
    expect_identical(out$table$feature_id, c("ruv_W1", "ruv_W2"))
    v <- t(as.matrix(out$table[hand$runs]))
    expect_equal(unname(crossprod(v)), diag(2), tolerance = 1e-12)
    expect_lt(max(abs(colSums(v))), 1e-12)
    expect_true(all(vapply(1:2, function(j) v[which.max(abs(v[, j])), j] > 0, logical(1))))
    centred <- scale(hand$w, scale = FALSE)
    expect_lt(max(abs(centred - v %*% crossprod(v, centred))), 1e-10)
    # The correction's declared number of factors: one for the omic, or one per data set.
    expect_identical(integrationInputDeclaredK(list(k = 2L), "a"), 2L)
    expect_identical(integrationInputDeclaredK(list(k = 2L, views = list(a = list(k = 3L))), "a"), 3L)
    expect_true(is.na(integrationInputDeclaredK(list(method = "x"), "a")))
    expect_true(is.na(integrationInputDeclaredK(list(k = 0L), "a")))
    expect_true(is.na(integrationInputDeclaredK(NULL, "a")))
})

test_that("KIT-184 K2: a third direction is approximate with its share; too little to go on is not determined", {
    hand <- intRfHand()
    runs <- hand$runs
    removed <- function(normalised, corrected, k = 2L) integrationInputRemovedFactors(intRfTable(normalised),
        intRfTable(corrected), k, runs)
    # (a) three of the complete features given a third direction.
    third <- c(1, 1, -1, -1, 2, -2, 0.5, -0.5)
    corrected <- hand$corrected
    corrected[10:12, ] <- corrected[10:12, ] - outer(c(1, 2, 3), third)
    out <- removed(hand$normalised, corrected)
    expect_identical(out$record$status, "approximate")
    expect_identical(out$record$k, 2L)
    expect_identical(out$record$rank, 3L)
    expect_identical(nrow(out$table), 2L)
    complete <- stats::complete.cases(hand$normalised)
    difference <- (hand$normalised - corrected)[complete, ]
    d <- svd(difference - rowMeans(difference))$d
    expect_equal(out$record$share, sum(d[1:2]^2) / sum(d^2), tolerance = 1e-12)
    expect_lt(out$record$share, 1)
    # (b) too few complete features, no declared number of factors, equal matrices: not determined, no table.
    few <- removed(hand$normalised[1:8, ], hand$corrected[1:8, ])
    none <- removed(hand$normalised, hand$corrected, NA_integer_)
    equal <- removed(hand$normalised, hand$normalised)
    # (c) a difference that is a constant per feature.
    constant <- removed(hand$normalised, hand$normalised - rep(c(0.5, 1), 15L))
    for (case in list(few, none, equal, constant)) {
        expect_identical(case$record$status, "not_determined")
        expect_null(case$table)
        expect_identical(case$record$k, 0L)
        expect_true(is.character(case$record$reason) && nzchar(case$record$reason))
    }
    expect_identical(few$record$features, 2L)
    expect_match(constant$record$reason, "constant", fixed = TRUE)
    expect_false(identical(few$record$reason, none$record$reason))
})

# --- K3: the contract ---

test_that("KIT-184 K3: the removed factors are written, listed and read back; a malformed file is refused", {
    hand <- intRfHand()
    runs <- hand$runs
    correction <- list(state = "ruv_corrected", method = "kit184-test", k = 2L)
    dir <- intRfWrite(withr::local_tempdir(), hand$normalised, hand$corrected, correction)
    expect_true(file.exists(file.path(dir, "removed_proteomics.parquet")))
    input <- readIntegrationInput(dir)
    entry <- input$manifest$views[[1L]]
    expect_identical(entry$removed_factors_file, "removed_proteomics.parquet")
    expect_identical(entry$removed_factors[c("method", "k_declared", "rank", "k", "features", "status")],
        list(method = "difference_svd", k_declared = 2L, rank = 2L, k = 2L, features = 24L, status = "exact"))
    sha <- digest::digest(file = file.path(dir, "removed_proteomics.parquet"), algo = "sha256")
    expect_identical(input$manifest$files[["removed_proteomics.parquet"]], sha)
    expect_identical(input$views$proteomics$removed$sha256, sha)
    expect_identical(input$views$proteomics$removed$k, 2L)
    expect_identical(input$views$proteomics$removed$status, "exact")
    factors <- readIntegrationInputRemovedFactors(dir, "proteomics", input)
    direct <- integrationInputRemovedFactors(intRfTable(hand$normalised), intRfTable(hand$corrected), 2L, runs)
    expect_identical(dimnames(factors), list(runs, c("ruv_W1", "ruv_W2")))
    expect_identical(unname(factors), unname(t(as.matrix(direct$table[runs]))))
    table <- direct$table
    refused <- function(edit = identity, table = NULL) {
        copy <- file.path(withr::local_tempdir(), "proteomics")
        dir.create(copy)
        stopifnot(all(file.copy(list.files(dir, full.names = TRUE), copy)))
        intRfTamper(copy, edit, if (!is.null(table)) "removed_proteomics.parquet", table)
        expect_error(readIntegrationInput(copy), class = "integration_input_invalid")
    }
    # A wrong row count, a missing value, another column order, a file without a corrected matrix.
    refused(table = table[1L, , drop = FALSE])
    missing <- table
    missing[1L, 3L] <- NA_real_
    refused(table = missing)
    refused(table = table[c("feature_id", rev(runs))])
    refused(function(manifest) {
        manifest$views[[1L]]$corrected_matrix_file <- NULL
        manifest$views[[1L]]$corrected_row_groups <- NULL
        manifest
    })
    # Not determined: a classed warning naming the data set, and the input written without the file.
    expect_warning(undetermined <- intRfWrite(withr::local_tempdir(), hand$normalised, hand$normalised - 0.5, correction),
        class = intRfUndetermined, regexp = "proteomics")
    expect_false(file.exists(file.path(undetermined, "removed_proteomics.parquet")))
    read <- readIntegrationInput(undetermined)
    expect_null(read$views$proteomics$removed)
    expect_null(read$manifest$views[[1L]]$removed_factors_file)
    expect_identical(read$manifest$views[[1L]]$removed_factors$status, "not_determined")
    # An input without a corrected matrix: the files and the data set's manifest fields the base writes.
    plain <- intRfWrite(withr::local_tempdir(), hand$normalised)
    expect_setequal(list.files(plain), c("integration_input.json", "matrix_proteomics.parquet",
        "features_proteomics.parquet", "samples.parquet"))
    plain_input <- readIntegrationInput(plain)
    expect_identical(names(plain_input$manifest$views[[1L]]), c("name", "matrix_file", "features_file", "rows", "columns",
        "matrix_stage", "row_groups"))
    expect_null(plain_input$views$proteomics$removed)
    expect_identical(names(plain_input$views$proteomics), c("features", "rows", "columns", "row_groups", "corrected",
        "corrected_row_groups"))
})

test_that("KIT-184 K3: the import carries the removed factors; an input without them gives the entry the base gives", {
    paths <- intRfProject()
    read <- integrationReadInputs(list.dirs(paths$integration_dir, recursive = FALSE))
    base_names <- c("omic", "path", "manifest_sha256", "views", "manifest", "samples", "biological_sample_column")
    expect_identical(names(read$record$lipidomics), base_names)
    expect_identical(names(read$record$transcriptomics), base_names)
    expect_identical(names(read$record$proteomics), c(base_names, "removed_factors"))
    held <- read$record$metabolomics$removed_factors
    expect_setequal(names(held), c("metabolomics_LCMS_Pos", "metabolomics_GCMS"))
    gcms <- held$metabolomics_GCMS
    expect_identical(gcms$status, "exact")
    expect_identical(gcms$k, 2L)
    dir <- file.path(paths$integration_dir, "metabolomics")
    input <- readIntegrationInput(dir)
    expect_identical(gcms$sha256, input$manifest$files[["removed_metabolomics_GCMS.parquet"]])
    expect_identical(gcms$factors, readIntegrationInputRemovedFactors(dir, "metabolomics_GCMS", input))
    expect_identical(rownames(gcms$factors), input$samples$Run)
    # The written factors span the two known directions the fixture removed.
    expect_lt(intRfResidual(intRfDirections(input$samples), gcms$factors), 1e-8)
    # The state holds them, and an input without them adds nothing to its omic's entry.
    fixture <- intRfModel(paths)
    object <- integrationImportedObject(fixture$manager, fields = "inputs")
    expect_identical(object@inputs$metabolomics$removed_factors$metabolomics_GCMS$factors, gcms$factors)
    expect_identical(names(object@inputs$lipidomics), base_names)
})

# K3's base parts (the base-store job only): the inputs the base commit's code wrote.
intRfBaseStore <- function(part) {
    root <- Sys.getenv("K184_BASE_STORE", unset = "")
    testthat::skip_if_not(nzchar(root) && dir.exists(file.path(root, part)),
        "the KIT-184 base store is built only inside base-store-job.R's runner job")
    file.path(root, part)
}
intRfBaseManager <- function(part, .local_envir = parent.frame()) {
    source_store <- intRfBaseStore(file.path(part, "store"))
    testthat::skip_if_not(workflowOwnerCapable(), "the integration store needs qs2, DuckDB and Arrow")
    root <- withr::local_tempdir(.local_envir = .local_envir)
    stopifnot(file.copy(source_store, root, recursive = TRUE))
    root <- file.path(root, basename(source_store))
    identity <- readSavedStoreIdentity(file.path(root, ".workflow-entities", "store.json"))
    manager <- newWorkflowStateOwner(newWorkflowEntityStore(root, identity$project_id, identity$workflow_id), reopen = TRUE)
    withr::defer(closeWorkflowStateOwner(manager), envir = .local_envir)
    manager
}

test_that("KIT-184 K3 (base store): an input written at base reads as before; one without a correction has base's bytes", {
    corrected_dirs <- list.dirs(intRfBaseStore(file.path("corrected", "integration")), recursive = FALSE)
    read <- integrationReadInputs(corrected_dirs)
    base_names <- c("omic", "path", "manifest_sha256", "views", "manifest", "samples", "biological_sample_column")
    for (omic in names(read$record)) expect_identical(names(read$record[[omic]]), base_names, info = omic)
    expect_true(read$views$proteomics$corrected)
    expect_null(read$inputs$proteomics$views$proteomics$removed)
    # The lipidomics input of the corpus rewritten by the base and by the head: the same files, fields and Parquet bytes.
    base_plain <- intRfBaseStore(file.path("plain", "lipidomics"))
    source_dir <- file.path(intRfData(), "corpus", "lipidomics")
    input <- readIntegrationInputWithMatrices(source_dir)
    manifest_fields <- input$manifest[c("omic_type", "workflow", "matrix_role", "scale", "sample_id_column", "group_column",
        "technical_replicate_column", "biological_sample_column")]
    manifest_fields$source_state <- list(name = "kit184-plain")
    head_plain <- file.path(withr::local_tempdir(), "lipidomics")
    writeIntegrationInput(head_plain, manifest_fields, input$samples, lapply(input$views, `[[`, "matrix"),
        lapply(input$views, `[[`, "features"))
    expect_identical(sort(list.files(head_plain)), sort(list.files(base_plain)))
    base_manifest <- jsonlite::read_json(file.path(base_plain, INTEGRATION_INPUT_MANIFEST), simplifyVector = FALSE)
    head_manifest <- jsonlite::read_json(file.path(head_plain, INTEGRATION_INPUT_MANIFEST), simplifyVector = FALSE)
    expect_identical(names(head_manifest), names(base_manifest))
    expect_identical(lapply(head_manifest$views, names), lapply(base_manifest$views, names))
    expect_identical(head_manifest$files, base_manifest$files)
})

# --- K5: the covariates at test time ---

test_that("KIT-184 K5: the covariates are the removed factors averaged as the values are, matched by name, as a basis", {
    biological <- c("s1", "s1", "s2", "s2", "s3", "s3", "s4", "s4", "s5", "s5", "s6", "s6")
    samples <- c("s3", "s1", "s6", "s2", "s5", "s4")
    between <- c(1, 1, -2, -2, 0.5, 0.5, 3, 3, -1, -1, -1.5, -1.5)
    within <- rep(c(1, -1), 6L)
    # A sample match that moves run r05 from s3 to s1: s1 then has three runs.
    match <- list(omic = "metabolomics", run = "r05", biological_sample = "s1")
    moved <- intRfUnit(biological, samples, factors = cbind(between, within + 0.3 * between), match = match)
    out <- integrationEnrichmentRemovedFactors(moved$object, moved$value, moved$model_parameters, "metab")
    expect_identical(rownames(out$covariates), samples)
    expect_identical(colnames(out$covariates), paste0("ruv_W", seq_len(ncol(out$covariates))))
    held <- moved$object@inputs$metabolomics$removed_factors$metab$factors
    averaged <- t(integrationAverageRuns(t(held[moved$runs, , drop = FALSE]), moved$key))[samples, , drop = FALSE]
    expect_identical(ncol(out$covariates), 2L)
    expect_lt(intRfResidual(averaged, out$covariates), 1e-12)
    expect_lt(intRfResidual(out$covariates, averaged), 1e-12)
    expect_equal(unname(crossprod(out$covariates)), diag(2), tolerance = 1e-12)
    expect_lt(max(abs(colSums(out$covariates))), 1e-12)
    # With the saved biological samples in place of the match, the span would be another one.
    saved <- t(integrationAverageRuns(t(held[moved$runs, , drop = FALSE]), biological))[samples, , drop = FALSE]
    expect_gt(intRfResidual(saved, out$covariates), 1e-6)
    detail <- out$detail
    expect_identical(detail[c("source", "k_declared", "k_input", "k_design", "status", "share", "features", "input_sha256")],
        list(source = "ruv_removed", k_declared = 2L, k_input = 2L, k_design = 2L, status = "exact", share = 1,
            features = 30L, input_sha256 = strrep("a", 64L)))
    expect_match(detail$covariates_sha256, "^[0-9a-f]{64}$")
    # A factor that only separates the two runs of each sample averages to zero and costs no degree of freedom.
    replicated <- intRfUnit(biological, samples, factors = cbind(between, within))
    one <- integrationEnrichmentRemovedFactors(replicated$object, replicated$value, replicated$model_parameters, "metab")
    expect_identical(ncol(one$covariates), 1L)
    expect_identical(one$detail$k_design, 1L)
    expect_identical(one$detail$k_input, 2L)
    expect_identical(colnames(one$covariates), "ruv_W1")
    # One run per biological sample: both directions stay.
    single <- intRfUnit(factors = intRfUnitFactors())
    both <- integrationEnrichmentRemovedFactors(single$object, single$value, single$model_parameters, "metab")
    expect_identical(dim(both$covariates), c(12L, 2L))
    expect_identical(rownames(both$covariates), single$samples)
    # Everything removed lay within biological samples: no covariate, still a run with the removed factors.
    zero <- intRfUnit(biological, samples, factors = cbind(within, c(1, -1, 0, 0, 2, -2, 0, 0, 1, -1, 0, 0)))
    none <- integrationEnrichmentRemovedFactors(zero$object, zero$value, zero$model_parameters, "metab")
    expect_identical(dim(none$covariates), c(6L, 0L))
    expect_identical(none$detail$k_design, 0L)
    resolved <- intRfUnitResolve(zero, adjust_ruv_factors = TRUE)
    expect_true(resolved$adjust_ruv_factors)
    expect_identical(resolved$design_covariates, paste0("ruv_removed:k=0:input=", strrep("a", 64L)))
    product <- intRfUnitProduct(zero, adjust_ruv_factors = TRUE)
    expect_identical(product$design_columns, c("intercept", "factor"))
    expect_identical(product$design_detail$k_design, 0L)
})

# --- K6: exactness ---

test_that("KIT-184 K6: with covariates the kernel is limma's camera on [intercept, factor, removed factors]; the same model on the corrected values", {
    unit <- intRfUnit(factors = intRfUnitFactors())
    product <- intRfUnitProduct(unit, adjust_ruv_factors = TRUE)
    covariates <- integrationEnrichmentRemovedFactors(unit$object, unit$value, unit$model_parameters, "metab")$covariates
    values <- intRfUnitValues(unit, covariates)
    y <- values$y
    f <- values$f
    x <- cbind(intercept = 1, factor = f, covariates)
    expect_identical(product$design_columns, c("intercept", "factor", "ruv_W1", "ruv_W2"))
    expect_identical(product$design_covariates, paste0("ruv_removed:k=2:input=", strrep("a", 64L)))
    expect_identical(product$design_detail$residual_df, 8L)
    expect_identical(product$design_detail$n_eq, 10L)
    expect_equal(unname(values$design), unname(x))
    # (i) limma's camera with the same design and the contrast by name.
    sets <- split(sub("^cpd:", "", unit$mapping$table$compound), unit$mapping$table$pathway)
    index <- lapply(sets, function(set) set[set %in% rownames(y)])
    direct <- limma::camera(y, index, design = x, contrast = which(colnames(x) == "factor"), inter.gene.cor = NA,
        allow.neg.cor = TRUE, trend.var = TRUE, use.ranks = FALSE, sort = FALSE)
    terms <- product$terms[match(names(index), product$terms$ID), , drop = FALSE]
    expect_equal(terms$pvalue, direct$PValue, tolerance = 1e-12)
    expect_identical(terms$Direction, direct$Direction)
    expect_equal(terms$correlation, direct$Correlation, tolerance = 1e-12)
    expect_identical(terms$setSize, as.integer(direct$NGenes))
    fit <- limma::eBayes(limma::lmFit(y, x), trend = TRUE)
    z <- limma::zscoreT(fit$t[, "factor"], df = fit$df.total, approx = TRUE, method = "hill")
    delta <- vapply(index, function(set) mean(z[rownames(y) %in% set]) - mean(z[!rownames(y) %in% set]), numeric(1))
    expect_equal(terms$effect, unname(delta), tolerance = 1e-12)
    # limma's default contrast (the last column) is another test.
    last <- limma::camera(y, index, design = x, inter.gene.cor = NA, allow.neg.cor = TRUE, trend.var = TRUE,
        use.ranks = FALSE, sort = FALSE)
    expect_gt(max(abs(last$PValue - direct$PValue)), 1e-6)
    # (ii) one feature against stats::lm.
    plain <- limma::lmFit(y, x)
    reference <- summary(stats::lm(y[3L, ] ~ f + covariates[, 1L] + covariates[, 2L]))
    expect_equal(unname(plain$coefficients[3L, "factor"]), reference$coefficients[2L, 1L], tolerance = 1e-10)
    expect_equal(plain$sigma[[3L]], reference$sigma, tolerance = 1e-10)
    expect_equal(unname(plain$coefficients[3L, "factor"] / (plain$stdev.unscaled[3L, "factor"] * plain$sigma[[3L]])),
        reference$coefficients[2L, 3L], tolerance = 1e-10)
    expect_equal(plain$df.residual[[3L]], 8)
    # (iii) corrected values y* = y - W a - c: the same model gives the same test (levels kept: every number; another
    # level: the correlation still).
    a <- withr::with_seed(20261843L, matrix(stats::rnorm(nrow(y) * 2L), nrow(y), 2L))
    removed <- a %*% t(covariates)
    test <- function(y_values, design) {
        out <- integrationEnrichmentPathwayTest(NULL, data.frame(term = rep(names(index), lengths(index)),
            gene = unlist(index, use.names = FALSE)), data.frame(term = names(index), name = names(index)), 1L,
            intRfSettings(), test = "camera", values = list(y = y_values, f = f, design = design, counts = values$counts,
                design_covariates = "x"))$terms
        out[match(names(index), out$ID), , drop = FALSE]
    }
    kept <- test(y - removed, x)
    expect_equal(kept$pvalue, terms$pvalue, tolerance = 1e-10)
    expect_identical(kept$Direction, terms$Direction)
    expect_equal(kept$correlation, terms$correlation, tolerance = 1e-10)
    shifted <- test(y - removed - seq_len(nrow(y)) / 10, x)
    expect_equal(shifted$correlation, terms$correlation, tolerance = 1e-10)
    # (iv) option 1, the corrected values on [intercept, factor], is another answer.
    option1 <- limma::camera(y - removed, index, design = cbind(1, f), contrast = 2, inter.gene.cor = NA,
        allow.neg.cor = TRUE, trend.var = TRUE, use.ranks = FALSE, sort = FALSE)
    expect_gt(max(abs(option1$PValue - terms$pvalue)), 1e-6)
    # (v) cameraPR on the z of the same design.
    pr <- integrationEnrichmentPathwayTest(NULL, data.frame(term = rep(names(index), lengths(index)),
        gene = unlist(index, use.names = FALSE)), data.frame(term = names(index), name = names(index)), 1L, intRfSettings(),
        test = "camera_pr", values = values)$terms
    direct_pr <- limma::cameraPR(stats::setNames(z, rownames(y)), index, inter.gene.cor = 0.01, sort = FALSE)
    expect_equal(pr$pvalue[match(names(index), pr$ID)], direct_pr$PValue, tolerance = 1e-12)
    # (vi) camera_floor with the same covariates: limma's floor, never below camera's p-value (KIT-191 8.2 item 5).
    floor <- intRfUnitProduct(unit, test = "camera_floor", adjust_ruv_factors = TRUE)
    floor_terms <- floor$terms[match(names(index), floor$terms$ID), , drop = FALSE]
    direct_floor <- limma::camera(y, index, design = x, contrast = which(colnames(x) == "factor"), inter.gene.cor = NA,
        allow.neg.cor = FALSE, trend.var = TRUE, use.ranks = FALSE, sort = FALSE)
    expect_equal(floor_terms$pvalue, direct_floor$PValue, tolerance = 1e-12)
    expect_identical(floor$design_columns, product$design_columns)
    expect_true(all(floor_terms$pvalue >= terms$pvalue - 1e-15))
    # Covariates in another column order give the same result.
    swapped <- intRfUnitValues(unit, covariates[, 2:1])
    expect_equal(test(y, swapped$design)$pvalue, terms$pvalue, tolerance = 1e-12)
})

# --- K8: refusals; K14: counts; K15: the test of record by n_eq (the small data set) ---

test_that("KIT-184 K8: too few samples for the design, a factor that is a combination of the removed factors, an approximate input", {
    six <- rep(sprintf("s%d", 1:6), each = 2L)
    between <- cbind(c(1, 1, -2, -2, 0.5, 0.5, 3, 3, -1, -1, -1.5, -1.5), c(1, 1, 2, 2, -3, -3, 1, 1, 0, 0, -1, -1))
    # (c) 6 biological samples and 2 covariates: 7 are needed; 7 samples run.
    unit <- intRfUnit(six, factors = between)
    condition <- expect_error(intRfUnitProduct(unit, adjust_ruv_factors = TRUE),
        class = "integration_enrichment_too_few_samples_for_test")
    expect_match(conditionMessage(condition), "With 2 removed unwanted factors in its model", fixed = TRUE)
    expect_match(conditionMessage(condition), "at least 7 biological samples in this data set; it has 6", fixed = TRUE)
    expect_s3_class(condition, "integration_refused")
    seven <- c(six[1:10], "s6", "s7")
    more <- intRfUnit(seven, factors = cbind(between[, 1L], c(between[1:11, 2L], 4)))
    expect_identical(intRfUnitProduct(more, adjust_ruv_factors = TRUE)$design_detail$residual_df, 3L)
    # Without covariates KIT-183's floor stands: 6 samples run without the adjustment, which is the default (VAL-160).
    off <- intRfUnitProduct(unit)
    expect_identical(off$design_columns, c("intercept", "factor"))
    expect_null(off$design_detail)
    expect_identical(intRfUnitProduct(unit, adjust_ruv_factors = FALSE)$terms, off$terms)
    # (d) the factor is a combination of the removed factors.
    plain <- intRfUnit()
    collinear <- intRfUnit(factors = cbind(1 + 2 * unname(plain$f), intRfUnitFactors()[, 1L]))
    condition <- expect_error(intRfUnitProduct(collinear, adjust_ruv_factors = TRUE),
        class = "integration_enrichment_factor_collinear")
    expect_match(conditionMessage(condition), "nothing is left to test", fixed = TRUE)
    expect_match(conditionMessage(condition), "metab", fixed = TRUE)
    # (e) an approximate input is not refused: the product carries the status and the share.
    record <- list(method = "difference_svd", k_declared = 2L, rank = 5L, k = 2L, share = 0.9, features = 30L,
        tolerance = 1e-8, status = "approximate", reason = NULL)
    approximate <- intRfUnitProduct(intRfUnit(factors = intRfUnitFactors(), record = record), adjust_ruv_factors = TRUE)
    expect_identical(approximate$design_detail[c("status", "share")], list(status = "approximate", share = 0.9))
    expect_identical(approximate$design_columns, c("intercept", "factor", "ruv_W1", "ruv_W2"))
    # Refusals of the setting itself.
    unit <- intRfUnit(factors = intRfUnitFactors())
    for (value in list("none", NA, "yes", 1L, c(TRUE, FALSE))) {
        expect_error(intRfUnitResolve(unit, adjust_ruv_factors = value), class = "integration_enrichment_parameters_invalid")
    }
    condition <- expect_error(intRfUnitResolve(unit, test = "gsea", adjust_ruv_factors = TRUE),
        class = "integration_enrichment_parameters_invalid")
    expect_match(conditionMessage(condition), "the pre-ranked test has no model that the removed factors could enter",
        fixed = TRUE)
    condition <- expect_error(intRfUnitResolve(unit, test = "permutation", adjust_ruv_factors = TRUE),
        class = "integration_enrichment_parameters_invalid")
    expect_match(conditionMessage(condition), "the permutation test as offered has no model the removed factors could enter",
        fixed = TRUE)
    # The pre-ranked and permutation tests keep nothing of the setting; FALSE is what they do.
    for (test in c("gsea", "permutation")) {
        resolved <- intRfUnitResolve(unit, test = test, adjust_ruv_factors = FALSE)
        expect_null(resolved$adjust_ruv_factors)
        expect_null(resolved$covariates)
        expect_identical(resolved$design_covariates, "none")
        expect_null(intRfUnitResolve(unit, test = test)$adjust_ruv_factors)
    }
    expect_identical(INTEGRATION_ENRICHMENT_DESIGN_TESTS, c("camera", "camera_floor", "camera_pr"))
})

test_that("KIT-184 K14: covariates change no count; run_complete and the factor's share by hand", {
    # Eight biological samples: four with two runs, four with one. m05 lacks run r02 (of s1, two runs): complete over the
    # samples, not over the runs.
    biological <- c("s1", "s1", "s2", "s2", "s3", "s3", "s4", "s4", "s5", "s6", "s7", "s8")
    unit <- intRfUnit(biological, factors = cbind(c(1, 1, -2, -2, 0.5, 0.5, 3, 3, -1, -1.5, 0.2, 0.9),
        c(-1, -1, 0, 0, 2, 2, 1, 1, -3, 0.5, 1.5, -1)))
    covariates <- integrationEnrichmentRemovedFactors(unit$object, unit$value, unit$model_parameters, "metab")$covariates
    with <- intRfUnitValues(unit, covariates)
    without <- intRfUnitValues(unit)
    expect_identical(with$counts, without$counts)
    expect_identical(with$y, without$y)
    expect_identical(with$feature_id, without$feature_id)
    expect_identical(with$run_complete, with$counts$universe - 1L)
    expect_null(without$run_complete)
    expect_null(without$factor_r2)
    reference <- summary(stats::lm(with$f ~ covariates[, 1L] + covariates[, 2L]))
    expect_equal(with$factor_r2, reference$r.squared, tolerance = 1e-12)
    product <- intRfUnitProduct(unit, adjust_ruv_factors = TRUE)
    detail <- product$design_detail
    expect_identical(names(detail), c("source", "k_declared", "k_input", "k_design", "status", "share", "features",
        "input_sha256", "covariates_sha256", "residual_df", "n_eq", "run_complete", "factor_r2"))
    expect_equal(detail$factor_r2, reference$r.squared, tolerance = 1e-12)
    expect_identical(detail$run_complete, with$run_complete)
    expect_identical(detail$residual_df, 4L)
    expect_identical(detail$n_eq, 6L)
    expect_identical(product$test_counts, without$counts)
    expect_equal(integrationEnrichmentAdjustedShare(detail$factor_r2, 8L, 2L), reference$adj.r.squared, tolerance = 1e-12)
    # Without the adjustment (the default) the product is KIT-183's: no detail, no count of complete runs.
    off <- intRfUnitProduct(unit)
    expect_null(off$design_detail)
    expect_false("run_complete" %in% names(off))
    expect_identical(off$design_covariates, "none")
})

test_that("KIT-184 K15: the test of record is asked with the samples left after one per removed factor", {
    asked <- integer()
    local_mocked_bindings(integrationEnrichmentTestOfRecord = function(n_samples) {
        asked <<- c(asked, n_samples)
        "camera"
    })
    two <- intRfUnit(factors = intRfUnitFactors(2L))
    expect_true(intRfUnitResolve(two, adjust_ruv_factors = TRUE)$adjust_ruv_factors)
    # The default is off on either model (VAL-160): the test of record is then asked with every sample.
    expect_false(intRfUnitResolve(two)$adjust_ruv_factors)
    expect_false(intRfUnitResolve(two, model_parameters = list(input_matrix = "uncorrected"))$adjust_ruv_factors)
    expect_null(intRfUnitResolve(two)$covariates)
    expect_identical(intRfUnitResolve(two)$design_covariates, "none")
    four <- intRfUnit(factors = intRfUnitFactors(4L))
    intRfUnitResolve(four, adjust_ruv_factors = TRUE)
    expect_identical(as.integer(asked), c(10L, 12L, 12L, 12L, 12L, 8L))
    texts <- integrationGuiTexts()
    entry <- list(database = "kegg", input_matrix = "ruv_corrected")
    few <- sprintf(texts[["enrich_test_few_df"]], 4L, 6L, 7L)
    expect_match(few, "With 4 removed factors in the model this test has 6 residual degrees of freedom", fixed = TRUE)
    expect_match(few, "(7)", fixed = TRUE)
    expect_true(few %in% integrationEnrichmentSummary(intRfUnitProduct(four, adjust_ruv_factors = TRUE), texts, entry))
    expect_false(any(grepl("residual degrees of freedom, fewer than", integrationEnrichmentSummary(intRfUnitProduct(two,
        adjust_ruv_factors = TRUE), texts, entry), fixed = TRUE)))
})

# --- K11: the sentences of a result, and the tab ---

test_that("KIT-184 K11: a result's sentences: the form's label, the note for the model's input, counts, share, overlap", {
    texts <- integrationGuiTexts()
    unit <- intRfUnit(factors = intRfUnitFactors())
    product <- intRfUnitProduct(unit, adjust_ruv_factors = TRUE)
    detail <- product$design_detail
    trained <- list(database = "kegg", input_matrix = "ruv_corrected")
    default <- list(database = "kegg", input_matrix = "uncorrected")
    summary <- integrationEnrichmentSummary(product, texts, trained)
    percent <- function(x) integrationEnrichmentPercent(100 * x)
    expect_true(sprintf(texts[["enrich_adjust_label_yes"]], 2L) %in% summary)
    expect_identical(sprintf(texts[["enrich_adjust_label_yes"]], 2L), "Adjusted for removed RUV factors: yes (2 factors).")
    expect_true(sprintf(texts[["enrich_values_ruv_covariates_note"]], 2L, 2L, 8L) %in% summary)
    expect_match(sprintf(texts[["enrich_values_ruv_covariates_note"]], 2L, 2L, 8L),
        "For features measured in every run this gives exactly the result of fitting the same model", fixed = TRUE)
    counts <- sprintf(texts[["enrich_test_covariates_counts"]], 2L, 12L, 8L, detail$run_complete,
        product$test_counts$universe, percent(detail$factor_r2), percent(2 / 11))
    expect_true(counts %in% summary)
    expect_match(counts, "were measured in every run; for the others the correction worked feature by feature", fixed = TRUE)
    expect_match(counts, "(about 18.2 % would be expected for an unrelated factor)", fixed = TRUE)
    expect_false(texts[["enrich_values_ruv_note"]] %in% summary)
    expect_false(texts[["enrich_values_ruv_unadjusted_note"]] %in% summary)
    expect_false(texts[["enrich_factor_removed_overlap"]] %in% summary)
    expect_false(any(grepl("lies within them", summary, fixed = TRUE)))
    # The overlap sentence follows the adjusted share: 0.55 at 12 samples and 2 factors is 0.45 adjusted (none); 0.62 is
    # 0.536 (shown).
    expect_equal(integrationEnrichmentAdjustedShare(0.55, 12L, 2L), 0.45, tolerance = 1e-12)
    below <- product
    below$design_detail$factor_r2 <- 0.55
    expect_false(texts[["enrich_factor_removed_overlap"]] %in% integrationEnrichmentSummary(below, texts, trained))
    above <- product
    above$design_detail$factor_r2 <- 0.62
    expect_true(texts[["enrich_factor_removed_overlap"]] %in% integrationEnrichmentSummary(above, texts, trained))
    # An approximate input: no sentence claims an exact match.
    approximate <- product
    approximate$design_detail$status <- "approximate"
    approximate$design_detail$share <- 0.9
    summary <- integrationEnrichmentSummary(approximate, texts, trained)
    expect_true(sprintf(texts[["enrich_removed_approximate"]], 2L, "90.0 %") %in% summary)
    expect_match(sprintf(texts[["enrich_removed_approximate"]], 2L, "90.0 %"),
        "so the match with the corrected values is approximate for every feature", fixed = TRUE)
    expect_true(any(grepl("The model has 8 residual degrees of freedom.", summary, fixed = TRUE)))
    expect_false(any(grepl("this gives exactly", summary, fixed = TRUE)))
    expect_false(any(grepl("were measured in every run", summary, fixed = TRUE)))
    # The adjustment switched on for a model trained on the values before RUV.
    summary <- integrationEnrichmentSummary(product, texts, default)
    expect_true(sprintf(texts[["enrich_values_adjusted_default_note"]], 2L) %in% summary)
    expect_true(sprintf(texts[["enrich_adjust_label_yes"]], 2L) %in% summary)
    expect_false(any(grepl("This model was trained on values corrected", summary, fixed = TRUE)))
    # Without the adjustment, the default (VAL-160): the label; on a model trained on corrected values the note of that
    # form, which names it as the default, not as a sensitivity analysis, and carries the simulation's numbers.
    off <- intRfUnitProduct(unit)
    summary <- integrationEnrichmentSummary(off, texts, trained)
    expect_true(texts[["enrich_adjust_label_no"]] %in% summary)
    expect_identical(texts[["enrich_adjust_label_no"]], "Adjusted for removed RUV factors: no.")
    expect_true(texts[["enrich_values_ruv_unadjusted_note"]] %in% summary)
    unadjusted_note <- texts[["enrich_values_ruv_unadjusted_note"]]
    expect_false(grepl("sensitivity analysis", unadjusted_note, fixed = TRUE))
    expect_false(grepl("can bias a set", unadjusted_note, fixed = TRUE))
    expect_match(unadjusted_note, "which is the default (the adjustment is optional)", fixed = TRUE)
    expect_match(unadjusted_note, paste("were called on the model's main factor in 2.0 to 11.2 % of cases in this form and in",
        "8.1 to 29.6 % with the adjustment"), fixed = TRUE)
    # The limit of the adjustment, with every result that has it and beside the control (results audit K184-R12 (d)).
    limit <- paste("Because the unwanted factors are estimated, sets of features that share them can still be called when",
        "that variation follows the factor (see Methods).")
    for (key in c("enrich_values_ruv_covariates_note", "enrich_values_ruv_covariates_note_approximate",
        "enrich_values_adjusted_default_note", "enrich_adjust_ruv_help")) {
        expect_match(texts[[key]], limit, fixed = TRUE, info = key)
    }
    help <- texts[["enrich_adjust_ruv_help"]]
    expect_match(help, "Off unless you switch it on.", fixed = TRUE)
    expect_false(grepl("This is the default for a model trained on", help, fixed = TRUE))
    expect_false(grepl("sensitivity analysis", help, fixed = TRUE))
    expect_match(help, "at most 5.4 %", fixed = TRUE)
    summary <- integrationEnrichmentSummary(off, texts, default)
    expect_true(texts[["enrich_adjust_label_no"]] %in% summary)
    expect_false(any(grepl("RUV-corrected values", summary, fixed = TRUE)))
    # camera_floor with the factors: its own sentence beside the note.
    floor <- integrationEnrichmentSummary(intRfUnitProduct(unit, test = "camera_floor", adjust_ruv_factors = TRUE), texts,
        trained)
    expect_true(texts[["enrich_floor_covariates_note"]] %in% floor)
    expect_match(texts[["enrich_floor_covariates_note"]],
        "False calls are at most those of the correlation-adjusted test with the same unwanted factors", fixed = TRUE)
    expect_match(texts[["enrich_floor_covariates_note"]], paste("were found in 13 to 32 % of cases with these factors in the",
        "model (the correlation-adjusted test: 30 to 48 %)"), fixed = TRUE)
    expect_false(grepl("has not been measured", texts[["enrich_floor_covariates_note"]], fixed = TRUE))
    expect_false(texts[["enrich_floor_covariates_note"]] %in% integrationEnrichmentSummary(product, texts, trained))
    # The pre-ranked and permutation tests of a model trained on corrected values: the limitation, no label.
    expect_match(texts[["enrich_values_ruv_note"]], paste("but this result uses the normalised values before RUV correction",
        "without the unwanted factors the correction removed"), fixed = TRUE)
    permutation <- suppressMessages(intRfUnitProduct(unit, test = "permutation"))
    summary <- integrationEnrichmentSummary(permutation, texts, trained)
    expect_true(texts[["enrich_values_ruv_note"]] %in% summary)
    expect_false(any(grepl("Adjusted for removed RUV factors", summary, fixed = TRUE)))
    expect_false(texts[["enrich_values_ruv_note"]] %in% integrationEnrichmentSummary(permutation, texts, default))
    # Beside the selector: the limitation for the tests without a design only.
    expect_true(texts[["enrich_values_ruv_note"]] %in% integrationEnrichmentTestNotes("gsea", 12L, "ruv_corrected",
        "metabolomics", texts))
    expect_true(texts[["enrich_values_ruv_note"]] %in% integrationEnrichmentTestNotes("permutation", 12L, "ruv_corrected",
        "metabolomics", texts))
    for (test in c("camera", "camera_floor")) {
        expect_false(texts[["enrich_values_ruv_note"]] %in% integrationEnrichmentTestNotes(test, 12L, "ruv_corrected",
            "metabolomics", texts))
    }
    # The request: the setting for a test that has a design only.
    request <- function(...) integrationEnrichmentRequest(list(factor = "Factor1", view = "v", database = "kegg",
        id_column = "kegg_id", ...))
    expect_true(request(test = "camera", adjust_ruv_factors = TRUE)$adjust_ruv_factors)
    expect_false(request(test = "camera_floor", adjust_ruv_factors = FALSE)$adjust_ruv_factors)
    expect_null(request(test = "gsea", adjust_ruv_factors = TRUE)$adjust_ruv_factors)
    expect_null(request(test = "permutation", adjust_ruv_factors = TRUE)$adjust_ruv_factors)
    expect_null(integrationEnrichmentRequest(list(factor = "Factor1", view = "v", database = "string", test = "camera",
        adjust_ruv_factors = TRUE))$adjust_ruv_factors)
    expect_identical(texts[["enrich_adjust_ruv"]], "Adjust for removed RUV factors")
})

test_that("KIT-184 K11: the tab: the control is off on either model; switched on, the adjusted form has its own key", {
    skip_if_not_installed("clusterProfiler")
    paths <- intRfProject()
    texts <- integrationGuiTexts()
    escaped <- function(text) as.character(htmltools::htmlEscape(text))
    withr::local_options(MultiScholaR.integration_gui = TRUE)
    no_progress <- function(message, value, expr) force(expr)
    inputs <- list(id_column = "kegg_id", correlation_method = "pearson", min_samples = 5L, seed = 20261121L,
        string_species = 9606L, string_fdr = 0.05, kegg_organism = "hsa", reactome_organism = "Homo sapiens",
        min_set_size = 5L, max_set_size = 500L)
    serve <- function(fixture, check) {
        sets <- intRfSets(fixture$manager, fixture$run$value)
        links <- paste(unlist(lapply(seq_along(sets), function(i) sprintf("path:map%05d\tcpd:%s", i, sets[[i]]))),
            collapse = "\n")
        names <- paste(sprintf("hsa%05d\t%s", seq_along(sets), names(sets)), collapse = "\n")
        replies <- list("https://rest.kegg.jp/info/kegg" = "kegg   Kyoto Encyclopedia\n                 Release fixture, Nov 25\n",
            "https://rest.kegg.jp/link/compound/pathway" = paste0(links, "\n"),
            "https://rest.kegg.jp/list/pathway/hsa" = paste0(names, "\n"))
        download <- function(url, timeout_seconds = 300) {
            if (is.null(replies[[url]])) stop("no network in this test: ", url)
            charToRaw(replies[[url]])
        }
        bus <- do.call(shiny::reactiveValues, omicWorkflowInitialBus(integrationGuiConfig()$shell))
        bus$state_manager <- fixture$manager
        bus$integration_model <- list(reference = fixture$run$reference, status = fixture$run$status)
        shiny::testServer(mod_integration_enrichment_server, args = list(workflow_data = bus, paths = paths,
            download = download, withProgressFn = no_progress), {
            session$flushReact()
            session$setInputs(factor = "Factor1", view = "metabolomics_GCMS", database = "kegg")
            check(session, output, state, bus)
        })
    }
    run <- function(session, ...) {
        do.call(session$setInputs, c(list(view = "metabolomics_GCMS", database = "kegg"), list(...), inputs))
        session$setInputs(download = 1)
        session$setInputs(run = as.integer(Sys.time()) %% 100000L + stats::runif(1))
    }
    # A model trained on corrected values: the control is shown for the tests with a design, and it is off (VAL-160).
    serve(intRfModel(paths, corrected = TRUE), function(session, output, state, bus) {
        session$setInputs(test = "camera")
        control <- paste(as.character(output$adjust_ui$html), collapse = "")
        expect_match(control, sprintf('id="%s"', session$ns("adjust_ruv_factors")), fixed = TRUE)
        expect_match(control, escaped(texts[["enrich_adjust_ruv"]]), fixed = TRUE)
        expect_match(control, escaped(texts[["enrich_adjust_ruv_help"]]), fixed = TRUE)
        expect_false(grepl('checked="checked"', control, fixed = TRUE))
        session$setInputs(test = "camera_floor")
        expect_match(paste(as.character(output$adjust_ui$html), collapse = ""), session$ns("adjust_ruv_factors"), fixed = TRUE)
        for (test in c("gsea", "permutation")) {
            session$setInputs(test = test)
            expect_identical(paste(as.character(output$adjust_ui$html), collapse = ""), "")
        }
        session$setInputs(view = "proteomics", database = "string")
        expect_identical(paste(as.character(output$adjust_ui$html), collapse = ""), "")
        run(session, test = "camera", adjust_ruv_factors = TRUE)
        expect_identical(bus$integration_enrichment$entry, "Factor1|metabolomics_GCMS|kegg|camera|adjusted")
        result <- as.character(output$result$html)
        expect_match(result, escaped(sprintf(texts[["enrich_adjust_label_yes"]], 2L)), fixed = TRUE)
        detail <- state$run$run$value$design_detail
        expect_match(result, escaped(sprintf(texts[["enrich_values_ruv_covariates_note"]], 2L, 2L, detail$residual_df)),
            fixed = TRUE)
        expect_false(grepl(escaped(texts[["enrich_values_ruv_note"]]), result, fixed = TRUE))
        run(session, test = "camera", adjust_ruv_factors = FALSE)
        expect_identical(bus$integration_enrichment$entry, "Factor1|metabolomics_GCMS|kegg|camera")
        result <- as.character(output$result$html)
        expect_match(result, escaped(texts[["enrich_adjust_label_no"]]), fixed = TRUE)
        expect_match(result, escaped(texts[["enrich_values_ruv_unadjusted_note"]]), fixed = TRUE)
        run(session, test = "gsea", adjust_ruv_factors = TRUE)
        expect_identical(bus$integration_enrichment$entry, "Factor1|metabolomics_GCMS|kegg")
        result <- as.character(output$result$html)
        expect_match(result, escaped(texts[["enrich_values_ruv_note"]]), fixed = TRUE)
        expect_false(grepl("Adjusted for removed RUV factors", result, fixed = TRUE))
    })
    # A model trained on the values before RUV: off as well; switched on, the adjusted form under its own key.
    serve(intRfModel(paths), function(session, output, state, bus) {
        session$setInputs(test = "camera")
        control <- paste(as.character(output$adjust_ui$html), collapse = "")
        expect_match(control, sprintf('id="%s"', session$ns("adjust_ruv_factors")), fixed = TRUE)
        expect_false(grepl('checked="checked"', control, fixed = TRUE))
        run(session, test = "camera", adjust_ruv_factors = FALSE)
        expect_identical(bus$integration_enrichment$entry, "Factor1|metabolomics_GCMS|kegg|camera")
        result <- as.character(output$result$html)
        expect_match(result, escaped(texts[["enrich_adjust_label_no"]]), fixed = TRUE)
        expect_false(grepl("RUV-corrected values", result, fixed = TRUE))
        run(session, test = "camera", adjust_ruv_factors = TRUE)
        expect_identical(bus$integration_enrichment$entry, "Factor1|metabolomics_GCMS|kegg|camera|adjusted")
        result <- as.character(output$result$html)
        expect_match(result, escaped(sprintf(texts[["enrich_adjust_label_yes"]], 2L)), fixed = TRUE)
        expect_match(result, escaped(sprintf(texts[["enrich_values_adjusted_default_note"]], 2L)), fixed = TRUE)
    })
})

# --- K7, K8, K9, K16: through the store ---

test_that("KIT-184 K7 (i) (a): the declarations are those of base; the setting is declared in none", {
    correlated <- integrationEnrichmentCorrelatedOperation()
    expect_identical(correlated[c("name", "input_class", "fields", "unit", "product_type", "factors", "parameters")],
        list(name = "integration_enrichment_correlated", input_class = "IntegrationViewsData",
            fields = c("features", "design_matrix", "inputs"), unit = "whole_matrix", product_type = "integration.enrichment",
            factors = list(product_schema = "multischolar.integration_enrichment_correlated.v1"),
            parameters = c("ranking", "database", "id_column", "mapping", "organism", "min_set_size", "max_set_size", "test",
                "test_settings", "values_slot", "design_covariates")))
    expect_identical(integrationEnrichmentOperation()$parameters, c("ranking", "database", "id_column", "service_response",
        "mapping", "organism", "string_species", "string_fdr", "string_rank_direction", "min_set_size", "max_set_size", "seed"))
    expect_identical(integrationEnrichmentRankingOperation()$parameters, c("model", "factor", "view", "correlation_method",
        "min_samples", "direction_rule"))
    expect_identical(unname(INTEGRATION_ENRICHMENT_TEST_SETTINGS[["camera"]]), paste0(
        "inter.gene.cor=NA;allow.neg.cor=TRUE;use.ranks=FALSE;trend.var=TRUE;directional=TRUE;",
        "statistic=zscoreT(moderated t, trend, approx=TRUE, method=hill);complete=all_samples;",
        "collapse=largest_mean_among_complete"))
})

test_that("KIT-184 K7 (ii), K16: the four combinations of the setting; camera_floor; the permutation; the corpus check", {
    skip_if_not_installed("clusterProfiler")
    paths <- intRfProject()
    trained <- intRfModel(paths, corrected = TRUE)
    default <- intRfModel(paths)
    view <- "metabolomics_GCMS"
    texts <- integrationGuiTexts()
    mapping <- intRfMapping(intRfSets(default$manager, default$run$value))
    run <- function(fixture, ..., view = "metabolomics_GCMS") suppressMessages(runIntegrationEnrichment(fixture$manager,
        intRfRequest(view, ...), mapping = mapping))
    samples <- intRfIndependent(trained$manager, trained$run$value, view)$samples
    n <- length(samples)
    expect_identical(n, 19L)
    labelled <- function(result, label) label %in% integrationEnrichmentSummary(result$run$value, texts, result$entry)
    yes <- sprintf(texts[["enrich_adjust_label_yes"]], 2L)
    no <- texts[["enrich_adjust_label_no"]]
    # A model trained on corrected values, the setting absent: KIT-183's test of the values before RUV, no covariates
    # (the default is off, VAL-160), under the key without a suffix.
    unadjusted <- run(trained)
    expect_identical(unadjusted$run$value$design_columns, c("intercept", "factor"))
    expect_identical(unadjusted$run$value$design_covariates, "none")
    expect_null(unadjusted$run$value$design_detail)
    expect_false(unadjusted$entry$settings$adjust_ruv_factors)
    expect_true(labelled(unadjusted, no))
    expect_true(texts[["enrich_values_ruv_unadjusted_note"]] %in% integrationEnrichmentSummary(unadjusted$run$value, texts,
        unadjusted$entry))
    expect_identical(names(intRfEntries(trained$manager)), "Factor1|metabolomics_GCMS|kegg|camera")
    expect_identical(run(trained, adjust_ruv_factors = FALSE)$run$status, "hit")
    # The same with the setting on: the removed factors in the design, under its own key.
    adjusted <- run(trained, adjust_ruv_factors = TRUE)
    value <- adjusted$run$value
    expect_identical(value$design_columns, c("intercept", "factor", "ruv_W1", "ruv_W2"))
    sha <- integrationImportedObject(trained$manager, fields = "inputs")@inputs$metabolomics$removed_factors[[view]]$sha256
    expect_identical(value$design_covariates, paste0("ruv_removed:k=2:input=", sha))
    expect_identical(value$design_detail$residual_df, n - 4L)
    expect_identical(value$design_detail$residual_df, 15L)
    expect_identical(value$design_detail[c("source", "k_declared", "k_input", "k_design", "status")],
        list(source = "ruv_removed", k_declared = 2L, k_input = 2L, k_design = 2L, status = "exact"))
    expect_identical(value$design_detail$input_sha256, sha)
    expect_true(adjusted$entry$settings$adjust_ruv_factors)
    expect_true(labelled(adjusted, yes))
    expect_setequal(names(intRfEntries(trained$manager)), c("Factor1|metabolomics_GCMS|kegg|camera",
        "Factor1|metabolomics_GCMS|kegg|camera|adjusted"))
    expect_identical(run(trained, adjust_ruv_factors = TRUE)$run$status, "hit")
    # The corpus check (6.6): the values tested are the values before RUV averaged per biological sample, never the
    # corrected ones; with the factors, limma's camera on them with the design's columns.
    export <- trained$manager$.context$snapshot$exports$integration_mofa_model
    object <- integrationImportedObject(trained$manager, fields = c("features", "design_matrix", "inputs"))
    before <- intRfIndependent(trained$manager, trained$run$value, view)
    after <- intRfIndependent(trained$manager, trained$run$value, view, slot = "corrected_views")
    mapped <- integrationEnrichmentMappedList(adjusted$ranking$value, list(database = "kegg", omic = "metabolomics",
        id_type = "kegg", id_column = "kegg_id"), object@features[[view]])
    covariates <- integrationEnrichmentRemovedFactors(object, trained$run$value, export$parameters, view)$covariates
    values <- integrationEnrichmentTestValues(object, trained$manager, export$source, trained$run$value, export$parameters,
        adjusted$ranking$value, mapped, view, covariates = covariates)
    y <- before$x[values$feature_id, samples, drop = FALSE]
    expect_equal(unname(values$y), unname(y), tolerance = 1e-12)
    rownames(y) <- rownames(values$y)
    sets <- split(sub("^cpd:", "", mapping$table$compound), mapping$table$pathway)
    index <- lapply(sets, function(set) set[set %in% rownames(y)])
    f <- values$f
    direct <- function(y_values, design) limma::camera(y_values, index, design = design,
        contrast = which(colnames(design) == "factor"), inter.gene.cor = NA, allow.neg.cor = TRUE, trend.var = TRUE,
        use.ranks = FALSE, sort = FALSE)
    design <- cbind(intercept = 1, factor = f, covariates[samples, , drop = FALSE])
    at <- function(terms) terms[match(names(index), terms$ID), , drop = FALSE]
    expect_equal(at(value$terms)$pvalue, direct(y, design)$PValue, tolerance = 1e-10)
    expect_equal(at(unadjusted$run$value$terms)$pvalue, direct(y, design[, 1:2])$PValue, tolerance = 1e-10)
    corrected <- after$x[values$feature_id, samples, drop = FALSE]
    rownames(corrected) <- rownames(y)
    on_corrected <- direct(corrected, design)
    # The same model on the corrected values: the same correlations (the removed part lies in the design's span); the
    # shifted levels move the trend covariate, so the p-values differ.
    expect_equal(at(value$terms)$correlation, on_corrected$Correlation, tolerance = 1e-8)
    expect_gt(max(abs(at(value$terms)$pvalue - on_corrected$PValue)), 1e-12)
    expect_gt(max(abs(at(unadjusted$run$value$terms)$pvalue - direct(corrected, design[, 1:2])$PValue)), 1e-8)
    # The covariates span the two known directions the fixture removed, averaged per biological sample.
    input <- readIntegrationInput(file.path(paths$integration_dir, "metabolomics"))
    known <- intRfDirections(input$samples)
    averaged <- t(integrationAverageRuns(t(known[before$runs, , drop = FALSE]), before$key))[samples, , drop = FALSE]
    expect_lt(intRfResidual(averaged, covariates[samples, , drop = FALSE]), 1e-8)
    # A model trained on the values before RUV, the setting absent: KIT-183's product as well (K7 (ii)).
    plain <- run(default)
    expect_identical(plain$run$value$design_columns, c("intercept", "factor"))
    expect_identical(plain$run$value$design_covariates, "none")
    expect_null(plain$run$value$design_detail)
    expect_false("run_complete" %in% names(plain$run$value))
    expect_identical(plain$run$value$test_counts$samples, n)
    expect_false(plain$entry$settings$adjust_ruv_factors)
    expect_true(labelled(plain, no))
    expect_false(any(grepl("RUV-corrected values", integrationEnrichmentSummary(plain$run$value, texts, plain$entry),
        fixed = TRUE)))
    expect_identical(run(default, adjust_ruv_factors = FALSE)$run$status, "hit")
    # The same model with the setting on: the removed factors, under its own key; the two forms side by side.
    on <- run(default, adjust_ruv_factors = TRUE)
    expect_identical(on$run$value$design_columns, c("intercept", "factor", "ruv_W1", "ruv_W2"))
    expect_true(labelled(on, yes))
    expect_true(sprintf(texts[["enrich_values_adjusted_default_note"]], 2L) %in% integrationEnrichmentSummary(on$run$value,
        texts, on$entry))
    expect_setequal(names(intRfEntries(default$manager)), c("Factor1|metabolomics_GCMS|kegg|camera",
        "Factor1|metabolomics_GCMS|kegg|camera|adjusted"))
    # K8 (f): switched on for a data set to which the omic applied no RUV.
    lipid <- names(integrationViewOmics(integrationImportedObject(default$manager, fields = "inputs")))[
        integrationViewOmics(integrationImportedObject(default$manager, fields = "inputs")) == "lipidomics"][[1L]]
    condition <- expect_error(runIntegrationEnrichment(default$manager, list(factor = "Factor1", view = lipid,
        database = "kegg", id_column = "display_name", adjust_ruv_factors = TRUE), mapping = mapping),
        class = "integration_enrichment_removed_factors_missing")
    expect_match(conditionMessage(condition), sprintf("No correction for unwanted variation (RUV) was applied to %s", lipid),
        fixed = TRUE)
    expect_match(conditionMessage(condition), "Switch 'Adjust for removed RUV factors' off.", fixed = TRUE)
    # A value that is neither on nor off; on with the tests that have no design.
    for (bad in list("none", NA, "yes")) expect_error(run(trained, adjust_ruv_factors = bad),
        class = "integration_enrichment_parameters_invalid")
    expect_error(run(trained, test = "gsea", adjust_ruv_factors = TRUE), class = "integration_enrichment_parameters_invalid")
    expect_error(run(trained, test = "permutation", adjust_ruv_factors = TRUE),
        class = "integration_enrichment_parameters_invalid")
    # camera_floor takes the setting as camera does (KIT-191 8.2 item 2).
    floor <- run(trained, test = "camera_floor", adjust_ruv_factors = TRUE)
    expect_identical(floor$run$value$design_columns, c("intercept", "factor", "ruv_W1", "ruv_W2"))
    expect_true(all(at(floor$run$value$terms)$pvalue >= at(value$terms)$pvalue - 1e-15))
    floor_off <- run(trained, test = "camera_floor")
    expect_identical(floor_off$run$value$design_columns, c("intercept", "factor"))
    floor_default <- run(default, test = "camera_floor")
    floor_on <- run(default, test = "camera_floor", adjust_ruv_factors = TRUE)
    expect_identical(floor_on$run$value$design_columns, c("intercept", "factor", "ruv_W1", "ruv_W2"))
    expect_identical(floor_default$run$value$design_covariates, "none")
    # The permutation has one form on either model: the values before RUV, no suffix, nothing kept, the limitation.
    permuted <- run(trained, test = "permutation")
    expect_null(permuted$entry$settings$adjust_ruv_factors)
    expect_identical(permuted$run$value$design_covariates, "none")
    summary <- integrationEnrichmentSummary(permuted$run$value, texts, permuted$entry)
    expect_true(texts[["enrich_values_ruv_note"]] %in% summary)
    expect_false(any(grepl("Adjusted for removed RUV factors", summary, fixed = TRUE)))
    gsea <- run(trained, test = "gsea")
    expect_null(gsea$entry$settings$adjust_ruv_factors)
    expect_true(texts[["enrich_values_ruv_note"]] %in% integrationEnrichmentSummary(gsea$run$value, texts, gsea$entry))
    expect_setequal(names(intRfEntries(trained$manager)), c("Factor1|metabolomics_GCMS|kegg|camera",
        "Factor1|metabolomics_GCMS|kegg|camera|adjusted", "Factor1|metabolomics_GCMS|kegg|camera_floor",
        "Factor1|metabolomics_GCMS|kegg|camera_floor|adjusted", "Factor1|metabolomics_GCMS|kegg|permutation",
        "Factor1|metabolomics_GCMS|kegg"))
    expect_true(all(c("Factor1|metabolomics_GCMS|kegg|camera_floor", "Factor1|metabolomics_GCMS|kegg|camera_floor|adjusted")
        %in% names(intRfEntries(default$manager))))
    # The key: the adjusted form has the suffix, whatever the model was trained on; the form without the adjustment,
    # which is the default and what an entry stored before KIT-184 is, has none.
    expect_identical(integrationEnrichmentKey("F", "v", "kegg", "camera", TRUE), "F|v|kegg|camera|adjusted")
    expect_identical(integrationEnrichmentKey("F", "v", "kegg", "camera", FALSE), "F|v|kegg|camera")
    expect_identical(integrationEnrichmentKey("F", "v", "kegg", "camera_floor", TRUE), "F|v|kegg|camera_floor|adjusted")
    expect_identical(integrationEnrichmentKey("F", "v", "kegg", "permutation", TRUE), "F|v|kegg|permutation")
    expect_identical(integrationEnrichmentKey("F", "v", "kegg", "gsea", FALSE), "F|v|kegg")
    expect_identical(integrationEnrichmentKey("F", "v", "string", "camera", TRUE), "F|v|string")
    expect_identical(integrationEnrichmentKey("F", "v", "kegg", "camera"), "F|v|kegg|camera")
    expect_false("input_matrix" %in% names(formals(integrationEnrichmentKey)))
    # K9 (b): the replay repeats each entry in the form it was saved in.
    for (fixture in list(trained, default)) {
        recorded <- integrationRecordedEnrichment(fixture$manager)
        comparison <- suppressMessages(integrationReplayEnrichment(fixture$manager, recorded))
        expect_true(all(comparison$terms_identical))
        expect_setequal(comparison$enrichment, names(recorded))
        repeated <- integrationRecordedEnrichment(fixture$manager)
        expect_setequal(names(repeated), names(recorded))
        for (key in names(recorded)) {
            expect_identical(repeated[[key]]$value$design_columns, recorded[[key]]$value$design_columns, info = key)
            expect_identical(repeated[[key]]$entry$settings$adjust_ruv_factors,
                recorded[[key]]$entry$settings$adjust_ruv_factors, info = key)
        }
    }
    # An entry without the setting (stored before KIT-184) replays without the adjustment, whatever the model.
    recorded <- integrationRecordedEnrichment(trained$manager)["Factor1|metabolomics_GCMS|kegg|camera"]
    expect_false(recorded[[1L]]$entry$settings$adjust_ruv_factors)
    recorded[[1L]]$entry$settings$adjust_ruv_factors <- NULL
    comparison <- integrationReplayEnrichment(trained$manager, recorded)
    expect_true(comparison$terms_identical)
})

test_that("KIT-184 K8 (a), (b), (g): an input without usable factors is refused in the words for its case", {
    skip_if_not_installed("clusterProfiler")
    paths <- intRfProject("mixed")
    trained <- intRfModel(paths, corrected = TRUE)
    default <- intRfModel(paths)
    texts <- integrationGuiTexts()
    mapping <- intRfMapping(intRfSets(default$manager, default$run$value))
    run <- function(fixture, view, ...) suppressMessages(runIntegrationEnrichment(fixture$manager, intRfRequest(view, ...),
        mapping = mapping))
    missing <- "integration_enrichment_removed_factors_missing"
    # (a) saved by an earlier version: "on" is refused on either model, in the same words, which say what to do; no
    # words depend on what the model was trained on, since the adjustment is nobody's default (VAL-160).
    for (fixture in list(trained, default)) {
        for (extra in list(list(adjust_ruv_factors = TRUE), list(test = "camera_floor", adjust_ruv_factors = TRUE))) {
            condition <- expect_error(do.call(run, c(list(fixture, "metabolomics_GCMS"), extra)), class = missing)
            message <- conditionMessage(condition)
            expect_match(message, "'Adjust for removed RUV factors' was switched on", fixed = TRUE)
            expect_match(message, "metabolomics_GCMS", fixed = TRUE)
            expect_match(message, "Save the metabolomics Summary again", fixed = TRUE)
            expect_match(message, "Switched off (the default), the test runs as it is.", fixed = TRUE)
            expect_false(grepl("trained on RUV-corrected values", message, fixed = TRUE))
            expect_s3_class(condition, "integration_enrichment_refused")
        }
        expect_length(intRfEntries(fixture$manager), 0L)
    }
    # The default runs on such an input, labelled as the form without the adjustment; so do the tests without a design.
    off <- run(trained, "metabolomics_GCMS")
    expect_identical(off$run$value$design_columns, c("intercept", "factor"))
    expect_false(off$entry$settings$adjust_ruv_factors)
    expect_true(texts[["enrich_adjust_label_no"]] %in% integrationEnrichmentSummary(off$run$value, texts, off$entry))
    expect_identical(run(trained, "metabolomics_GCMS", adjust_ruv_factors = FALSE)$run$status, "hit")
    expect_identical(run(trained, "metabolomics_GCMS", test = "gsea")$run$value$service$method, "gsea")
    expect_identical(run(trained, "metabolomics_GCMS", test = "permutation")$run$value$test, "permutation")
    # (b) not determined: the reason, and no advice to save again; (g) the same words under either model.
    for (fixture in list(trained, default)) {
        condition <- expect_error(run(fixture, "metabolomics_LCMS_Pos", adjust_ruv_factors = TRUE), class = missing)
        message <- conditionMessage(condition)
        expect_match(message, "'Adjust for removed RUV factors' was switched on", fixed = TRUE)
        expect_match(message, "could not be determined", fixed = TRUE)
        expect_match(message, "constant", fixed = TRUE)
        expect_false(grepl("Save", message, fixed = TRUE))
        expect_match(message, "Switched off (the default), the test runs as it is.", fixed = TRUE)
    }
    # The sentences of the two former forms for a model trained on the values before RUV are gone with the default.
    expect_false(any(grepl("removed_factors_missing_default", names(INTEGRATION_ENRICHMENT_REFUSALS), fixed = TRUE)))
    # The default of a model trained on the values before RUV is untouched too.
    expect_identical(run(default, "metabolomics_GCMS")$run$value$design_covariates, "none")
})

# --- K10: the record; K12: the workbook ---

test_that("KIT-184 K10: the record and its Methods name each entry's form; the inputs sheet names the removed factors", {
    skip_if_not_installed("clusterProfiler")
    paths <- intRfProject()
    trained <- intRfModel(paths, corrected = TRUE)
    view <- "metabolomics_GCMS"
    texts <- integrationRecordTexts()
    mapping <- intRfMapping(intRfSets(trained$manager, trained$run$value))
    run <- function(fixture, ...) suppressMessages(runIntegrationEnrichment(fixture$manager, intRfRequest(view, ...),
        mapping = mapping))
    adjusted <- run(trained, adjust_ruv_factors = TRUE)
    run(trained)
    run(trained, test = "gsea")
    record <- buildWorkflowAnalysisRecord(trained$manager)
    block <- record$analyses$integration_mofa_model$integration
    tables <- integrationRecordEnrichmentTables(block, texts)
    ranking <- tables$summary[[texts[["col_enrich_ranking"]]]]
    detail <- adjusted$run$value$design_detail
    counts <- adjusted$run$value$test_counts
    yes <- "Adjusted for removed RUV factors: yes (2 factors)."
    no <- "Adjusted for removed RUV factors: no."
    expect_setequal(ranking, c(
        paste0(sprintf(texts[["enrich_ranking_value_camera_ruv"]], 2L, counts$samples, detail$residual_df, counts$universe),
            ". ", yes),
        paste0(sprintf(texts[["enrich_ranking_value_camera"]], counts$samples, counts$universe), ". ", no),
        sprintf(texts[["enrich_ranking_value"]], "pearson", 5L)))
    expect_match(sprintf(texts[["enrich_ranking_value_camera_ruv"]], 2L, 19L, 15L, 100L),
        "regression of each feature on the factor and 2 removed unwanted factors over 19 biological samples (15 residual",
        fixed = TRUE)
    expect_identical(ncol(tables$summary), 11L)
    methods <- integrationMethodsEnrichment(block, texts)
    fragment <- sprintf(texts[["methods_test_covariates_ruv"]], 2L, detail$features, detail$input_sha256, detail$residual_df)
    expect_match(fragment, "an intercept and the 2 unwanted factors that the RUV correction removed from this data set",
        fixed = TRUE)
    camera_sentences <- methods[grepl("pathways were tested with a competitive gene-set test", methods, fixed = TRUE)]
    expect_length(camera_sentences, 2L)
    expect_identical(sum(grepl(fragment, camera_sentences, fixed = TRUE)), 1L)
    expect_identical(sum(grepl(texts[["methods_test_covariates_ruv_off"]], camera_sentences, fixed = TRUE)), 1L)
    expect_false(any(grepl("an intercept and no other covariate,", camera_sentences, fixed = TRUE)))
    expect_false(any(grepl("%", gsub("5 %|95 %|[0-9.]+ %", "", camera_sentences), fixed = TRUE)))
    limits <- methods[startsWith(methods, texts[["methods_enrich_limits_tests"]])]
    expect_length(limits, 1L)
    expect_match(limits, texts[["methods_enrich_limits_ruv"]], fixed = TRUE)
    # The clause as the owner decided it (VAL-160; the results audit's K184-R12): the adjustment optional and off unless
    # chosen, what it kept, and its limit with the simulation's numbers.
    clause <- texts[["methods_enrich_limits_ruv"]]
    expect_false(grepl("by default add to their model", clause, fixed = TRUE))
    expect_false(grepl("sensitivity analysis", clause, fixed = TRUE))
    for (phrase in c("a setting that is off unless chosen",
        "in at most 5.4 % (95 % interval of that largest figure 4.5 to 6.4 %) of cases",
        "in at most 4.1 % of cases with the adjustment when that variation was unrelated to the groups",
        "when it followed the groups (correlation 0.6), in 8.1 to 29.6 % with the adjustment",
        "29.6 %, 95 % interval 27.8 to 31.3 %", "11.2 to 13.3 %", "2.0 to 11.2 % without it", "9.8 to 12.6 %",
        "1.6 to 2.8 % had the true unwanted factors been known")) {
        expect_match(clause, phrase, fixed = TRUE, info = phrase)
    }
    expect_false(grepl("sensitivity analysis", texts[["methods_test_covariates_ruv_off"]], fixed = TRUE))
    expect_false(grepl("variation the correction removed stays in the values tested. P-values", limits, fixed = TRUE))
    # The inputs sheet: six columns; the correction cell of a data set in the model names its removed factors.
    sheets <- integrationRecordTables(record, texts)
    inputs <- sheets$integration_inputs
    expect_identical(ncol(inputs), 6L)
    cell <- inputs[[texts[["col_correction"]]]][inputs[[texts[["col_omic"]]]] == "metabolomics"]
    sha <- integrationImportedObject(trained$manager, fields = "inputs")@inputs$metabolomics$removed_factors[[view]]$sha256
    expect_match(cell, sprintf("2 unwanted factors recorded as removed from %s (exact; SHA-256 %s)", view, sha), fixed = TRUE)
    held <- block$inputs[[which(vapply(block$inputs, function(input) identical(input$omic, "metabolomics"), logical(1)))]]
    expect_identical(held$removed_factors[[view]][c("k", "status", "sha256")], list(k = 2L, status = "exact", sha256 = sha))
    # The record's block keeps each entry's form from its product.
    forms <- vapply(block$enrichment, function(entry) entry$design_covariates %||% NA_character_, character(1))
    expect_identical(sum(startsWith(forms, "ruv_removed:k=2:input="), na.rm = TRUE), 1L)
    expect_identical(sum(forms %in% "none"), 1L)
    # The summary of each entry as the tab shows it: the note of the form without the adjustment for that entry only.
    gui <- integrationGuiTexts()
    entries <- intRfEntries(trained$manager)
    recorded <- integrationRecordedEnrichment(trained$manager)
    shown <- lapply(recorded, function(item) integrationEnrichmentSummary(item$value, gui, item$entry))
    expect_true(gui[["enrich_values_ruv_unadjusted_note"]] %in% shown[["Factor1|metabolomics_GCMS|kegg|camera"]])
    expect_false(gui[["enrich_values_ruv_unadjusted_note"]] %in% shown[["Factor1|metabolomics_GCMS|kegg|camera|adjusted"]])
    expect_true(gui[["enrich_values_ruv_note"]] %in% shown[["Factor1|metabolomics_GCMS|kegg"]])
    expect_false(gui[["enrich_values_ruv_note"]] %in% shown[["Factor1|metabolomics_GCMS|kegg|camera|adjusted"]])
    expect_length(entries, 3L)
    # A model trained on the values before RUV: the adjustment switched on is said so; without it, no RUV clause.
    default <- intRfModel(paths)
    plain <- run(default)
    on <- run(default, adjust_ruv_factors = TRUE)
    block <- buildWorkflowAnalysisRecord(default$manager)$analyses$integration_mofa_model$integration
    methods <- integrationMethodsEnrichment(block, texts)
    on_detail <- on$run$value$design_detail
    expect_identical(sum(grepl(paste0(sprintf(texts[["methods_test_covariates_ruv"]], 2L, on_detail$features,
        on_detail$input_sha256, on_detail$residual_df), texts[["methods_test_covariates_model_uncorrected"]]), methods,
        fixed = TRUE)), 1L)
    plain_counts <- plain$run$value$test_counts
    expect_true(sprintf(texts[["methods_enrich_camera"]], "Factor1", view, "KEGG", texts[["methods_test_name_camera"]],
        sprintf(texts[["methods_test_settings_camera"]], plain_counts$samples, plain_counts$universe,
            plain$run$value$mapping_counts$mapped), "KEGG Release fixture", "hsa", 5L, 500L,
        plain$run$value$service$pathways_testable) %in% methods)
    limits <- methods[startsWith(methods, texts[["methods_enrich_limits_tests"]])]
    expect_false(grepl("RUV-corrected values", limits, fixed = TRUE))
    expect_match(limits, "the values the model is trained on by default. P-values were adjusted", fixed = TRUE)
    ranking <- integrationRecordEnrichmentTables(block, texts)$summary[[texts[["col_enrich_ranking"]]]]
    expect_true(paste0(sprintf(texts[["enrich_ranking_value_camera"]], plain_counts$samples, plain_counts$universe), ". ",
        no) %in% ranking)
    cell <- integrationRecordTables(buildWorkflowAnalysisRecord(default$manager), texts)$integration_inputs[[
        texts[["col_correction"]]]]
    expect_false(any(grepl("recorded as removed", cell, fixed = TRUE)))
})

test_that("KIT-184 K12: the workbook: the default is without the adjustment; the override switches it on; STRING takes none", {
    skip_if_not_installed("clusterProfiler")
    paths <- intRfProject()
    trained <- intRfModel(paths, corrected = TRUE)
    texts <- integrationGuiTexts()
    mapping <- intRfMapping(intRfSets(trained$manager, trained$run$value))
    local_mocked_bindings(integrationDownloadKeggMapping = function(organism, download = NULL) mapping,
        integrationDownloadReactomeMapping = function(download = NULL) list(database = "Reactome", release = "Reactome 0",
            sha256 = strrep("f", 64L), url = "fixture", downloaded_at = "x", table = data.frame(chebi_id = "1",
                reactome_id = "R-HSA-1", reactome_term = "x", organism = "Homo sapiens", stringsAsFactors = FALSE)))
    projectOf <- function(fixture) list(manager = fixture$manager,
        store_root = dirname(fixture$manager$.context$store$directory), paths = paths, views = "metabolomics_GCMS")
    results <- list(results = list(views = "metabolomics_GCMS"), factor = "Factor1")
    id_columns <- list(metabolomics_GCMS = "kegg_id")
    project <- projectOf(trained)
    # Nothing named: the form without the adjustment (VAL-160), with its note for a model trained on corrected values.
    first <- integrationEnrichFactor(project, results, id_columns = id_columns)
    expect_false("Factor1|metabolomics_GCMS|kegg|camera|adjusted" %in% names(intRfEntries(trained$manager)))
    entry <- intRfEntries(trained$manager)[["Factor1|metabolomics_GCMS|kegg|camera"]]
    expect_false(entry$settings$adjust_ruv_factors)
    summary <- first$enrichments[["metabolomics_GCMS kegg"]]$summary
    expect_true(texts[["enrich_adjust_label_no"]] %in% summary)
    expect_true(texts[["enrich_values_ruv_unadjusted_note"]] %in% summary)
    off <- integrationEnrichFactor(project, results, id_columns = id_columns, overrides = list(adjust_ruv_factors = FALSE))
    expect_identical(off$enrichments[["metabolomics_GCMS kegg"]]$terms, first$enrichments[["metabolomics_GCMS kegg"]]$terms)
    # Switched on in `overrides`: the adjusted form, under its own key.
    on_trained <- integrationEnrichFactor(project, results, id_columns = id_columns,
        overrides = list(adjust_ruv_factors = TRUE))
    expect_true(sprintf(texts[["enrich_adjust_label_yes"]], 2L) %in%
        on_trained$enrichments[["metabolomics_GCMS kegg"]]$summary)
    adjusted <- intRfEntries(trained$manager)[["Factor1|metabolomics_GCMS|kegg|camera|adjusted"]]
    expect_true(adjusted$settings$adjust_ruv_factors)
    direct <- suppressMessages(runIntegrationEnrichment(trained$manager, intRfRequest(adjust_ruv_factors = TRUE),
        mapping = mapping))
    expect_identical(direct$run$status, "hit")
    expect_identical(on_trained$enrichments[["metabolomics_GCMS kegg"]]$terms, as.data.frame(direct$run$value$terms))
    # A name that is not a setting of the table stops the call before any request.
    before <- names(intRfEntries(trained$manager))
    condition <- expect_error(integrationEnrichFactor(project, results, id_columns = id_columns,
        overrides = list(design_covariates = "none")))
    expect_match(conditionMessage(condition), "these are not settings of the table: design_covariates", fixed = TRUE)
    expect_identical(names(intRfEntries(trained$manager)), before)
    settings <- integrationSettings()$enrichment
    expect_identical(settings$default[settings$setting == "adjust_ruv_factors"],
        "no (optional: yes adds the unwanted factors the RUV correction removed to the test's model)")
    expect_false(integrationEnrichmentParameterSpec()$adjust_ruv_factors$default)
    # A model trained on the values before RUV, switched on: the data sets with a correction are adjusted, the others
    # reported as refused in the words for a data set without RUV.
    default <- intRfModel(paths)
    object <- integrationImportedObject(default$manager, fields = "inputs")
    lipid <- names(integrationViewOmics(object))[integrationViewOmics(object) == "lipidomics"][[1L]]
    on <- integrationEnrichFactor(projectOf(default), list(results = list(views = c("metabolomics_GCMS", lipid)),
        factor = "Factor1"), id_columns = stats::setNames(list("kegg_id", "display_name"), c("metabolomics_GCMS", lipid)),
        overrides = list(adjust_ruv_factors = TRUE))
    expect_true(sprintf(texts[["enrich_adjust_label_yes"]], 2L) %in% on$enrichments[["metabolomics_GCMS kegg"]]$summary)
    expect_true(any(grepl(sprintf("No correction for unwanted variation (RUV) was applied to %s", lipid), on$skipped,
        fixed = TRUE)))
    expect_true("Factor1|metabolomics_GCMS|kegg|camera|adjusted" %in% names(intRfEntries(default$manager)))
    # The requests sent: the setting only when given, only with a test that has a design, never to STRING.
    sent <- list()
    local_mocked_bindings(runIntegrationEnrichment = function(manager, parameters, ...) {
        sent[[length(sent) + 1L]] <<- parameters
        integrationEnrichmentRefuse("integration_enrichment_parameters_invalid", "stand-in")
    })
    both <- list(results = list(views = c("proteomics", "metabolomics_GCMS")), factor = "Factor1")
    integrationEnrichFactor(project, both, api_key = "k", id_columns = id_columns, overrides = list(adjust_ruv_factors = FALSE))
    by_database <- stats::setNames(sent, vapply(sent, `[[`, character(1), "database"))
    expect_false("adjust_ruv_factors" %in% names(by_database$string))
    expect_false(by_database$kegg$adjust_ruv_factors)
    sent <- list()
    integrationEnrichFactor(project, both, api_key = "k", id_columns = id_columns)
    expect_false(any(vapply(sent, function(request) "adjust_ruv_factors" %in% names(request), logical(1))))
    sent <- list()
    integrationEnrichFactor(project, results, id_columns = id_columns, overrides = list(test = "gsea",
        adjust_ruv_factors = TRUE))
    expect_false(any(vapply(sent, function(request) "adjust_ruv_factors" %in% names(request), logical(1))))
})

# --- K7 (i) (b), (c), (iii) and K9 (a): the base store (the base-store job only) ---

intRfIdentities <- function(manager, entry, extra = list()) {
    export <- manager$.context$snapshot$exports$integration_mofa_model
    s <- entry$settings
    declaration <- integrationEnrichmentCorrelatedOperation()
    parameters <- c(list(ranking = entry$ranking$entity_id, database = s$database, id_column = s$id_column,
        mapping = entry$stored$sha256, organism = s$kegg_organism, min_set_size = s$min_set_size,
        max_set_size = s$max_set_size, test = s$test, test_settings = INTEGRATION_ENRICHMENT_TEST_SETTINGS[[s$test]],
        values_slot = "views", design_covariates = "none"), extra)
    head <- withWorkflowStateScope(manager, function() workflowComponentScientificIdentity(manager, export$source,
        declaration$name, parameters[declaration$parameters]))
    base <- withWorkflowStateScope(manager, function() describeSavedWorkflowEntity(manager$.context$store,
        entry$reference$entity_id)$identity)
    list(head = head, base = base)
}

test_that("KIT-184 K7 (i) (b), (c): a camera entry made at base on a default model has the head's identity apart from the build", {
    skip_if_not_installed("clusterProfiler")
    manager <- intRfBaseManager("default")
    recorded <- integrationRecordedEnrichment(manager)[["Factor1|metabolomics_GCMS|kegg|camera"]]
    entry <- recorded$entry
    expect_null(entry$settings$adjust_ruv_factors)
    keys <- c("input", "design", "parameters", "selectors")
    identities <- intRfIdentities(manager, entry)
    expect_identical(identities$head[keys], identities$base[keys])
    expect_false(identical(identities$head$software, identities$base$software))
    s <- entry$settings
    request <- list(factor = s$factor, view = s$view, database = s$database, id_column = s$id_column, test = "camera")
    mapping <- c(lapply(entry$mapping, identity), list(table = recorded$stored))
    first <- runIntegrationEnrichment(manager, request, mapping = mapping)
    second <- runIntegrationEnrichment(manager, request, mapping = mapping)
    expect_identical(second$run$status, "hit")
    expect_identical(second$run$reference$entity_id, first$run$reference$entity_id)
    expect_identical(first$run$value$design_covariates, "none")
    expect_false(first$entry$settings$adjust_ruv_factors)
    expect_equal(as.data.frame(first$run$value$terms), as.data.frame(recorded$value$terms), tolerance = 1e-12)
})

test_that("KIT-184 K7 (iii), K9 (a): an entry made at base on a model trained on corrected values is the default form, without the adjustment", {
    skip_if_not_installed("clusterProfiler")
    manager <- intRfBaseManager("corrected")
    texts <- integrationGuiTexts()
    recorded <- integrationRecordedEnrichment(manager)
    expect_setequal(names(recorded), c("Factor1|metabolomics_GCMS|kegg|camera", "Factor1|metabolomics_GCMS|kegg"))
    item <- recorded[["Factor1|metabolomics_GCMS|kegg|camera"]]
    entry <- item$entry
    expect_identical(entry$input_matrix, "ruv_corrected")
    expect_null(entry$settings$adjust_ruv_factors)
    # (iii) the same computation has the same identity apart from the build.
    keys <- c("input", "design", "parameters", "selectors")
    identities <- intRfIdentities(manager, entry)
    expect_identical(identities$head[keys], identities$base[keys])
    expect_false(identical(identities$head$software, identities$base$software))
    # The stored entry is shown as what it is: the form without the adjustment, with its note.
    summary <- integrationEnrichmentSummary(item$value, texts, entry)
    expect_true(texts[["enrich_adjust_label_no"]] %in% summary)
    expect_true(texts[["enrich_values_ruv_unadjusted_note"]] %in% summary)
    s <- entry$settings
    request <- list(factor = s$factor, view = s$view, database = s$database, id_column = s$id_column, test = "camera")
    mapping <- c(lapply(entry$mapping, identity), list(table = item$stored))
    # Switched on, the adjustment needs the factors, which an input saved at base does not hold.
    expect_error(runIntegrationEnrichment(manager, c(request, list(adjust_ruv_factors = TRUE)), mapping = mapping),
        class = "integration_enrichment_removed_factors_missing")
    # The default of this model is what base computed: no covariates, under the key the base entry has (VAL-160).
    first <- runIntegrationEnrichment(manager, request, mapping = mapping)
    second <- runIntegrationEnrichment(manager, c(request, list(adjust_ruv_factors = FALSE)), mapping = mapping)
    expect_identical(second$run$status, "hit")
    expect_identical(second$run$reference$entity_id, first$run$reference$entity_id)
    expect_identical(first$run$value$design_covariates, "none")
    expect_equal(as.data.frame(first$run$value$terms), as.data.frame(item$value$terms), tolerance = 1e-12)
    # K9 (a): the replay repeats both entries of the base as they were saved.
    comparison <- integrationReplayEnrichment(manager, recorded)
    expect_true(all(comparison$terms_identical))
    expect_setequal(names(intRfEntries(manager)), c("Factor1|metabolomics_GCMS|kegg|camera", "Factor1|metabolomics_GCMS|kegg"))
    replayed <- intRfEntries(manager)[["Factor1|metabolomics_GCMS|kegg|camera"]]
    expect_false(replayed$settings$adjust_ruv_factors)
    value <- integrationRecordedEnrichment(manager)[["Factor1|metabolomics_GCMS|kegg|camera"]]$value
    expect_identical(value$design_columns, c("intercept", "factor"))
    # K10, second case: a record holding the entry stored before this ticket names its form.
    block <- buildWorkflowAnalysisRecord(manager)$analyses$integration_mofa_model$integration
    ranking <- integrationRecordEnrichmentTables(block)$summary[[integrationRecordTexts()[["col_enrich_ranking"]]]]
    expect_true(any(endsWith(ranking, "Adjusted for removed RUV factors: no.")))
    expect_false(any(grepl("yes (", ranking, fixed = TRUE)))
})

# --- No placeholder of the simulation's numbers in any shipped text (VAL-160; the numbers are written in) ---

test_that("KIT-184: no placeholder of the simulation's numbers is left in a shipped text", {
    placeholder <- "\\[[RIP][A-Za-z0-9]{1,3}\\]"
    expect_true(grepl(placeholder, "called in [R3] of cases"))
    expect_true(grepl(placeholder, "found in [PF] of cases"))
    expect_true(grepl(placeholder, "(95 % interval [I3])"))
    shipped <- c(unlist(integrationGuiTexts()), unlist(integrationRecordTexts()), INTEGRATION_ENRICHMENT_REFUSALS)
    expect_gt(length(shipped), 100L)
    expect_identical(names(shipped)[grepl(placeholder, shipped)], character())
    expect_identical(names(shipped)[grepl("has not been measured yet", shipped, fixed = TRUE)], character())
    starter <- system.file("workbooks", "integration", "integration_workflow_starter.rmd", package = "MultiScholaR")
    expect_true(nzchar(starter))
    expect_identical(grep(placeholder, readLines(starter, warn = FALSE), value = TRUE), character())
})

# --- The cross-platform record (ticket 7.5): the removed factors of K1's fixture and the terms of K6's ---

test_that("KIT-184 cross-platform: K1's removed factors and K6's terms agree with the Linux reference", {
    hand <- intRfHand()
    factors <- integrationInputRemovedFactors(intRfTable(hand$normalised), intRfTable(hand$corrected), 2L, hand$runs)$table
    unit <- intRfUnit(factors = intRfUnitFactors())
    # K6's terms are those of the adjusted test, which the reference holds: asked for by name, the default being off.
    terms <- intRfUnitProduct(unit, adjust_ruv_factors = TRUE)$terms
    terms <- terms[order(terms$ID), c("ID", "Direction", "pvalue", "p.adjust", "correlation", "effect")]
    text <- c(utils::capture.output(utils::write.table(format(factors, digits = 17), sep = "\t", quote = FALSE,
        row.names = FALSE)), utils::capture.output(utils::write.table(format(terms, digits = 17), sep = "\t", quote = FALSE,
        row.names = FALSE)))
    cat("\nKIT-184 cross-platform factors and terms (sha256 ", digest::digest(text, algo = "sha256"), "):\n", sep = "")
    cat(text, sep = "\n")
    reference <- file.path(intRfData(), "kit184-cross-platform.tsv")
    testthat::skip_if_not(file.exists(reference), "the cross-platform reference is not recorded yet")
    lines <- readLines(reference)
    split_at <- which(startsWith(lines, "ID\t"))
    ref_factors <- utils::read.delim(text = lines[seq_len(split_at - 1L)], stringsAsFactors = FALSE, check.names = FALSE)
    ref_terms <- utils::read.delim(text = lines[split_at:length(lines)], stringsAsFactors = FALSE, check.names = FALSE)
    here <- t(as.matrix(factors[hand$runs]))
    there <- t(as.matrix(ref_factors[hand$runs]))
    # The largest principal angle between the two spans (the sine, for small angles the angle itself).
    angle <- max(abs(here - there %*% crossprod(there, here)))
    expect_lte(angle, 1e-8, label = sprintf("largest principal angle to the reference's factors (%g)", angle))
    expect_identical(terms$ID, ref_terms$ID)
    expect_identical(terms$Direction, ref_terms$Direction)
    for (column in c("pvalue", "p.adjust", "correlation", "effect")) {
        difference <- max(abs(terms[[column]] - ref_terms[[column]]))
        expect_lte(difference, 1e-6, label = sprintf("largest absolute difference of %s (%g)", column, difference))
    }
})
