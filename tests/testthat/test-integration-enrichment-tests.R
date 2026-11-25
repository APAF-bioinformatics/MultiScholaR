# KIT-183 (VAL-51, VAL-53): the integration enrichment's test of record, limma's camera on each feature's regression on
# the oriented factor (the correlation between a set's features estimated), the pre-ranked test (GSEA) selectable, the
# labelled fallback cameraPR implemented for the simulation. Declared before the run (ticket section 10; each test's
# failing value on the base is stated beside it):
#   T1  the entry point's camera and cameraPR equal direct limma calls; `effect` is camera's own contrast (camera's
#       p-value recomputed from it); BH within the call (base: no `test` argument / no design function);
#   T2  a correlated null set with a small common shift: camera p > 0.05 with its correlation >= 0.4, cameraPR p < 0.01;
#   T3  sets counted and tested on U_cam, not the ranked list; absent identifiers change nothing; a complete feature
#       the model did not select is in U_cam;
#   T4  every factor score and weight multiplied by -1 gives an identical terms table;
#   T5  the base store's KEGG entry (made at base, `kit183-base-store`) replays with gsea, terms identical;
#   T6  collapse among complete features by the largest mean; counts by hand;
#   T7  refusals; STRING ignores the test (key, record note, never refused);
#   T8  restated (VAL-85 (3)): (a) the declarations are base's; (b) the base store's identities equal the head's apart
#       from the build; (c) a hit within one build;
#   T9  camera declares no seed;
#   T10 record and Methods name each entry's test; the limits sentence once; factor strength;
#   T11 the tab: selector, default, notes, counts, caution, relabelled seed, key with the test;
#   T12 the workbook: gsea override gives base's terms, the default is camera;
#   T13 corpus check (6.6), with the size of the change written to the output;
#   T15 integrationEnrichmentGsea unchanged (formals, body, value).
# KIT-191 restated, where its texts changed by design (VAL-85 (2), VAL-97): T10 (the score column, the limits paragraph
# per test and kind, the strength sentence per test and kind), T11 (the seed label, the notes per kind, the
# few-samples sentences per test and kind), T15 under D9 (a) (the body differs from base only by the three RNG kinds).

# --- Fixtures ---

intEtData <- function() normalizePath(testthat::test_path("..", "testdata", "integration"), mustWork = FALSE)
intEtBaseStore <- function() file.path(intEtData(), "kit183-base-store")
intEtStored <- function(name) utils::read.delim(file.path(intEtData(), "mofa_reference", paste0(name, ".tsv")),
    check.names = FALSE, stringsAsFactors = FALSE)
intEtSettings <- function(min_size = 5L, max_size = 500L) {
    integrationGseaSettings(p_adjust_method = "BH", min_size = min_size, max_size = max_size, p_cutoff = 1, tie_jitter = 0)
}

# The stored reference model as the trainer's value (as test-integration-enrichment.R); `sign` multiplies every factor
# score and weight.
intEtTrainer <- function(sign = 1) {
    function(views, options) {
        factors <- intEtStored("factors")
        f <- sign * as.matrix(factors[-1L])
        rownames(f) <- factors$biological_sample
        f <- f[colnames(views[[1L]]), , drop = FALSE]
        w <- intEtStored("weights")
        weights <- lapply(stats::setNames(names(views), names(views)), function(view) {
            rows <- w[w$view == view, , drop = FALSE]
            m <- sign * as.matrix(rows[-(1:2)])
            rownames(m) <- paste(view, rows$feature_id, sep = ":")
            m[rownames(views[[view]]), , drop = FALSE]
        })
        v <- intEtStored("variance_explained")
        r2 <- matrix(v$r2_percent, nrow = length(unique(v$factor)), dimnames = list(unique(v$factor), unique(v$view)))
        r2 <- r2[, names(views), drop = FALSE]
        list(factors = f, weights = weights, r2_per_factor = r2, r2_total = colSums(r2), elbo = c(-2, -1), iterations = 186L,
            options = options, versions = list(MOFA2 = "stored reference"), messages = character(),
            model_file = list(bytes = as.raw(1:3), sha256 = strrep("d", 64L), size = 3))
    }
}

# Invented KEGG compound identifiers for metabolomics_GCMS (the same rule as capture-base-store.R): one per feature,
# the first 20 features in pairs on one identifier.
intEtKeggIds <- function(n) {
    ids <- sprintf("C%05d", 10000L + seq_len(n))
    pairs <- seq(2L, 20L, by = 2L)
    ids[pairs] <- ids[pairs - 1L]
    ids
}

# A project holding the corpus, metabolomics_GCMS given a `kegg_id` column (rewritten with the package's writer).
intEtProject <- function(omics = c("lipidomics", "metabolomics", "proteomics", "transcriptomics"),
    .local_envir = parent.frame()) {
    testthat::skip_if_not(dir.exists(file.path(intEtData(), "corpus")), "the tier-1 corpus is not present")
    testthat::skip_if_not(workflowOwnerCapable(), "the integration store needs qs2, DuckDB and Arrow")
    root <- withr::local_tempdir(.local_envir = .local_envir)
    results_dir <- file.path(root, "results", "integration")
    paths <- list(base_dir = root, integration_dir = file.path(root, "integration"),
        source_dir = file.path(root, "scripts", "integration"), results_dir = results_dir,
        mofa_tables_dir = file.path(results_dir, "mofa", "tables"), mofa_plots_dir = file.path(results_dir, "mofa", "plots"),
        study_report_dir = file.path(root, "results_summary", "integration", "Study_report"))
    dir.create(paths$source_dir, recursive = TRUE)
    for (omic in omics) {
        dir <- file.path(paths$integration_dir, omic)
        dir.create(dir, recursive = TRUE)
        stopifnot(all(file.copy(list.files(file.path(intEtData(), "corpus", omic), full.names = TRUE), dir)))
    }
    dir <- file.path(paths$integration_dir, "metabolomics")
    input <- readIntegrationInputWithMatrices(dir)
    read <- input$manifest
    manifest <- list(omic_type = read$omic_type, workflow = read$workflow, matrix_role = read$matrix_role, scale = read$scale,
        correction = NULL, source_state = list(name = "kit183-kegg"), sample_id_column = read$sample_id_column,
        group_column = read$group_column, technical_replicate_column = read$technical_replicate_column,
        biological_sample_column = read$biological_sample_column)
    features <- lapply(input$views, `[[`, "features")
    features$metabolomics_GCMS$kegg_id <- intEtKeggIds(nrow(features$metabolomics_GCMS))
    writeIntegrationInput(dir, manifest, input$samples, lapply(input$views, `[[`, "matrix"), features)
    paths
}

intEtModel <- function(paths, parameters = list(), trainer = intEtTrainer(), .local_envir = parent.frame()) {
    manager <- newIntegrationStateOwner(file.path(paths$source_dir, paste0("store-", basename(tempfile()))))
    withr::defer(closeWorkflowStateOwner(manager), envir = .local_envir)
    scan <- integrationScanProject(paths$integration_dir)
    integrationCommitImport(manager, integrationImportRequest(scan, scan$views$view))
    run <- runIntegrationModel(manager, parameters, trainer = trainer)
    list(manager = manager, run = run)
}

# A KEGG mapping in `integrationDownloadKeggMapping()`'s form, from a named list of sets of compound identifiers.
intEtMapping <- function(sets) {
    pathway <- paste0("path:hsa", sprintf("%05d", seq_along(sets)))
    list(database = "KEGG", release = "KEGG Release fixture", sha256 = strrep("e", 64L), url = "fixture",
        downloaded_at = "2026-11-23T00:00:00+0000",
        table = data.frame(pathway = rep(pathway, lengths(sets)), compound = paste0("cpd:", unlist(sets, use.names = FALSE)),
            name = rep(names(sets), lengths(sets)), stringsAsFactors = FALSE))
}

# The data set averaged per biological sample over S_d, computed on the whole matrix (independently of the pages).
intEtIndependentValues <- function(manager, value, view) {
    object <- integrationImportedObject(manager)
    table <- object@views[[view]]
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
    list(x = x, samples = samples, features = object@features[[view]])
}

# The determinism fixture of test-integration-enrichment.R (`pathwayFixture`), as capture-base-store.R builds it.
intEtGseaFixture <- function() {
    ids <- sprintf("G%02d", 1:60)
    ranked <- data.frame(identifier = ids, value = seq(3, -3, length.out = 60), stringsAsFactors = FALSE)
    term2gene <- data.frame(term = rep(c("top", "bottom", "spread", "middle"), c(10, 10, 10, 15)),
        gene = c(ids[1:10], ids[51:60], ids[seq(1, 60, by = 6)], ids[21:35]), stringsAsFactors = FALSE)
    term2name <- data.frame(term = c("top", "bottom", "spread", "middle"), name = c("Top set", "Bottom set",
        "Spread set", "Middle set"), stringsAsFactors = FALSE)
    list(ranked = ranked, term2gene = term2gene, term2name = term2name)
}

# A small data set in memory (no store): 40 metabolites in 12 biological samples, one run each, with a factor.
#   C00001: m01 (mean 18, weakly related) and m02 (mean 15, |r| larger), both complete: camera keeps m01.
#   C00002: m03 (|r| the largest, missing in s05) and m04 (complete): camera tests C00002 through m04.
#   m05, m06 incomplete. Counts by hand: 40 mapped features, 37 complete, 3 incomplete, 1 collapsed, 36 identifiers.
intEtUnit <- function() {
    samples <- sprintf("s%02d", 1:12)
    runs <- sprintf("r%02d", 1:12)
    f <- c(-1.6, -1.2, -0.9, -0.6, -0.4, -0.1, 0.1, 0.4, 0.6, 0.9, 1.2, 1.6)
    ids <- sprintf("m%02d", 1:40)
    kegg <- c("C00001", "C00001", "C00002", "C00002", sprintf("C%05d", 5:40))
    noise <- withr::with_seed(20261183L, matrix(stats::rnorm(40 * 12, sd = 0.3), 40, 12))
    mu <- c(18, 15, 16, 17, withr::with_seed(7L, stats::runif(36, 14, 20)))
    x <- mu + noise
    x[2L, ] <- x[2L, ] + 0.9 * f
    x[3L, ] <- x[3L, ] + 1.2 * f
    x[11:15, ] <- x[11:15, ] + 0.5 * matrix(f, 5, 12, byrow = TRUE)
    x[3L, 5L] <- NA
    x[5L, 2L] <- NA
    x[6L, 9L] <- NA
    dimnames(x) <- list(ids, runs)
    table <- data.frame(feature_id = ids, x, check.names = FALSE, stringsAsFactors = FALSE)
    features <- data.frame(feature_id = ids, display_name = ids, kegg = kegg, stringsAsFactors = FALSE)
    object <- methods::new("IntegrationViewsData", views = list(metab = table[0L, , drop = FALSE]),
        features = list(metab = features),
        design_matrix = data.frame(omic = "metabolomics", run = runs, group = rep(c("A", "B"), each = 6),
            biological_sample = samples, stringsAsFactors = FALSE),
        inputs = list(metabolomics = list(omic = "metabolomics", views = "metab")))
    value <- list(factors = data.frame(biological_sample = samples, Factor1 = f, stringsAsFactors = FALSE),
        weights = list(metab = "unused"))
    correlation <- vapply(seq_len(40), function(i) {
        ok <- !is.na(x[i, ])
        stats::cor(x[i, ok], f[ok])
    }, numeric(1))
    ranking <- list(factor = "Factor1", view = "metab", min_samples = 5L,
        direction = list(rule = "largest_absolute_weight_positive", sign = 1L, view = "metab", feature_id = "m03", weight = 1),
        counts = list(features = 40L, ranked = 40L, too_few_samples = 0L, no_variation = 0L, internal_standard = 0L),
        ranking = data.frame(feature_id = ids, n_samples = as.integer(rowSums(!is.na(x))), correlation = correlation,
            status = "ranked", stringsAsFactors = FALSE))
    p <- list(factor = "Factor1", view = "metab", omic = "metabolomics", database = "kegg", id_type = "kegg",
        id_column = "kegg", correlation_method = "pearson", min_samples = 5L, kegg_organism = "hsa", min_set_size = 5L,
        max_set_size = 500L, seed = 20261121L, test = "camera")
    pages <- function(view, columns, kernel) kernel(table[columns], 0L)
    list(object = object, value = value, ranking = ranking, p = p, pages = pages, x = x, f = f, kegg = kegg)
}

intEtUnitValues <- function(unit, covariates = NULL) {
    mapped <- integrationEnrichmentMappedList(unit$ranking, unit$p, unit$object@features$metab)
    integrationEnrichmentTestValues(unit$object, NULL, NULL, unit$value, list(), unit$ranking, mapped, "metab",
        pages = unit$pages, covariates = covariates)
}

intEtUnitKegg <- function(extra = list()) {
    sets <- c(list(P6 = sprintf("C%05d", 5:10), PA = sprintf("C%05d", 11:20), PB = sprintf("C%05d", 21:30),
        PC = c("C00001", "C00002", sprintf("C%05d", 31:35))), extra)
    intEtMapping(sets)
}

# --- T1, T2: exactness on a matrix fixture ---

test_that("KIT-183 T1: camera and cameraPR at the entry point equal direct limma calls; effect is camera's contrast", {
    n <- 12L
    g <- 40L
    f <- c(-1.5, -1.3, -1.0, -0.7, -0.4, -0.2, 0.1, 0.3, 0.6, 0.9, 1.2, 1.6)
    mu <- withr::with_seed(11L, stats::runif(g, 14, 22))
    sd <- 0.2 + 0.6 * exp(-(mu - 14) / 2)
    y <- mu + sd * withr::with_seed(12L, matrix(stats::rnorm(g * n), g, n))
    y[1:8, ] <- y[1:8, ] + 0.8 * matrix(f, 8, n, byrow = TRUE)
    y[9:16, ] <- y[9:16, ] - 0.6 * matrix(f, 8, n, byrow = TRUE)
    rownames(y) <- sprintf("C%05d", seq_len(g))
    ids <- rownames(y)
    sets <- list(down = ids[9:16], mixed = ids[c(1:4, 20:25)], null = ids[26:35], up = ids[1:8])
    term2gene <- data.frame(term = rep(names(sets), lengths(sets)), gene = unlist(sets, use.names = FALSE),
        stringsAsFactors = FALSE)
    term2name <- data.frame(term = names(sets), name = paste(names(sets), "set"), stringsAsFactors = FALSE)
    design <- integrationEnrichmentTestDesign(f)
    expect_identical(colnames(design), c("intercept", "factor"))
    expect_equal(unname(design), unname(cbind(1, f)))
    values <- list(y = y, design = design, counts = list(universe = g), design_covariates = "none")
    ranked <- data.frame(identifier = ids, value = stats::cor(t(y), f)[, 1L], stringsAsFactors = FALSE)
    out <- integrationEnrichmentPathwayTest(ranked, term2gene, term2name, 1L, intEtSettings(), test = "camera",
        values = values)
    index <- split(term2gene$gene, term2gene$term)
    direct <- limma::camera(y, index, cbind(1, f), contrast = 2, inter.gene.cor = NA, allow.neg.cor = TRUE,
        trend.var = TRUE, use.ranks = FALSE, sort = FALSE)
    terms <- out$terms[match(names(index), out$terms$ID), , drop = FALSE]
    expect_identical(out$method, "camera")
    expect_identical(terms$ID, names(index))
    expect_equal(terms$pvalue, direct$PValue)
    expect_identical(terms$Direction, direct$Direction)
    expect_equal(terms$correlation, direct$Correlation)
    expect_identical(terms$setSize, as.integer(direct$NGenes))
    expect_identical(terms$Description, paste(names(index), "set"))
    expect_identical(unique(terms$test), "camera")
    expect_equal(out$terms$p.adjust, stats::p.adjust(out$terms$pvalue, "BH"))
    # effect is camera's own delta: its two-sample t, from effect, the correlation and var(z), gives camera's p-value.
    fit <- limma::eBayes(limma::lmFit(y, cbind(1, f)), trend = TRUE)
    z <- limma::zscoreT(fit$t[, 2L], df = fit$df.total, approx = TRUE, method = "hill")
    m <- terms$setSize
    vif <- 1 + (m - 1) * terms$correlation
    pooled <- ((g - 1) * stats::var(z) - terms$effect^2 * m * (g - m) / g) / (g - 2)
    t <- terms$effect / sqrt(pooled * (vif / m + 1 / (g - m)))
    expect_equal(2 * stats::pt(-abs(t), df = min(n - 2L, g - 2L)), terms$pvalue, tolerance = 1e-6)
    expect_identical(terms$Direction, ifelse(terms$effect > 0, "Up", "Down"))
    expect_identical(terms$Direction[terms$ID %in% c("up", "down")], c("Down", "Up"))
    # cameraPR on camera's own statistic, with the correlation assumed to be 0.01.
    pr <- integrationEnrichmentPathwayTest(ranked, term2gene, term2name, 1L, intEtSettings(), test = "camera_pr",
        values = values)
    direct_pr <- limma::cameraPR(stats::setNames(z, ids), index, inter.gene.cor = 0.01, sort = FALSE)
    pr_terms <- pr$terms[match(names(index), pr$terms$ID), , drop = FALSE]
    expect_equal(pr_terms$pvalue, direct_pr$PValue)
    expect_identical(pr_terms$Direction, direct_pr$Direction)
    expect_identical(unique(pr_terms$correlation), 0.01)
    expect_identical(unique(pr_terms$test), "camera_pr")
    # The sensitivity arm's floor (allow.neg.cor = FALSE) through the settings only.
    floor <- integrationEnrichmentPathwayTest(ranked, term2gene, term2name, 1L, c(intEtSettings(), list(allow_neg_cor = FALSE)),
        test = "camera", values = values)
    direct_floor <- limma::camera(y, index, cbind(1, f), contrast = 2, inter.gene.cor = NA, allow.neg.cor = FALSE,
        trend.var = TRUE, use.ranks = FALSE, sort = FALSE)
    expect_equal(floor$terms$pvalue[match(names(index), floor$terms$ID)], direct_floor$PValue)
    # The pre-ranked test is unchanged at the entry point.
    gsea <- integrationEnrichmentPathwayTest(ranked, term2gene, term2name, 1L, intEtSettings())
    expect_identical(gsea$method, "gsea")
    # An added covariate column never moves the tested coefficient: camera with the contrast by name.
    w <- matrix(withr::with_seed(13L, stats::rnorm(n)), ncol = 1L, dimnames = list(NULL, "ruv_W1"))
    covariate_design <- integrationEnrichmentTestDesign(f, w)
    expect_identical(colnames(covariate_design), c("intercept", "factor", "ruv_W1"))
    with_w <- integrationEnrichmentPathwayTest(ranked, term2gene, term2name, 1L, intEtSettings(), test = "camera",
        values = list(y = y, design = covariate_design, counts = list(), design_covariates = "x"))
    direct_w <- limma::camera(y, index, cbind(1, w, f), contrast = 3, inter.gene.cor = NA, allow.neg.cor = TRUE,
        trend.var = TRUE, use.ranks = FALSE, sort = FALSE)
    expect_equal(with_w$terms$pvalue[match(names(index), with_w$terms$ID)], direct_w$PValue)
    expect_error(integrationEnrichmentTestDesign(f, matrix(1, n, 1L, dimnames = list(NULL, "factor"))),
        class = "integration_enrichment_parameters_invalid")
    expect_error(integrationEnrichmentTestDesign(f, cbind(a = f)), class = "integration_enrichment_parameters_invalid")
})

test_that("KIT-183 T2: a correlated null set with a small common shift: camera does not call it, cameraPR does", {
    n <- 12L
    g <- 400L
    f <- c(-1.5, -1.3, -1.0, -0.7, -0.4, -0.2, 0.1, 0.3, 0.6, 0.9, 1.2, 1.6)
    fs <- as.numeric(scale(f))
    w <- stats::residuals(stats::lm(withr::with_seed(21L, stats::rnorm(n)) ~ f))
    w <- as.numeric(scale(w))
    u <- 0.4 * fs + sqrt(0.84) * w
    e <- withr::with_seed(22L, matrix(stats::rnorm(g * n), g, n))
    y <- 16 + 0.5 * e
    y[1:20, ] <- 16 + 0.5 * (sqrt(0.6) * matrix(u, 20, n, byrow = TRUE) + sqrt(0.4) * e[1:20, ])
    ids <- sprintf("C%05d", seq_len(g))
    rownames(y) <- ids
    term2gene <- data.frame(term = "block", gene = ids[1:20], stringsAsFactors = FALSE)
    term2name <- data.frame(term = "block", name = "Correlated block", stringsAsFactors = FALSE)
    values <- list(y = y, design = integrationEnrichmentTestDesign(f), counts = list(), design_covariates = "none")
    ranked <- data.frame(identifier = ids, value = stats::cor(t(y), f)[, 1L], stringsAsFactors = FALSE)
    camera <- integrationEnrichmentPathwayTest(ranked, term2gene, term2name, 1L, intEtSettings(), test = "camera",
        values = values)$terms
    pr <- integrationEnrichmentPathwayTest(ranked, term2gene, term2name, 1L, intEtSettings(), test = "camera_pr",
        values = values)$terms
    expect_gt(camera$pvalue, 0.05)
    expect_gte(camera$correlation, 0.4)
    expect_lt(pr$pvalue, 0.01)
})

# --- T3, T6, T7: universe, collapse, counts and refusals on the small data set ---

test_that("KIT-183 T6: among complete features on one identifier, camera keeps the larger mean; counts by hand", {
    unit <- intEtUnit()
    values <- intEtUnitValues(unit)
    expect_identical(values$counts[c("samples", "mapped", "complete", "incomplete", "collapsed", "universe")],
        list(samples = 12L, mapped = 40L, complete = 37L, incomplete = 3L, collapsed = 1L, universe = 36L))
    expect_identical(values$feature_id[match("C00001", rownames(values$y))], "m01")
    expect_identical(values$feature_id[match("C00002", rownames(values$y))], "m04")
    expect_false(any(c("m02", "m03", "m05", "m06") %in% values$feature_id))
    expect_identical(rownames(values$y), unique(unit$kegg[-c(3L, 5L, 6L)]))
    expect_equal(unname(values$y["C00001", ]), unname(unit$x["m01", ]))
    expect_identical(values$samples, sprintf("s%02d", 1:12))
    expect_identical(values$design_columns, c("intercept", "factor"))
    expect_identical(values$design_covariates, "none")
    # The pre-ranked test keeps the larger |r| (INT-19, unchanged).
    mapped <- integrationEnrichmentMappedList(unit$ranking, unit$p, unit$object@features$metab)
    expect_identical(mapped$table$feature_id[mapped$table$identifier == "C00001"], "m02")
    expect_identical(mapped$table$feature_id[mapped$table$identifier == "C00002"], "m03")
    # Covariates keyed by biological sample, a superset allowed; a missing sample refused.
    w <- matrix(seq_len(13), ncol = 1L, dimnames = list(c(sprintf("s%02d", 13:1)), "ruv_W1"))
    with_w <- intEtUnitValues(unit, w)
    expect_identical(with_w$design_columns, c("intercept", "factor", "ruv_W1"))
    expect_equal(unname(with_w$design[, "ruv_W1"]), as.numeric(13:2))
    expect_identical(with_w$y, values$y)
    expect_error(intEtUnitValues(unit, w[rownames(w) != "s05", , drop = FALSE]),
        class = "integration_enrichment_parameters_invalid")
})

test_that("KIT-183 T3: sets are counted and tested on the complete identifiers; absent identifiers change nothing", {
    skip_if_not_installed("clusterProfiler")
    unit <- intEtUnit()
    values <- intEtUnitValues(unit)
    meta <- intEtUnitKegg()[c("database", "release", "sha256", "url", "downloaded_at")]
    kernel <- function(test, mapping) {
        p <- unit$p
        p$test <- test
        integrationEnrichmentKernel(unit$object, list(p = p, ranking_value = unit$ranking, mapping_table = mapping$table,
            mapping_meta = meta, values = if (test != "gsea") values))
    }
    camera <- kernel("camera", intEtUnitKegg())
    gsea <- kernel("gsea", intEtUnitKegg())
    # P6: 6 ranked members, 4 of them complete: tested by gsea, not by camera.
    expect_identical(camera$service$pathways_testable, 3L)
    expect_setequal(camera$terms$ID, c("path:hsa00002", "path:hsa00003", "path:hsa00004"))
    expect_identical(gsea$service$pathways_testable, 4L)
    expect_true("path:hsa00001" %in% gsea$terms$ID)
    expect_identical(camera$terms$setSize[match("path:hsa00004", camera$terms$ID)], 7L)
    expect_identical(camera$test, "camera")
    expect_identical(camera$test_counts$universe, 36L)
    expect_null(camera$service$seed)
    expect_identical(camera$service$test_settings, INTEGRATION_ENRICHMENT_TEST_SETTINGS[["camera"]])
    expect_null(gsea[["test"]])
    # 300 identifiers absent from the data, in a set of their own and added to a tested set: nothing changes.
    absent <- sprintf("C%05d", 50000L + seq_len(300L))
    padded <- intEtUnitKegg(list(PX = absent))
    padded$table <- rbind(padded$table, data.frame(pathway = "path:hsa00002", compound = paste0("cpd:", absent),
        name = "PA", stringsAsFactors = FALSE))
    again <- kernel("camera", padded)
    expect_identical(again$terms, camera$terms)
    expect_identical(again$service$pathways_testable, camera$service$pathways_testable)
})

test_that("KIT-183 T7: too few complete identifiers or samples are refused; no other test runs in its place", {
    unit <- intEtUnit()
    few <- unit
    # Every feature but nine missing in one sample: 9 complete identifiers.
    table <- data.frame(feature_id = sprintf("m%02d", 1:40), unit$x, check.names = FALSE, stringsAsFactors = FALSE)
    table[-c(1, 7:14), 2L] <- NA
    few$pages <- function(view, columns, kernel) kernel(table[columns], 0L)
    condition <- expect_error(intEtUnitValues(few), class = "integration_enrichment_too_few_complete")
    expect_s3_class(condition, "integration_refused")
    expect_lt(condition$counts$universe, 10L)
    expect_match(conditionMessage(condition), "pre-ranked test can use features with missing values", fixed = TRUE)
    # Three biological samples.
    small <- unit
    small$object@design_matrix$biological_sample <- rep(sprintf("s%02d", 1:3), each = 4L)
    small$value$factors <- data.frame(biological_sample = sprintf("s%02d", 1:3), Factor1 = c(-1, 0, 1))
    condition <- expect_error(intEtUnitValues(small), class = "integration_enrichment_too_few_samples_for_test")
    expect_match(conditionMessage(condition), "needs at least 5 biological samples", fixed = TRUE)
    # The floor counts the design's columns: 12 samples, 9 covariates need 13.
    w <- matrix(withr::with_seed(5L, stats::rnorm(12 * 9)), 12, 9,
        dimnames = list(sprintf("s%02d", 1:12), paste0("ruv_W", 1:9)))
    expect_error(intEtUnitValues(unit, w), class = "integration_enrichment_too_few_samples_for_test")
    # A limma error is refused; no other test is run in its place.
    values <- intEtUnitValues(unit)
    local_mocked_bindings(camera = function(...) stop("the engine broke"), .package = "limma")
    gsea_calls <- 0L
    local_mocked_bindings(integrationEnrichmentGsea = function(...) {
        gsea_calls <<- gsea_calls + 1L
        stop("the pre-ranked test must not run")
    })
    term2gene <- data.frame(term = "PA", gene = sprintf("C%05d", 11:20), stringsAsFactors = FALSE)
    condition <- expect_error(integrationEnrichmentPathwayTest(data.frame(identifier = "C00011", value = 1), term2gene,
        data.frame(term = "PA", name = "PA"), 1L, intEtSettings(), test = "camera", values = values),
        class = "integration_enrichment_method_failed")
    expect_match(condition$engine_message, "the engine broke", fixed = TRUE)
    expect_identical(gsea_calls, 0L)
})

test_that("KIT-183 T7: the test setting: resolved, refused when not offered, ignored by STRING (key, record note)", {
    paths <- intEtProject()
    fixture <- intEtModel(paths)
    manager <- fixture$manager
    value <- fixture$run$value
    object <- integrationImportedObject(manager, fields = c("features", "design_matrix", "inputs"))
    resolve <- function(...) integrationResolveEnrichmentParameters(list(factor = "Factor1", ...), value, object, 20L)
    kegg <- resolve(view = "metabolomics_GCMS", database = "kegg", id_column = "kegg_id")
    expect_identical(kegg$test, "camera")
    expect_false(kegg$test_ignored)
    expect_identical(resolve(view = "metabolomics_GCMS", database = "kegg", id_column = "kegg_id", test = "gsea")$test, "gsea")
    expect_error(resolve(view = "metabolomics_GCMS", database = "kegg", id_column = "kegg_id", test = "camera_pr"),
        class = "integration_enrichment_test_not_offered")
    expect_error(resolve(view = "metabolomics_GCMS", database = "kegg", id_column = "kegg_id", test = "fry"),
        class = "integration_enrichment_parameters_invalid")
    string <- resolve(view = "proteomics", database = "string", test = "camera")
    expect_true(is.na(string$test))
    expect_true(string$test_ignored)
    expect_false(resolve(view = "proteomics", database = "string")$test_ignored)
    expect_identical(integrationEnrichmentKey("Factor1", "proteomics", "string", "camera"), "Factor1|proteomics|string")
    expect_identical(integrationEnrichmentKey("Factor1", "m", "kegg", "camera"), "Factor1|m|kegg|camera")
    expect_identical(integrationEnrichmentKey("Factor1", "m", "kegg", "gsea"), "Factor1|m|kegg")
    texts <- integrationRecordTexts()
    entry <- list(database = "string", settings = list(correlation_method = "pearson", min_samples = 5L, test_ignored = TRUE))
    expect_identical(integrationRecordEnrichmentRanking(entry, texts), paste0(sprintf(texts[["enrich_ranking_value"]],
        "pearson", 5L), "; ", texts[["enrich_string_test_note"]]))
    entry$settings$test_ignored <- FALSE
    expect_identical(integrationRecordEnrichmentRanking(entry, texts), sprintf(texts[["enrich_ranking_value"]], "pearson", 5L))
})

# --- Store level on the corpus: T3 (U_cam on the store), T4, T9, T10, T13 ---

test_that("KIT-183 T3, T9, T10: camera through the store: U_cam from the views, no seed, record and Methods per test", {
    skip_if_not_installed("clusterProfiler")
    paths <- intEtProject()
    fixture <- intEtModel(paths, list(features_per_view = 50L))
    manager <- fixture$manager
    value <- fixture$run$value
    view <- "metabolomics_GCMS"
    independent <- intEtIndependentValues(manager, value, view)
    complete <- rownames(independent$x)[stats::complete.cases(independent$x)]
    ids <- stats::setNames(independent$features$kegg_id, independent$features$feature_id)
    single <- names(ids)[!duplicated(ids) & !duplicated(ids, fromLast = TRUE)]
    pool <- intersect(single, complete)
    sets <- list(first = ids[pool[1:10]], second = ids[pool[11:20]], third = ids[pool[21:30]])
    mapping <- intEtMapping(sets)
    request <- list(factor = "Factor1", view = view, database = "kegg", id_column = "kegg_id")
    camera <- runIntegrationEnrichment(manager, request, mapping = mapping)
    expect_identical(camera$run$value$test, "camera")
    expect_identical(camera$entry$settings$test, "camera")
    expect_identical(camera$run$value$design_columns, c("intercept", "factor"))
    expect_identical(camera$run$value$design_covariates, "none")
    # T1 store variant: the values read through the store equal an independent read of `views` averaged per sample.
    export <- manager$.context$snapshot$exports$integration_mofa_model
    object <- integrationImportedObject(manager, fields = c("features", "design_matrix", "inputs"))
    mapped <- integrationEnrichmentMappedList(camera$ranking$value, list(database = "kegg", omic = "metabolomics",
        id_type = "kegg", id_column = "kegg_id"), object@features[[view]])
    values <- integrationEnrichmentTestValues(object, manager, export$source, value, export$parameters,
        camera$ranking$value, mapped, view)
    expect_equal(unname(values$y), unname(independent$x[values$feature_id, values$samples, drop = FALSE]), tolerance = 1e-12)
    expect_identical(values$samples, independent$samples)
    # A complete feature that the model did not select is in U_cam and counted.
    weights <- integrationModelWeights(manager, value)
    selected <- weights$feature_id[weights$view == view]
    expect_true(any(!values$feature_id %in% selected))
    expect_identical(camera$run$value$test_counts$complete, sum(stats::complete.cases(independent$x[names(ids)[
        names(ids) %in% mapped$mapping$feature_id[mapped$mapping$status %in% c("mapped", "duplicated")]], , drop = FALSE])))
    # T9: camera declares no seed: another seed is the same product.
    reseeded <- runIntegrationEnrichment(manager, c(request, list(seed = 7L)), mapping = mapping)
    expect_identical(reseeded$run$status, "hit")
    expect_identical(reseeded$run$reference$entity_id, camera$run$reference$entity_id)
    expect_identical(reseeded$run$value$terms, camera$run$value$terms)
    # The pre-ranked test beside it, under its own key.
    gsea <- runIntegrationEnrichment(manager, c(request, list(test = "gsea")), mapping = mapping)
    expect_identical(gsea$run$value$service$method, "gsea")
    expect_false(identical(gsea$run$reference$entity_id, camera$run$reference$entity_id))
    keys <- names(manager$.context$snapshot$exports$integration_mofa_model$enrichment)
    expect_setequal(keys, c("Factor1|metabolomics_GCMS|kegg|camera", "Factor1|metabolomics_GCMS|kegg"))
    # T10: the record names each entry's test.
    texts <- integrationRecordTexts()
    record <- buildWorkflowAnalysisRecord(manager)
    block <- record$analyses$integration_mofa_model$integration
    tables <- integrationRecordEnrichmentTables(block, texts)
    ranking <- tables$summary[[texts[["col_enrich_ranking"]]]]
    counts <- camera$run$value$test_counts
    # KIT-184 (D04): every result of a test that has a design names its form in its row.
    expect_setequal(ranking, c(paste0(sprintf(texts[["enrich_ranking_value_camera"]], counts$samples, counts$universe),
        ". Adjusted for removed RUV factors: no."), sprintf(texts[["enrich_ranking_value"]], "pearson", 5L)))
    strength <- integrationEnrichmentFactorStrength(value, "Factor1", view)
    # KIT-191: each entry's own test's figure (metabolites: camera 6.7 %, gsea 34.4 %).
    strengthOf <- function(figure, phrase) sprintf(texts[["enrich_factor_strength"]],
        integrationEnrichmentPercent(strength$factor_r2_percent), integrationEnrichmentPercent(strength$largest_r2_percent),
        "metabolite data", texts[["enrich_strength_terms_metabolites"]], figure, phrase)
    expect_setequal(unique(tables$summary[[texts[["col_enrich_strength"]]]]),
        c(strengthOf("6.7 %", "correlation-adjusted test"), strengthOf("34.4 %", "pre-ranked test")))
    expect_identical(texts[["col_enrich_score"]], paste("Enrichment score (NES), or effect: the difference in mean z-score",
        "(correlation-adjusted tests) or in mean Fisher-transformed correlation (permutation) between the term and the other",
        "features"))
    methods <- integrationMethodsEnrichment(block, texts)
    camera_sentence <- sprintf(texts[["methods_enrich_camera"]], "Factor1", view, "KEGG", texts[["methods_test_name_camera"]],
        sprintf(texts[["methods_test_settings_camera"]], counts$samples, counts$universe,
            camera$run$value$mapping_counts$mapped), "KEGG Release fixture", "hsa", 5L, 500L,
        camera$run$value$service$pathways_testable)
    expect_true(camera_sentence %in% methods)
    expect_true(any(startsWith(methods, sprintf("Factor1 in %s: KEGG pathways were tested by ranked gene-set", view))))
    expect_identical(sum(startsWith(methods, texts[["methods_enrich_limits_tests"]])), 1L)
    expect_identical(sum(startsWith(methods, "For enrichment, every measured feature")), 1L)
    expect_false(any(grepl("clusterProfiler", camera_sentence, fixed = TRUE)))
    # Only the camera entry: the ranking sentence is still printed; only STRING entries: no limits sentence.
    only_camera <- block
    only_camera$enrichment <- Filter(function(entry) identical(entry$test, "camera"), block$enrichment)
    expect_identical(sum(startsWith(integrationMethodsEnrichment(only_camera, texts), "For enrichment, every measured")), 1L)
    string <- list(complete = TRUE, database = "string", factor = "Factor1", view = "proteomics",
        settings = list(correlation_method = "pearson", min_samples = 5L), mapping_counts = list(mapped = 10L),
        service = list(version = list(string_version = "12.5"), species = "9606", queried_at = "x"))
    expect_false(any(startsWith(integrationMethodsEnrichment(list(enrichment = list(string)), texts),
        texts[["methods_enrich_limits_tests"]])))
})

test_that("KIT-183 T4: every factor score and weight multiplied by -1 gives an identical terms table", {
    skip_if_not_installed("clusterProfiler")
    view <- "metabolomics_GCMS"
    run <- function(sign) {
        paths <- intEtProject(.local_envir = parent.frame())
        fixture <- intEtModel(paths, trainer = intEtTrainer(sign), .local_envir = parent.frame())
        independent <- intEtIndependentValues(fixture$manager, fixture$run$value, view)
        ids <- stats::setNames(independent$features$kegg_id, independent$features$feature_id)
        complete <- rownames(independent$x)[stats::complete.cases(independent$x)]
        pool <- intersect(names(ids)[!duplicated(ids) & !duplicated(ids, fromLast = TRUE)], complete)
        mapping <- intEtMapping(list(a = ids[pool[1:8]], b = ids[pool[9:20]], c = ids[pool[21:30]]))
        out <- runIntegrationEnrichment(fixture$manager, list(factor = "Factor1", view = view, database = "kegg",
            id_column = "kegg_id"), mapping = mapping)
        out$run$value
    }
    plus <- run(1)
    minus <- run(-1)
    expect_identical(plus$test, "camera")
    expect_identical(plus$direction$sign, -minus$direction$sign)
    expect_identical(plus$terms, minus$terms)
})

test_that("KIT-183 T13: corpus check: planted sets found in the expected direction, a decoy not; the size of the change", {
    skip_if_not_installed("clusterProfiler")
    paths <- intEtProject()
    fixture <- intEtModel(paths)
    manager <- fixture$manager
    value <- fixture$run$value
    export <- manager$.context$snapshot$exports$integration_mofa_model
    truth <- utils::read.delim(file.path(intEtData(), "corpus", "truth_scores.tsv"), stringsAsFactors = FALSE)
    loadings <- utils::read.delim(file.path(intEtData(), "corpus", "truth_loadings.tsv"), stringsAsFactors = FALSE)
    factors <- setdiff(names(value$factors), "biological_sample")
    oriented <- vapply(factors, function(factor) {
        sign <- integrationEnrichmentFactorDirection(manager, value, factor)$sign
        stats::cor(sign * value$factors[[factor]], truth$Z1[match(value$factors$biological_sample, truth$biological_sample)])
    }, numeric(1))
    factor <- factors[[which.max(abs(oriented))]]
    expected_up <- if (oriented[[factor]] > 0) "Up" else "Down"
    expected_down <- setdiff(c("Up", "Down"), expected_up)
    object <- integrationImportedObject(manager, fields = c("features", "design_matrix", "inputs"))
    size_rows <- list()
    for (view in c("proteomics", "metabolomics_GCMS")) {
        omic <- integrationViewOmics(object)[[view]]
        p <- if (omic == "proteomics") list(database = "string", omic = omic, id_type = "uniprot") else
            list(database = "kegg", omic = omic, id_type = "kegg", id_column = "kegg_id")
        ranking <- integrationEnrichmentRank(object, manager, export$source, value, export$parameters, factor, view,
            "pearson", 5L)
        mapped <- integrationEnrichmentMappedList(ranking, p, object@features[[view]])
        values <- integrationEnrichmentTestValues(object, manager, export$source, value, export$parameters, ranking,
            mapped, view)
        # The counts recomputed independently from the whole matrix.
        independent <- intEtIndependentValues(manager, value, view)
        candidates <- mapped$mapping$feature_id[mapped$mapping$status %in% c("mapped", "duplicated")]
        identifier <- stats::setNames(mapped$mapping$identifier, mapped$mapping$feature_id)
        whole <- independent$x[candidates, , drop = FALSE]
        is_complete <- stats::complete.cases(whole)
        expect_identical(values$counts$complete, sum(is_complete), info = view)
        expect_identical(values$counts$incomplete, sum(!is_complete), info = view)
        expect_identical(values$counts$collapsed, sum(is_complete) - length(unique(identifier[candidates[is_complete]])),
            info = view)
        # Sets drawn once, seed 20261184, from complete features on single-feature identifiers.
        rows <- loadings[loadings$view == view, , drop = FALSE]
        kept <- values$feature_id
        id_of <- stats::setNames(rownames(values$y), kept)
        pool <- function(rule) intersect(rows$feature_id[rule], kept)
        size <- if (omic == "proteomics") 20L else 10L
        up_pool <- pool(rows$Z1 > 0)
        down_pool <- pool(rows$Z1 < 0)
        decoy_pool <- pool(rows$Z1 == 0 & rows$Z2 == 0 & (is.na(rows$Z3) | rows$Z3 == 0))
        expect_gte(length(up_pool), 5L)
        expect_gte(length(down_pool), 5L)
        draw <- function(x, k) withr::with_seed(20261184L, sample(x, min(k, length(x))))
        sets <- list(up = id_of[draw(up_pool, size)], down = id_of[draw(down_pool, size)],
            decoy = id_of[draw(decoy_pool, size)])
        term2gene <- data.frame(term = rep(names(sets), lengths(sets)), gene = unlist(sets, use.names = FALSE),
            stringsAsFactors = FALSE)
        term2name <- data.frame(term = names(sets), name = names(sets), stringsAsFactors = FALSE)
        ranked <- mapped$table[c("identifier", "value")]
        camera <- integrationEnrichmentPathwayTest(ranked, term2gene, term2name, 1L, intEtSettings(), "camera",
            values)$terms
        gsea <- integrationEnrichmentPathwayTest(ranked, term2gene, term2name, 20261121L, intEtSettings())$terms
        row <- function(terms, id) terms[terms$ID == id, , drop = FALSE]
        expect_lte(row(camera, "up")$p.adjust, 0.05)
        expect_identical(row(camera, "up")$Direction, expected_up, info = view)
        expect_lte(row(camera, "down")$p.adjust, 0.05)
        expect_identical(row(camera, "down")$Direction, expected_down, info = view)
        expect_gt(row(camera, "decoy")$p.adjust, 0.05)
        expect_lte(row(gsea, "up")$p.adjust, 0.05)
        expect_identical(sign(row(gsea, "up")$NES), if (expected_up == "Up") 1 else -1)
        expect_lte(row(gsea, "down")$p.adjust, 0.05)
        expect_identical(sign(row(gsea, "down")$NES), if (expected_up == "Up") -1 else 1)
        expect_gt(row(gsea, "decoy")$p.adjust, 0.05)
        # Descriptive: 40 invented pathways of 5 to 20 complete identifiers (seed 20261186); calls at BH <= 0.05.
        universe <- rownames(values$y)
        invented <- withr::with_seed(20261186L, lapply(seq_len(40L), function(i) sample(universe, sample(5:20, 1L))))
        names(invented) <- sprintf("P%02d", seq_along(invented))
        t2g <- data.frame(term = rep(names(invented), lengths(invented)), gene = unlist(invented, use.names = FALSE),
            stringsAsFactors = FALSE)
        t2n <- data.frame(term = names(invented), name = names(invented), stringsAsFactors = FALSE)
        called_camera <- integrationEnrichmentPathwayTest(ranked, t2g, t2n, 1L, intEtSettings(), "camera", values)$terms
        called_gsea <- integrationEnrichmentPathwayTest(ranked, t2g, t2n, 20261121L, intEtSettings())$terms
        a <- called_camera$ID[called_camera$p.adjust <= 0.05]
        b <- called_gsea$ID[called_gsea$p.adjust <= 0.05]
        size_rows[[view]] <- data.frame(view = view, factor = factor, pathways = 40L, called_camera = length(a),
            called_gsea = length(b), called_both = length(intersect(a, b)), universe_camera = values$counts$universe,
            ranked_gsea = nrow(ranked), stringsAsFactors = FALSE)
    }
    size <- do.call(rbind, size_rows)
    cat("\nKIT-183 T13 size of the change (pathways called at BH <= 0.05):\n")
    utils::write.table(size, stdout(), sep = "\t", quote = FALSE, row.names = FALSE)
    out <- Sys.getenv("MULTISCHOLAR_K183_OUT", unset = "")
    if (nzchar(out)) {
        dir.create(out, recursive = TRUE, showWarnings = FALSE)
        utils::write.table(size, file.path(out, "t13-size-of-change.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
    }
    expect_identical(nrow(size), 2L)
})

# --- T5, T8, T15: guards against the base (the store made at base, kit183-base-store) ---

# The base store is built by the base commit's code inside the runner job that runs these tests
# (docs/qa/workflow-execution/kit/kit-183/base-store-job.R, VAL-64 (2)), which names it in K183_BASE_STORE; it is never
# kept in the repository. Without it T5 and T8 skip.
intEtBaseStoreCopy <- function(.local_envir = parent.frame()) {
    source_store <- Sys.getenv("K183_BASE_STORE", unset = "")
    testthat::skip_if_not(nzchar(source_store) && dir.exists(file.path(source_store, ".workflow-entities")),
        "the KIT-183 base store is built only inside base-store-job.R's runner job")
    testthat::skip_if_not(workflowOwnerCapable(), "the integration store needs qs2, DuckDB and Arrow")
    root <- withr::local_tempdir(.local_envir = .local_envir)
    stopifnot(file.copy(source_store, root, recursive = TRUE))
    root <- file.path(root, basename(source_store))
    stopifnot(dir.exists(file.path(root, ".workflow-entities")))
    identity <- readSavedStoreIdentity(file.path(root, ".workflow-entities", "store.json"))
    manager <- newWorkflowStateOwner(newWorkflowEntityStore(root, identity$project_id, identity$workflow_id), reopen = TRUE)
    withr::defer(closeWorkflowStateOwner(manager), envir = .local_envir)
    manager
}

test_that("KIT-183 T5: a KEGG entry made at base replays with the pre-ranked test, terms identical", {
    skip_if_not_installed("clusterProfiler")
    manager <- intEtBaseStoreCopy()
    recorded <- integrationRecordedEnrichment(manager)
    expect_identical(names(recorded), "Factor1|metabolomics_GCMS|kegg")
    expect_null(recorded[[1L]]$entry$settings$test)
    comparison <- integrationReplayEnrichment(manager, recorded)
    expect_true(comparison$terms_identical)
    expect_true(comparison$ranking_equal)
    entry <- manager$.context$snapshot$exports$integration_mofa_model$enrichment[["Factor1|metabolomics_GCMS|kegg"]]
    expect_identical(entry$settings$test, "gsea")
})

# T8 restated (AUDIT-KIT-183-impl I03, VAL-85 (3)): the pre-ranked test's identity is unchanged by this ticket apart
# from the build. A cache hit across builds is not expected and not tested: every saved product's identity includes the
# build's digest (R/utils_workflow_dia_products.R:273-277, R/utils_workflow_dia_state.R:1655-1665).
test_that("KIT-183 T8: (a) the ranking and pre-ranked declarations are those of base", {
    ranking <- integrationEnrichmentRankingOperation()
    enrichment <- integrationEnrichmentOperation()
    expect_identical(ranking[c("name", "input_class", "fields", "unit", "product_type", "factors", "parameters")],
        list(name = "integration_enrichment_ranking", input_class = "IntegrationViewsData",
            fields = c("features", "design_matrix", "inputs"), unit = "whole_matrix",
            product_type = "integration.enrichment_ranking",
            factors = list(product_schema = "multischolar.integration_enrichment_ranking.v1"),
            parameters = c("model", "factor", "view", "correlation_method", "min_samples", "direction_rule")))
    expect_identical(enrichment[c("name", "input_class", "fields", "unit", "product_type", "factors", "parameters")],
        list(name = "integration_enrichment", input_class = "IntegrationViewsData", fields = c("features", "inputs"),
            unit = "whole_matrix", product_type = "integration.enrichment",
            factors = list(product_schema = "multischolar.integration_enrichment.v1"),
            parameters = c("ranking", "database", "id_column", "service_response", "mapping", "organism", "string_species",
                "string_fdr", "string_rank_direction", "min_set_size", "max_set_size", "seed")))
})

test_that("KIT-183 T8: (b) the base store's ranking and gsea identities equal the head's apart from the build; (c) a hit in one build", {
    skip_if_not_installed("clusterProfiler")
    manager <- intEtBaseStoreCopy()
    export <- manager$.context$snapshot$exports$integration_mofa_model
    recorded <- integrationRecordedEnrichment(manager)[[1L]]
    entry <- recorded$entry
    s <- entry$settings
    keys <- c("input", "design", "parameters", "selectors")
    identityAtHead <- function(declaration, parameters) withWorkflowStateScope(manager, function()
        workflowComponentScientificIdentity(manager, export$source, declaration$name, parameters[declaration$parameters]))
    recordedIdentity <- function(entity_id) withWorkflowStateScope(manager, function()
        describeSavedWorkflowEntity(manager$.context$store, entity_id)$identity)
    # (b) the ranking.
    ranking <- integrationEnrichmentRankingOperation()
    head_ranking <- identityAtHead(ranking, list(model = export$reference$entity_id, factor = s$factor, view = s$view,
        correlation_method = s$correlation_method, min_samples = s$min_samples,
        direction_rule = INTEGRATION_ENRICHMENT_DIRECTION_RULE))
    base_ranking <- recordedIdentity(entry$ranking$entity_id)
    expect_identical(head_ranking[keys], base_ranking[keys])
    expect_false(identical(head_ranking$software, base_ranking$software))
    # (b) the KEGG gsea product, with the enrichment's `ranking` set to the recorded ranking's entity id.
    enrichment <- integrationEnrichmentOperation()
    head_gsea <- identityAtHead(enrichment, list(ranking = entry$ranking$entity_id, database = s$database,
        id_column = s$id_column, service_response = "none", mapping = entry$stored$sha256, organism = s$kegg_organism,
        string_species = s$string_species, string_fdr = s$string_fdr, string_rank_direction = s$string_rank_direction,
        min_set_size = s$min_set_size, max_set_size = s$max_set_size, seed = s$seed, test = "gsea"))
    base_gsea <- recordedIdentity(entry$reference$entity_id)
    expect_identical(head_gsea[keys], base_gsea[keys])
    expect_false(identical(head_gsea$software, base_gsea$software))
    # (c) same build: the recorded request run twice at head; the second run is a hit with the first run's entity id.
    request <- list(factor = s$factor, view = s$view, database = s$database, id_column = s$id_column, test = "gsea")
    mapping <- c(lapply(entry$mapping, identity), list(table = recorded$stored))
    first <- runIntegrationEnrichment(manager, request, mapping = mapping)
    second <- runIntegrationEnrichment(manager, request, mapping = mapping)
    expect_identical(second$run$status, "hit")
    expect_identical(second$run$reference$entity_id, first$run$reference$entity_id)
    expect_identical(second$ranking$status, "hit")
})

test_that("KIT-183 T15 (KIT-191 D9 (a)): integrationEnrichmentGsea: formals and value as at base, body but the RNG kinds", {
    skip_if_not_installed("clusterProfiler")
    path <- file.path(intEtBaseStore(), "gsea-deparse.txt")
    testthat::skip_if_not(file.exists(path), "the KIT-183 base store is not present")
    gsea <- utils::removeSource(integrationEnrichmentGsea)
    # The body differs from base only by the three RNG-kind arguments of the seeded call.
    code <- body(gsea)
    at <- which(vapply(as.list(code), function(e) startsWith(paste(deparse(e), collapse = ""),
        "outcome <- withr::with_seed("), logical(1)))
    expect_length(at, 1L)
    call <- code[[at]][[3L]]
    kinds <- c(".rng_kind", ".rng_normal_kind", ".rng_sample_kind")
    expect_identical(lapply(as.list(call)[kinds], identity), list(.rng_kind = "Mersenne-Twister",
        .rng_normal_kind = "Inversion", .rng_sample_kind = "Rejection"))
    parts <- as.list(call)
    code[[at]][[3L]] <- as.call(parts[!names(parts) %in% kinds])
    expect_identical(c(deparse(formals(gsea)), deparse(code)), readLines(path))
    fixture <- intEtGseaFixture()
    value <- integrationEnrichmentGsea(fixture$ranked, fixture$term2gene, fixture$term2name, seed = 11L)
    expect_identical(value, readRDS(file.path(intEtBaseStore(), "gsea-fixture-value.rds")))
})

# --- T11: the tab; T12: the workbook ---

intEtKeggDownload <- function(sets) {
    links <- paste(unlist(lapply(seq_along(sets), function(i) sprintf("path:map%05d\tcpd:%s", i, sets[[i]]))), collapse = "\n")
    names <- paste(sprintf("hsa%05d\t%s", seq_along(sets), names(sets)), collapse = "\n")
    replies <- list("https://rest.kegg.jp/info/kegg" = "kegg   Kyoto Encyclopedia\n                 Release fixture, Nov 25\n",
        "https://rest.kegg.jp/link/compound/pathway" = paste0(links, "\n"),
        "https://rest.kegg.jp/list/pathway/hsa" = paste0(names, "\n"))
    function(url, timeout_seconds = 300) {
        if (is.null(replies[[url]])) stop("no network in this test: ", url)
        charToRaw(replies[[url]])
    }
}

intEtCorpusSets <- function(manager, value, view = "metabolomics_GCMS") {
    independent <- intEtIndependentValues(manager, value, view)
    ids <- stats::setNames(independent$features$kegg_id, independent$features$feature_id)
    complete <- rownames(independent$x)[stats::complete.cases(independent$x)]
    pool <- intersect(names(ids)[!duplicated(ids) & !duplicated(ids, fromLast = TRUE)], complete)
    list(first = unname(ids[pool[1:10]]), second = unname(ids[pool[11:20]]), third = unname(ids[pool[21:30]]))
}

test_that("KIT-183 T11: the tab: the test selector, its default, notes, counts, caution, seed label and key", {
    skip_if_not_installed("clusterProfiler")
    paths <- intEtProject()
    fixture <- intEtModel(paths)
    manager <- fixture$manager
    texts <- integrationGuiTexts()
    escaped <- function(text) as.character(htmltools::htmlEscape(text))
    withr::local_options(MultiScholaR.integration_gui = TRUE)
    html <- as.character(htmltools::renderTags(mod_integration_enrichment_ui("integration_enrichment"))$html)
    expect_match(html, escaped(texts[["enrich_seed"]]), fixed = TRUE)
    expect_identical(texts[["enrich_seed"]], "Seed of the pre-ranked and permutation tests (the correlation-adjusted tests use none)")
    expect_match(html, 'data-testid="integration-enrichment-test"', fixed = TRUE)
    sets <- intEtCorpusSets(manager, fixture$run$value)
    bus <- do.call(shiny::reactiveValues, omicWorkflowInitialBus(integrationGuiConfig()$shell))
    bus$state_manager <- manager
    bus$integration_model <- list(reference = fixture$run$reference, status = fixture$run$status)
    no_progress <- function(message, value, expr) force(expr)
    shiny::testServer(mod_integration_enrichment_server, args = list(workflow_data = bus, paths = paths,
        download = intEtKeggDownload(sets), withProgressFn = no_progress), {
        session$flushReact()
        session$setInputs(factor = "Factor1", view = "metabolomics_GCMS", database = "kegg")
        selector <- as.character(output$test_ui$html)
        expect_match(selector, sprintf('id="%s"', session$ns("test")), fixed = TRUE)
        expect_match(selector, '<option value="camera" selected>', fixed = TRUE)
        expect_match(selector, escaped(texts[["enrich_test_camera"]]), fixed = TRUE)
        expect_match(selector, escaped(texts[["enrich_test_gsea"]]), fixed = TRUE)
        expect_false(grepl("camera_pr", selector, fixed = TRUE))
        session$setInputs(test = "camera")
        expect_match(as.character(output$test_note$html), escaped(texts[["enrich_test_camera_note_metabolites"]]), fixed = TRUE)
        session$setInputs(test = "gsea")
        expect_match(as.character(output$test_note$html), escaped(texts[["enrich_test_gsea_note_metabolites"]]), fixed = TRUE)
        expect_false(grepl(escaped(texts[["enrich_test_camera_note_metabolites"]]), as.character(output$test_note$html),
            fixed = TRUE))
        session$setInputs(database = "reactome")
        expect_match(as.character(output$test_ui$html), sprintf('id="%s"', session$ns("test")), fixed = TRUE)
        session$setInputs(database = "kegg", test = "camera", id_column = "kegg_id", correlation_method = "pearson",
            min_samples = 5L, seed = 20261121L, string_species = 9606L, string_fdr = 0.05, kegg_organism = "hsa",
            reactome_organism = "Homo sapiens", min_set_size = 5L, max_set_size = 500L)
        session$setInputs(download = 1)
        session$setInputs(run = 1)
        expect_identical(bus$integration_enrichment$entry, "Factor1|metabolomics_GCMS|kegg|camera")
        result <- as.character(output$result$html)
        value <- state$run$run$value
        counts <- value$test_counts
        expect_match(result, escaped(sprintf(texts[["enrich_test_counts"]], counts$universe, counts$samples,
            counts$incomplete, counts$collapsed)), fixed = TRUE)
        expect_match(result, escaped(texts[["enrich_caution"]]), fixed = TRUE)
        expect_match(result, escaped(texts[["enrich_direction_camera"]]), fixed = TRUE)
        strength <- integrationEnrichmentFactorStrength(fixture$run$value, "Factor1", "metabolomics_GCMS")
        expect_match(result, escaped(sprintf(texts[["enrich_factor_strength"]],
            integrationEnrichmentPercent(strength$factor_r2_percent),
            integrationEnrichmentPercent(strength$largest_r2_percent), "metabolite data",
            texts[["enrich_strength_terms_metabolites"]], "6.7 %", "correlation-adjusted test")), fixed = TRUE)
        expect_false(grepl("Terms at the top of the list", result, fixed = TRUE))
        # STRING: no selector, no note.
        session$setInputs(view = "proteomics", database = "string")
        hidden <- paste(as.character(output$test_ui$html), collapse = "")
        expect_false(grepl(sprintf('id="%s"', session$ns("test")), hidden, fixed = TRUE))
    })
    # The summary of a STRING result carries neither the caution nor the factor-strength sentence.
    string <- list(database = "string", factor = "Factor1", view = "proteomics", service = list(version = list(string_version = "12.5")),
        ranking_counts = list(features = 1L, ranked = 1L, too_few_samples = 0L, no_variation = 0L, internal_standard = 0L),
        mapping_counts = list(mapped = 1L, unmapped = 0L, duplicated = 0L),
        direction = list(feature_id = "P1", view = "proteomics", sign = 1L), tests = 1L)
    summary <- integrationEnrichmentSummary(string, texts, list(database = "string", factor_r2_percent = 10,
        largest_r2_percent = 20))
    expect_false(texts[["enrich_caution"]] %in% summary)
    expect_false(any(grepl("this factor explains", summary, fixed = TRUE)))
    # The notes: the few-samples sentence below 11 samples. KIT-184: the RUV note for a model trained on RUV-corrected
    # values stays with the tests that have no design (the pre-ranked and permutation tests); camera has the setting
    # "Adjust for removed RUV factors" and its own sentences with each result.
    notes <- integrationEnrichmentTestNotes("camera", 8L, "ruv_corrected", "metabolomics", texts)
    expect_identical(notes, c(texts[["enrich_test_camera_note_metabolites"]],
        sprintf(texts[["enrich_test_few_samples_camera_metabolites"]], 8L)))
    expect_identical(integrationEnrichmentTestNotes("gsea", 8L, "ruv_corrected", "metabolomics", texts),
        c(texts[["enrich_test_gsea_note_metabolites"]], sprintf(texts[["enrich_test_few_samples_gsea_metabolites"]], 8L),
            texts[["enrich_values_ruv_note"]]))
    expect_identical(integrationEnrichmentTestNotes("gsea", 12L, "uncorrected", "metabolomics", texts),
        texts[["enrich_test_gsea_note_metabolites"]])
    expect_identical(integrationEnrichmentPercent(0.04), "less than 0.1 %")
    expect_identical(integrationEnrichmentPercent(12.345), "12.3 %")
})

test_that("KIT-183 T12: the workbook: the gsea override gives base's terms; the default is the test of record", {
    skip_if_not_installed("clusterProfiler")
    base_terms <- file.path(intEtBaseStore(), "base-terms.tsv")
    paths <- intEtProject()
    fixture <- intEtModel(paths)
    manager <- fixture$manager
    sets <- intEtCorpusSets(manager, fixture$run$value)
    mapping <- intEtMapping(sets)
    mapping$database <- "KEGG"
    local_mocked_bindings(
        integrationDownloadKeggMapping = function(organism, download = NULL) mapping,
        integrationDownloadReactomeMapping = function(download = NULL) list(database = "Reactome", release = "Reactome 0",
            sha256 = strrep("f", 64L), url = "fixture", downloaded_at = "x", table = data.frame(chebi_id = "1",
                reactome_id = "R-HSA-1", reactome_term = "x", organism = "Homo sapiens", stringsAsFactors = FALSE)))
    project <- list(manager = manager, store_root = dirname(manager$.context$store$directory), paths = paths,
        views = c("metabolomics_GCMS", "metabolomics_LCMS_Pos"))
    results <- list(results = list(views = "metabolomics_GCMS"), factor = "Factor1")
    default <- integrationEnrichFactor(project, results, id_columns = list(metabolomics_GCMS = "kegg_id"))
    entry <- manager$.context$snapshot$exports$integration_mofa_model$enrichment[["Factor1|metabolomics_GCMS|kegg|camera"]]
    expect_identical(entry$settings$test, "camera")
    expect_true("camera" %in% default$enrichments[["metabolomics_GCMS kegg"]]$terms$test)
    earlier <- integrationEnrichFactor(project, results, id_columns = list(metabolomics_GCMS = "kegg_id"),
        overrides = list(test = "gsea"))
    gsea <- earlier$enrichments[["metabolomics_GCMS kegg"]]
    texts <- integrationGuiTexts()
    expect_identical(sum(gsea$summary == texts[["enrich_caution"]]), 1L)
    expect_identical(sum(grepl("this factor explains", gsea$summary, fixed = TRUE)), 1L)
    testthat::skip_if_not(file.exists(base_terms), "the KIT-183 base store is not present")
    base <- utils::read.delim(base_terms, stringsAsFactors = FALSE, check.names = FALSE)
    head_terms <- gsea$terms
    head_terms[] <- lapply(head_terms, function(x) if (is.numeric(x)) signif(x, 12) else as.character(x))
    base[] <- lapply(base, function(x) if (is.numeric(x)) signif(x, 12) else as.character(x))
    expect_identical(names(head_terms), names(base))
    expect_equal(head_terms, base, tolerance = 1e-10, ignore_attr = TRUE)
    base_summary <- readLines(file.path(intEtBaseStore(), "base-summary.txt"))
    added <- c(texts[["enrich_caution"]], gsea$summary[grepl("this factor explains", gsea$summary, fixed = TRUE)])
    expect_identical(setdiff(gsea$summary, added), base_summary)
})

# --- The cross-platform record (ticket 7): T1's fixture with camera and gsea; the terms printed at full precision (the
# reference, from a Linux run-tests job, is tests/testdata/integration/kit183-cross-platform.tsv), and on another
# platform compared by the rule: SHA-256 equal, or the same IDs and directions and every number within 1e-6 ---

intEtPlatformTerms <- function() {
    n <- 12L
    g <- 40L
    f <- c(-1.5, -1.3, -1.0, -0.7, -0.4, -0.2, 0.1, 0.3, 0.6, 0.9, 1.2, 1.6)
    mu <- withr::with_seed(11L, stats::runif(g, 14, 22))
    sd <- 0.2 + 0.6 * exp(-(mu - 14) / 2)
    y <- mu + sd * withr::with_seed(12L, matrix(stats::rnorm(g * n), g, n))
    y[1:8, ] <- y[1:8, ] + 0.8 * matrix(f, 8, n, byrow = TRUE)
    y[9:16, ] <- y[9:16, ] - 0.6 * matrix(f, 8, n, byrow = TRUE)
    ids <- sprintf("C%05d", seq_len(g))
    rownames(y) <- ids
    sets <- list(down = ids[9:16], mixed = ids[c(1:4, 20:25)], null = ids[26:35], up = ids[1:8])
    term2gene <- data.frame(term = rep(names(sets), lengths(sets)), gene = unlist(sets, use.names = FALSE))
    term2name <- data.frame(term = names(sets), name = names(sets))
    values <- list(y = y, design = integrationEnrichmentTestDesign(f), counts = list(), design_covariates = "none")
    ranked <- data.frame(identifier = ids, value = stats::cor(t(y), f)[, 1L])
    camera <- integrationEnrichmentPathwayTest(ranked, term2gene, term2name, 1L, intEtSettings(), "camera", values)$terms
    gsea <- integrationEnrichmentPathwayTest(ranked, term2gene, term2name, 20261121L, intEtSettings())$terms
    rbind(data.frame(test = "camera", ID = camera$ID, direction = camera$Direction, pvalue = camera$pvalue,
            p.adjust = camera$p.adjust, correlation = camera$correlation, effect = camera$effect, NES = NA_real_),
        data.frame(test = "gsea", ID = gsea$ID, direction = ifelse(gsea$NES > 0, "Up", "Down"), pvalue = gsea$pvalue,
            p.adjust = gsea$p.adjust, correlation = NA_real_, effect = NA_real_, NES = gsea$NES))
}

test_that("KIT-183 cross-platform: camera and gsea on T1's fixture agree with the Linux reference within 1e-6", {
    skip_if_not_installed("clusterProfiler")
    terms <- intEtPlatformTerms()
    terms <- terms[order(terms$test, terms$ID), ]
    rownames(terms) <- NULL
    text <- utils::capture.output(utils::write.table(format(terms, digits = 17), sep = "\t", quote = FALSE,
        row.names = FALSE))
    cat("\nKIT-183 cross-platform terms (sha256 ", digest::digest(text, algo = "sha256"), "):\n", sep = "")
    cat(text, sep = "\n")
    reference <- file.path(intEtData(), "kit183-cross-platform.tsv")
    testthat::skip_if_not(file.exists(reference), "the cross-platform reference is not recorded yet")
    ref <- utils::read.delim(reference, stringsAsFactors = FALSE)
    ref <- ref[order(ref$test, ref$ID), ]
    expect_identical(terms$ID, ref$ID)
    expect_identical(terms$direction, ref$direction)
    for (column in c("pvalue", "p.adjust", "correlation", "effect", "NES")) {
        difference <- max(abs(terms[[column]] - ref[[column]]), na.rm = TRUE)
        expect_lte(difference, 1e-6, label = sprintf("largest absolute difference of %s (%g)", column, difference))
    }
})
