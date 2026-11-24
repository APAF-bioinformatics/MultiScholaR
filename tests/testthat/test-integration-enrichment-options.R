# KIT-191 (VAL-85 (2), VAL-97): the integration enrichment's tests offered as a list of options with their measured
# trade-offs: `camera` (the default, unchanged), `camera_floor` (negative estimated correlations counted as zero),
# `permutation` (PCGSE's permutation of the factor's values as KIT-183's simulation measured it) and `gsea`. Declared
# before the fix (ticket section 10; each test's failing value on the tests-only commit is that the behaviour is absent):
#   U1  camera_floor equals limma::camera(allow.neg.cor = FALSE); p_floor >= p_camera, equal where the estimate >= 0,
#       larger where it is negative; with a covariate column the same;
#   U2  permutation equals a reference written as loops from 5.2's formulas (same reorderings); a 4-sample fixture whose
#       reorderings hold the identity: each counted, p >= (1 + d) / (1 + B), d recorded;
#   U4  the permutation's (and, D9 (a), gsea's) result does not depend on the session's RNGkind(); .Random.seed kept;
#   U5  seeds: different seeds, different products; same seed, a hit; camera ignores the seed;
#   U6  resolver: four tests, camera_pr refused, unknown refused naming the four, STRING ignores the test, the sample
#       floor, permutations 1000 or 10000;
#   U7  (a) camera's canonical text and the declarations' parameters as at base;
#   U8  keys, entries (seed, permutations), replay of camera_floor and permutation entries; a permutation reproduced in
#       a second store;
#   U9  record: evidence sheet, Methods per test (no format code left), limits per test and kind, strength per test and
#       kind, several tests on one factor, no sheet without a pathway entry;
#   U10 the tab: four choices, camera pre-selected, notes and evidence of the kind, the lipid line, the seed and the
#       reorderings shown only for their tests, few-samples sentences, the floor sentence, the permutation's counts;
#   U11 the workbook: permutation and camera_floor overrides run and record; STRING's request carries no test;
#   U12 the cost guard: the prediction, refusal above the tab's limit, a message in scripts, L1's coefficients;
#   U13 the corpus control for camera_floor and permutation.
# U3, U7 (b), U8's branch-point store and U14 run as scripts in docs/qa/workflow-execution/kit/kit-191/.

# --- Fixtures (KIT-183's, test-integration-enrichment-tests.R) ---

intEoData <- function() normalizePath(testthat::test_path("..", "testdata", "integration"), mustWork = FALSE)
intEoStored <- function(name) utils::read.delim(file.path(intEoData(), "mofa_reference", paste0(name, ".tsv")),
    check.names = FALSE, stringsAsFactors = FALSE)
intEoSettings <- function(min_size = 5L, max_size = 500L, ...) {
    c(integrationGseaSettings(p_adjust_method = "BH", min_size = min_size, max_size = max_size, p_cutoff = 1,
        tie_jitter = 0), list(...))
}

intEoTrainer <- function() {
    function(views, options) {
        factors <- intEoStored("factors")
        f <- as.matrix(factors[-1L])
        rownames(f) <- factors$biological_sample
        f <- f[colnames(views[[1L]]), , drop = FALSE]
        w <- intEoStored("weights")
        weights <- lapply(stats::setNames(names(views), names(views)), function(view) {
            rows <- w[w$view == view, , drop = FALSE]
            m <- as.matrix(rows[-(1:2)])
            rownames(m) <- paste(view, rows$feature_id, sep = ":")
            m[rownames(views[[view]]), , drop = FALSE]
        })
        v <- intEoStored("variance_explained")
        r2 <- matrix(v$r2_percent, nrow = length(unique(v$factor)), dimnames = list(unique(v$factor), unique(v$view)))
        r2 <- r2[, names(views), drop = FALSE]
        list(factors = f, weights = weights, r2_per_factor = r2, r2_total = colSums(r2), elbo = c(-2, -1), iterations = 186L,
            options = options, versions = list(MOFA2 = "stored reference"), messages = character(),
            model_file = list(bytes = as.raw(1:3), sha256 = strrep("d", 64L), size = 3))
    }
}

intEoKeggIds <- function(n) {
    ids <- sprintf("C%05d", 10000L + seq_len(n))
    pairs <- seq(2L, 20L, by = 2L)
    ids[pairs] <- ids[pairs - 1L]
    ids
}

intEoProject <- function(.local_envir = parent.frame()) {
    testthat::skip_if_not(dir.exists(file.path(intEoData(), "corpus")), "the tier-1 corpus is not present")
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
        stopifnot(all(file.copy(list.files(file.path(intEoData(), "corpus", omic), full.names = TRUE), dir)))
    }
    dir <- file.path(paths$integration_dir, "metabolomics")
    input <- readIntegrationInputWithMatrices(dir)
    read <- input$manifest
    manifest <- list(omic_type = read$omic_type, workflow = read$workflow, matrix_role = read$matrix_role, scale = read$scale,
        correction = NULL, source_state = list(name = "kit183-kegg"), sample_id_column = read$sample_id_column,
        group_column = read$group_column, technical_replicate_column = read$technical_replicate_column,
        biological_sample_column = read$biological_sample_column)
    features <- lapply(input$views, `[[`, "features")
    features$metabolomics_GCMS$kegg_id <- intEoKeggIds(nrow(features$metabolomics_GCMS))
    writeIntegrationInput(dir, manifest, input$samples, lapply(input$views, `[[`, "matrix"), features)
    paths
}

intEoModel <- function(paths, .local_envir = parent.frame()) {
    manager <- newIntegrationStateOwner(file.path(paths$source_dir, paste0("store-", basename(tempfile()))))
    withr::defer(closeWorkflowStateOwner(manager), envir = .local_envir)
    scan <- integrationScanProject(paths$integration_dir)
    integrationCommitImport(manager, integrationImportRequest(scan, scan$views$view))
    run <- runIntegrationModel(manager, list(), trainer = intEoTrainer())
    list(manager = manager, run = run)
}

intEoMapping <- function(sets) {
    pathway <- paste0("path:hsa", sprintf("%05d", seq_along(sets)))
    list(database = "KEGG", release = "KEGG Release fixture", sha256 = strrep("e", 64L), url = "fixture",
        downloaded_at = "2026-11-24T00:00:00+0000",
        table = data.frame(pathway = rep(pathway, lengths(sets)), compound = paste0("cpd:", unlist(sets, use.names = FALSE)),
            name = rep(names(sets), lengths(sets)), stringsAsFactors = FALSE))
}

intEoIndependentValues <- function(manager, value, view) {
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

intEoCorpusSets <- function(manager, value, view = "metabolomics_GCMS") {
    independent <- intEoIndependentValues(manager, value, view)
    ids <- stats::setNames(independent$features$kegg_id, independent$features$feature_id)
    complete <- rownames(independent$x)[stats::complete.cases(independent$x)]
    pool <- intersect(names(ids)[!duplicated(ids) & !duplicated(ids, fromLast = TRUE)], complete)
    list(first = unname(ids[pool[1:10]]), second = unname(ids[pool[11:20]]), third = unname(ids[pool[21:30]]))
}

intEoKeggDownload <- function(sets) {
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

# KIT-183 T1's matrix fixture (40 x 12), with features 36 to 40 given a common residual of opposite signs (36 to 38
# against 39 and 40), so that the set of those five has a negative estimated correlation.
intEoMatrix <- function(covariate = FALSE) {
    n <- 12L
    g <- 40L
    f <- c(-1.5, -1.3, -1.0, -0.7, -0.4, -0.2, 0.1, 0.3, 0.6, 0.9, 1.2, 1.6)
    mu <- withr::with_seed(11L, stats::runif(g, 14, 22))
    sd <- 0.2 + 0.6 * exp(-(mu - 14) / 2)
    y <- mu + sd * withr::with_seed(12L, matrix(stats::rnorm(g * n), g, n))
    y[1:8, ] <- y[1:8, ] + 0.8 * matrix(f, 8, n, byrow = TRUE)
    y[9:16, ] <- y[9:16, ] - 0.6 * matrix(f, 8, n, byrow = TRUE)
    u <- as.numeric(scale(stats::residuals(stats::lm(withr::with_seed(31L, stats::rnorm(n)) ~ f))))
    small <- 0.05 * withr::with_seed(32L, matrix(stats::rnorm(5L * n), 5L, n))
    y[36:40, ] <- 17 + 0.5 * matrix(u, 5L, n, byrow = TRUE) * c(1, 1, 1, -1, -1) + small
    ids <- sprintf("C%05d", seq_len(g))
    rownames(y) <- ids
    sets <- list(down = ids[9:16], mixed = ids[c(1:4, 20:25)], negative = ids[36:40], up = ids[1:8])
    w <- if (covariate) matrix(withr::with_seed(13L, stats::rnorm(n)), ncol = 1L, dimnames = list(NULL, "ruv_W1"))
    design <- integrationEnrichmentTestDesign(f, w)
    list(y = y, f = f, ids = ids, sets = sets, design = design,
        values = list(y = y, f = f, design = design, counts = list(universe = g, samples = n), design_covariates = "none"),
        term2gene = data.frame(term = rep(names(sets), lengths(sets)), gene = unlist(sets, use.names = FALSE),
            stringsAsFactors = FALSE),
        term2name = data.frame(term = names(sets), name = paste(names(sets), "set"), stringsAsFactors = FALSE),
        ranked = data.frame(identifier = ids, value = stats::cor(t(y), f)[, 1L], stringsAsFactors = FALSE))
}

intEoTest <- function(fixture, test, settings = intEoSettings(), seed = 7L) {
    suppressMessages(integrationEnrichmentPathwayTest(fixture$ranked, fixture$term2gene, fixture$term2name, seed, settings,
        test = test, values = fixture$values))
}

# A data set of `n` samples and `g` features in memory, with a factor and five sets; for the permutation's reference.
intEoPermutationFixture <- function(n = 12L, g = 60L) {
    f <- withr::with_seed(41L, stats::rnorm(n))
    y <- 10 + withr::with_seed(42L, matrix(stats::rnorm(g * n), g, n))
    y[1:6, ] <- y[1:6, ] + 0.9 * matrix(f, 6L, n, byrow = TRUE)
    y[7:12, ] <- y[7:12, ] - 0.7 * matrix(f, 6L, n, byrow = TRUE)
    ids <- sprintf("C%05d", seq_len(g))
    rownames(y) <- ids
    sets <- list(a_up = ids[1:6], b_down = ids[7:12], c_mixed = ids[c(1:3, 30:36)], d_null = ids[40:54], e_tail = ids[55:60])
    list(y = y, f = f, ids = ids, sets = sets,
        values = list(y = y, f = f, design = integrationEnrichmentTestDesign(f), counts = list(universe = g, samples = n),
            design_covariates = "none"),
        term2gene = data.frame(term = rep(names(sets), lengths(sets)), gene = unlist(sets, use.names = FALSE),
            stringsAsFactors = FALSE),
        term2name = data.frame(term = names(sets), name = names(sets), stringsAsFactors = FALSE),
        ranked = data.frame(identifier = ids, value = 0, stringsAsFactors = FALSE))
}

# 5.2 written as loops: per permutation the factor reordered and each feature's correlation recomputed from the
# row-standardised values; per set the two-sample t from the means and sums of squares in and out; a reordering that
# leaves the factor unchanged counted without comparing; |t*| within 1e-9 below |t| counted (AUDIT-KIT-193 U2).
intEoPermutationReference <- function(y, f, sets, seed, permutations) {
    n <- ncol(y)
    fisher <- function(r) atanh(pmax(pmin(r, 1 - 1e-12), -1 + 1e-12))
    setT <- function(z, members) {
        inside <- z[members]
        outside <- z[-members]
        m1 <- length(inside)
        m2 <- length(outside)
        pooled <- sqrt((sum((inside - mean(inside))^2) + sum((outside - mean(outside))^2)) / (m1 + m2 - 2))
        c(t = (mean(inside) - mean(outside)) / (pooled * sqrt(1 / m1 + 1 / m2)), effect = mean(inside) - mean(outside))
    }
    members <- lapply(sets, function(set) match(set, rownames(y)))
    z <- fisher(as.numeric(stats::cor(t(y), f)))
    observed <- t(vapply(members, function(m) setT(z, m), numeric(2)))
    perms <- withr::with_seed(seed, vapply(seq_len(permutations), function(i) sample(n), integer(n)),
        .rng_kind = "Mersenne-Twister", .rng_normal_kind = "Inversion", .rng_sample_kind = "Rejection")
    ys <- t(scale(t(y)))
    count <- numeric(length(sets))
    identity_draws <- 0L
    for (b in seq_len(permutations)) {
        fb <- f[perms[, b]]
        if (identical(fb, f)) {
            identity_draws <- identity_draws + 1L
            count <- count + 1
            next
        }
        fs <- as.numeric(scale(fb))
        zb <- fisher(vapply(seq_len(nrow(y)), function(i) sum(ys[i, ] * fs) / (n - 1), numeric(1)))
        for (k in seq_along(members)) {
            if (abs(setT(zb, members[[k]])[["t"]]) >= abs(observed[k, "t"]) - 1e-9) count[[k]] <- count[[k]] + 1
        }
    }
    pvalue <- (1 + count) / (1 + permutations)
    data.frame(ID = names(sets), t = observed[, "t"], effect = observed[, "effect"], count = count, pvalue = pvalue,
        Direction = ifelse(observed[, "t"] > 0, "Up", "Down"), identity_draws = identity_draws, stringsAsFactors = FALSE)
}

# --- U1, U2: exactness ---

test_that("KIT-191 U1: camera_floor equals limma camera with allow.neg.cor = FALSE; never below camera's p-value", {
    fixture <- intEoMatrix()
    index <- split(fixture$term2gene$gene, fixture$term2gene$term)
    camera <- intEoTest(fixture, "camera")$terms
    camera <- camera[match(names(index), camera$ID), , drop = FALSE]
    # The guard on the fixture: the negative set's estimated correlation is below 0.
    expect_lt(camera$correlation[camera$ID == "negative"], 0)
    floor <- intEoTest(fixture, "camera_floor")
    expect_identical(floor$method, "camera_floor")
    terms <- floor$terms[match(names(index), floor$terms$ID), , drop = FALSE]
    direct <- limma::camera(fixture$y, index, fixture$design, contrast = which(colnames(fixture$design) == "factor"),
        inter.gene.cor = NA, allow.neg.cor = FALSE, trend.var = TRUE, use.ranks = FALSE, sort = FALSE)
    expect_identical(terms$ID, rownames(direct))
    expect_identical(terms$setSize, as.integer(direct$NGenes))
    expect_identical(terms$Direction, as.character(direct$Direction))
    expect_identical(terms$correlation, as.numeric(direct$Correlation))
    expect_identical(terms$pvalue, as.numeric(direct$PValue))
    expect_identical(unique(terms$test), "camera_floor")
    expect_equal(floor$terms$p.adjust, stats::p.adjust(floor$terms$pvalue, "BH"))
    expect_identical(terms$effect, camera$effect)
    expect_identical(floor$settings$test_settings, INTEGRATION_ENRICHMENT_TEST_SETTINGS[["camera_floor"]])
    # p_floor >= p_camera on every set; equal where the estimate is >= 0, larger where it is negative.
    nonnegative <- camera$correlation >= 0
    expect_true(all(terms$pvalue >= camera$pvalue))
    expect_identical(terms$pvalue[nonnegative], camera$pvalue[nonnegative])
    expect_true(all(terms$pvalue[!nonnegative] > camera$pvalue[!nonnegative]))
    expect_true("negative" %in% camera$ID[!nonnegative])
    # The internal hook keeps its meaning for camera only: the floor is set by the test.
    hooked <- intEoTest(fixture, "camera_floor", intEoSettings(allow_neg_cor = TRUE))$terms
    expect_identical(hooked, floor$terms)
    expect_identical(intEoTest(fixture, "camera", intEoSettings(allow_neg_cor = FALSE))$terms$pvalue, floor$terms$pvalue)
    # With a covariate column in the design (KIT-184's case): the same, the contrast by name.
    with_w <- intEoMatrix(covariate = TRUE)
    floor_w <- intEoTest(with_w, "camera_floor")$terms
    floor_w <- floor_w[match(names(index), floor_w$ID), , drop = FALSE]
    camera_w <- intEoTest(with_w, "camera")$terms
    camera_w <- camera_w[match(names(index), camera_w$ID), , drop = FALSE]
    direct_w <- limma::camera(with_w$y, index, with_w$design, contrast = which(colnames(with_w$design) == "factor"),
        inter.gene.cor = NA, allow.neg.cor = FALSE, trend.var = TRUE, use.ranks = FALSE, sort = FALSE)
    expect_identical(floor_w$pvalue, as.numeric(direct_w$PValue))
    expect_true(all(floor_w$pvalue >= camera_w$pvalue))
})

test_that("KIT-191 U2: the permutation test equals the loops of 5.2; reorderings that leave the factor unchanged are counted", {
    fixture <- intEoPermutationFixture()
    out <- intEoTest(fixture, "permutation", intEoSettings(permutations = 1000L), seed = 7L)
    expect_identical(out$method, "permutation")
    reference <- intEoPermutationReference(fixture$y, fixture$f, fixture$sets, 7L, 1000L)
    terms <- out$terms[match(reference$ID, out$terms$ID), , drop = FALSE]
    t_product <- integrationEnrichmentSetT(integrationEnrichmentFisherZ(as.numeric(stats::cor(t(fixture$y), fixture$f))),
        lapply(fixture$sets, function(set) sort(match(set, fixture$ids))))$t[, 1L]
    expect_equal(unname(t_product), unname(reference$t), tolerance = 1e-12)
    expect_identical(terms$pvalue, reference$pvalue)
    expect_identical(terms$Direction, reference$Direction)
    expect_equal(terms$effect, unname(reference$effect), tolerance = 1e-12)
    expect_identical(out$terms$p.adjust, stats::p.adjust(out$terms$pvalue, "BH"))
    expect_true(all(is.na(out$terms$correlation)))
    expect_identical(unique(out$terms$test), "permutation")
    expect_identical(out$terms$setSize[match(reference$ID, out$terms$ID)], unname(lengths(fixture$sets)))
    expect_identical(out$permutation[c("seed", "permutations", "identity_draws", "samples")],
        list(seed = 7L, permutations = 1000L, identity_draws = 0L, samples = 12L))
    expect_identical(out$settings$test_settings, sprintf(INTEGRATION_ENRICHMENT_PERMUTATION_SETTINGS, 1000L, 7L))
    # Four samples: 24 orderings, so about 42 of 1,000 reorderings are the identity; each is counted.
    small <- intEoPermutationFixture(n = 4L, g = 30L)
    small$sets <- list(a_up = small$ids[1:6], d_null = small$ids[13:20], e_tail = small$ids[21:30])
    small$term2gene <- data.frame(term = rep(names(small$sets), lengths(small$sets)),
        gene = unlist(small$sets, use.names = FALSE), stringsAsFactors = FALSE)
    small$term2name <- data.frame(term = names(small$sets), name = names(small$sets), stringsAsFactors = FALSE)
    four <- intEoTest(small, "permutation", intEoSettings(permutations = 1000L), seed = 11L)
    reference_four <- intEoPermutationReference(small$y, small$f, small$sets, 11L, 1000L)
    d <- reference_four$identity_draws[[1L]]
    expect_gt(d, 0L)
    expect_identical(four$permutation$identity_draws, as.integer(d))
    expect_identical(four$permutation$min_p, (1 + d) / 1001)
    terms_four <- four$terms[match(reference_four$ID, four$terms$ID), , drop = FALSE]
    expect_identical(terms_four$pvalue, reference_four$pvalue)
    expect_true(all(four$terms$pvalue >= (1 + d) / 1001))
    # A factor symmetric about its mean (AUDIT-KIT-193 U2): the reversal of (1, 2, 3, 4) gives t* = -t for every set, an
    # exact tie of |t| that the product counts whatever the rounding; with d unchanged draws and r reversals every
    # p-value is at least (1 + d + r) / 1001, and equals the loops'.
    symmetric <- small
    symmetric$f <- c(1, 2, 3, 4)
    symmetric$values$f <- symmetric$f
    symmetric$values$design <- integrationEnrichmentTestDesign(symmetric$f)
    reflected <- intEoTest(symmetric, "permutation", intEoSettings(permutations = 1000L), seed = 11L)
    perms <- withr::with_seed(11L, vapply(seq_len(1000L), function(i) sample(4L), integer(4L)),
        .rng_kind = "Mersenne-Twister", .rng_normal_kind = "Inversion", .rng_sample_kind = "Rejection")
    d_sym <- sum(apply(perms, 2L, function(x) identical(x, 1:4)))
    r_sym <- sum(apply(perms, 2L, function(x) identical(x, 4:1)))
    expect_gt(r_sym, 0L)
    expect_identical(reflected$permutation$identity_draws, as.integer(d_sym))
    expect_true(all(reflected$terms$pvalue >= (1 + d_sym + r_sym) / 1001))
    reference_sym <- intEoPermutationReference(symmetric$y, symmetric$f, symmetric$sets, 11L, 1000L)
    expect_identical(reflected$terms$pvalue[match(reference_sym$ID, reflected$terms$ID)], reference_sym$pvalue)
    expect_identical(reflected$settings$test_settings, sprintf(INTEGRATION_ENRICHMENT_PERMUTATION_SETTINGS, 1000L, 11L))
    expect_match(reflected$settings$test_settings, "ties=abs_within_1e-9", fixed = TRUE)
})

# --- U4, U5, U6: randomness, identity, resolver ---

test_that("KIT-191 U4: the permutation's (and gsea's) result does not depend on RNGkind(); the session's state is kept", {
    fixture <- intEoPermutationFixture()
    reference <- intEoTest(fixture, "permutation", intEoSettings(permutations = 1000L), seed = 5L)$terms
    gsea_fixture <- list(ranked = data.frame(identifier = sprintf("G%02d", 1:60), value = seq(3, -3, length.out = 60)),
        term2gene = data.frame(term = rep(c("top", "bottom", "spread", "middle"), c(10, 10, 10, 15)),
            gene = sprintf("G%02d", c(1:10, 51:60, seq(1, 60, by = 6), 21:35))),
        term2name = data.frame(term = c("top", "bottom", "spread", "middle"), name = c("Top", "Bottom", "Spread", "Middle")))
    gsea <- function() {
        skip_if_not_installed("clusterProfiler")
        integrationEnrichmentGsea(gsea_fixture$ranked, gsea_fixture$term2gene, gsea_fixture$term2name, seed = 11L)$terms
    }
    has_gsea <- requireNamespace("clusterProfiler", quietly = TRUE)
    gsea_default <- if (has_gsea) gsea()
    kinds <- RNGkind()
    withr::defer(RNGkind(kinds[[1L]], kinds[[2L]], kinds[[3L]]))
    for (setting in list(list("L'Ecuyer-CMRG"), list(sample.kind = "Rounding"))) {
        suppressWarnings(do.call(RNGkind, setting))
        set.seed(99L)
        before <- .Random.seed
        again <- intEoTest(fixture, "permutation", intEoSettings(permutations = 1000L), seed = 5L)$terms
        expect_identical(.Random.seed, before)
        expect_identical(again, reference)
        if (has_gsea) {
            expect_identical(gsea(), gsea_default)
            expect_identical(.Random.seed, before)
        }
        RNGkind(kinds[[1L]], kinds[[2L]], kinds[[3L]])
    }
})

test_that("KIT-191 U5, U8: seeds and identities; keys, entries and the replay of the new tests", {
    paths <- intEoProject()
    fixture <- intEoModel(paths)
    manager <- fixture$manager
    view <- "metabolomics_GCMS"
    mapping <- intEoMapping(intEoCorpusSets(manager, fixture$run$value))
    request <- list(factor = "Factor1", view = view, database = "kegg", id_column = "kegg_id")
    run <- function(...) suppressMessages(runIntegrationEnrichment(manager, c(request, list(...)), mapping = mapping))
    first <- run(test = "permutation", seed = 1L)
    other <- run(test = "permutation", seed = 2L)
    again <- run(test = "permutation", seed = 1L)
    # U5: the seed is in the permutation's identity; the same seed is the same product.
    expect_false(identical(first$run$reference$entity_id, other$run$reference$entity_id))
    expect_identical(again$run$status, "hit")
    expect_identical(again$run$reference$entity_id, first$run$reference$entity_id)
    ten <- run(test = "permutation", seed = 1L, permutations = 10000L)
    expect_false(identical(ten$run$reference$entity_id, first$run$reference$entity_id))
    # camera ignores the seed (KIT-183 T9 unchanged).
    camera_a <- run(test = "camera", seed = 1L)
    camera_b <- run(test = "camera", seed = 2L)
    expect_identical(camera_b$run$status, "hit")
    expect_identical(camera_b$run$reference$entity_id, camera_a$run$reference$entity_id)
    floor <- run(test = "camera_floor")
    expect_false(identical(floor$run$reference$entity_id, camera_a$run$reference$entity_id))
    # U8: keys and entries.
    expect_identical(integrationEnrichmentKey("Factor1", view, "kegg", "camera_floor"), "Factor1|metabolomics_GCMS|kegg|camera_floor")
    expect_identical(integrationEnrichmentKey("Factor1", view, "kegg", "permutation"), "Factor1|metabolomics_GCMS|kegg|permutation")
    entries <- manager$.context$snapshot$exports$integration_mofa_model$enrichment
    expect_setequal(names(entries), paste0("Factor1|metabolomics_GCMS|kegg|", c("permutation", "camera", "camera_floor")))
    permutation_entry <- entries[["Factor1|metabolomics_GCMS|kegg|permutation"]]
    expect_identical(permutation_entry$settings$test, "permutation")
    expect_identical(permutation_entry$settings$seed, 1L)
    expect_identical(permutation_entry$settings$permutations, 10000L)
    expect_null(entries[["Factor1|metabolomics_GCMS|kegg|camera"]]$settings$permutations)
    value <- ten$run$value
    expect_identical(value$test, "permutation")
    expect_identical(value$service$seed, 1L)
    expect_identical(value$service$test_settings, sprintf(INTEGRATION_ENRICHMENT_PERMUTATION_SETTINGS, 10000L, 1L))
    expect_identical(value$permutation$permutations, 10000L)
    # The replay restores each entry's test, seed and number of reorderings.
    recorded <- integrationRecordedEnrichment(manager)
    comparison <- suppressMessages(integrationReplayEnrichment(manager, recorded))
    expect_true(all(comparison$terms_identical))
    expect_setequal(comparison$enrichment, names(entries))
    # A second store, the same request: the permutation's terms are reproduced.
    second <- intEoModel(intEoProject())
    reproduced <- suppressMessages(runIntegrationEnrichment(second$manager, c(request, list(test = "permutation", seed = 1L,
        permutations = 10000L)), mapping = mapping))
    expect_identical(reproduced$run$value$terms, value$terms)
})

test_that("KIT-191 U6: the resolver offers four tests; camera_pr and unknown tests refused; the number of reorderings", {
    paths <- intEoProject()
    fixture <- intEoModel(paths)
    object <- integrationImportedObject(fixture$manager, fields = c("features", "design_matrix", "inputs"))
    resolve <- function(...) integrationResolveEnrichmentParameters(list(factor = "Factor1", view = "metabolomics_GCMS",
        database = "kegg", id_column = "kegg_id", ...), fixture$run$value, object, 20L)
    expect_identical(INTEGRATION_ENRICHMENT_TESTS_OFFERED, c("camera", "camera_floor", "permutation", "gsea"))
    for (test in INTEGRATION_ENRICHMENT_TESTS_OFFERED) expect_identical(resolve(test = test)$test, test)
    expect_identical(resolve()$test, "camera")
    expect_identical(resolve()$permutations, 1000L)
    expect_error(resolve(test = "camera_pr"), class = "integration_enrichment_test_not_offered")
    condition <- expect_error(resolve(test = "fry"), class = "integration_enrichment_parameters_invalid")
    for (test in INTEGRATION_ENRICHMENT_TESTS_OFFERED) expect_match(conditionMessage(condition), test, fixed = TRUE)
    expect_identical(resolve(test = "permutation", permutations = 10000)$permutations, 10000L)
    expect_error(resolve(test = "permutation", permutations = 500L), class = "integration_enrichment_parameters_invalid")
    expect_error(resolve(test = "permutation", permutations = 0L), class = "integration_enrichment_parameters_invalid")
    string <- integrationResolveEnrichmentParameters(list(factor = "Factor1", view = "proteomics", database = "string",
        test = "permutation"), fixture$run$value, object, 20L)
    expect_true(is.na(string$test))
    expect_true(string$test_ignored)
    # The new tests need the data set's values, and read them through the same sample floor as camera.
    matrix_fixture <- intEoMatrix()
    for (test in c("camera_floor", "permutation")) {
        expect_error(integrationEnrichmentPathwayTest(matrix_fixture$ranked, matrix_fixture$term2gene,
            matrix_fixture$term2name, 1L, intEoSettings(), test = test, values = NULL),
            class = "integration_enrichment_parameters_invalid")
    }
    expect_error(integrationEnrichmentPermutation(matrix_fixture$values, matrix_fixture$term2gene, matrix_fixture$term2name,
        intEoSettings(permutations = 0L), 1L), class = "integration_enrichment_parameters_invalid")
    # AUDIT-KIT-191-impl I05: the refusals raised for every test on the values name the permutation test too.
    expect_match(INTEGRATION_ENRICHMENT_REFUSALS[["integration_enrichment_too_few_complete"]], "permutation", fixed = TRUE)
    expect_match(INTEGRATION_ENRICHMENT_REFUSALS[["integration_enrichment_too_few_samples_for_test"]], "permutation",
        fixed = TRUE)
    # AUDIT-KIT-191-impl I04: a pre-ranked product (it has `tests`, no `test`) is not read as a correlated one.
    gsea_value <- list(database = "kegg", factor = "Factor1", view = "m", tests = 2L, service = list(release = "r"),
        ranking_counts = list(features = 1L, ranked = 1L, too_few_samples = 0L, no_variation = 0L, internal_standard = 0L),
        mapping_counts = list(mapped = 1L, unmapped = 0L, duplicated = 0L),
        direction = list(feature_id = "x", view = "m", sign = 1L))
    summary <- integrationEnrichmentSummary(gsea_value, integrationGuiTexts())
    expect_true(any(startsWith(summary, "Direction: the factor is oriented")))
    expect_false(any(grepl("Used by the", summary, fixed = TRUE)))
})

test_that("KIT-191 U7 (a): camera's canonical text and the declarations are those of base", {
    expect_identical(INTEGRATION_ENRICHMENT_TEST_SETTINGS[["camera"]], paste0(
        "inter.gene.cor=NA;allow.neg.cor=TRUE;use.ranks=FALSE;trend.var=TRUE;directional=TRUE;",
        "statistic=zscoreT(moderated t, trend, approx=TRUE, method=hill);complete=all_samples;",
        "collapse=largest_mean_among_complete"))
    expect_identical(integrationEnrichmentTestSettingsText("camera", list(seed = 1L, permutations = 1000L)),
        INTEGRATION_ENRICHMENT_TEST_SETTINGS[["camera"]])
    expect_identical(integrationEnrichmentCorrelatedOperation()$parameters, c("ranking", "database", "id_column", "mapping",
        "organism", "min_set_size", "max_set_size", "test", "test_settings", "values_slot", "design_covariates"))
    expect_identical(integrationEnrichmentOperation()$parameters, c("ranking", "database", "id_column", "service_response",
        "mapping", "organism", "string_species", "string_fdr", "string_rank_direction", "min_set_size", "max_set_size", "seed"))
    expect_identical(integrationEnrichmentRankingOperation()$parameters, c("model", "factor", "view", "correlation_method",
        "min_samples", "direction_rule"))
})

# --- U9, U10, U11: record, tab, workbook ---

test_that("KIT-191 U9: the record: evidence sheet, Methods per test, limits per test and kind, strength, several tests", {
    paths <- intEoProject()
    fixture <- intEoModel(paths)
    manager <- fixture$manager
    view <- "metabolomics_GCMS"
    mapping <- intEoMapping(intEoCorpusSets(manager, fixture$run$value))
    request <- list(factor = "Factor1", view = view, database = "kegg", id_column = "kegg_id")
    floor <- runIntegrationEnrichment(manager, c(request, list(test = "camera_floor")), mapping = mapping)
    permutation <- suppressMessages(runIntegrationEnrichment(manager, c(request, list(test = "permutation", seed = 3L)),
        mapping = mapping))
    camera <- runIntegrationEnrichment(manager, request, mapping = mapping)
    texts <- integrationRecordTexts()
    record <- buildWorkflowAnalysisRecord(manager)
    block <- record$analyses$integration_mofa_model$integration
    tables <- integrationRecordTables(record, texts)
    sheet <- tables[[INTEGRATION_RECORD_EVIDENCE_SHEET]]
    expect_identical(nrow(sheet), nrow(INTEGRATION_ENRICHMENT_TEST_EVIDENCE))
    expect_identical(sheet[[texts[["col_evidence_value"]]]], INTEGRATION_ENRICHMENT_TEST_EVIDENCE$value_text)
    expect_true(all(startsWith(sheet[[texts[["col_evidence_source"]]]], "KIT-183 simulation, run 45225854502")))
    expect_setequal(unique(sheet[[texts[["col_evidence_kind"]]]]), c("metabolite", "gene and protein"))
    # Methods: each entry's sentence, its fragments formatted, no format code left.
    methods <- integrationMethodsEnrichment(block, texts)
    expect_false(any(grepl("%[sd]", methods)))
    counts <- permutation$run$value$test_counts
    permutation_sentence <- sprintf(texts[["methods_enrich_permutation"]], "Factor1", view, "KEGG",
        sprintf(texts[["methods_test_name_permutation"]], 1000L, 1000L),
        sprintf(texts[["methods_test_settings_permutation"]], 3L, counts$samples, counts$universe,
            permutation$run$value$mapping_counts$mapped), "KEGG Release fixture", "hsa", 5L, 500L,
        permutation$run$value$service$pathways_testable)
    expect_true(permutation_sentence %in% methods)
    floor_sentence <- sprintf(texts[["methods_enrich_camera_floor"]], "Factor1", view, "KEGG",
        texts[["methods_test_name_camera_floor"]], sprintf(texts[["methods_test_settings_camera_floor"]], counts$samples,
            counts$universe, floor$run$value$mapping_counts$mapped), "KEGG Release fixture", "hsa", 5L, 500L,
        floor$run$value$service$pathways_testable)
    expect_true(floor_sentence %in% methods)
    # The limits: one paragraph, a clause per test used and kind, with that test's metabolite figures.
    limits <- methods[startsWith(methods, texts[["methods_enrich_limits_tests"]])]
    expect_length(limits, 1L)
    expect_match(limits, "the permutation test, on metabolite data, called terms with no link to the factor in 4.2 to 5.2 %",
        fixed = TRUE)
    expect_match(limits, "in up to 8.6 %", fixed = TRUE)
    expect_match(limits, "negative correlations counted as zero, on metabolite data, called terms with no link to the factor in 1.7 to 3.6 %",
        fixed = TRUE)
    expect_false(grepl("pre-ranked test, on", limits, fixed = TRUE))
    expect_false(grepl("gene and protein", limits, fixed = TRUE))
    # Several tests on one factor, data set and database.
    expect_true(sprintf(texts[["methods_enrich_several_tests"]], "Factor1", view, "KEGG", 3L,
        "correlation-adjusted test, correlation-adjusted test with negative correlations counted as zero, permutation test") %in%
        methods)
    # The strength sentence carries the figure of the entry's test and kind.
    strength <- tables$integration_enrichment[[texts[["col_enrich_strength"]]]]
    expect_true(any(grepl("up to 8.6 % of cases by the permutation test", strength, fixed = TRUE)))
    expect_true(any(grepl("up to 6.7 % of cases by the correlation-adjusted test with negative", strength, fixed = TRUE)))
    entry <- function(omic, test) list(database = "kegg", omic = omic, test = test, factor_r2_percent = 5,
        largest_r2_percent = 20)
    expect_match(integrationEnrichmentStrengthText(entry("proteomics", "permutation"), texts), "up to 12.7 % of cases by the permutation test",
        fixed = TRUE)
    expect_match(integrationEnrichmentStrengthText(entry("proteomics", "gsea"), texts), "up to 60.1 % of cases by the pre-ranked test",
        fixed = TRUE)
    expect_match(integrationEnrichmentStrengthText(entry("metabolomics", "gsea"), texts), "up to 34.4 % of cases",
        fixed = TRUE)
    expect_match(integrationEnrichmentStrengthText(entry("lipidomics", "camera"), texts), "no lipid data were simulated",
        fixed = TRUE)
    # AUDIT-KIT-191-impl I02: gsea's gene-and-protein weaker-factor figure (60.1 %) is of block members the model did not
    # use, so the limits clause names them with the kind's phrase.
    gp <- list(complete = TRUE, database = "kegg", omic = "proteomics", test = "gsea", settings = list(test = "gsea"))
    expect_match(integrationMethodsEnrichmentLimits(list(gp), texts), paste("terms made of features the model used, or of",
        "features that vary together with them in up to 60.1 %"), fixed = TRUE)
    # AUDIT-KIT-191-impl I03: a block made by KIT-183's merged code holds the pre-ranked product's count of tests in `test`;
    # the declared test is read first, and a value that is no offered test is not taken for camera.
    old <- list(database = "kegg", omic = "metabolomics", test = 12L, settings = list(test = "gsea", omic = "metabolomics"),
        factor_r2_percent = 5, largest_r2_percent = 20)
    expect_identical(integrationEnrichmentEntryTest(old), "gsea")
    expect_match(integrationEnrichmentStrengthText(old, texts), "up to 34.4 % of cases by the pre-ranked test", fixed = TRUE)
    unknown <- list(database = "kegg", omic = "metabolomics", test = "fry", settings = list(test = "fry"),
        factor_r2_percent = 5, largest_r2_percent = 20)
    expect_true(is.na(integrationEnrichmentEntryTest(unknown)))
    expect_identical(integrationEnrichmentStrengthText(unknown, texts), "")
    # Only STRING entries: no evidence sheet.
    string_only <- record
    string_only$analyses$integration_mofa_model$integration$enrichment <- list(list(complete = TRUE, database = "string",
        factor = "Factor1", view = "proteomics", settings = list(correlation_method = "pearson", min_samples = 5L),
        direction = list(feature_id = "P1", view = "proteomics", sign = 1L), ranking_counts = list(features = 1L, ranked = 1L),
        mapping_counts = list(mapped = 10L), service = list(version = list(string_version = "12.5"), host = "x",
            species = "9606", queried_at = "x"), tests = 0L, terms = list()))
    expect_null(integrationRecordTables(string_only, texts)[[INTEGRATION_RECORD_EVIDENCE_SHEET]])
    expect_false(any(startsWith(integrationMethodsEnrichment(string_only$analyses$integration_mofa_model$integration, texts),
        texts[["methods_enrich_limits_tests"]])))
})

test_that("KIT-191 U10: the tab: four tests, camera pre-selected, the kind's notes and evidence, seed, floor, counts", {
    skip_if_not_installed("clusterProfiler")
    paths <- intEoProject()
    fixture <- intEoModel(paths)
    manager <- fixture$manager
    texts <- integrationGuiTexts()
    escaped <- function(text) as.character(htmltools::htmlEscape(text))
    attribute <- function(text) as.character(htmltools::htmlEscape(text, attribute = TRUE))
    withr::local_options(MultiScholaR.integration_gui = TRUE)
    html <- as.character(htmltools::renderTags(mod_integration_enrichment_ui("integration_enrichment"))$html)
    # The seed is shown for the pre-ranked and permutation tests only; the reorderings for the permutation test only.
    expect_match(html, attribute(paste("(input['integration_enrichment-database'] == 'kegg' ||",
        "input['integration_enrichment-database'] == 'reactome') && (input['integration_enrichment-test'] == 'gsea' ||",
        "input['integration_enrichment-test'] == 'permutation')")), fixed = TRUE)
    expect_match(html, attribute("input['integration_enrichment-test'] == 'permutation'"), fixed = TRUE)
    expect_match(html, 'data-testid="integration-enrichment-seed"', fixed = TRUE)
    expect_match(html, 'data-testid="integration-enrichment-permutations"', fixed = TRUE)
    expect_match(html, escaped(texts[["enrich_permutations_10000"]]), fixed = TRUE)
    sets <- intEoCorpusSets(manager, fixture$run$value)
    bus <- do.call(shiny::reactiveValues, omicWorkflowInitialBus(integrationGuiConfig()$shell))
    bus$state_manager <- manager
    bus$integration_model <- list(reference = fixture$run$reference, status = fixture$run$status)
    no_progress <- function(message, value, expr) force(expr)
    shiny::testServer(mod_integration_enrichment_server, args = list(workflow_data = bus, paths = paths,
        download = intEoKeggDownload(sets), withProgressFn = no_progress), {
        session$flushReact()
        session$setInputs(factor = "Factor1", view = "metabolomics_GCMS", database = "kegg")
        selector <- as.character(output$test_ui$html)
        expect_match(selector, '<option value="camera" selected>', fixed = TRUE)
        for (test in INTEGRATION_ENRICHMENT_TESTS_OFFERED) {
            expect_match(selector, sprintf('<option value="%s"', test), fixed = TRUE)
            expect_match(selector, escaped(texts[[paste0("enrich_test_", test)]]), fixed = TRUE)
        }
        expect_false(grepl("camera_pr", selector, fixed = TRUE))
        # The metabolite rows and notes.
        evidence <- as.character(output$evidence$html)
        for (measure in c("A", "B", "D", "C3", "I")) expect_match(evidence, escaped(texts[[paste0("enrich_evidence_", measure)]]),
            fixed = TRUE)
        expect_match(evidence, escaped(texts[["enrich_evidence_E_metabolites"]]), fixed = TRUE)
        expect_false(grepl(escaped(texts[["enrich_evidence_G"]]), evidence, fixed = TRUE))
        expect_match(evidence, escaped(integrationEnrichmentEvidenceValue("metabolites", "permutation", "E")), fixed = TRUE)
        expect_match(evidence, escaped(sprintf(texts[["enrich_evidence_heading"]], "metabolite")), fixed = TRUE)
        expect_false(grepl(escaped(texts[["enrich_evidence_lipids"]]), evidence, fixed = TRUE))
        for (test in INTEGRATION_ENRICHMENT_TESTS_OFFERED) {
            session$setInputs(test = test)
            expect_match(as.character(output$test_note$html), escaped(texts[[paste0("enrich_test_", test, "_note_metabolites")]]),
                fixed = TRUE)
        }
        # A lipid data set: the metabolite figures, with the line saying no lipid data were simulated.
        lipid <- names(state$model$omics)[state$model$omics == "lipidomics"][[1L]]
        session$setInputs(view = lipid, database = "kegg")
        expect_match(as.character(output$evidence$html), escaped(texts[["enrich_evidence_lipids"]]), fixed = TRUE)
        # A permutation run: its counts and its floor.
        session$setInputs(view = "metabolomics_GCMS", database = "kegg", test = "permutation", permutations = "1000",
            id_column = "kegg_id", correlation_method = "pearson", min_samples = 5L, seed = 20261121L, string_species = 9606L,
            string_fdr = 0.05, kegg_organism = "hsa", reactome_organism = "Homo sapiens", min_set_size = 5L,
            max_set_size = 500L)
        session$setInputs(download = 1)
        session$setInputs(run = 1)
        expect_identical(bus$integration_enrichment$entry, "Factor1|metabolomics_GCMS|kegg|permutation")
        result <- as.character(output$result$html)
        value <- state$run$run$value
        counts <- value$test_counts
        expect_match(result, escaped(sprintf(texts[["enrich_test_counts_permutation"]], counts$universe, counts$samples,
            counts$incomplete, counts$collapsed)), fixed = TRUE)
        expect_match(result, escaped(integrationEnrichmentPermutationFloorText(value$permutation, as.integer(value$tests),
            texts)), fixed = TRUE)
        expect_match(result, "up to 8.6 % of cases by the permutation test", fixed = TRUE)
        # STRING: no selector, no evidence.
        session$setInputs(view = "proteomics", database = "string")
        expect_false(grepl(sprintf('id="%s"', session$ns("test")), paste(as.character(output$test_ui$html), collapse = ""),
            fixed = TRUE))
        expect_identical(paste(as.character(output$evidence$html), collapse = ""), "")
    })
    # Gene and protein data sets (Gene Ontology, KIT-181): rows A, B, D, E, G, H.
    genes <- as.character(htmltools::renderTags(integrationEnrichmentEvidenceUi("proteomics", texts))$html)
    for (measure in c("A", "B", "D", "G", "H")) expect_match(genes, escaped(texts[[paste0("enrich_evidence_", measure)]]),
        fixed = TRUE)
    expect_false(grepl(escaped(texts[["enrich_evidence_I"]]), genes, fixed = TRUE))
    expect_match(genes, escaped(integrationEnrichmentEvidenceValue("genes_proteins", "camera", "B")), fixed = TRUE)
    # The few-samples sentence per test and kind, with the most variable features' figure.
    notes <- integrationEnrichmentTestNotes("camera", 8L, "uncorrected", "proteomics", texts)
    expect_identical(notes, c(texts[["enrich_test_camera_note_genes_proteins"]],
        sprintf(texts[["enrich_test_few_samples_camera_genes_proteins"]], 8L)))
    expect_match(notes[[2L]], "up to 9.9 %", fixed = TRUE)
    for (test in INTEGRATION_ENRICHMENT_TESTS_OFFERED) for (omic in c("metabolomics", "proteomics", "lipidomics")) {
        few <- integrationEnrichmentTestNotes(test, 6L, "uncorrected", omic, texts)
        expect_length(few, 2L)
        expect_match(few[[2L]], "This data set has 6 biological samples", fixed = TRUE)
    }
    expect_length(integrationEnrichmentTestNotes("permutation", 12L, "uncorrected", "metabolomics", texts), 1L)
    # The floor (corrected, AUDIT-KIT-193 A05): 2,000 terms, 1,000 reorderings, none unchanged: none is called unless at
    # least 40 are called together.
    # The floor printed exactly as (1 + d)/(1 + B) (AUDIT-KIT-191-impl I01: a three-figure decimal rounded to nearest can
    # exceed it), and k_min exact.
    floorOf <- function(B, d, n, m) integrationEnrichmentPermutationFloorText(list(permutations = B, identity_draws = d,
        min_p = (1 + d) / (1 + B), samples = n), m, texts)
    expect_identical(integrationEnrichmentFloorFraction(0L, 1000L), "1/1,001")
    floor <- floorOf(1000L, 0L, 12L, 2000L)
    expect_identical(floor, sprintf(texts[["enrich_test_permutation_floor"]], 1000L, 12L, "1/1,001", 2000L, 40L))
    expect_match(floor, "none can be called at a 5 % false discovery rate unless at least 40 are called together",
        fixed = TRUE)
    expect_identical(floorOf(1000L, 8L, 5L, 30L), sprintf(texts[["enrich_test_permutation_floor"]], 1000L, 5L, "9/1,001",
        30L, 6L))
    # Cases a rounded decimal printed above the floor: 1/10,001 (was "0.0001") and 2/1,001 (was "0.002").
    expect_identical(floorOf(10000L, 0L, 12L, 2000L), sprintf(texts[["enrich_test_permutation_floor"]], 10000L, 12L,
        "1/10,001", 2000L, 4L))
    expect_identical(floorOf(1000L, 1L, 6L, 100L), sprintf(texts[["enrich_test_permutation_floor"]], 1000L, 6L, "2/1,001",
        100L, 4L))
    # A whole-number bound: B = 1,000, d = 10, m = 273: m (1 + d) / (0.05 (1 + B)) = 60 exactly, so k_min = 60, not 61.
    expect_identical(floorOf(1000L, 10L, 5L, 273L), sprintf(texts[["enrich_test_permutation_floor"]], 1000L, 5L, "11/1,001",
        273L, 60L))
    # p_min above 0.05 (4 samples, 50 unchanged draws): no term can be called; the _none sentence.
    expect_identical(floorOf(1000L, 50L, 4L, 10L), sprintf(texts[["enrich_test_permutation_floor_none"]], 1000L, 4L,
        "51/1,001"))
    # The sentence is true of A05's counterexample: 100 terms, four at 2/1001, the rest at 1: BH calls the four, with
    # no p-value at or below 1/1000, and four is at least the sentence's ceiling(100 (1/1001) / 0.05) = 2.
    p <- c(rep(2 / 1001, 4L), rep(1, 96L))
    expect_identical(sum(stats::p.adjust(p, "BH") <= 0.05), 4L)
    expect_identical(sum(p <= 0.001), 0L)
    expect_gte(sum(stats::p.adjust(p, "BH") <= 0.05), ceiling(100 * (1 / 1001) / 0.05))
})

test_that("KIT-191 U11: the workbook: the permutation and camera_floor overrides run and record; STRING has no test", {
    paths <- intEoProject()
    fixture <- intEoModel(paths)
    manager <- fixture$manager
    mapping <- intEoMapping(intEoCorpusSets(manager, fixture$run$value))
    local_mocked_bindings(integrationDownloadKeggMapping = function(organism, download = NULL) mapping,
        integrationDownloadReactomeMapping = function(download = NULL) list(database = "Reactome", release = "Reactome 0",
            sha256 = strrep("f", 64L), url = "fixture", downloaded_at = "x", table = data.frame(chebi_id = "1",
                reactome_id = "R-HSA-1", reactome_term = "x", organism = "Homo sapiens", stringsAsFactors = FALSE)))
    project <- list(manager = manager, store_root = dirname(manager$.context$store$directory), paths = paths,
        views = c("metabolomics_GCMS", "metabolomics_LCMS_Pos"))
    results <- list(results = list(views = "metabolomics_GCMS"), factor = "Factor1")
    id_columns <- list(metabolomics_GCMS = "kegg_id")
    expect_message(permuted <- integrationEnrichFactor(project, results, id_columns = id_columns,
        overrides = list(test = "permutation", seed = 3L)), "predicted to take about")
    terms <- permuted$enrichments[["metabolomics_GCMS kegg"]]$terms
    expect_true(all(terms$test == "permutation"))
    entry <- manager$.context$snapshot$exports$integration_mofa_model$enrichment[["Factor1|metabolomics_GCMS|kegg|permutation"]]
    expect_identical(entry$settings$seed, 3L)
    expect_identical(entry$settings$permutations, 1000L)
    floor <- integrationEnrichFactor(project, results, id_columns = id_columns, overrides = list(test = "camera_floor"))
    expect_true(all(floor$enrichments[["metabolomics_GCMS kegg"]]$terms$test == "camera_floor"))
    expect_true("permutations" %in% integrationSettings()$enrichment$setting)
    # STRING with a test given: the request sent to STRING carries none (STRING ignores it).
    sent <- list()
    local_mocked_bindings(runIntegrationEnrichment = function(manager, parameters, ...) {
        sent[[length(sent) + 1L]] <<- parameters
        integrationEnrichmentRefuse("integration_enrichment_parameters_invalid", "stand-in")
    })
    project$views <- "proteomics"
    integrationEnrichFactor(project, list(results = list(views = "proteomics"), factor = "Factor1"), api_key = "k",
        overrides = list(test = "permutation"))
    expect_length(sent, 1L)
    expect_identical(sent[[1L]]$database, "string")
    expect_null(sent[[1L]]$test)
})

# --- U12: the cost guard ---

test_that("KIT-191 U12: the guard: the prediction, refused above the tab's limit, a message in scripts; L1's coefficients", {
    cost <- list(a = 1e-9, c = 2e-8, s = 3e-8)
    expect_equal(integrationEnrichmentPermutationSeconds(12, 4000, 1000, 15000, cost),
        1e-9 * 12 * 4000 * 1000 + 2e-8 * 4000 * 1000 + 3e-8 * 15000 * 1000)
    fixture <- intEoMatrix()
    texts <- integrationGuiTexts()
    local_mocked_bindings(integrationEnrichmentPermutationSeconds = function(...) 119)
    runs <- integrationEnrichmentPathwayTest(fixture$ranked, fixture$term2gene, fixture$term2name, 1L,
        intEoSettings(permutations = 1000L, time_limit = INTEGRATION_ENRICHMENT_TAB_SECONDS), test = "permutation",
        values = fixture$values)
    expect_identical(nrow(runs$terms), 4L)
    local_mocked_bindings(integrationEnrichmentPermutationSeconds = function(...) 121)
    refused <- expect_error(integrationEnrichmentPathwayTest(fixture$ranked, fixture$term2gene, fixture$term2name, 1L,
        intEoSettings(permutations = 1000L, time_limit = INTEGRATION_ENRICHMENT_TAB_SECONDS), test = "permutation",
        values = fixture$values), class = "integration_enrichment_test_too_large")
    expect_identical(conditionMessage(refused), sprintf(INTEGRATION_ENRICHMENT_REFUSALS[["integration_enrichment_test_too_large"]],
        "2.0", 40L, 12L, 4L, 1000L))
    expect_identical(refused$predicted_seconds, 121)
    # Scripts (no limit): the message, then the run.
    expect_message(scripted <- integrationEnrichmentPathwayTest(fixture$ranked, fixture$term2gene, fixture$term2name, 1L,
        intEoSettings(permutations = 1000L), test = "permutation", values = fixture$values),
        sprintf(texts[["enrich_test_permutation_time"]], "121", 40L, 12L, 4L, 1000L), fixed = TRUE)
    expect_identical(scripted$terms, runs$terms)
    expect_identical(INTEGRATION_ENRICHMENT_TAB_SECONDS, 120)
    # The coefficients are L1's (docs/qa/workflow-execution/kit/kit-191/lifecycle.md).
    expect_identical(INTEGRATION_ENRICHMENT_PERMUTATION_COST, list(a = 5.12e-9, c = 2.13e-7, s = 5.48e-8, run = "45285106025"))
})

# --- U13: the corpus control ---

test_that("KIT-191 U13: corpus control: camera_floor and permutation find the planted sets in their direction, not the decoy", {
    paths <- intEoProject()
    fixture <- intEoModel(paths)
    manager <- fixture$manager
    value <- fixture$run$value
    export <- manager$.context$snapshot$exports$integration_mofa_model
    truth <- utils::read.delim(file.path(intEoData(), "corpus", "truth_scores.tsv"), stringsAsFactors = FALSE)
    loadings <- utils::read.delim(file.path(intEoData(), "corpus", "truth_loadings.tsv"), stringsAsFactors = FALSE)
    factors <- setdiff(names(value$factors), "biological_sample")
    oriented <- vapply(factors, function(factor) {
        sign <- integrationEnrichmentFactorDirection(manager, value, factor)$sign
        stats::cor(sign * value$factors[[factor]], truth$Z1[match(value$factors$biological_sample, truth$biological_sample)])
    }, numeric(1))
    factor <- factors[[which.max(abs(oriented))]]
    expected_up <- if (oriented[[factor]] > 0) "Up" else "Down"
    expected_down <- setdiff(c("Up", "Down"), expected_up)
    object <- integrationImportedObject(manager, fields = c("features", "design_matrix", "inputs"))
    for (view in c("proteomics", "metabolomics_GCMS")) {
        omic <- integrationViewOmics(object)[[view]]
        p <- if (omic == "proteomics") list(database = "string", omic = omic, id_type = "uniprot") else
            list(database = "kegg", omic = omic, id_type = "kegg", id_column = "kegg_id")
        ranking <- integrationEnrichmentRank(object, manager, export$source, value, export$parameters, factor, view,
            "pearson", 5L)
        mapped <- integrationEnrichmentMappedList(ranking, p, object@features[[view]])
        values <- integrationEnrichmentTestValues(object, manager, export$source, value, export$parameters, ranking,
            mapped, view)
        rows <- loadings[loadings$view == view, , drop = FALSE]
        kept <- values$feature_id
        id_of <- stats::setNames(rownames(values$y), kept)
        pool <- function(rule) intersect(rows$feature_id[rule], kept)
        size <- if (omic == "proteomics") 20L else 10L
        draw <- function(x, k) withr::with_seed(20261184L, sample(x, min(k, length(x))))
        sets <- list(up = id_of[draw(pool(rows$Z1 > 0), size)], down = id_of[draw(pool(rows$Z1 < 0), size)],
            decoy = id_of[draw(pool(rows$Z1 == 0 & rows$Z2 == 0 & (is.na(rows$Z3) | rows$Z3 == 0)), size)])
        term2gene <- data.frame(term = rep(names(sets), lengths(sets)), gene = unlist(sets, use.names = FALSE),
            stringsAsFactors = FALSE)
        term2name <- data.frame(term = names(sets), name = names(sets), stringsAsFactors = FALSE)
        ranked <- mapped$table[c("identifier", "value")]
        for (test in c("camera_floor", "permutation")) {
            terms <- suppressMessages(integrationEnrichmentPathwayTest(ranked, term2gene, term2name, 20261121L,
                intEoSettings(permutations = 1000L), test, values))$terms
            row <- function(id) terms[terms$ID == id, , drop = FALSE]
            info <- paste(view, test)
            expect_lte(row("up")$p.adjust, 0.05)
            expect_identical(row("up")$Direction, expected_up, info = info)
            expect_lte(row("down")$p.adjust, 0.05)
            expect_identical(row("down")$Direction, expected_down, info = info)
            expect_gt(row("decoy")$p.adjust, 0.05)
        }
    }
})

# --- The cross-platform record (ticket 6): camera_floor on U1's fixture and the permutation on U2's, the terms printed
# at full precision (the reference, from a Linux run-tests job, is tests/testdata/integration/kit191-cross-platform.tsv),
# and on another platform compared by the rule: the same IDs and directions and every number within 1e-6 ---

intEoPlatformTerms <- function() {
    floor <- intEoTest(intEoMatrix(), "camera_floor")$terms
    permutation <- intEoTest(intEoPermutationFixture(), "permutation", intEoSettings(permutations = 1000L), seed = 7L)$terms
    rbind(data.frame(test = "camera_floor", ID = floor$ID, direction = floor$Direction, pvalue = floor$pvalue,
            p.adjust = floor$p.adjust, effect = floor$effect),
        data.frame(test = "permutation", ID = permutation$ID, direction = permutation$Direction, pvalue = permutation$pvalue,
            p.adjust = permutation$p.adjust, effect = permutation$effect))
}

test_that("KIT-191 cross-platform: camera_floor and the permutation agree with the Linux reference within 1e-6", {
    terms <- intEoPlatformTerms()
    terms <- terms[order(terms$test, terms$ID), ]
    rownames(terms) <- NULL
    text <- utils::capture.output(utils::write.table(format(terms, digits = 17), sep = "\t", quote = FALSE,
        row.names = FALSE))
    cat("\nKIT-191 cross-platform terms (sha256 ", digest::digest(text, algo = "sha256"), "):\n", sep = "")
    cat(text, sep = "\n")
    reference <- file.path(intEoData(), "kit191-cross-platform.tsv")
    testthat::skip_if_not(file.exists(reference), "the cross-platform reference is not recorded yet")
    ref <- utils::read.delim(reference, stringsAsFactors = FALSE)
    ref <- ref[order(ref$test, ref$ID), ]
    expect_identical(terms$ID, ref$ID)
    expect_identical(terms$direction, ref$direction)
    for (column in c("pvalue", "p.adjust", "effect")) {
        difference <- max(abs(terms[[column]] - ref[[column]]), na.rm = TRUE)
        expect_lte(difference, 1e-6, label = sprintf("largest absolute difference of %s (%g)", column, difference))
    }
})
