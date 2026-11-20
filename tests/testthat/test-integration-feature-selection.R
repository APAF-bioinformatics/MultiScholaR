# KIT-151 (INT-12): the integration model's feature selection on a fixture built here where the selection acts (audit
# SCI-13: the tier-1 corpus has at most 400 features per data set, so a cap of 1,000 never acts on it). Declared before
# the first run:
#   - ENG-20: on a data set made by the selection evidence's generator (KIT-151 run b, `k151b.R` on branch
#     run-mofa-evidence-k151b, sha256 cc313336..., scenario `full`: a mean-variance trend and abundance-dependent missing
#     values), the package's alignment keeps the same features, in the same order, as the evidence script's own
#     selection (its `prepare_views`, written out below as it is there: observed in at least half the samples, then the
#     variance minus the lowess trend at the mean, the 1,000 largest, ties by position); also for raw variance;
#   - the alignment in row pages keeps the same features (identical value);
#   - SCI-03: a planted feature with by far the largest variance, observed in 3 of the 20 samples, is left out by the
#     default selection (counted as sparse) and would be kept by raw variance if the observed share were not required.

intFsWriteInput <- function(root, omic, view, values) {
    biological <- colnames(values)
    runs <- sprintf("%s%02d", substr(omic, 1L, 1L), seq_along(biological))
    colnames(values) <- runs
    matrix <- data.frame(feature_id = rownames(values), values, check.names = FALSE, stringsAsFactors = FALSE)
    features <- data.frame(feature_id = rownames(values), display_name = tolower(rownames(values)), stringsAsFactors = FALSE)
    samples <- data.frame(Run = runs, group = sub("_.*$", "", biological), replicates = 1L, tech_rep_group = biological,
        stringsAsFactors = FALSE)
    manifest <- list(omic_type = omic, workflow = "TEST", matrix_role = "normalised", scale = "log2", correction = NULL,
        source_state = list(name = "kit151-test"), sample_id_column = "Run", group_column = "group",
        technical_replicate_column = "tech_rep_group", biological_sample_column = "tech_rep_group")
    dir <- file.path(root, "integration", omic)
    writeIntegrationInput(dir, manifest, samples, stats::setNames(list(matrix), view), stats::setNames(list(features), view))
    dir
}

# The evidence generator's data set D1 in scenario `full` (k151b.R `generate`, its constants), 3,000 features x 20
# samples, without the planted factors (the selection does not use them), plus one planted sparse feature (SCI-03).
intFsFixture <- function(.local_envir = parent.frame()) {
    root <- withr::local_tempdir(.local_envir = .local_envir)
    biological <- sprintf("%s_%d", rep(c("A", "B"), each = 10L), rep(seq_len(10L), 2L))
    d <- 3000L
    n <- length(biological)
    x <- withr::with_seed(151L, {
        offset <- stats::rnorm(d, 20, 2)
        sd_j <- 1.5 * 2^(-(offset - 20) / 4)
        x <- matrix(stats::rnorm(d * n), d, n) * sd_j + offset
        p <- 1 - (1 - stats::plogis((17.5 - x) / 0.5)) * (1 - 0.02)
        x[matrix(stats::runif(length(x)), nrow(x)) < p] <- NA_real_
        x
    })
    # SCI-03: by far the largest variance, observed in 3 of the 20 samples.
    x[d, ] <- NA_real_
    x[d, c(1L, 7L, 15L)] <- c(5, 35, 20)
    dimnames(x) <- list(sprintf("PRO%05d", seq_len(d)), biological)
    small <- matrix(withr::with_seed(152L, stats::rnorm(50L * n, 25, 2)), 50L, n,
        dimnames = list(sprintf("LIP%03d", seq_len(50L)), biological))
    dirs <- c(intFsWriteInput(root, "proteomics", "proteomics", x), intFsWriteInput(root, "lipidomics", "lipidomics_LCMS_Pos", small))
    list(dirs = dirs, values = x, sparse = rownames(x)[[d]])
}

# The evidence script's selection (k151b.R, `fitted_trend` and `prepare_views` with eligibility `half`, no z-score),
# as written there.
intFsEvidenceSelection <- function(m, select, cap) {
    TREND_SPAN <- 0.3
    HALF_SHARE <- 0.5
    fitted_trend <- function(abundance, y) {
        fit <- stats::lowess(abundance, y, f = TREND_SPAN)
        stats::approx(fit$x, fit$y, xout = abundance, ties = base::mean, rule = 2)$y
    }
    need <- ceiling(HALF_SHARE * ncol(m))
    m <- m[rowSums(!is.na(m)) >= need, , drop = FALSE]
    if (nrow(m) > cap) {
        variance <- apply(m, 1L, stats::var, na.rm = TRUE)
        score <- switch(select,
            raw = variance,
            trend_diff = variance - fitted_trend(rowMeans(m, na.rm = TRUE), variance))
        keep <- sort(order(-score, seq_along(score))[seq_len(cap)])
        m <- m[keep, , drop = FALSE]
    }
    rownames(m)
}

test_that("KIT-151 (ENG-20): on the evidence generator's data the package keeps the features the evidence script keeps", {
    fixture <- intFsFixture()
    object <- integrationViewsFromInputs(fixture$dirs)$object
    for (case in list(c("variance_above_trend", "trend_diff"), c("top_variance", "raw"))) {
        parameters <- integrationResolveMofaParameters(object, list(feature_selection = case[[1L]]))
        expect_identical(parameters$features_per_view, 1000L)
        expect_identical(parameters$min_observed_fraction, 0.5)
        aligned <- integrationAlignViews(object, parameters)
        kept <- rownames(aligned$views$proteomics)
        expected <- intFsEvidenceSelection(fixture$values[, colnames(aligned$views$proteomics)], case[[2L]], 1000L)
        expect_identical(kept, expected, info = case[[1L]])
        expect_length(kept, 1000L)
        # The small data set keeps every feature.
        expect_identical(nrow(aligned$views$lipidomics_LCMS_Pos), 50L, info = case[[1L]])
        # The same in row pages.
        expect_identical(integrationAlignViewsPartitioned(object, parameters, integrationObjectViewPages(object)), aligned,
            info = case[[1L]])
    }
    # The two rules choose differently here, and the default is the trend's.
    trend <- rownames(integrationAlignViews(object, integrationResolveMofaParameters(object))$views$proteomics)
    raw <- rownames(integrationAlignViews(object, integrationResolveMofaParameters(object,
        list(feature_selection = "top_variance")))$views$proteomics)
    expect_identical(trend, intFsEvidenceSelection(fixture$values[, sort(colnames(fixture$values), method = "radix")],
        "trend_diff", 1000L))
    expect_false(identical(trend, raw))
})

test_that("KIT-151 (SCI-03): a feature seen in few samples is not selected for its chance variance", {
    fixture <- intFsFixture()
    object <- integrationViewsFromInputs(fixture$dirs)$object
    aligned <- integrationAlignViews(object, integrationResolveMofaParameters(object))
    summary <- aligned$view_summary
    expect_false(fixture$sparse %in% rownames(aligned$views$proteomics))
    observed <- rowSums(!is.na(fixture$values))
    expect_identical(summary$sparse_dropped[summary$view == "proteomics"], sum(observed < 10L))
    # Without the observed share (any share up to 0.1 asks for two samples) raw variance would take it first.
    loose <- integrationAlignViews(object, integrationResolveMofaParameters(object, list(feature_selection = "top_variance",
        min_observed_fraction = 0.1)))
    expect_identical(rownames(loose$views$proteomics)[[1000L]], fixture$sparse)
    expect_identical(integrationMinObserved(list(feature_selection = "top_variance", min_observed_fraction = 0.3), 10L), 3L)
    expect_identical(integrationMinObserved(list(feature_selection = "none", min_observed_fraction = NA_real_), 20L), 2L)
})
