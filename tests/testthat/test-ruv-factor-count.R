# KIT-129 W6: the number of unwanted factors (k). The estimated count is `sva::num.sv(method = "be", B = 20, seed)` on
# the features complete in every sample, as KIT-124's simulators call it (arms 7, 8, 9, 10, 15, 16); `sva` is an
# optional dependency. The diagnostic's quantities (observed and permuted shares of residual variance, p-values) are
# not returned by `num.sv`, so they are computed with its algorithm on its random stream; the first test holds the
# two together. The separation rule (`findBestKElbow`, `trxFindBestK`) is not touched.

ruvCountFixture <- function(seed, features = 200L, samples = 12L, factors = 1L, groups = 2L, missing = 0) {
    set.seed(seed)
    group <- rep(LETTERS[seq_len(groups)], length.out = samples)
    x <- matrix(stats::rnorm(features * samples), features, samples)
    if (factors > 0L) {
        x <- x + matrix(stats::rnorm(features * factors), features) %*% t(matrix(stats::rnorm(samples * factors), samples))
    }
    if (missing > 0) x[sample(length(x), round(missing * length(x)))] <- NA_real_
    dimnames(x) <- list(sprintf("F%04d", seq_len(features)), sprintf("s%02d", seq_len(samples)))
    list(x = x, group = group)
}

test_that("the estimated count is sva::num.sv's on seeded data sets, with the quantities of the diagnostic", {
    skip_if_not_installed("sva")
    cells <- list(
        list(seed = 1L, features = 200L, samples = 8L, factors = 1L),
        list(seed = 2L, features = 1000L, samples = 12L, factors = 2L),
        list(seed = 3L, features = 40L, samples = 24L, factors = 3L),
        list(seed = 4L, features = 200L, samples = 9L, factors = 0L, groups = 3L),
        list(seed = 5L, features = 1000L, samples = 6L, factors = 1L),
        list(seed = 6L, features = 300L, samples = 12L, factors = 2L, missing = 0.02),
        list(seed = 7L, features = 300L, samples = 12L, factors = 0L),
        list(seed = 8L, features = 150L, samples = 16L, factors = 3L, groups = 4L)
    )
    counts <- integer()
    for (cell in cells) {
        fixture <- do.call(ruvCountFixture, cell)
        complete <- fixture$x[rowSums(!is.finite(fixture$x)) == 0, , drop = FALSE]
        expected <- sva::num.sv(complete, stats::model.matrix(~ group, data.frame(group = fixture$group)), method = "be",
            B = 20, seed = 900L + cell$seed)
        estimate <- estimateRuvFactorCount(fixture$x, fixture$group, seed = 900L + cell$seed)
        info <- paste("seed", cell$seed)
        expect_identical(estimate$k, as.integer(expected), info = info)
        expect_identical(estimate$rule, "num_sv_be", info = info)
        expect_identical(estimate$seed, 900L + cell$seed, info = info)
        expect_identical(estimate$permutations, 20L, info = info)
        expect_identical(estimate$features_used, nrow(complete), info = info)
        expect_identical(estimate$features_total, nrow(fixture$x), info = info)
        # the diagnostic: one row per residual component, the count read from the same p-values
        components <- estimate$components
        residual <- min(nrow(complete), ncol(complete)) - length(unique(fixture$group))
        expect_identical(components$component, seq_len(residual), info = info)
        expect_equal(sum(components$observed_share), 1, tolerance = 1e-12, info = info)
        expect_identical(dim(estimate$null_shares), c(20L, residual), info = info)
        expect_equal(rowSums(estimate$null_shares), rep(1, 20L), tolerance = 1e-12, info = info)
        expect_equal(components$p_value, vapply(seq_len(residual), function(i) {
            mean(estimate$null_shares[, i] >= components$observed_share[[i]])
        }, numeric(1)), info = info)
        expect_identical(components$p_value_monotone, cummax(components$p_value), info = info)
        expect_identical(sum(components$p_value_monotone <= estimate$p_threshold), estimate$k, info = info)
        expect_identical(components$counted, components$component <= estimate$k, info = info)
        counts <- c(counts, estimate$k)
    }
    # the data sets cover no factor, one and several
    expect_true(all(c(0L, 1L, 2L) %in% counts))
})

test_that("the estimate is reproducible, leaves the caller's random stream alone and has a recorded default seed", {
    skip_if_not_installed("sva")
    fixture <- ruvCountFixture(11L)
    set.seed(42L)
    before <- .Random.seed
    first <- estimateRuvFactorCount(fixture$x, fixture$group)
    expect_identical(.Random.seed, before)
    expect_identical(estimateRuvFactorCount(fixture$x, fixture$group), first)
    expect_identical(first$seed, RUV_FACTOR_COUNT_SEED)
    # the groups may be a factor; the count is that of the same labels as text
    expect_identical(estimateRuvFactorCount(fixture$x, factor(fixture$group, levels = c("B", "A", "unused")))$k, first$k)
})

test_that("the estimate is refused with a classed condition and its reason", {
    fixture <- ruvCountFixture(12L, features = 60L)
    testthat::local_mocked_bindings(ruvFactorCountSvaAvailable = function() FALSE)
    absent <- expect_error(estimateRuvFactorCount(fixture$x, fixture$group), class = "ruv_factor_count_needs_sva")
    expect_s3_class(absent, "ruv_factor_count_unavailable")
    expect_match(conditionMessage(absent), "sva", fixed = TRUE)
    testthat::local_mocked_bindings(ruvFactorCountSvaAvailable = function() TRUE)
    # too few features complete in every sample (19 of 60)
    sparse <- fixture$x
    sparse[20:60, 1L] <- NA_real_
    few <- expect_error(estimateRuvFactorCount(sparse, fixture$group), class = "ruv_factor_count_too_few_features")
    expect_s3_class(few, "ruv_factor_count_unavailable")
    expect_match(conditionMessage(few), "19", fixed = TRUE)
    # three samples in two groups leave one residual component: nothing to compare
    small <- expect_error(estimateRuvFactorCount(fixture$x[, 1:3], fixture$group[1:3]), class = "ruv_factor_count_too_few_samples")
    expect_s3_class(small, "ruv_factor_count_unavailable")
    # wrong arguments are not refusals of the estimate
    expect_error(estimateRuvFactorCount(fixture$x, fixture$group[-1L]))
    expect_error(estimateRuvFactorCount(fixture$x, rep("A", ncol(fixture$x))))
})

test_that("the diagnostic plot shows the observed shares against the permutation null and marks the count", {
    skip_if_not_installed("sva")
    fixture <- ruvCountFixture(2L, features = 1000L, samples = 12L, factors = 2L)
    estimate <- estimateRuvFactorCount(fixture$x, fixture$group, seed = 902L)
    expect_identical(estimate$k, 2L)
    plot <- plotRuvFactorCount(estimate)
    expect_s3_class(plot, "ggplot")
    built <- ggplot2::ggplot_build(plot)
    layers <- unname(vapply(plot$layers, function(layer) class(layer$geom)[[1L]], character(1)))
    expect_identical(layers, c("GeomLine", "GeomLine", "GeomPoint", "GeomVline"))
    shown <- nrow(estimate$components)
    # the permutation null: one line per permutation over the components shown
    expect_identical(nrow(built$data[[1L]]), 20L * shown)
    expect_equal(sort(built$data[[1L]]$y), sort(as.vector(estimate$null_shares)))
    # the observed shares, and the count marked after its last component
    expect_equal(built$data[[2L]]$y, estimate$components$observed_share)
    expect_equal(built$data[[3L]]$y, estimate$components$observed_share)
    expect_equal(built$data[[4L]]$xintercept, 2.5)
    expect_match(plot$labels$subtitle, "2", fixed = TRUE)
    expect_match(plot$labels$caption, "902", fixed = TRUE)
    expect_identical(plot$labels$y, "Share of residual variance")
    # no factor: nothing is marked
    none <- ruvCountFixture(7L, features = 300L, factors = 0L)
    zero <- estimateRuvFactorCount(none$x, none$group, seed = 907L)
    expect_identical(zero$k, 0L)
    zero_layers <- vapply(plotRuvFactorCount(zero)$layers, function(layer) class(layer$geom)[[1L]], character(1))
    expect_false("GeomVline" %in% zero_layers)
    # many components: the first ones are shown, the counted ones always
    expect_identical(nrow(ggplot2::ggplot_build(plotRuvFactorCount(estimate, max_components = 4L))$data[[2L]]), 4L)
    # saved as today's k diagnostic is (the canonical-correlation plot: a grob in the product, a PNG by ggsave)
    expect_s3_class(ggplot2::ggplotGrob(stripProtNormQcPlotOwners(plot)), "gtable")
    path <- file.path(withr::local_tempdir(), "ruv_factor_count.png")
    ggplot2::ggsave(path, plot, width = 8, height = 6, dpi = 150)
    expect_gt(file.size(path), 1000)
})

test_that("the largest k a design can carry is the smallest of its three limits", {
    limit <- ruvFactorCountLimit(samples = 12L, replicate_sets = 2L, groups = 2L, controls = 500L)
    expect_identical(limit$k_max, 9L)
    expect_match(limit$binding, "residual degree of freedom")
    # declared technical replicates as replicate sets: 16 samples in 8 units and 2 groups
    expect_identical(ruvFactorCountLimit(16L, 8L, 2L, 500L)$k_max, 8L)
    expect_match(ruvFactorCountLimit(16L, 8L, 2L, 500L)$binding, "replicate sets")
    expect_identical(ruvFactorCountLimit(12L, 2L, 2L, 3L)$k_max, 2L)
    expect_match(ruvFactorCountLimit(12L, 2L, 2L, 3L)$binding, "controls")
    expect_identical(ruvFactorCountLimit(3L, 2L, 2L, 500L)$k_max, 0L)
})

test_that("the resolver returns k with its recorded source for the count, the separation rule and a number", {
    limit <- ruvFactorCountLimit(12L, 2L, 2L, 500L)
    estimate <- list(k = 2L, rule = "num_sv_be", seed = 1L)
    count <- function() estimate
    never <- function() stop("not to be called")
    # the estimated count, defined for every control rule
    for (control_rule in c("all_features", "supplied_list", "first_pass_drop_top_25")) {
        resolved <- resolveRuvFactorCount("num_sv_be", control_rule, limit, estimate_count = count, separation_rule = never)
        expect_identical(resolved[c("k", "k_rule", "k_estimate", "status")],
            list(k = 2L, k_rule = "num_sv_be", k_estimate = 2L, status = "ok"), info = control_rule)
        expect_identical(resolved$estimate, estimate, info = control_rule)
        expect_identical(resolved$control_rule, control_rule, info = control_rule)
    }
    # the separation rule, where a separate control set exists
    for (control_rule in c("supplied_list", "first_pass_drop_top_25", "sweep_least_significant")) {
        resolved <- resolveRuvFactorCount("sweep_bounded", control_rule, limit, estimate_count = never, separation_rule = function() 3)
        expect_identical(resolved[c("k", "k_rule", "k_estimate", "status")],
            list(k = 3L, k_rule = "sweep_bounded", k_estimate = 3L, status = "ok"), info = control_rule)
        expect_null(resolved$estimate)
    }
    # ... and never with all features as controls
    undefined <- expect_error(resolveRuvFactorCount("sweep_bounded", "all_features", limit, estimate_count = never, separation_rule = never),
        class = "ruv_k_separation_undefined")
    expect_s3_class(undefined, "ruv_k_refused")
    expect_match(conditionMessage(undefined), "all features as controls", fixed = TRUE)
    expect_error(resolveRuvFactorCount("sweep_bounded", "supplied_list", limit, separation_rule = function() NA_integer_),
        class = "ruv_k_not_estimated")
    # a number: the user's k is used and the estimate is recorded beside it
    chosen <- resolveRuvFactorCount("user", "all_features", limit, k_user = 4, estimate_count = count, separation_rule = never)
    expect_identical(chosen[c("k", "k_rule", "k_estimate", "status")], list(k = 4L, k_rule = "user", k_estimate = 2L, status = "ok"))
    expect_identical(chosen$estimate, estimate)
    expect_error(resolveRuvFactorCount("user", "all_features", limit, estimate_count = count))
    expect_error(resolveRuvFactorCount("user", "all_features", limit, k_user = 0, estimate_count = count))
    expect_error(resolveRuvFactorCount("user", "all_features", limit, k_user = 1.5, estimate_count = count))
    expect_error(resolveRuvFactorCount("elbow", "all_features", limit, estimate_count = count))
    expect_error(resolveRuvFactorCount("num_sv_be", "some", limit, estimate_count = count))
})

test_that("the resolver refuses a k the design cannot carry, reports no factor, and keeps a user's k without sva", {
    limit <- ruvFactorCountLimit(6L, 2L, 2L, 500L)
    expect_identical(limit$k_max, 3L)
    for (rule in c("num_sv_be", "user")) {
        above <- expect_error(resolveRuvFactorCount(rule, "all_features", limit, k_user = 4L, estimate_count = function() list(k = 4L)),
            class = "ruv_k_above_limit")
        expect_s3_class(above, "ruv_k_refused")
        expect_match(conditionMessage(above), "at most 3", fixed = TRUE)
    }
    expect_error(resolveRuvFactorCount("sweep_bounded", "supplied_list", limit, separation_rule = function() 4L), class = "ruv_k_above_limit")
    # an estimate of 0 is a result, not a refusal: no RUV for the run, with a notice
    zero <- NULL
    expect_message(zero <- resolveRuvFactorCount("num_sv_be", "all_features", limit, estimate_count = function() list(k = 0L)),
        class = "ruv_k_zero")
    expect_identical(zero[c("k", "k_rule", "k_estimate", "status")], list(k = 0L, k_rule = "num_sv_be", k_estimate = 0L, status = "none"))
    expect_match(zero$reason, "no unwanted factor", fixed = TRUE)
    # the estimate's refusals pass through under the automatic rule ...
    testthat::local_mocked_bindings(ruvFactorCountSvaAvailable = function() FALSE)
    fixture <- ruvCountFixture(13L, features = 60L, samples = 6L)
    count <- function() estimateRuvFactorCount(fixture$x, fixture$group)
    expect_error(resolveRuvFactorCount("num_sv_be", "all_features", limit, estimate_count = count), class = "ruv_factor_count_needs_sva")
    # ... and are recorded, not fatal, beside a user's k
    chosen <- resolveRuvFactorCount("user", "all_features", limit, k_user = 1L, estimate_count = count)
    expect_identical(chosen[c("k", "k_rule", "k_estimate", "status")], list(k = 1L, k_rule = "user", k_estimate = NA_integer_, status = "ok"))
    expect_match(chosen$k_estimate_reason, "sva", fixed = TRUE)
    expect_null(chosen$estimate)
})

# ---- Audit of 2026-11-18 (KIT-129, register F11, F12, F14, and the flag for F1) ----

test_that("a missing group and a single group are refusals of the estimate, with their class and reason (audit F11)", {
    fixture <- ruvCountFixture(14L, features = 60L)
    missing_group <- fixture$group
    missing_group[[1L]] <- NA
    missing <- expect_error(estimateRuvFactorCount(fixture$x, missing_group), class = "ruv_factor_count_group_missing")
    expect_s3_class(missing, "ruv_factor_count_unavailable")
    expect_match(conditionMessage(missing), "1 of 12", fixed = TRUE)
    single <- expect_error(estimateRuvFactorCount(fixture$x, rep("A", ncol(fixture$x))), class = "ruv_factor_count_one_group")
    expect_s3_class(single, "ruv_factor_count_unavailable")
    # a factor with one level in use is one group
    expect_error(estimateRuvFactorCount(fixture$x, factor(rep("A", ncol(fixture$x)), levels = c("A", "B"))),
        class = "ruv_factor_count_one_group")
})

test_that("a user's k stands whatever error the estimate raises; the reason is recorded (audit F11)", {
    limit <- ruvFactorCountLimit(12L, 2L, 2L, 500L)
    for (failure in list(\() stop("a plain error of the estimate"),
        \() cli::cli_abort("the diagnostic disagrees", class = "ruv_factor_count_inconsistent"),
        \() estimateRuvFactorCount(matrix(0, 30L, 4L), c("A", NA, "B", "B")))) {
        chosen <- resolveRuvFactorCount("user", "all_features", limit, k_user = 2L, estimate_count = failure)
        expect_identical(chosen[c("k", "k_rule", "k_estimate", "status")],
            list(k = 2L, k_rule = "user", k_estimate = NA_integer_, status = "ok"))
        expect_true(nzchar(chosen$k_estimate_reason))
        expect_null(chosen$estimate)
    }
    # under the automatic rule the same errors are not swallowed
    expect_error(resolveRuvFactorCount("num_sv_be", "all_features", limit, estimate_count = \() stop("a plain error of the estimate")),
        "a plain error")
})

test_that("the estimate draws under a fixed and recorded RNG kind, whatever the caller's (audit F12)", {
    skip_if_not_installed("sva")
    fixture <- ruvCountFixture(2L, features = 300L, samples = 12L, factors = 2L)
    reference <- estimateRuvFactorCount(fixture$x, fixture$group, seed = 902L)
    expect_identical(reference[c("rng_kind", "rng_normal_kind", "rng_sample_kind")],
        list(rng_kind = "Mersenne-Twister", rng_normal_kind = "Inversion", rng_sample_kind = "Rejection"))
    callers_kind <- RNGkind()
    withr::defer(do.call(RNGkind, as.list(callers_kind)))
    suppressWarnings(RNGkind("Wichmann-Hill", "Box-Muller", "Rounding"))
    set.seed(7L)
    before <- .Random.seed
    other <- estimateRuvFactorCount(fixture$x, fixture$group, seed = 902L)
    # the caller's kind and stream are as they were
    expect_identical(RNGkind(), c("Wichmann-Hill", "Box-Muller", "Rounding"))
    expect_identical(.Random.seed, before)
    expect_identical(other, reference)
})

test_that("the estimate says when the design declares technical replicates and the count removed the groups only (audit F1, flag only)", {
    skip_if_not_installed("sva")
    fixture <- ruvCountFixture(11L)
    plain <- estimateRuvFactorCount(fixture$x, fixture$group)
    expect_identical(plain[c("removed", "technical_replicates_ignored")], list(removed = "groups", technical_replicates_ignored = FALSE))
    flagged <- estimateRuvFactorCount(fixture$x, fixture$group, technical_replicates = TRUE)
    expect_identical(flagged[c("removed", "technical_replicates_ignored")], list(removed = "groups", technical_replicates_ignored = TRUE))
    # what is counted does not change
    expect_identical(flagged[c("k", "components", "null_shares")], plain[c("k", "components", "null_shares")])
})

test_that("the plot's legend: the permutation null is a line only, the counted marks in the observed colour (audit F14)", {
    skip_if_not_installed("sva")
    fixture <- ruvCountFixture(2L, features = 300L, samples = 12L, factors = 2L)
    plot <- plotRuvFactorCount(estimateRuvFactorCount(fixture$x, fixture$group, seed = 902L))
    built <- ggplot2::ggplot_build(plot)
    colour <- built$plot$guides$get_params("colour")
    shape <- built$plot$guides$get_params("shape")
    # the colour key: a point on Observed only
    expect_identical(colour$override.aes$shape, c(16, NA))
    # the shape key: the observed series' colour, no line through the marks
    expect_identical(shape$override.aes$colour, "#0072B2")
    expect_identical(shape$override.aes$linetype, 0)
})
