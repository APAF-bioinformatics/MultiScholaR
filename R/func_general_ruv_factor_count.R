# KIT-129 W6: the number of unwanted factors (k) of an RUV step, shared by every omic's new path (the user's
# direction of 2026-11-18: two rules, each offered where it is defined).
#   The estimated count: `sva::num.sv(method = "be", B = 200, seed)` on the features complete in every sample, as
#   KIT-124's exactness round calls it (arms 31 to 33; decision M30: 200 permutations, where arms 7 to 16 used 20). `sva` is an optional dependency; without it the estimate
#   is refused, never replaced.
#   The separation rule: today's estimate on the canonical-correlation plot (`findBestKElbow`, `trxFindBestK`, the
#   sweeps' `best_k`). It is not changed here; an omic hands it to the resolver as a function.
#   A number: the user's k.
# `resolveRuvFactorCount` returns k with its recorded source. Nothing here is called by today's default paths.

#' The seed of the estimated count when the caller gives none: fixed, and recorded with every estimate.
#' @noRd
RUV_FACTOR_COUNT_SEED <- 20261118L

#' The RNG kind the estimated count draws under (R's defaults since 3.6.0), whatever the caller's: fixed, and
#' recorded with every estimate.
#' @noRd
RUV_FACTOR_COUNT_RNG <- c(kind = "Mersenne-Twister", normal_kind = "Inversion", sample_kind = "Rejection")

#' The k rules and the control rules, by the names the analysis record uses (KIT-129 section 1e).
#' @noRd
RUV_K_RULES <- c("num_sv_be", "sweep_bounded", "user")

#' @noRd
RUV_CONTROL_RULES <- c("all_features", "supplied_list", "first_pass_drop_top_25", "sweep_least_significant")

#' Whether `sva` can be loaded (a function of its own so a test can stand in for a library without it).
#' @noRd
ruvFactorCountSvaAvailable <- function() requireNamespace("sva", quietly = TRUE)

#' The shares of residual variance per component, observed and under row-wise permutation: the quantities
#' `sva::num.sv(method = "be")` computes and does not return, by its own lines (sva 3.60.0) so that the same seed
#' draws the same permutations.
#' @param dat Complete feature-by-sample matrix.
#' @param mod The model matrix whose column space is removed.
#' @param permutations,seed As `sva::num.sv`'s `B` and `seed`.
#' @return `list(observed, null)`: a vector over the residual components and a permutation-by-component matrix.
#' @noRd
permuteRuvResidualShares <- function(dat, mod, permutations, seed) {
    set.seed(seed)
    hat <- mod %*% solve(t(mod) %*% mod) %*% t(mod)
    res <- dat - t(hat %*% t(dat))
    components <- seq_len(min(dim(dat)) - ceiling(sum(diag(hat))))
    shares <- function(residuals) {
        d <- svd(residuals)$d[components]
        d^2 / sum(d^2)
    }
    observed <- shares(res)
    null <- matrix(0, nrow = permutations, ncol = length(components))
    for (i in seq_len(permutations)) {
        res0 <- t(apply(res, 1, sample, replace = FALSE))
        null[i, ] <- shares(res0 - t(hat %*% t(res0)))
    }
    list(observed = observed, null = null)
}

#' The estimated number of unwanted factors (KIT-129 W6)
#'
#' The permutation estimate of Buja and Eyuboglu (1992) as `sva::num.sv(dat, model.matrix(~ group), method = "be",
#' B = 200, seed)` computes it, on the features observed in every sample: the groups are removed, and a component of
#' the residuals is counted while at most 10% of the row-wise permutations reach its share of the residual
#' variance. The count is `sva`'s own return value; the diagnostic's quantities are computed beside it on the same
#' random stream and must give the same count. The draws are made under a fixed RNG kind (`RUV_FACTOR_COUNT_RNG`,
#' recorded in the result), so the count does not depend on the caller's `RNGkind()`; the caller's kind and random
#' stream are left as they were.
#'
#' Refused with a condition of class `ruv_factor_count_unavailable` and one of `ruv_factor_count_group_missing`
#' (a sample without a group), `ruv_factor_count_one_group` (fewer than two groups), `ruv_factor_count_needs_sva`
#' (`sva` is not installed), `ruv_factor_count_too_few_features` (fewer than `min_features` features observed in
#' every sample) or `ruv_factor_count_too_few_samples` (fewer than two residual components beside the groups).
#'
#' @param values Numeric feature-by-sample matrix; missing values allowed.
#' @param groups The group of each column of `values` (the levels removed before the count).
#' @param seed Whole number; recorded in the result.
#' @param permutations Number of permutations (`num.sv`'s `B`): 200 (decision M30; the exactness round's arms 31 to 33).
#' @param min_features The fewest complete features the estimate is made on (20, the rule of simulator S4).
#' @param technical_replicates Whether the design declares technical replicates. Without `units` the count removes
#'   the groups whatever its value, and it is recorded as `technical_replicates_ignored` so that a caller can warn.
#' @param units The biological unit of each column (the declared technical-replicate id) on a design with technical
#'   replicates, or NULL. Given, the count removes the units, not only the groups (decision M28; the exactness
#'   round's arm 29: `sva::num.sv(dat, model.matrix(~ unit), method = "be", B, seed)`), and `removed` says so.
#' @return `list(k, rule = "num_sv_be", seed, permutations, p_threshold, features_used, features_total, samples,
#'   groups, components, null_shares, sva_version, rng_kind, rng_normal_kind, rng_sample_kind, removed,
#'   technical_replicates_ignored)`: `components` has one row per residual component (`component`,
#'   `observed_share`, `null_mean`, `p_value`, `p_value_monotone`, `counted`); `null_shares` is the
#'   permutation-by-component matrix of shares; `removed` names what the count removed before counting (`"groups"`,
#'   or `"biological units"`).
#' @noRd
estimateRuvFactorCount <- function(values, groups, seed = RUV_FACTOR_COUNT_SEED, permutations = 200L,
                                   min_features = 20L, technical_replicates = FALSE, units = NULL) {
    checkmate::assert_matrix(values, mode = "numeric", min.cols = 2L)
    checkmate::assert_atomic_vector(groups, len = ncol(values))
    seed <- checkmate::assert_int(seed, coerce = TRUE)
    permutations <- checkmate::assert_count(permutations, positive = TRUE, coerce = TRUE)
    checkmate::assert_count(min_features, positive = TRUE)
    checkmate::assert_flag(technical_replicates)
    if (!is.null(units)) checkmate::assert_atomic_vector(units, len = ncol(values), any.missing = FALSE)
    refuse <- function(reason, class) {
        cli::cli_abort("the number of unwanted factors cannot be estimated: {reason}",
            class = c(class, "ruv_factor_count_unavailable"), call = NULL)
    }
    if (anyNA(groups)) {
        refuse(sprintf("%d of %d samples have no group", sum(is.na(groups)), length(groups)), "ruv_factor_count_group_missing")
    }
    group <- if (is.factor(groups)) droplevels(groups) else as.character(groups)
    n_groups <- length(unique(group))
    if (n_groups < 2L) {
        refuse(sprintf("the samples are in %d group; the count removes the differences between two groups or more", n_groups),
            "ruv_factor_count_one_group")
    }
    if (!ruvFactorCountSvaAvailable()) {
        refuse("the package sva is not installed (the estimate is sva::num.sv); install it, or give the number of factors",
            "ruv_factor_count_needs_sva")
    }
    dat <- unname(values[rowSums(!is.finite(values)) == 0L, , drop = FALSE])
    if (nrow(dat) < min_features) {
        refuse(sprintf("%d of %d features are observed in every sample; at least %d are needed", nrow(dat), nrow(values),
            min_features), "ruv_factor_count_too_few_features")
    }
    removed_columns <- if (is.null(units)) n_groups else length(unique(as.character(units)))
    if (min(dim(dat)) - removed_columns < 2L) {
        refuse(sprintf("%d samples in %d %s leave fewer than two residual components to compare", ncol(dat), removed_columns,
            if (is.null(units)) "groups" else "biological units"), "ruv_factor_count_too_few_samples")
    }
    mod <- if (is.null(units)) stats::model.matrix(~ group, data.frame(group = group)) else
        stats::model.matrix(~ unit, data.frame(unit = as.character(units)))
    # Both calls seed themselves with `seed`; the kind they draw under is fixed here, and withr restores the
    # caller's kind and stream.
    counted <- withr::with_seed(seed, list(
        k = as.integer(sva::num.sv(dat, mod, method = "be", B = permutations, seed = seed)),
        shares = permuteRuvResidualShares(dat, mod, permutations, seed)
    ), .rng_kind = RUV_FACTOR_COUNT_RNG[["kind"]], .rng_normal_kind = RUV_FACTOR_COUNT_RNG[["normal_kind"]],
        .rng_sample_kind = RUV_FACTOR_COUNT_RNG[["sample_kind"]])
    k <- counted$k
    observed <- counted$shares$observed
    null <- counted$shares$null
    p_value <- vapply(seq_along(observed), \(i) mean(null[, i] >= observed[[i]]), numeric(1))
    p_threshold <- 0.1 # fixed inside `sva::num.sv`
    if (!identical(sum(cummax(p_value) <= p_threshold), k)) {
        cli::cli_abort(c("the diagnostic of the estimated count disagrees with sva::num.sv ({sum(cummax(p_value) <= p_threshold)} against {k})",
            i = "sva {utils::packageVersion('sva')} may compute the estimate differently from the version this follows (3.60.0)"),
            class = "ruv_factor_count_inconsistent", call = NULL)
    }
    components <- data.frame(component = seq_along(observed), observed_share = observed, null_mean = colMeans(null),
        p_value = p_value, p_value_monotone = cummax(p_value), counted = seq_along(observed) <= k)
    list(k = k, rule = "num_sv_be", seed = seed, permutations = permutations, p_threshold = p_threshold,
        features_used = nrow(dat), features_total = nrow(values), samples = ncol(dat), groups = n_groups,
        components = components, null_shares = null, sva_version = as.character(utils::packageVersion("sva")),
        rng_kind = RUV_FACTOR_COUNT_RNG[["kind"]], rng_normal_kind = RUV_FACTOR_COUNT_RNG[["normal_kind"]],
        rng_sample_kind = RUV_FACTOR_COUNT_RNG[["sample_kind"]], removed = if (is.null(units)) "groups" else "biological units",
        technical_replicates_ignored = technical_replicates && is.null(units))
}

#' The diagnostic plot of the estimated count (KIT-129 W6; user requirement of 2026-11-18)
#'
#' The observed share of residual variance per component against the shares under permutation, with the count
#' marked: the counterpart, for the estimated count, of the canonical-correlation plot of the separation rule.
#' @param estimate The value of `estimateRuvFactorCount()`.
#' @param max_components The components shown; the counted ones and the two after them are always shown.
#' @return A ggplot object.
#' @noRd
plotRuvFactorCount <- function(estimate, max_components = 15L) {
    checkmate::assert_list(estimate)
    checkmate::assert_names(names(estimate), must.include = c("k", "seed", "permutations", "p_threshold",
        "features_used", "components", "null_shares"))
    checkmate::assert_count(max_components, positive = TRUE)
    shown <- seq_len(min(nrow(estimate$components), max(max_components, estimate$k + 2L)))
    observed <- estimate$components[shown, , drop = FALSE]
    observed$status <- factor(ifelse(observed$counted, "Counted", "Not counted"), levels = c("Counted", "Not counted"))
    null <- data.frame(component = rep(shown, each = estimate$permutations),
        permutation = rep(seq_len(estimate$permutations), times = length(shown)),
        share = as.vector(estimate$null_shares[, shown, drop = FALSE]))
    # The series names are literals in the mappings, so the plot holds no reference to this function's frame.
    plot <- ggplot2::ggplot() +
        ggplot2::geom_line(data = null, ggplot2::aes(x = .data$component, y = .data$share, group = .data$permutation,
            colour = "Permuted (null)"), linewidth = 0.3, alpha = 0.6) +
        ggplot2::geom_line(data = observed, ggplot2::aes(x = .data$component, y = .data$observed_share,
            colour = "Observed"), linewidth = 0.8) +
        ggplot2::geom_point(data = observed, ggplot2::aes(x = .data$component, y = .data$observed_share,
            colour = "Observed", shape = .data$status), size = 2.5, show.legend = TRUE)
    # The count is marked after its last component; with no factor counted there is nothing to mark.
    if (estimate$k > 0L) {
        plot <- plot + ggplot2::geom_vline(xintercept = estimate$k + 0.5, linetype = "dashed", colour = "grey30")
    }
    plot +
        ggplot2::scale_colour_manual(name = NULL, values = c(Observed = "#0072B2", `Permuted (null)` = "#D55E00")) +
        ggplot2::scale_shape_manual(name = NULL, values = c(Counted = 16, `Not counted` = 1), drop = FALSE) +
        # Audit F14: the permutation null's key is a line without a point; the counted marks' key is in the
        # observed series' colour, without a line through it.
        ggplot2::guides(colour = ggplot2::guide_legend(order = 1L, override.aes = list(shape = c(16, NA))),
            shape = ggplot2::guide_legend(order = 2L, override.aes = list(colour = "#0072B2", linetype = 0))) +
        ggplot2::scale_x_continuous(breaks = shown) +
        ggplot2::labs(x = "Component of the residuals (groups removed)", y = "Share of residual variance",
            title = "Number of unwanted factors: estimated count",
            subtitle = sprintf("Estimated number of factors: %d", estimate$k),
            caption = sprintf(paste0("sva::num.sv, method \"be\": a component is counted while at most %s%% of %d ",
                "permutations reach its share.\n%d features observed in every sample; seed %d."),
                format(100 * estimate$p_threshold), estimate$permutations, estimate$features_used, estimate$seed)) +
        ggplot2::theme_bw() +
        ggplot2::theme(legend.position = "bottom")
}

#' The estimated count recorded beside a k the user chose (KIT-129 integration item 3; section 1e's `ruv_k_estimate`,
#' "what the estimator returned when the user chose"). The run path calls it once per RUV step under the k rule
#' `user`; it never changes k. An error of the estimate (no `sva`, too few complete features or samples) is kept as
#' its reason.
#' @param values Feature-by-sample matrix the count is computed on (the step's own: the normalised matrix, the TMM
#'   logCPM, the assay's values).
#' @param groups The groups, one per column of `values`.
#' @param seed The count's seed (recorded with it).
#' @param units The biological units on a design with declared technical replicates (decision M28), else NULL.
#' @return `list(k_estimate, k_estimate_reason, estimate, seed)`: the count or `NA` with the reason, its record for
#'   `plotRuvFactorCount` (NULL on an error) and the seed.
#' @noRd
recordRuvFactorCountBesideUserK <- function(values, groups, seed = RUV_FACTOR_COUNT_SEED, units = NULL) {
    estimate <- tryCatch(estimateRuvFactorCount(values, groups, seed = seed, technical_replicates = !is.null(units), units = units),
        error = \(e) e)
    if (inherits(estimate, "condition")) {
        return(list(k_estimate = NA_integer_, k_estimate_reason = conditionMessage(estimate), estimate = NULL, seed = as.integer(seed)))
    }
    list(k_estimate = as.integer(estimate$k), k_estimate_reason = NA_character_, estimate = estimate, seed = as.integer(seed))
}

#' The largest k a design can carry (KIT-129 section 2, "Details"): k cannot exceed the samples minus the replicate
#' sets (the rank RUV-III's replicates leave), must leave a residual degree of freedom in a model of the groups and
#' the factors, and must be below the number of controls observed in every sample.
#' On a design with declared technical replicates (decision M28) the units, not the injections, carry the test: k must
#' also leave a residual degree of freedom at the level of the biological units, `units - between_unit_rank - k >= 1`.
#' @param samples,replicate_sets,groups,controls Counts.
#' @param units The number of biological units when the design declares technical replicates; NULL otherwise (the
#'   limit is then as before).
#' @param between_unit_rank The rank of the design between units (the groups and any covariate between units).
#' @return `list(k_max, binding, samples, replicate_sets, groups, controls, units, between_unit_rank, k_max_units)`;
#'   `binding` names the limit that binds; `k_max_units` the units' limit alone (NA without technical replicates).
#' @noRd
ruvFactorCountLimit <- function(samples, replicate_sets, groups, controls, units = NULL, between_unit_rank = groups) {
    limits <- c(`samples minus replicate sets` = samples - replicate_sets,
        `a residual degree of freedom beside the groups` = samples - groups - 1L,
        `fewer factors than controls observed in every sample` = controls - 1L)
    k_max_units <- NA_integer_
    if (!is.null(units)) {
        k_max_units <- as.integer(max(units - between_unit_rank - 1L, 0L))
        limits <- c(limits, `a residual degree of freedom between the biological units (the design has technical replicates)` =
            units - between_unit_rank - 1L)
    }
    list(k_max = as.integer(max(min(limits), 0L)), binding = names(limits)[[which.min(limits)]],
        samples = as.integer(samples), replicate_sets = as.integer(replicate_sets), groups = as.integer(groups),
        controls = as.integer(controls), units = if (is.null(units)) NA_integer_ else as.integer(units),
        between_unit_rank = as.integer(between_unit_rank), k_max_units = k_max_units)
}

#' Whether RUV-III's replicate sets are the biological units of a design with declared technical replicates (KIT-129
#' integration, technical-replicate design; arms 15 and 16 of the simulators): the two vectors partition the samples
#' the same way.
#' @param replicate_sets The replicate set of each sample.
#' @param units The biological unit of each sample, or NULL (no declared technical replicates).
#' @return TRUE or FALSE.
#' @noRd
ruvReplicateSetsAreUnits <- function(replicate_sets, units) {
    if (is.null(units)) return(FALSE)
    sets <- as.character(replicate_sets)
    units <- as.character(units)
    n_units <- length(unique(units))
    length(unique(sets)) == n_units && length(unique(paste(sets, units, sep = "\r"))) == n_units
}

#' The residual degrees of freedom left between biological units by a model (decision M28): the units minus the rank
#' of the model's unit means. NA when the design declares no technical replicates.
#' @param model A model matrix over the samples (rows), the formula's terms and any W columns.
#' @param units The biological unit of each row, or NULL.
#' @return An integer, or NA.
#' @noRd
ruvUnitResidualDf <- function(model, units) {
    if (is.null(units)) return(NA_integer_)
    units <- as.character(units)
    totals <- rowsum(model, units)
    unit_means <- totals / as.vector(table(units)[rownames(totals)])
    as.integer(nrow(unit_means) - qr(unit_means)$rank)
}

#' The number of unwanted factors and its recorded source (KIT-129 W6)
#'
#' One resolver for every omic: given the control rule of the chosen design and the user's choice of k rule, it
#' returns k. The estimated count (`"num_sv_be"`) is defined for every control rule. The separation rule
#' (`"sweep_bounded"`, today's estimate) reads the gap between the canonical correlations of all features and of
#' the controls, so it is refused when the controls are all features. A number (`"user"`) is used as given, with
#' the estimated count recorded beside it when it can be computed.
#'
#' Refused with a condition of class `ruv_k_refused` and one of `ruv_k_separation_undefined`, `ruv_k_not_estimated`
#' (the separation rule returned none) or `ruv_k_above_limit` (never reduced); the refusals of
#' `estimateRuvFactorCount()` pass through under `"num_sv_be"`; under `"user"` no error of the estimate stops the
#' run (it is recorded in `k_estimate_reason`). A k of 0 is a result: no RUV for the run, with a
#' message of class `ruv_k_zero`.
#'
#' @param k_rule One of `RUV_K_RULES`.
#' @param control_rule One of `RUV_CONTROL_RULES`.
#' @param limit The value of `ruvFactorCountLimit()`.
#' @param k_user The user's k (a whole number of 1 or more), for `"user"`.
#' @param estimate_count Function of no argument returning `estimateRuvFactorCount()`'s value; required for
#'   `"num_sv_be"`, optional for `"user"`.
#' @param separation_rule Function of no argument returning the omic's separation-rule k; for `"sweep_bounded"`.
#' @return `list(k, k_rule, k_estimate, k_estimate_reason, k_max, control_rule, status, reason, estimate)`:
#'   `status` is `"ok"` or `"none"` (k = 0, with `reason`); `k_estimate` what the rule returned (under `"user"`,
#'   the estimated count, or `NA` with `k_estimate_reason`); `estimate` the count's record for the diagnostic.
#' @noRd
resolveRuvFactorCount <- function(k_rule, control_rule, limit, k_user = NULL, estimate_count = NULL,
                                  separation_rule = NULL) {
    checkmate::assert_choice(k_rule, RUV_K_RULES)
    checkmate::assert_choice(control_rule, RUV_CONTROL_RULES)
    checkmate::assert_list(limit)
    checkmate::assert_int(limit$k_max, lower = 0)
    checkmate::assert_function(estimate_count, null.ok = !identical(k_rule, "num_sv_be"))
    checkmate::assert_function(separation_rule, null.ok = !identical(k_rule, "sweep_bounded"))
    refuse <- function(reason, class) {
        cli::cli_abort("the number of unwanted factors (k) is refused: {reason}", class = c(class, "ruv_k_refused"),
            call = NULL)
    }
    estimate <- NULL
    k_estimate <- NA_integer_
    k_estimate_reason <- NA_character_
    if (identical(k_rule, "sweep_bounded")) {
        if (identical(control_rule, "all_features")) {
            refuse(paste("the separation rule reads the gap between the canonical correlations of all features and of",
                "the controls, and with all features as controls the two curves are the same;",
                "choose the estimated count or give a number"), "ruv_k_separation_undefined")
        }
        k <- separation_rule()
        if (!checkmate::test_int(k, lower = 0)) refuse("the separation rule returned no k", "ruv_k_not_estimated")
        k <- k_estimate <- as.integer(k)
    } else if (identical(k_rule, "num_sv_be")) {
        estimate <- estimate_count()
        k <- k_estimate <- as.integer(estimate$k)
    } else {
        k <- checkmate::assert_int(k_user, lower = 1, coerce = TRUE)
        if (!is.null(estimate_count)) {
            # The user's k stands whether or not the count can be estimated; the reason is kept with the record.
            # Audit F11: an error of any kind, classed as a refusal of the estimate or not.
            estimate <- tryCatch(estimate_count(), error = \(e) e)
            if (inherits(estimate, "condition")) {
                k_estimate_reason <- conditionMessage(estimate)
                estimate <- NULL
            } else {
                k_estimate <- as.integer(estimate$k)
            }
        }
    }
    if (k > limit$k_max) {
        refuse(sprintf("k = %d and at most %d can be estimated (%s: %d samples, %d replicate sets, %d groups, %d controls)",
            k, limit$k_max, limit$binding, limit$samples, limit$replicate_sets, limit$groups, limit$controls),
            "ruv_k_above_limit")
    }
    reason <- NA_character_
    if (k == 0L) {
        reason <- if (identical(k_rule, "num_sv_be")) {
            "the estimated count found no unwanted factor (no component above the permutation null)"
        } else {
            "the separation rule returned no unwanted factor"
        }
        cli::cli_inform("RUV is not applied: {reason}", class = "ruv_k_zero")
    }
    list(k = k, k_rule = k_rule, k_estimate = k_estimate, k_estimate_reason = k_estimate_reason, k_max = limit$k_max,
        control_rule = control_rule, status = if (k > 0L) "ok" else "none", reason = reason, estimate = estimate)
}

#' Each comparison fitted as a coefficient: its exact standard error per feature (decision M29)
#'
#' `limma::contrasts.fit` forms a comparison's standard error from one correlation matrix of the coefficients for
#' every feature. With RUV factors in the model the design is not orthogonal, so for a feature with missing values (or
#' per-observation weights) that is an approximation. Here each comparison is fitted directly, as the exactness round's
#' arm 25 did: the model is reparameterised so that the comparison is its first coefficient
#' (`limma::contrastAsCoef`), refitted, and that coefficient and its unscaled standard deviation replace
#' contrasts.fit's. The residual variances and their degrees of freedom do not depend on the parameterisation and are
#' contrasts.fit's (the fit's own). A comparison contrasts.fit leaves missing for a feature (not estimable on the
#' feature's observed samples) stays missing.
#' @param fit The uncontrasted `MArrayLM`.
#' @param contrasts A limma contrast matrix over the columns of `design`.
#' @param design The model matrix `fit` was fitted with.
#' @param refit Function of a model matrix returning the `MArrayLM` of the same data with the same block,
#'   correlation and weights.
#' @return The value of `limma::contrasts.fit(fit, contrasts)` with the exact coefficients and unscaled standard
#'   deviations.
#' @noRd
ruvExactContrastsFit <- function(fit, contrasts, design, refit) {
    contrasted <- limma::contrasts.fit(fit, contrasts)
    for (j in seq_len(ncol(contrasts))) {
        exact <- refit(limma::contrastAsCoef(design, contrasts[, j, drop = FALSE], first = TRUE)$design)
        if (!identical(dimnames(exact$coefficients)[[1L]], dimnames(contrasted$coefficients)[[1L]]) ||
            nrow(exact$coefficients) != nrow(contrasted$coefficients)) {
            cli::cli_abort("the refitted model does not hold the features of the fit", class = "ruv_exact_contrasts_invalid")
        }
        estimable <- !is.na(contrasted$coefficients[, j])
        contrasted$coefficients[estimable, j] <- exact$coefficients[estimable, 1L]
        contrasted$stdev.unscaled[estimable, j] <- exact$stdev.unscaled[estimable, 1L]
    }
    contrasted
}

# ------------------------------------------------------------------------------------------------
# KIT-129 fourth pass item 4: the warning that an estimated unwanted factor follows the groups (decisions M11, M23,
# M24, M27). The statistic is KIT-129 section 3's (the multiple correlation of each estimated factor with the groups,
# as `docs/qa/workflow-execution/kit/kit-129/r0/association.R` on branch `sim-k129-r0b` computes it); the thresholds
# are the calibration's, validated on runs never used to set them, for the design and k rule each was measured with.
# ------------------------------------------------------------------------------------------------

#' The validated thresholds of the factor-tied-to-groups warning (decision M27), one row per omic, design and k rule.
#' A combination without a row has no validated threshold: its value is recorded and shown without a verdict. The
#' threshold fires when the largest association is strictly above it. `false_alarms` and `detection` are the shares
#' fired on the validation runs at correlation 0 (C0, C3, C4) and at 0.3 or more (C5, S11, S3).
#' @return A data frame: `omic` (`proteomics` is classic proteomics), `ruv_design`, `ruv_k_rule`, `threshold`,
#'   `simulator`, `arm`, `false_alarms`, `detection`, `rule` (`per_design`, M23; `follow_up`, M24), `source`,
#'   `setting_run`, `validation_run`.
#' @noRd
ruvGroupAssociationThresholds <- function() {
    r0 <- "docs/qa/workflow-execution/kit/kit-129/r0/"
    row <- function(omic, ruv_design, threshold, simulator, arm, false_alarms, detection, rule, source, validation_run) {
        data.frame(omic = omic, ruv_design = ruv_design, ruv_k_rule = if (identical(arm, 5L)) "sweep_bounded" else "num_sv_be",
            threshold = threshold, simulator = simulator, arm = arm, false_alarms = false_alarms, detection = detection, rule = rule,
            source = paste0(r0, source), setting_run = "44203265622", validation_run = validation_run, stringsAsFactors = FALSE)
    }
    rbind(
        row("proteomics", "groups_all_features", 0.058077853, "S1", 7L, 0.051, 1.000, "follow_up", "threshold-followup-S1.tsv", "44238222543"),
        row("proteomics", "supplied_controls", 0.070068989, "S1", 8L, 0.043, 1.000, "follow_up", "threshold-followup-S1.tsv", "44238222543"),
        row("proteomics", "first_pass_nulls", 0.16376377, "S1", 5L, 0.047, 0.989, "follow_up", "threshold-followup-S1.tsv", "44238222543"),
        row("transcriptomics", "groups_all_features", 0.038965628, "S3", 7L, 0.091, 1.000, "per_design", "threshold-by-design-S3.tsv",
            "44227129836"),
        row("lipidomics", "groups_all_features", 0.088611623, "S4", 7L, 0.095, 0.929, "per_design", "threshold-by-design-S4.tsv",
            "44227129836"),
        row("metabolomics", "groups_all_features", 0.088611623, "S4", 7L, 0.095, 0.929, "per_design", "threshold-by-design-S4.tsv",
            "44227129836"),
        row("lipidomics", "supplied_controls", 0.19494193, "S4", 8L, 0.060, 0.915, "follow_up", "threshold-followup-S4.tsv", "44238222543"),
        row("metabolomics", "supplied_controls", 0.19494193, "S4", 8L, 0.060, 0.915, "follow_up", "threshold-followup-S4.tsv", "44238222543"))
}

#' The association of each estimated factor with the groups: the square root of R squared of `lm(W_j ~ group)`.
#' @param W The factors, samples in rows (a matrix or a data frame).
#' @param group The group of each row.
#' @return One number per factor.
#' @noRd
ruvFactorGroupAssociations <- function(W, group) {
    W <- as.matrix(W)
    stopifnot(nrow(W) == length(group), !anyNA(W))
    group <- factor(group)
    vapply(seq_len(ncol(W)), function(j) sqrt(summary(stats::lm(W[, j] ~ group))$r.squared), numeric(1L))
}

#' Each approach's measured consequences of an unwanted factor associated with the groups (decision M36): the mean
#' false discovery proportion at a 5% target (`mean_fdp`, subset `full`) in KIT-124's cells with association 0 (C0),
#' 0.3 (C5) and 0.6 (stress S3), the approach's own arm, and without correction (arm 1) at 0.6. One row per omic,
#' design and k rule: the arm each shipped combination reproduces (7 recommended, 8 supplied list, with the estimated
#' number; first-pass nulls with the separation rule: arm 5, in transcriptomics arm 24). NA where the cell was not
#' measured for that arm (never another approach's figure).
#' @return A data frame: `omic`, `ruv_design`, `ruv_k_rule`, `simulator`, `arm`, `fdp_0`, `fdp_03`, `fdp_06`,
#'   `uncorrected_06`, `source`.
#' @noRd
ruvGroupAssociationConsequences <- function() {
    scored <- "docs/qa/workflow-execution/kit/kit-124/scored/"
    row <- function(omic, ruv_design, ruv_k_rule, simulator, arm, fdp_0, fdp_03, fdp_06, uncorrected_06, source) {
        data.frame(omic = omic, ruv_design = ruv_design, ruv_k_rule = ruv_k_rule, simulator = simulator, arm = arm, fdp_0 = fdp_0,
            fdp_03 = fdp_03, fdp_06 = fdp_06, uncorrected_06 = uncorrected_06,
            source = paste0(scored, source, if (!identical(source, sprintf("final-20261118/%s/k124-cells.tsv", tolower(simulator))))
                sprintf("; arm 1: %sfinal-20261118/%s/k124-cells.tsv", scored, tolower(simulator))), stringsAsFactors = FALSE)
    }
    final <- function(simulator) sprintf("final-20261118/%s/k124-cells.tsv", tolower(simulator))
    lipid_metab <- function(omic) rbind(
        row(omic, "groups_all_features", "num_sv_be", "S4", 7L, 0.0574542, 0.0723922, 0.130459, 0.641749, final("S4")),
        row(omic, "supplied_controls", "num_sv_be", "S4", 8L, 0.0500134, 0.0676851, 0.120356, 0.641749, final("S4")),
        row(omic, "first_pass_nulls", "sweep_bounded", "S4", 5L, 0.0974275, 0.238564, 0.537242, 0.641749, final("S4")))
    rbind(
        row("proteomics", "groups_all_features", "num_sv_be", "S1", 7L, 0.0535118, 0.0560489, 0.0616284, 0.593446, final("S1")),
        row("proteomics", "supplied_controls", "num_sv_be", "S1", 8L, 0.0506339, 0.0530153, 0.0584843, 0.593446, final("S1")),
        row("proteomics", "first_pass_nulls", "sweep_bounded", "S1", 5L, 0.0515094, 0.077624, 0.272569, 0.593446, final("S1")),
        row("dia_limpa", "groups_all_features", "num_sv_be", "S2", 7L, 0.0514045, 0.0525753, 0.0651367, 0.532737, final("S2")),
        row("dia_limpa", "supplied_controls", "num_sv_be", "S2", 8L, 0.0514643, 0.0542155, 0.0727574, 0.532737, final("S2")),
        row("dia_limpa", "first_pass_nulls", "sweep_bounded", "S2", 5L, 0.0935895, 0.149448, 0.300619, 0.532737, final("S2")),
        row("transcriptomics", "groups_all_features", "num_sv_be", "S3", 7L, 0.048062, 0.0480989, 0.0524544, 0.450032, final("S3")),
        row("transcriptomics", "supplied_controls", "num_sv_be", "S3", 8L, 0.0484711, 0.0475697, 0.0479324, 0.450032, final("S3")),
        # arm 24 (the separation rule on the first pass's own controls) was not run in stress S3
        row("transcriptomics", "first_pass_nulls", "sweep_bounded", "S3", 24L, 0.0467308, 0.0643636, NA_real_, 0.450032,
            "newarms-gha/k124-cells.tsv"),
        lipid_metab("lipidomics"), lipid_metab("metabolomics"))
}

#' The check of one run's factors (decisions M27, M35, M36): the largest association, and against the detection limit
#' validated for the omic, design and k rule an outcome (`detectable` above it, `not_detectable` at or below it), or
#' `no_limit` where none was validated. With a typed k (M35) the estimated number's limit applies only when the typed k
#' equals the estimated number recorded beside it (the calibrated calculation); otherwise there is no limit. With an
#' outcome against a limit the approach's measured consequences (`ruvGroupAssociationConsequences`) are recorded.
#' Section 1e's record fields.
#' @param W The run's estimated factors (samples in rows), or NULL (no factor).
#' @param group The group of each row of `W`.
#' @param omic `proteomics` (classic), `dia_limpa`, `transcriptomics`, `lipidomics` or `metabolomics`.
#' @param ruv_design,ruv_k_rule The run's design and k rule.
#' @param k,k_estimate The k applied and, under a typed k, the estimated number recorded beside it (NA when none).
#' @return `list(ruv_group_association, ruv_group_association_factors, ruv_group_association_threshold,
#'   ruv_group_association_flag, ruv_group_association_outcome, ruv_group_association_source,
#'   ruv_group_association_basis, ruv_group_association_fdp_0, _fdp_03, _fdp_06, _uncorrected_06)`; the outcome is
#'   `detectable`, `not_detectable`, `no_limit` (the value without a verdict; the flag NA), `no_factors`, or
#'   `not_computed` (the factors' samples not found in the design). The basis names the k whose limit applied
#'   (`estimated_number`, `separation_rule`, `typed_k_equals_estimated_number`) or, with no limit under a typed k,
#'   `typed_k_differs_from_estimated_number` / `typed_k_without_estimated_number`.
#' @noRd
ruvGroupAssociationCheck <- function(W, group, omic, ruv_design, ruv_k_rule, k = NA_integer_, k_estimate = NA_integer_) {
    out <- list(ruv_group_association = NA_real_, ruv_group_association_factors = numeric(), ruv_group_association_threshold = NA_real_,
        ruv_group_association_flag = NA, ruv_group_association_outcome = "no_factors", ruv_group_association_source = NA_character_,
        ruv_group_association_basis = NA_character_, ruv_group_association_fdp_0 = NA_real_, ruv_group_association_fdp_03 = NA_real_,
        ruv_group_association_fdp_06 = NA_real_, ruv_group_association_uncorrected_06 = NA_real_)
    if (is.null(W) || !NCOL(W)) return(out)
    if (length(group) != NROW(W) || anyNA(group)) {
        out$ruv_group_association_outcome <- "not_computed"
        return(out)
    }
    values <- ruvFactorGroupAssociations(W, group)
    out$ruv_group_association <- max(values)
    out$ruv_group_association_factors <- values
    rule <- ruv_k_rule %||% ""
    basis <- switch(rule, num_sv_be = "estimated_number", sweep_bounded = "separation_rule", NA_character_)
    if (identical(rule, "user")) {
        typed <- suppressWarnings(as.numeric(k %||% NA))
        estimated <- suppressWarnings(as.numeric(k_estimate %||% NA))
        if (length(typed) != 1L || length(estimated) != 1L || is.na(estimated)) {
            basis <- "typed_k_without_estimated_number"
        } else if (!is.na(typed) && typed == estimated) {
            rule <- "num_sv_be"
            basis <- "typed_k_equals_estimated_number"
        } else {
            basis <- "typed_k_differs_from_estimated_number"
        }
    }
    out$ruv_group_association_basis <- basis
    thresholds <- ruvGroupAssociationThresholds()
    row <- thresholds[thresholds$omic %in% omic & thresholds$ruv_design %in% ruv_design & thresholds$ruv_k_rule %in% rule, , drop = FALSE]
    if (!nrow(row)) {
        out$ruv_group_association_outcome <- "no_limit"
        return(out)
    }
    out$ruv_group_association_threshold <- row$threshold
    out$ruv_group_association_flag <- out$ruv_group_association > row$threshold
    out$ruv_group_association_outcome <- if (out$ruv_group_association_flag) "detectable" else "not_detectable"
    out$ruv_group_association_source <- sprintf("%s (simulator %s, arm %d; set on run %s, validated on run %s)", row$source,
        row$simulator, row$arm, row$setting_run, row$validation_run)
    consequences <- ruvGroupAssociationConsequences()
    figures <- consequences[consequences$omic %in% omic & consequences$ruv_design %in% ruv_design & consequences$ruv_k_rule %in% rule, ,
        drop = FALSE]
    if (nrow(figures) == 1L) {
        out[c("ruv_group_association_fdp_0", "ruv_group_association_fdp_03", "ruv_group_association_fdp_06",
            "ruv_group_association_uncorrected_06")] <- as.list(unlist(figures[c("fdp_0", "fdp_03", "fdp_06", "uncorrected_06")]))
    }
    out
}

#' The group of each row of a factor matrix, from the design (rows named by sample).
#' @noRd
ruvFactorGroups <- function(W, design, sample_id, group_id) {
    if (is.null(W)) return(NULL)
    design <- as.data.frame(design)
    as.character(design[[group_id]][match(rownames(W), as.character(design[[sample_id]]))])
}

#' The fields of the check a step's record keeps (KIT-129 section 1e: `ruv_group_association`, its flag; with the
#' detection limit (`_threshold`), the outcome and the limit's source; M35's basis and M36's measured consequences).
#' @noRd
ruvGroupAssociationRecordFields <- c("ruv_group_association", "ruv_group_association_threshold", "ruv_group_association_flag",
    "ruv_group_association_outcome", "ruv_group_association_source", "ruv_group_association_basis", "ruv_group_association_fdp_0",
    "ruv_group_association_fdp_03", "ruv_group_association_fdp_06", "ruv_group_association_uncorrected_06")

#' The new designs' record fields, one set in every omic (KIT-129 sixth pass item 4, W8, ticket section 1e): the design,
#' the replicate structure and its column, the k rule, the typed k, the applied k, the estimated number with its seed and
#' permutations, the control rule, a supplied list's digest, count and found count, the controls used, the first pass's
#' test, and the DA's standard-error and test methods. A step of today's design records none of them.
#' @noRd
RUV_DESIGN_RECORD_FIELDS <- c("ruv_design", "ruv_replicate_structure", "ruv_grouping_variable", "ruv_k_rule", "ruv_k_typed", "ruv_k",
    "ruv_k_estimate", "ruv_k_seed", "ruv_k_permutations", "ruv_control_rule", "ruv_control_list_sha256", "ruv_control_list_count",
    "ruv_control_list_found", "ruv_control_count", "ruv_first_pass_rule", "ruv_da_method", "ruv_se_method", "ruv_test_method")

#' One set of `RUV_DESIGN_RECORD_FIELDS`: the values given (a NULL or empty value reads NA), NA for the others; the
#' standard-error and test methods follow from `ruv_da_method` (`ruvDesignDaMethods`).
#' @param ... Named values of `RUV_DESIGN_RECORD_FIELDS` (one each).
#' @return A named list in the order of `RUV_DESIGN_RECORD_FIELDS`.
#' @noRd
ruvDesignRecordFields <- function(...) {
    given <- list(...)
    stopifnot(all(names(given) %in% RUV_DESIGN_RECORD_FIELDS))
    out <- stats::setNames(rep(list(NA), length(RUV_DESIGN_RECORD_FIELDS)), RUV_DESIGN_RECORD_FIELDS)
    for (name in names(given)) {
        value <- given[[name]]
        if (length(value) == 1L && !is.list(value)) out[[name]] <- if (is.factor(value)) as.character(value) else value
    }
    methods <- ruvDesignDaMethods(out$ruv_da_method)
    out$ruv_se_method <- methods[["se"]]
    out$ruv_test_method <- methods[["test"]]
    out
}

#' The standard-error and test methods of a new design's DA, by the method its state names (`workflowRuvDaMethodOfArgs`;
#' transcriptomics `edger_qlf_factors_in_model`, `edger_qlf_none`). Decisions M29 (each comparison's standard error
#' exact, `limma::contrastAsCoef`) and M30 (classic proteomics: treat() with the trend and robust settings).
#' @param method The DA method code, or NA.
#' @return `c(se, test)`, NA for an unknown code.
#' @noRd
ruvDesignDaMethods <- function(method) {
    table <- list(
        ruv_factors_exact_se_trend_robust = c(se = "exact_per_comparison", test = "limma_treat_trend_robust"),
        ruv_none_trend_robust = c(se = "contrasts_fit", test = "limma_treat_trend_robust"),
        dpcde_ruv_factors_exact_se = c(se = "exact_per_comparison", test = "limpa_dpcde"),
        dpcde_ruv_none = c(se = "contrasts_fit", test = "limpa_dpcde"),
        limma_ruv_factors_exact_se = c(se = "exact_per_comparison", test = "limma_da_settings"),
        edger_qlf_factors_in_model = c(se = "edger_contrast", test = "edger_qlf"),
        edger_qlf_none = c(se = "edger_contrast", test = "edger_qlf"))
    if (!checkmate::test_string(method) || is.null(table[[method]])) return(c(se = NA_character_, test = NA_character_))
    table[[method]]
}

#' The request parameters of a new design that a step's record keeps beside the registry's own names, so that the saved
#' declaration reproduces the request (KIT-129 sixth pass item 5): the design, the rules, the typed k and seed, and the
#' supplied list itself (`ruv_control_ids`; the identity of classic proteomics and lipid/metab holds only its digest and
#' count, so the list is kept in the step's provenance and written to the analysis record). A legacy request names none.
#' @noRd
RUV_DESIGN_REQUEST_PARAMETERS <- c("ruv_design", "ruv_k_rule", "ruv_k_seed", "ruv_control_rule", "ruv_control_ids",
    "ruv_control_list_sha256", "ruv_control_count", "ruv_k")

#' Check a replayed request's supplied list against its recorded digest (KIT-129 sixth pass item 5): the list the record
#' keeps must hash to the recorded `ruv_control_list_sha256` and count, else the replay is refused.
#' @param parameters The replayed parameters. @param step The record's step (for the message).
#' @return `parameters`, invisibly.
#' @noRd
assertReplayRuvControlList <- function(parameters, step) {
    recorded <- parameters[["ruv_control_list_sha256", exact = TRUE]]
    if (is.null(recorded) || all(is.na(recorded))) return(invisible(parameters))
    ids <- parameters[["ruv_control_ids", exact = TRUE]]
    listed <- if (is.character(ids) && length(ids)) trxCanonicalControlIds(ids) else character()
    if (!length(listed) || !identical(trxControlListSha256(listed), recorded) ||
        !identical(as.integer(parameters[["ruv_control_count", exact = TRUE]]), length(listed))) {
        cli::cli_abort("step {.val {step$state}}: the record's supplied control list does not match its recorded digest",
            class = "workflow_replay_invalid")
    }
    invisible(parameters)
}

#' The identity parameters of a supplied control list (KIT-129 fifth pass item 3): a classic proteomics or lipid/metab
#' request under the design `supplied_controls` is identified by its list's sha256 (`trxControlListSha256`, section 1e:
#' the hash of the sorted ids) and the number of ids listed, not by the ids themselves; the list travels beside them as
#' `ruv_control_ids`, which the kernel reads and checks against the hash.
#' @noRd
RUV_SUPPLIED_LIST_PARAMETERS <- c("ruv_control_list_sha256", "ruv_control_count")

#' A request's parameters with the identity of its supplied list (`RUV_SUPPLIED_LIST_PARAMETERS`) set from its
#' `ruv_control_ids`: NA and 0 for an absent or empty list (the kernel refuses it).
#' @param parameters A request's parameters.
#' @return `parameters`.
#' @noRd
ruvSuppliedListIdentity <- function(parameters) {
    ids <- parameters[["ruv_control_ids", exact = TRUE]]
    listed <- if (length(ids)) trxCanonicalControlIds(ids) else character()
    parameters[RUV_SUPPLIED_LIST_PARAMETERS] <- list(if (length(listed)) trxControlListSha256(listed) else NA_character_,
        length(listed))
    parameters
}

#' The supplied list of a declared request, checked (KIT-129 fifth pass item 3): given under `supplied_controls` and
#' under no other design, one id or more, no missing id, and its hash the request's identity
#' (`ruvSuppliedListIdentity`). Refuses with the legacy refusal class of every omic (`workflow_execution_invalid`).
#' @param parameters The declaration's parameters.
#' @return The canonical list (`trxCanonicalControlIds`), or NULL for another design.
#' @noRd
ruvDeclaredSuppliedList <- function(parameters) {
    invalid <- function(message) cli::cli_abort(message, class = "workflow_execution_invalid", call = NULL)
    supplied <- identical(parameters$ruv_design, "supplied_controls")
    ids <- parameters[["ruv_control_ids", exact = TRUE]]
    if (!supplied) {
        if (!is.null(ids)) invalid("ruv_control_ids names the controls of the design supplied_controls, and is given for no other design")
        return(NULL)
    }
    if (!checkmate::test_character(ids, min.len = 1L, any.missing = FALSE, min.chars = 1L)) {
        invalid("the design supplied_controls needs a list of one control identifier or more (ruv_control_ids)")
    }
    listed <- trxCanonicalControlIds(ids)
    if (!identical(parameters[["ruv_control_list_sha256", exact = TRUE]], trxControlListSha256(listed)) ||
        !identical(as.integer(parameters[["ruv_control_count", exact = TRUE]]), length(listed))) {
        invalid("the request's ruv_control_list_sha256 and ruv_control_count are not those of its ruv_control_ids")
    }
    listed
}

#' Warn when a supplied list names identifiers the data does not hold (KIT-129 fifth pass item 3, as DIA with limpa's
#' `dia_ruv_control_ids_unmatched`): a warning of class `ruv_control_ids_unmatched` with the fields `matched`,
#' `listed` and `unmatched`.
#' @param listed The canonical list. @param present The ids the data holds.
#' @param what The features in words (`proteins`, `lipids`, `metabolites`).
#' @return The unmatched ids.
#' @noRd
warnRuvSuppliedListUnmatched <- function(listed, present, what) {
    unmatched <- setdiff(listed, present)
    if (length(unmatched)) {
        matched <- length(listed) - length(unmatched)
        cli::cli_warn(c("{matched} of the {length(listed)} listed control identifiers are {what} of this data set; the controls are those {matched}",
            i = "not found: {.val {utils::head(unmatched, 5L)}}{if (length(unmatched) > 5L) ' and others'}"),
            class = "ruv_control_ids_unmatched", matched = matched, listed = length(listed), unmatched = unmatched)
    }
    unmatched
}
