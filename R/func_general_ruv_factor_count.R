# KIT-129 W6: the number of unwanted factors (k) of an RUV step, shared by every omic's new path (the user's
# direction of 2026-11-18: two rules, each offered where it is defined).
#   The estimated count: `sva::num.sv(method = "be", B = 20, seed)` on the features complete in every sample, as
#   KIT-124's simulators call it (arms 7, 8, 9, 10, 15, 16). `sva` is an optional dependency; without it the estimate
#   is refused, never replaced.
#   The separation rule: today's estimate on the canonical-correlation plot (`findBestKElbow`, `trxFindBestK`, the
#   sweeps' `best_k`). It is not changed here; an omic hands it to the resolver as a function.
#   A number: the user's k.
# `resolveRuvFactorCount` returns k with its recorded source. Nothing here is called by today's default paths.

#' The seed of the estimated count when the caller gives none: fixed, and recorded with every estimate.
#' @noRd
RUV_FACTOR_COUNT_SEED <- 20261118L

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
#' B = 20, seed)` computes it, on the features observed in every sample: the groups are removed, and a component of
#' the residuals is counted while at most 10% of the row-wise permutations reach its share of the residual
#' variance. The count is `sva`'s own return value; the diagnostic's quantities are computed beside it on the same
#' random stream and must give the same count. The caller's random stream is left as it was.
#'
#' Refused with a condition of class `ruv_factor_count_unavailable` and one of `ruv_factor_count_needs_sva` (`sva`
#' is not installed), `ruv_factor_count_too_few_features` (fewer than `min_features` features observed in every
#' sample) or `ruv_factor_count_too_few_samples` (fewer than two residual components beside the groups).
#'
#' @param values Numeric feature-by-sample matrix; missing values allowed.
#' @param groups The group of each column of `values` (the levels removed before the count).
#' @param seed Whole number; recorded in the result.
#' @param permutations Number of permutations (`num.sv`'s `B`; 20 in the simulations).
#' @param min_features The fewest complete features the estimate is made on (20, the rule of simulator S4).
#' @return `list(k, rule = "num_sv_be", seed, permutations, p_threshold, features_used, features_total, samples,
#'   groups, components, null_shares, sva_version)`: `components` has one row per residual component
#'   (`component`, `observed_share`, `null_mean`, `p_value`, `p_value_monotone`, `counted`); `null_shares` is the
#'   permutation-by-component matrix of shares.
#' @noRd
estimateRuvFactorCount <- function(values, groups, seed = RUV_FACTOR_COUNT_SEED, permutations = 20L,
                                   min_features = 20L) {
    checkmate::assert_matrix(values, mode = "numeric", min.cols = 2L)
    checkmate::assert_atomic_vector(groups, any.missing = FALSE, len = ncol(values))
    seed <- checkmate::assert_int(seed, coerce = TRUE)
    permutations <- checkmate::assert_count(permutations, positive = TRUE, coerce = TRUE)
    checkmate::assert_count(min_features, positive = TRUE)
    group <- if (is.factor(groups)) droplevels(groups) else as.character(groups)
    n_groups <- length(unique(group))
    if (n_groups < 2L) cli::cli_abort("{.arg groups} must name at least two groups")
    refuse <- function(reason, class) {
        cli::cli_abort("the number of unwanted factors cannot be estimated: {reason}",
            class = c(class, "ruv_factor_count_unavailable"), call = NULL)
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
    if (min(dim(dat)) - n_groups < 2L) {
        refuse(sprintf("%d samples in %d groups leave fewer than two residual components to compare", ncol(dat), n_groups),
            "ruv_factor_count_too_few_samples")
    }
    mod <- stats::model.matrix(~ group, data.frame(group = group))
    counted <- withr::with_preserve_seed(list(
        k = as.integer(sva::num.sv(dat, mod, method = "be", B = permutations, seed = seed)),
        shares = permuteRuvResidualShares(dat, mod, permutations, seed)
    ))
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
        components = components, null_shares = null, sva_version = as.character(utils::packageVersion("sva")))
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

#' The largest k a design can carry (KIT-129 section 2, "Details"): k cannot exceed the samples minus the replicate
#' sets (the rank RUV-III's replicates leave), must leave a residual degree of freedom in a model of the groups and
#' the factors, and must be below the number of controls observed in every sample.
#' @param samples,replicate_sets,groups,controls Counts.
#' @return `list(k_max, binding, samples, replicate_sets, groups, controls)`; `binding` names the limit that binds.
#' @noRd
ruvFactorCountLimit <- function(samples, replicate_sets, groups, controls) {
    limits <- c(`samples minus replicate sets` = samples - replicate_sets,
        `a residual degree of freedom beside the groups` = samples - groups - 1L,
        `fewer factors than controls observed in every sample` = controls - 1L)
    list(k_max = as.integer(max(min(limits), 0L)), binding = names(limits)[[which.min(limits)]],
        samples = as.integer(samples), replicate_sets = as.integer(replicate_sets), groups = as.integer(groups),
        controls = as.integer(controls))
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
#' `estimateRuvFactorCount()` pass through under `"num_sv_be"`. A k of 0 is a result: no RUV for the run, with a
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
            estimate <- tryCatch(estimate_count(), ruv_factor_count_unavailable = function(e) e)
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
