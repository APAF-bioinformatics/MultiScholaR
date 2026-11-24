# Shared metrics M1 to M11 (KIT-159 section 4.2). Statistic classes (A22): `accuracy` for M1, M2,
# M3, M8, M9; `threshold_count` for M4, M5, M6, M10, M11; `rank` for M7.

METRIC_CLASSES <- c(
    valBias = "accuracy", valMedianAbsError = "accuracy", valCv = "accuracy",
    valConfusion = "threshold_count", valFdp = "threshold_count",
    valSensitivity = "threshold_count", valPartialAuc = "rank", valCorrelation = "accuracy",
    valQuartetSnr = "accuracy", valMcc = "threshold_count", valCoverage = "threshold_count"
)

#' Median of errors with the n_min rule, pooled and per class
#' @noRd
valSummariseErrors <- function(err, class, n_min) {
    medianOrNa <- function(e) {
        e_ok <- e[is.finite(e)]
        if (length(e_ok) < n_min) NA_real_ else stats::median(e_ok)
    }
    ok <- is.finite(err)
    n <- sum(ok)
    value <- medianOrNa(err)
    per_class <- tapply(err, class, medianOrNa)
    list(
        value = value, n = as.integer(n), n_missing = as.integer(sum(!ok)),
        status = if (is.na(value)) "too_few" else "ok",
        mean_of_classes = if (anyNA(per_class)) NA_real_ else mean(per_class),
        per_class = data.frame(class = names(per_class), value = as.numeric(per_class),
                               stringsAsFactors = FALSE)
    )
}

#' Rows with an expected value for the accuracy metrics M1 and M2
#' @noRd
valAccuracyRows <- function(scored, class) {
    s <- valScoredOnly(scored)
    if (!is.null(class)) s <- s[s$class %in% class, , drop = FALSE]
    s[is.finite(s$expected), , drop = FALSE]
}

#' M1: median bias of the estimates
#'
#' Median over scored keys of class `class` of (estimate - expected). Missing estimates are
#' excluded and counted. With `class = NULL` both the pooled median (`value`) and the mean of the
#' per-class medians (`mean_of_classes`) are returned.
#'
#' @param scored A scored table (one arm, unit and contrast).
#' @param class Class or classes, or `NULL` for all.
#' @param n_min Fewest keys for a value.
#' @return A list: `value`, `n`, `n_missing`, `status`, `mean_of_classes`, `per_class`.
#' @export
valBias <- function(scored, class = NULL, n_min = 10L) {
    s <- valAccuracyRows(scored, class)
    valSummariseErrors(s$estimate - s$expected, s$class, n_min)
}

#' M2: median absolute error of the estimates
#'
#' @inheritParams valBias
#' @return As [valBias()].
#' @export
valMedianAbsError <- function(scored, class = NULL, n_min = 10L) {
    s <- valAccuracyRows(scored, class)
    valSummariseErrors(abs(s$estimate - s$expected), s$class, n_min)
}

#' M3: coefficient of variation on the log2 matrix
#'
#' Per feature and group: linear CV = sd(2^x) / mean(2^x), geometric CV =
#' sqrt(exp((ln 2 sd(x))^2) - 1). A group with fewer than two finite values is skipped; a feature
#' with no group left is skipped. Value: median over features of the per-feature median over groups.
#'
#' @param matrix Numeric log2 matrix, features in rows.
#' @param groups Group of each column.
#' @param scale `linear` (primary, D10) or `geometric`.
#' @param n_min Fewest features for a value.
#' @return A list: `value`, `n`, `status`.
#' @export
valCv <- function(matrix, groups, scale = c("linear", "geometric"), n_min = 10L) {
    scale <- match.arg(scale)
    checkmate::assert_matrix(matrix, mode = "numeric")
    checkmate::assert_atomic_vector(groups, len = ncol(matrix))
    cvOne <- function(x) {
        x <- x[is.finite(x)]
        if (length(x) < 2L) return(NA_real_)
        if (scale == "linear") {
            stats::sd(2^x) / mean(2^x)
        } else {
            sqrt(exp((log(2) * stats::sd(x))^2) - 1)
        }
    }
    per_feature <- apply(matrix, 1L, function(row) {
        cv <- vapply(split(row, groups), cvOne, numeric(1))
        if (all(is.na(cv))) NA_real_ else stats::median(cv, na.rm = TRUE)
    })
    per_feature <- per_feature[is.finite(per_feature)]
    n <- length(per_feature)
    value <- if (n < n_min) NA_real_ else stats::median(per_feature)
    list(value = value, n = as.integer(n), status = if (is.na(value)) "too_few" else "ok")
}

#' Positives, negatives and the grey zone of a scored table
#' @noRd
valTruthSides <- function(s, positive_rule, delta) {
    if (positive_rule == "class") {
        list(pos = s$class %in% "changed", neg = s$class %in% "null",
             grey = rep(FALSE, nrow(s)))
    } else {
        e <- s$expected
        list(pos = is.finite(e) & abs(e) > delta, neg = is.finite(e) & e == 0,
             grey = is.finite(e) & e != 0 & abs(e) <= delta)
    }
}

#' M4: confusion table in the own and union views
#'
#' A call on a positive with the expected sign is a TP. A call on a positive with the wrong sign
#' is a false discovery and its key stays a miss: it adds 1 to FP, to `sign_errors` and to FN
#' (A03, D5). A call on a negative is FP; an uncalled negative is TN. The own view holds the keys
#' the arm tested; the union view every scored key given (an untested positive is a miss).
#'
#' @param scored A scored table (one arm, unit and contrast; restricted to a view beforehand).
#' @param alpha Call threshold.
#' @param call_rule `own` (`s_own`) or `bh` (`s_bh`).
#' @param positive_rule `expected` (expected != 0 / = 0) or `class` (`changed` / `null`).
#' @param delta Grey zone: abs(expected) <= delta, delta > 0, is neither P nor N.
#' @param lambda Fold-change threshold of the call rule (0: none).
#' @return A data frame with one row per view (`own`, `union`): `P`, `N`, `TP`, `FP`, `FN`, `TN`,
#'   `sign_errors`, `no_statistic`, `grey_zone`.
#' @export
valConfusion <- function(scored, alpha = 0.05, call_rule = c("own", "bh"),
                         positive_rule = c("expected", "class"), delta = 0, lambda = 0) {
    call_rule <- match.arg(call_rule)
    positive_rule <- match.arg(positive_rule)
    checkmate::assert_number(alpha, lower = 0, upper = 1)
    checkmate::assert_number(delta, lower = 0)
    s <- valScoredOnly(scored)
    stat <- if (call_rule == "own") s$s_own else s$s_bh
    testable <- is.finite(s$p)
    quantified <- is.finite(s$estimate)
    fc_ok <- if (lambda > 0) quantified & abs(s$estimate) >= lambda else rep(TRUE, nrow(s))
    called <- testable & is.finite(stat) & stat < alpha & fc_ok
    sides <- valTruthSides(s, positive_rule, delta)
    has_dir <- sides$pos & is.finite(s$expected) & s$expected != 0 & quantified
    right_sign <- !has_dir | sign(s$estimate) == sign(s$expected)
    tp <- called & sides$pos & right_sign
    sign_err <- called & sides$pos & !right_sign
    fp <- (called & sides$neg) | sign_err
    count <- function(rows) {
        data.frame(
            P = sum(sides$pos[rows]), N = sum(sides$neg[rows]), TP = sum(tp[rows]),
            FP = sum(fp[rows]), FN = sum((sides$pos & !tp)[rows]),
            TN = sum((sides$neg & !called)[rows]), sign_errors = sum(sign_err[rows])
        )
    }
    out <- rbind(count(testable), count(rep(TRUE, nrow(s))))
    out[] <- lapply(out, as.integer)
    out <- cbind(view = c("own", "union"), out, stringsAsFactors = FALSE)
    out$no_statistic <- as.integer(sum(quantified & !testable))
    out$grey_zone <- as.integer(sum(sides$grey))
    out$call_rule <- call_rule
    out$alpha <- alpha
    out
}

#' M5: false discovery proportion
#'
#' FDP = FP / (TP + FP); with no call FDP = 0 and `no_calls = TRUE` (the BH convention). An all-null
#' truth with any call gives FDP = 1.
#'
#' @param conf Output of [valConfusion()].
#' @return A list: `value`, `no_calls`, `status`.
#' @export
valFdp <- function(conf) {
    r <- conf[conf$view == "union", ]
    calls <- r$TP + r$FP
    if (calls == 0L) return(list(value = 0, no_calls = TRUE, status = "ok"))
    list(value = r$FP / calls, no_calls = FALSE, status = "ok")
}

#' M6: sensitivity in a view
#'
#' @param conf Output of [valConfusion()].
#' @param view `own` or `union`.
#' @return A list: `value`, `status` (`no_positives` when P = 0).
#' @export
valSensitivity <- function(conf, view = c("union", "own")) {
    view <- match.arg(view)
    r <- conf[conf$view == view, ]
    if (r$P == 0L) return(list(value = NA_real_, status = "no_positives"))
    list(value = r$TP / r$P, status = "ok")
}

#' M10: Matthews correlation from M4's table
#'
#' A wrong-sign key counts in FP and in FN, so the four cells sum to n + sign errors (A26).
#'
#' @param conf Output of [valConfusion()].
#' @param view `own` or `union`.
#' @return A list: `value`, `status` (`marginal_zero` when any marginal is zero).
#' @export
valMcc <- function(conf, view = c("own", "union")) {
    view <- match.arg(view)
    r <- conf[conf$view == view, ]
    tp <- as.numeric(r$TP)
    fp <- as.numeric(r$FP)
    fn <- as.numeric(r$FN)
    tn <- as.numeric(r$TN)
    den <- (tp + fp) * (tp + fn) * (tn + fp) * (tn + fn)
    if (den == 0) return(list(value = NA_real_, status = "marginal_zero"))
    list(value = (tp * tn - fp * fn) / sqrt(den), status = "ok")
}

#' The ranking p of M7
#'
#' One-sided p in the expected direction for a positive with an expected sign (p / 2 when the
#' estimate has that sign, 1 - p / 2 otherwise); the two-sided p for every other key; 1 for an
#' untested key (the union view's last block).
#' @noRd
valRankingP <- function(s, pos) {
    dir <- ifelse(pos & is.finite(s$expected), sign(s$expected), 0)
    tested <- is.finite(s$p)
    one_sided <- dir != 0 & tested & is.finite(s$estimate)
    rp <- s$p
    rp[one_sided] <- ifelse(sign(s$estimate[one_sided]) == dir[one_sided], s$p[one_sided] / 2,
                            1 - s$p[one_sided] / 2)
    rp[!tested] <- 1
    rp
}

#' Area under a ROC polyline up to x_max (trapezoid rule, linear interpolation at x_max)
#' @noRd
valAreaTo <- function(fpr, tpr, x_max) {
    area <- 0
    for (i in seq_len(length(fpr) - 1L)) {
        x0 <- fpr[i]
        x1 <- fpr[i + 1L]
        if (x0 >= x_max) break
        y0 <- tpr[i]
        y1 <- tpr[i + 1L]
        if (x1 > x_max) {
            y1 <- y0 + (y1 - y0) * (x_max - x0) / (x1 - x0)
            x1 <- x_max
        }
        area <- area + (x1 - x0) * (y0 + y1) / 2
    }
    area
}

#' M7: partial area under the ROC curve to `fpr_max`
#'
#' Keys ranked by the ranking p (see section 4.2, M7); the ROC curve starts at the origin; tied
#' ranking p form one block traversed as a straight segment; linear interpolation to `fpr_max`;
#' area by the trapezoid rule divided by `fpr_max` (1 = perfect, `fpr_max` / 2 = random).
#'
#' @param scored A scored table (one arm, unit and contrast; restricted to a view beforehand).
#' @param fpr_max Upper end of the false positive rate (no default; SCI-23).
#' @param min_nulls_in_region Fewest nulls floor(fpr_max * N) in the region.
#' @param min_events Fewest positives.
#' @param positive_rule `expected` or `class`.
#' @param delta Grey zone, as [valConfusion()].
#' @return A list: `value` (normalised), `raw`, `full_auc`, `n_pos`, `n_null`, `roc`, `ranking_p`,
#'   `status`.
#' @export
valPartialAuc <- function(scored, fpr_max, min_nulls_in_region = 10L, min_events = 5L,
                          positive_rule = c("expected", "class"), delta = 0) {
    positive_rule <- match.arg(positive_rule)
    if (missing(fpr_max) || !checkmate::test_number(fpr_max, lower = 0, upper = 1) ||
            fpr_max <= 0) {
        valAbort("{.arg fpr_max} must be declared, in (0, 1]", "val_spec_invalid")
    }
    s <- valScoredOnly(scored)
    sides <- valTruthSides(s, positive_rule, delta)
    keep <- sides$pos | sides$neg
    s <- s[keep, , drop = FALSE]
    pos <- sides$pos[keep]
    n_pos <- sum(pos)
    n_null <- sum(!pos)
    if (floor(fpr_max * n_null + 1e-9) < min_nulls_in_region) {
        valAbort("too few nulls in the region: floor({fpr_max} x {n_null}) < {min_nulls_in_region}",
                 "val_spec_invalid")
    }
    if (n_pos < min_events || n_pos == 0L || n_null == 0L) {
        valAbort("too few positives for the curve: {n_pos} < {min_events}", "val_spec_invalid")
    }
    rp <- valRankingP(s, pos)
    block <- signif(rp, 12L)
    pos_b <- as.vector(rowsum(as.numeric(pos), block, reorder = TRUE))
    neg_b <- as.vector(rowsum(as.numeric(!pos), block, reorder = TRUE))
    fpr <- c(0, cumsum(neg_b) / n_null)
    tpr <- c(0, cumsum(pos_b) / n_pos)
    raw <- valAreaTo(fpr, tpr, fpr_max)
    list(
        value = raw / fpr_max, raw = raw, full_auc = valAreaTo(fpr, tpr, 1),
        n_pos = as.integer(n_pos), n_null = as.integer(n_null),
        roc = data.frame(fpr = fpr, tpr = tpr), ranking_p = stats::setNames(rp, s$key),
        fpr_max = fpr_max, status = "ok"
    )
}

#' M8: correlation of the estimates with the expected values
#'
#' @param scored A scored table.
#' @param method `pearson` or `spearman`.
#' @param n_min Fewest keys with both values finite.
#' @return A list: `value`, `n`, `status`.
#' @export
valCorrelation <- function(scored, method = c("pearson", "spearman"), n_min = 10L) {
    method <- match.arg(method)
    s <- valScoredOnly(scored)
    ok <- is.finite(s$estimate) & is.finite(s$expected)
    n <- sum(ok)
    if (n < n_min || n < 2L) return(list(value = NA_real_, n = as.integer(n), status = "too_few"))
    list(value = stats::cor(s$estimate[ok], s$expected[ok], method = method),
         n = as.integer(n), status = "ok")
}

#' M9: Quartet signal-to-noise ratio from principal-component scores
#'
#' D(i, j) = w1 dPC1^2 + w2 dPC2^2 over ordered distinct pairs; SNR = 10 log10(mean inter-group D /
#' mean intra-group D).
#'
#' @param scores Matrix with the first two principal-component scores in its first two columns.
#' @param weights The two variance-explained proportions (not renormalised).
#' @param groups Group of each row.
#' @return The SNR (dB), or NA when a group has fewer than two samples.
#' @export
valSnrFromScores <- function(scores, weights, groups) {
    checkmate::assert_matrix(scores, min.cols = 2L)
    checkmate::assert_numeric(weights, len = 2L)
    checkmate::assert_atomic_vector(groups, len = nrow(scores))
    if (any(table(groups) < 2L)) return(NA_real_)
    d1 <- outer(scores[, 1L], scores[, 1L], "-")^2
    d2 <- outer(scores[, 2L], scores[, 2L], "-")^2
    d <- weights[[1L]] * d1 + weights[[2L]] * d2
    same <- outer(groups, groups, "==")
    off <- row(d) != col(d)
    10 * log10(mean(d[off & !same]) / mean(d[off & same]))
}

#' M9: Quartet signal-to-noise ratio from a matrix
#'
#' As Quartet's `calc_signoise_ratio`: per-feature z-score, `prcomp(t(z), scale. = FALSE)`, weights
#' the variance-explained proportions of PC1 and PC2 (`rounded5`: rounded to 5 decimals as
#' `summary(prcomp)$importance` gives them; `exact`: unrounded). Features with a non-finite value or
#' zero variance are dropped and counted.
#'
#' @param matrix Numeric matrix, features in rows (the declared scale and filter applied before).
#' @param groups Group of each column.
#' @param weights `rounded5` or `exact`.
#' @return A list: `value`, `weights`, `n_features`, `n_dropped`.
#' @export
valQuartetSnr <- function(matrix, groups, weights = c("rounded5", "exact")) {
    weights <- match.arg(weights)
    checkmate::assert_matrix(matrix, mode = "numeric")
    checkmate::assert_atomic_vector(groups, len = ncol(matrix))
    finite <- apply(matrix, 1L, function(r) all(is.finite(r)))
    vars <- apply(matrix, 1L, stats::var)
    keep <- finite & is.finite(vars) & vars > 0
    x <- matrix[keep, , drop = FALSE]
    if (any(table(groups) < 2L) || nrow(x) < 2L) {
        return(list(value = NA_real_, weights = c(NA_real_, NA_real_),
                    n_features = sum(keep), n_dropped = sum(!keep)))
    }
    z <- t(scale(t(x)))
    pc <- stats::prcomp(t(z), scale. = FALSE)
    w <- if (weights == "rounded5") {
        summary(pc)$importance[2L, 1:2]
    } else {
        (pc$sdev^2 / sum(pc$sdev^2))[1:2]
    }
    list(value = valSnrFromScores(pc$x[, 1:2], w, groups), weights = unname(w),
         n_features = as.integer(sum(keep)), n_dropped = as.integer(sum(!keep)))
}

#' M11: coverage per class
#'
#' @param scored A scored table (one arm, unit and contrast).
#' @param matrix The arm's matrix, for its completeness (share of finite cells), or `NULL`.
#' @param alpha Call threshold.
#' @return A data frame per class: `keys`, `quantified`, `testable`, `called_own`, `called_bh`,
#'   `completeness`.
#' @export
valCoverage <- function(scored, matrix = NULL, alpha = 0.05) {
    s <- valScoredOnly(scored)
    completeness <- if (is.null(matrix)) NA_real_ else mean(is.finite(matrix))
    testable <- is.finite(s$p)
    called_own <- testable & is.finite(s$s_own) & s$s_own < alpha
    called_bh <- testable & is.finite(s$s_bh) & s$s_bh < alpha
    classes <- sort(unique(s$class))
    do.call(rbind, lapply(classes, function(cl) {
        i <- s$class %in% cl
        data.frame(class = cl, keys = sum(i), quantified = sum(is.finite(s$estimate[i])),
                   testable = sum(testable[i]), called_own = sum(called_own[i]),
                   called_bh = sum(called_bh[i]), completeness = completeness,
                   stringsAsFactors = FALSE)
    }))
}
