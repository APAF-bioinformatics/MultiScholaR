# KIT-210: the reference RUV-III-C of the candidates C3, C4a and C7 (evidence code under docs/qa, not package code).
# A copy of `kit-185/ruviiic-ref.R` (SHA-256 019ee0a6a039f3bf1b842dbff6543e964dda08ce393908a933ade1ea40e39130) with
# the marked changes (3) and (4) and the guard; with `directions = "zero_filled"` or `"complete"` and `guard = NULL` it
# computes exactly what that file computes (check E1').
#
# Upstream: RUVIIIC 1.0.19, `R/RUVIII_C_Varying.R`, function `RUVIIIC_Varying_R` (lines 75-261), SHA-256
# 6809fe96aa00946c40a7413becd8ba911e70958987c6fb19daac9a103e9314ab (byte-identical to cran/RUVIIIC at c534151f, the
# commit `packaging/multischolar.lock.json` locks; KIT-210 step A, run 45351569441). MIT licence (upstream DESCRIPTION:
# "License: MIT + file LICENSE"); copyright Children's Medical Research Institute; authors Rohan Shah, Sean Peters,
# Qing Zhong. Permission is granted under the MIT licence to use, copy, modify and distribute the software; the notice
# above is kept with this copy.
#
# The algorithm is upstream's, per target variable: the rows where the target is observed, the potential controls
# observed on all of those rows, the replicate matrix reduced to the sets with a row, the residual dimension
# (rows - sets); fewer than k residual dimensions or controls: not corrected (NA); exactly k: the GLS formula on the
# controls only (upstream :170-191); otherwise the k leading eigenvectors of P S[rows, rows] P' (P the projection off the
# replicate sets on those rows), alpha = V' Y0[rows, ], W = Yc ac' (ac ac')^-1, corrected = y - W alpha (upstream
# :195-208). Missing rows of the target stay NA.
#
# MARKED CHANGES (the only changes to what is computed):
#   (1) [KIT-185] `directions = "complete"`: S is the cross-product over the features observed in every run,
#       S <- Yc %*% t(Yc), in place of upstream's S <- YWithoutNA %*% t(YWithoutNA) over every feature with missing
#       values set to 0 (upstream :110-111, :127). `directions = "zero_filled"` keeps upstream's S. (Candidate C3.)
#   (2) [KIT-185] Under every setting the condition `ncol(Msubset) < 2 || nrow(Msubset) < 2 ||` (upstream :159) is
#       removed, so that a feature observed in one replicate set is corrected as the compiled code corrects it
#       (`src/RUVIIIC_Varying.cpp:166-171`); the rule left is upstream's `min(residualDimensions, controls) < k`.
#   (3) [KIT-210] `directions = "subset"` (candidate C4a): for the target's observed rows, S is the cross-product over
#       F(rows) = the features observed on every one of those rows, Y[rows, F] %*% t(Y[rows, F]); no zero enters.
#   (4) [KIT-210] `directions = "nipals"` (candidate C7): S is the cross-product of `nipals(Y, ncomp)` (the package's
#       `imputeNipalsMatrix`, passed in as `nipals`), in place of the zero-filled Y; targets, controls, alpha and W are
#       taken from the observed values as upstream (they are observed on the target's rows by construction).
#   Guard (KIT-210 section 5, D7): with `guard = g`, `"complete"` stops with class `candidate_too_few_complete` when
#       fewer than g features are observed in every run, and `"subset"` leaves a target NA when |F(rows)| < g (counted in
#       the attribute `guard_targets`). `guard = NULL` (the default) applies no guard.
#
# Edits that change no computed value (KIT-185 7.2; verified by E1 (a) against the compiled code): the k leading
# eigenvectors only (eigen(symmetric = TRUE) of the m-by-m matrix; upstream takes min(rd, controls) of them with
# RSpectra and keeps the first k, the same subspace, and W alpha does not depend on the basis chosen in it); targets
# that share one pattern of observed rows share the rows, the controls, P, V and W, so they are corrected together;
# controls per pattern from a precomputed missingness matrix; no progress bar, no checkpoint file, no `try()`
# printing (an error in one pattern leaves its targets NA, as upstream's `try()` does, and is counted).

#' RUV-III-C, varying controls, upstream's R reference with the marked changes.
#' @param k Number of unwanted factors.
#' @param Y Runs-by-features matrix with NA for missing values; row names required.
#' @param M Replicate matrix (runs by replicate sets, no intercept).
#' @param toCorrect,potentialControls Feature names.
#' @param directions `"zero_filled"` (upstream), `"complete"` (1), `"subset"` (3) or `"nipals"` (4).
#' @param guard `NULL` (no guard) or the minimum number of complete features g.
#' @param nipals,ncomp For `"nipals"`: the imputation function (runs-by-features in and out) and its component count.
#' @return Runs-by-`toCorrect` matrix, NA where not corrected; attributes `residual_dimensions` (per target),
#'   `pattern_errors` (count), `gap_ratio` (per eigen-branch pattern: lambda_k / (lambda_k - lambda_k+1)) and
#'   `guard_targets` (targets left NA by the guard).
ruviiicRef <- function(k, Y, M, toCorrect, potentialControls, directions = c("zero_filled", "complete", "subset", "nipals"),
                       guard = NULL, nipals = NULL, ncomp = NULL) {
    directions <- match.arg(directions)
    stopifnot(is.matrix(Y), is.matrix(M), nrow(M) == nrow(Y), !is.null(rownames(Y)), k >= 1, k < nrow(Y),
        k <= length(potentialControls), !anyDuplicated(potentialControls), all(potentialControls %in% colnames(Y)),
        all(toCorrect %in% colnames(Y)))
    missing <- is.na(Y)
    Y0 <- Y
    Y0[missing] <- 0
    symmetrised <- if (identical(directions, "zero_filled")) {
        Y0 %*% t(Y0)
    } else if (identical(directions, "complete")) {
        # MARKED CHANGE (1): the directions from the features observed in every run only.
        complete <- colSums(missing) == 0L
        if (!is.null(guard) && sum(complete) < guard) {
            stop(structure(class = c("candidate_too_few_complete", "error", "condition"),
                list(message = sprintf("%d features observed in every run, fewer than the guard %d", sum(complete), guard), call = NULL)))
        }
        Yc <- Y[, complete, drop = FALSE]
        Yc %*% t(Yc)
    } else if (identical(directions, "nipals")) {
        # MARKED CHANGE (4): the cross-product of the NIPALS-imputed matrix.
        stopifnot(is.function(nipals), !is.null(ncomp))
        Yn <- nipals(Y, ncomp)
        Yn %*% t(Yn)
    } else NULL # MARKED CHANGE (3): formed per pattern below
    control_index <- match(potentialControls, colnames(Y))
    target_index <- match(toCorrect, colnames(Y))
    out <- matrix(NA_real_, nrow(Y), length(toCorrect), dimnames = list(rownames(Y), toCorrect))
    residual_dimensions <- stats::setNames(rep(-1L, length(toCorrect)), toCorrect)
    observed <- !missing[, target_index, drop = FALSE]
    keys <- apply(observed, 2L, function(v) paste(which(v), collapse = ","))
    errors <- 0L
    guard_targets <- 0L
    gap_ratio <- numeric()
    for (key in unique(keys)) {
        cols <- which(keys == key)
        indices <- which(observed[, cols[[1L]]])
        if (!length(indices)) next
        Msubset <- M[indices, , drop = FALSE]
        Msubset <- Msubset[, colSums(Msubset) > 0, drop = FALSE]
        m <- length(indices)
        rd <- m - ncol(Msubset)
        residual_dimensions[cols] <- rd
        controls <- control_index[colSums(missing[indices, control_index, drop = FALSE]) == 0L]
        # MARKED CHANGE (2): upstream's `ncol(Msubset) < 2 || nrow(Msubset) < 2 ||` is not part of the rule.
        if (min(rd, length(controls)) < k) next
        S_rows <- if (identical(directions, "subset")) {
            # MARKED CHANGE (3): the features observed on every row of this pattern.
            on_rows <- colSums(missing[indices, , drop = FALSE]) == 0L
            if (!is.null(guard) && sum(on_rows) < guard) {
                guard_targets <- guard_targets + length(cols)
                next
            }
            Yf <- Y[indices, on_rows, drop = FALSE]
            Yf %*% t(Yf)
        } else symmetrised[indices, indices, drop = FALSE]
        submatrix <- Y0[indices, , drop = FALSE]
        targets <- submatrix[, target_index[cols], drop = FALSE]
        controls_values <- submatrix[, controls, drop = FALSE]
        adjusted <- tryCatch({
            if (rd == k) {
                invMatrix <- solve(controls_values %*% t(controls_values))
                Msubset %*% solve(t(Msubset) %*% invMatrix %*% Msubset) %*% t(Msubset) %*% invMatrix %*% targets
            } else {
                projection <- diag(1, m) - Msubset %*% solve(t(Msubset) %*% Msubset) %*% t(Msubset)
                decomposition <- eigen(projection %*% S_rows %*% t(projection), symmetric = TRUE)
                values <- decomposition$values
                if (length(values) > k) gap_ratio <- c(gap_ratio, values[[k]] / (values[[k]] - values[[k + 1L]]))
                vectors <- decomposition$vectors[, seq_len(k), drop = FALSE]
                ac <- t(vectors) %*% controls_values
                W <- controls_values %*% t(ac) %*% solve(ac %*% t(ac))
                targets - W %*% (t(vectors) %*% targets)
            }
        }, error = function(e) NULL)
        if (is.null(adjusted)) {
            errors <- errors + 1L
            next
        }
        out[indices, cols] <- adjusted
    }
    attr(out, "residual_dimensions") <- residual_dimensions
    attr(out, "pattern_errors") <- errors
    attr(out, "gap_ratio") <- gap_ratio
    attr(out, "guard_targets") <- guard_targets
    out
}
