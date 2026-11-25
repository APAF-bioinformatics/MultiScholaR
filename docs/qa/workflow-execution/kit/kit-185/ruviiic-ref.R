# KIT-185: the comparator C and the reference Z of the measurement (evidence code under docs/qa, not package code).
#
# Upstream: RUVIIIC 1.0.19, `R/RUVIII_C_Varying.R`, function `RUVIIIC_Varying_R` (lines 75-261), SHA-256
# 6809fe96aa00946c40a7413becd8ba911e70958987c6fb19daac9a103e9314ab (byte-identical to cran/RUVIIIC at c534151f, the
# commit `packaging/multischolar.lock.json` locks). MIT licence (upstream DESCRIPTION: "License: MIT + file LICENSE");
# copyright Children's Medical Research Institute; authors Rohan Shah, Sean Peters, Qing Zhong. Permission is granted
# under the MIT licence to use, copy, modify and distribute the software; the notice above is kept with this copy.
#
# The algorithm is upstream's, per target variable: the rows where the target is observed, the potential controls
# observed on all of those rows, the replicate matrix reduced to the sets with a row, the residual dimension
# (rows - sets); fewer than k residual dimensions or controls: not corrected (NA); exactly k: the GLS formula on the
# controls only (upstream :170-191); otherwise the k leading eigenvectors of P S[rows, rows] P' (P the projection off the
# replicate sets on those rows), alpha = V' Y0[rows, ], W = Yc ac' (ac ac')^-1, corrected = y - W alpha (upstream
# :195-208). Missing rows of the target stay NA.
#
# MARKED CHANGES (the only changes to what is computed):
#   (1) `directions = "complete"`: S is the cross-product over the features observed in every run,
#       S <- Yc %*% t(Yc), in place of upstream's S <- YWithoutNA %*% t(YWithoutNA) over every feature with missing
#       values set to 0 (upstream :110-111, :127). `directions = "zero_filled"` keeps upstream's S.
#   (2) Under both settings the condition `ncol(Msubset) < 2 || nrow(Msubset) < 2 ||` (upstream :159) is removed, so
#       that a feature observed in one replicate set is corrected as the compiled code corrects it
#       (`src/RUVIIIC_Varying.cpp:166-171`); the rule left is upstream's `min(residualDimensions, controls) < k`.
#
# Edits that change no computed value (KIT-185 7.2; verified by E1 (a) against the compiled code): the k leading
# eigenvectors only (eigen(symmetric = TRUE) of the m-by-m matrix; upstream takes min(rd, controls) of them with
# RSpectra and keeps the first k, the same subspace, and W alpha does not depend on the basis chosen in it); targets
# that share one pattern of observed rows share the rows, the controls, P, V and W, so they are corrected together;
# controls per pattern from a precomputed missingness matrix; no progress bar, no checkpoint file, no `try()`
# printing (an error in one pattern leaves its targets NA, as upstream's `try()` does, and is counted).

#' RUV-III-C, varying controls, upstream's R reference with the two marked changes.
#' @param k Number of unwanted factors.
#' @param Y Runs-by-features matrix with NA for missing values; row names required.
#' @param M Replicate matrix (runs by replicate sets, no intercept).
#' @param toCorrect,potentialControls Feature names.
#' @param directions `"zero_filled"` (upstream) or `"complete"` (marked change 1).
#' @return Runs-by-`toCorrect` matrix, NA where not corrected; attributes `residual_dimensions` (per target),
#'   `pattern_errors` (count) and `gap_ratio` (per eigen-branch pattern: lambda_k / (lambda_k - lambda_k+1)).
ruviiicRef <- function(k, Y, M, toCorrect, potentialControls, directions = c("zero_filled", "complete")) {
    directions <- match.arg(directions)
    stopifnot(is.matrix(Y), is.matrix(M), nrow(M) == nrow(Y), !is.null(rownames(Y)), k >= 1, k < nrow(Y),
        k <= length(potentialControls), !anyDuplicated(potentialControls), all(potentialControls %in% colnames(Y)),
        all(toCorrect %in% colnames(Y)))
    missing <- is.na(Y)
    Y0 <- Y
    Y0[missing] <- 0
    symmetrised <- if (identical(directions, "zero_filled")) {
        Y0 %*% t(Y0)
    } else {
        # MARKED CHANGE (1): the directions from the features observed in every run only.
        Yc <- Y[, colSums(missing) == 0L, drop = FALSE]
        Yc %*% t(Yc)
    }
    control_index <- match(potentialControls, colnames(Y))
    target_index <- match(toCorrect, colnames(Y))
    out <- matrix(NA_real_, nrow(Y), length(toCorrect), dimnames = list(rownames(Y), toCorrect))
    residual_dimensions <- stats::setNames(rep(-1L, length(toCorrect)), toCorrect)
    observed <- !missing[, target_index, drop = FALSE]
    keys <- apply(observed, 2L, function(v) paste(which(v), collapse = ","))
    errors <- 0L
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
        submatrix <- Y0[indices, , drop = FALSE]
        targets <- submatrix[, target_index[cols], drop = FALSE]
        controls_values <- submatrix[, controls, drop = FALSE]
        adjusted <- tryCatch({
            if (rd == k) {
                invMatrix <- solve(controls_values %*% t(controls_values))
                Msubset %*% solve(t(Msubset) %*% invMatrix %*% Msubset) %*% t(Msubset) %*% invMatrix %*% targets
            } else {
                projection <- diag(1, m) - Msubset %*% solve(t(Msubset) %*% Msubset) %*% t(Msubset)
                decomposition <- eigen(projection %*% symmetrised[indices, indices, drop = FALSE] %*% t(projection),
                    symmetric = TRUE)
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
    out
}
