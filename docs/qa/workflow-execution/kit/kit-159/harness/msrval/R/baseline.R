# The baseline arm (KIT-159 section 5.5), fixed before any validation data are read.

#' Parse declared contrasts `"X-Y"` (X minus Y) into ids and group pairs
#' @noRd
valParseContrasts <- function(contrasts, groups) {
    checkmate::assert_character(contrasts, min.len = 1L, any.missing = FALSE)
    ids <- names(contrasts)
    if (is.null(ids)) ids <- contrasts
    ids[!nzchar(ids)] <- contrasts[!nzchar(ids)]
    parts <- strsplit(contrasts, "-", fixed = TRUE)
    ok <- vapply(parts, function(p) length(p) == 2L && all(trimws(p) %in% groups), logical(1))
    if (!all(ok)) {
        valAbort("a contrast is not of the form <group>-<group> over the design's groups",
                 "val_spec_invalid")
    }
    data.frame(id = ids, g1 = vapply(parts, function(p) trimws(p[[1L]]), character(1)),
               g2 = vapply(parts, function(p) trimws(p[[2L]]), character(1)),
               stringsAsFactors = FALSE)
}

#' The baseline arm: the upstream tool's quantities through a plain stated test
#'
#' Intensities: zeros to NA, log2, `normalise` (`none` for quantities the tool already normalised;
#' `median`: subtract each sample's median over finite values and add the median of the sample
#' medians), a per-contrast filter (a key enters a contrast with at least `min_per_group` finite
#' values in each of that contrast's two groups; the fit uses every key entering at least one
#' contrast), one `limma::lmFit(~ 0 + group [+ batch])`, `contrasts.fit`, `eBayes()` (no trend, no
#' robust), BH within contrast. Counts: `edgeR::filterByExpr` at its defaults, TMM with
#' `keep.lib.sizes = FALSE`, `limma::voom`, the same model. No imputation, RUV or fold-change
#' threshold.
#'
#' @param quant Numeric matrix of the tool's quantities (features in rows, samples in columns).
#' @param design Data frame with `sample`, `group` and the optional batch column.
#' @param contrasts Character `"X-Y"` (X minus Y), named by contrast id.
#' @param data_type `intensity` or `counts`.
#' @param normalise `none` or `median` (intensities).
#' @param batch Name of a design column added as a fixed effect, or `NULL`.
#' @param min_per_group Fewest finite values per group of a contrast.
#' @param alpha Call threshold.
#' @return A list: `da` (key, contrast, estimate, se, t, p, s_own, called_own, call_statistic),
#'   `filter` (contrast, n_keys, n_keys_r11), `normalised`, `sample_medians`, `design_columns`,
#'   `design_matrix`, `contrast_matrix`.
#' @export
valBaselineArm <- function(quant, design, contrasts, data_type = c("intensity", "counts"),
                           normalise = c("none", "median"), batch = NULL, min_per_group = 2L,
                           alpha = 0.05) {
    data_type <- match.arg(data_type)
    normalise <- match.arg(normalise)
    checkmate::assert_matrix(quant, mode = "numeric", row.names = "unique", col.names = "unique")
    checkmate::assert_data_frame(design)
    checkmate::assert_names(names(design), must.include = c("sample", "group", batch))
    if (!setequal(design$sample, colnames(quant))) {
        valAbort("the design's samples differ from the matrix columns", "val_spec_invalid")
    }
    q <- quant[, design$sample, drop = FALSE]
    groups <- as.character(design$group)
    ct <- valParseContrasts(contrasts, unique(groups))
    sizes <- table(groups)
    small <- sizes[c(ct$g1, ct$g2)] < min_per_group
    if (any(small)) {
        valAbort("a contrast group has fewer than {min_per_group} samples", "val_baseline_refused")
    }
    design_columns <- design[, c("group", batch), drop = FALSE]
    x_mat <- valBaselineDesign(groups, if (is.null(batch)) NULL else design[[batch]])
    cm <- matrix(0, nrow = ncol(x_mat), ncol = nrow(ct),
                 dimnames = list(colnames(x_mat), ct$id))
    for (j in seq_len(nrow(ct))) {
        cm[ct$g1[j], j] <- 1
        cm[ct$g2[j], j] <- -1
    }
    fitted <- if (data_type == "intensity") {
        valBaselineIntensity(q, groups, ct, x_mat, normalise, min_per_group)
    } else {
        valBaselineCounts(q, groups, ct, x_mat)
    }
    fit <- limma::lmFit(fitted$x, x_mat)
    fit <- limma::contrasts.fit(fit, cm)
    fit <- limma::eBayes(fit)
    s2_post <- stats::setNames(fit$s2.post, rownames(fit$coefficients))
    da <- do.call(rbind, lapply(seq_len(nrow(ct)), function(j) {
        keys <- rownames(fitted$enter)[fitted$enter[, j]]
        p <- fit$p.value[keys, j]
        data.frame(key = keys, contrast = ct$id[j], estimate = fit$coefficients[keys, j],
                   se = fit$stdev.unscaled[keys, j] * sqrt(s2_post[keys]),
                   t = fit$t[keys, j], p = p, s_own = stats::p.adjust(p, method = "BH"),
                   call_statistic = "BH", stringsAsFactors = FALSE, row.names = NULL)
    }))
    da$called_own <- is.finite(da$s_own) & da$s_own < alpha
    list(da = da, filter = fitted$filter, normalised = fitted$normalised,
         sample_medians = fitted$sample_medians, design_columns = design_columns,
         design_matrix = x_mat, contrast_matrix = cm, data_type = data_type,
         normalise = normalise)
}

#' Group-means design with an optional batch term; refuses a non-estimable design (A42)
#' @noRd
valBaselineDesign <- function(groups, batch) {
    grp <- factor(groups, levels = unique(groups))
    x_mat <- if (is.null(batch)) {
        stats::model.matrix(~ 0 + grp)
    } else {
        bt <- factor(batch)  # nolint
        stats::model.matrix(~ 0 + grp + bt)
    }
    colnames(x_mat)[seq_len(nlevels(grp))] <- levels(grp)
    if (qr(x_mat)$rank < ncol(x_mat)) {
        valAbort("the batch term leaves a contrast not estimable", "val_baseline_refused")
    }
    x_mat
}

#' Intensity path: zeros to NA, log2, normalisation, per-contrast filter
#' @noRd
valBaselineIntensity <- function(q, groups, ct, x_mat, normalise, min_per_group) {
    x <- q
    x[!is.na(x) & x == 0] <- NA
    x <- log2(x)
    sample_medians <- NULL
    if (normalise == "median") {
        sample_medians <- apply(x, 2L, stats::median, na.rm = TRUE)
        x <- sweep(x, 2L, sample_medians) + stats::median(sample_medians)
    }
    finite_in <- function(g) rowSums(is.finite(x[, groups == g, drop = FALSE]))
    counts <- vapply(unique(groups), finite_in, numeric(nrow(x)))
    counts <- matrix(counts, nrow = nrow(x), dimnames = list(rownames(x), unique(groups)))
    enter <- vapply(seq_len(nrow(ct)), function(j) {
        counts[, ct$g1[j]] >= min_per_group & counts[, ct$g2[j]] >= min_per_group
    }, logical(nrow(x)))
    enter <- matrix(enter, nrow = nrow(x), dimnames = list(rownames(x), ct$id))
    contrast_groups <- unique(c(ct$g1, ct$g2))
    r11 <- apply(counts[, contrast_groups, drop = FALSE] >= min_per_group, 1L, all)
    filter <- data.frame(contrast = ct$id, n_keys = as.integer(colSums(enter)),
                         n_keys_r11 = as.integer(colSums(enter & r11)), stringsAsFactors = FALSE)
    in_fit <- rowSums(enter) > 0
    list(x = x[in_fit, , drop = FALSE], enter = enter[in_fit, , drop = FALSE], filter = filter,
         normalised = x, sample_medians = sample_medians)
}

#' Counts path: filterByExpr, TMM, voom
#' @noRd
valBaselineCounts <- function(q, groups, ct, x_mat) {
    y <- edgeR::DGEList(counts = q, group = factor(groups))
    keep <- edgeR::filterByExpr(y, group = factor(groups))
    y <- y[keep, , keep.lib.sizes = FALSE]
    y <- edgeR::calcNormFactors(y, method = "TMM")
    v <- limma::voom(y, x_mat)
    enter <- matrix(TRUE, nrow = nrow(v$E), ncol = nrow(ct),
                    dimnames = list(rownames(v$E), ct$id))
    filter <- data.frame(contrast = ct$id, n_keys = nrow(v$E), n_keys_r11 = nrow(v$E),
                         stringsAsFactors = FALSE)
    list(x = v, enter = enter, filter = filter, normalised = v$E, sample_medians = NULL)
}
