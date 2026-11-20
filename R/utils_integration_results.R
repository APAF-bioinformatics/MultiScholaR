# MultiScholaR: Interactive Multi-Omics Analysis
# Copyright (C) 2024-2026 Ignatius Pang, William Klare, and APAF-bioinformatics
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

# KIT-063: the standard results of a trained integration model (MOFA2), read from the model's saved product
# (`integrationModelProduct`, KIT-062) and the imported inputs state; nothing here trains or changes the model. The
# results are pure reads of the stored model, so they are no declared operation; saving them to the project goes
# through the owner's action-output publication (`publishWorkflowActionOutputs`, the journal the omics' publications
# use), drawn in a short folder and every file checked (`withWorkflowShortStaging`, `checkWorkflowFigureFiles`;
# KIT-052). Every choice the user's integration workbook would make (how many features, which factors, the statistic,
# figure size and format) is a declared parameter with a provisional default in one table,
# `integrationResultsParameterSpec()`. Enrichment on the weights is KIT-141, the record and report KIT-142 (a saved
# publication is kept with its model for the record, `integrationRecordResultsPublication`).

#' The label of the results' publication on the owner's action-output journal.
#' @noRd
INTEGRATION_RESULTS_LABEL <- "integration-results"

#' The Results tab's parameters with their provisional defaults: the one table the user's integration workbook settles
#' (`settled_by`), in the form of `integrationMofaParameterSpec()`.
#' @return Named list by parameter: `default`, `choices` (NULL for a number or a free value), `settled_by`.
#' @noRd
integrationResultsParameterSpec <- function() {
    workbook <- "workbook (provisional)"
    list(
        # Which factors the figures and tables show: "all", or factor names (Factor1, ...).
        factors_shown = list(default = "all", choices = NULL, settled_by = workbook),
        # How many features per factor and data set are listed and drawn (the starter workbook's plot_weights: 10).
        top_features = list(default = 10L, choices = NULL, settled_by = workbook),
        # The weights' plot: scaled by the largest absolute weight of the factor in the data set (MOFA2's
        # plot_weights(scale = TRUE)), or as estimated.
        weight_scale = list(default = "scaled", choices = c("scaled", "raw"), settled_by = workbook),
        # The design column the factors are related to (the integration design holds the group only).
        design_covariate = list(default = "group", choices = "group", settled_by = workbook),
        # The factor-design statistic: one-way analysis of variance of each factor on the group, R^2 as the effect.
        association_test = list(default = "anova", choices = "anova", settled_by = workbook),
        # The p-values' adjustment across the factors tested.
        p_adjust_method = list(default = "BH", choices = c("BH", "holm", "bonferroni", "none"), settled_by = workbook),
        figure_format = list(default = "png", choices = c("png", "pdf"), settled_by = workbook),
        figure_width = list(default = 8, choices = NULL, settled_by = workbook),
        figure_height = list(default = 6, choices = NULL, settled_by = workbook),
        figure_dpi = list(default = 300L, choices = NULL, settled_by = workbook))
}

#' The Results parameter table (one row per parameter), for the ticket, the ledger and the hand-off.
#' @return `data.frame(parameter, default, choices, settled_by)`.
#' @noRd
integrationResultsParameterTable <- function() {
    spec <- integrationResultsParameterSpec()
    data.frame(parameter = names(spec),
        default = vapply(spec, function(entry) paste(format(entry$default), collapse = ", "), character(1)),
        choices = vapply(spec, function(entry) if (is.null(entry$choices)) "" else paste(entry$choices, collapse = ", "),
            character(1)),
        settled_by = vapply(spec, function(entry) entry$settled_by, character(1)),
        stringsAsFactors = FALSE, row.names = NULL)
}

#' Resolve the Results parameters (absent ones take the table's default) and check them. Refused in plain words
#' (`integration_parameters_invalid`): a value outside its choices, an unknown factor, a number that is not positive.
#' @param parameters The caller's parameters (named list).
#' @param factors The model's factor names.
#' @return The parameters, named in the table's order.
#' @noRd
integrationResolveResultsParameters <- function(parameters = list(), factors) {
    spec <- integrationResultsParameterSpec()
    invalid <- function(message) integrationRefuse(message, "integration_parameters_invalid")
    out <- list()
    for (name in names(spec)) {
        value <- parameters[[name, exact = TRUE]]
        out[name] <- list(if (is.null(value)) spec[[name]]$default else value)
        choices <- spec[[name]]$choices
        if (!is.null(choices) && !(is.character(out[[name]]) && length(out[[name]]) == 1L && isTRUE(out[[name]] %in% choices))) {
            invalid(sprintf("%s must be one of %s", name, paste(choices, collapse = ", ")))
        }
    }
    shown <- out$factors_shown
    if (!is.character(shown) || !length(shown) || anyNA(shown) ||
        !(identical(shown, "all") || all(shown %in% factors))) {
        invalid(sprintf("factors_shown must be all or factor names of the model (%s)", paste(factors, collapse = ", ")))
    }
    out$factors_shown <- if (identical(shown, "all")) factors else factors[factors %in% shown]
    positive <- function(name, whole) {
        value <- out[[name]]
        if (!(is.numeric(value) && length(value) == 1L && !is.na(value) && value > 0 && (!whole || value == round(value)))) {
            invalid(sprintf("%s must be a positive %s", name, if (whole) "whole number" else "number"))
        }
        if (whole) as.integer(value) else as.numeric(value)
    }
    out$top_features <- positive("top_features", TRUE)
    out$figure_width <- positive("figure_width", FALSE)
    out$figure_height <- positive("figure_height", FALSE)
    out$figure_dpi <- positive("figure_dpi", TRUE)
    out
}

#' Which data sets were corrected for unwanted variation with the groups as replicate sets (KIT-060 D10): a corrected
#' matrix whose correction names no technical-replicate structure (`replicate_structure` or
#' `ruv_replicate_structure`), except data sets the correction lists as not corrected. An unrecorded structure on a
#' corrected matrix counts as the groups (the cautious reading).
#' @param object The imported inputs state (`inputs` hydrated).
#' @param views The model's data sets.
#' @return Named logical by data set.
#' @noRd
integrationGroupCorrectedViews <- function(object, views) {
    out <- stats::setNames(rep(FALSE, length(views)), views)
    for (entry in object@inputs) {
        manifest <- entry$manifest
        mine <- intersect(entry$views, views)
        if (!length(mine) || !identical(manifest$matrix_role, "corrected") || !is.list(manifest$correction)) next
        correction <- manifest$correction
        structure <- c(integrationGuiScalar(correction$replicate_structure), integrationGuiScalar(correction$ruv_replicate_structure))
        if (any(structure %in% "technical_replicates")) next
        out[mine] <- TRUE
        if (is.list(correction$views)) {
            for (view in correction$views) {
                applied <- as.logical(integrationGuiScalar(view$applied))
                assay <- as.character(integrationGuiScalar(view$assay))
                named <- mine[mine == assay | mine == paste(entry$omic, assay, sep = "_")]
                if (length(named) && identical(applied, FALSE)) out[named] <- FALSE
            }
        }
    }
    out
}

#' The relation of each factor to the design's groups (KIT-063; ledger row "integration: factors against the design"):
#' per factor, over the biological samples with a group, a one-way analysis of variance of the factor's values on the
#' group; the share of the factor's variance between the groups (R^2 = between-group over total sum of squares, in
#' percent) is the effect size, the F test's p-value is reported raw and adjusted across the factors tested. Without a
#' test (`tested = FALSE`, D10) the p-values are NA; with fewer than two groups or no residual degree of freedom
#' nothing is computed.
#' @param values `data.frame(biological_sample, group, <factor columns>)`.
#' @param factors The factors to relate.
#' @param tested Whether p-values are given (FALSE under D10's claim limit).
#' @param adjust The `p.adjust` method.
#' @return `data.frame(factor, samples, groups, r2_percent, f_statistic, df_group, df_residual, p_value, p_adjusted,
#'   status)`; `status` is `tested`, `not_tested` (D10) or `cannot_test`.
#' @noRd
integrationFactorAssociation <- function(values, factors, tested = TRUE, adjust = "BH") {
    rows <- lapply(factors, function(name) {
        y <- values[[name]]
        group <- as.character(values$group)
        keep <- !is.na(y) & !is.na(group) & nzchar(group)
        y <- y[keep]
        group <- factor(group[keep])
        k <- nlevels(group)
        n <- length(y)
        if (k < 2L || n - k < 1L) {
            return(data.frame(factor = name, samples = n, groups = k, r2_percent = NA_real_, f_statistic = NA_real_,
                df_group = NA_integer_, df_residual = NA_integer_, p_value = NA_real_, status = "cannot_test",
                stringsAsFactors = FALSE))
        }
        table <- stats::anova(stats::lm(y ~ group))
        between <- table[["Sum Sq"]][[1L]]
        within <- table[["Sum Sq"]][[2L]]
        data.frame(factor = name, samples = n, groups = k, r2_percent = 100 * between / (between + within),
            f_statistic = table[["F value"]][[1L]], df_group = as.integer(table[["Df"]][[1L]]),
            df_residual = as.integer(table[["Df"]][[2L]]),
            p_value = if (isTRUE(tested)) table[["Pr(>F)"]][[1L]] else NA_real_,
            status = if (isTRUE(tested)) "tested" else "not_tested", stringsAsFactors = FALSE)
    })
    out <- do.call(rbind, rows)
    out$p_adjusted <- NA_real_
    done <- out$status == "tested" & !is.na(out$p_value)
    if (any(done)) out$p_adjusted[done] <- stats::p.adjust(out$p_value[done], method = adjust)
    rownames(out) <- NULL
    out[c("factor", "samples", "groups", "r2_percent", "f_statistic", "df_group", "df_residual", "p_value", "p_adjusted",
        "status")]
}

#' MOFA2's own warnings, kept with the model (`integrationTrainMofa`'s `messages`, prefix "WARNING:"), in plain words
#' (KIT-063). The warning that factors follow each sample's overall level ("strongly correlated with the average
#' expression of features", earlier "with the total number of expressed features") is named with its factors; any other warning is said to exist, its text as the detail.
#' @param value The model's product.
#' @param texts `integrationGuiTexts()`.
#' @return A list of `list(kind, sentence, detail)`.
#' @noRd
integrationModelAdvisories <- function(value, texts = integrationGuiTexts()) {
    raw <- as.character(value$messages %||% character())
    raw <- unique(trimws(sub("^WARNING:\\s*", "", raw[startsWith(raw, "WARNING:")])))
    raw <- raw[nzchar(raw)]
    lapply(raw, function(message) {
        # MOFA2 1.22.0's wording is "strongly correlated with the average expression of features" (MOFA2 lane run
        # 44616769378); earlier releases said "the total number of expressed features".
        if (grepl("correlated with the (average expression of features|total number of expressed features)", message)) {
            listed <- regmatches(message, regexpr("Factor\\(s\\)\\s*[0-9][0-9, ]*", message))
            numbers <- if (length(listed)) regmatches(listed, gregexpr("[0-9]+", listed))[[1L]] else character()
            factors <- if (length(numbers)) paste(paste0("Factor", numbers), collapse = ", ") else "?"
            return(list(kind = "overall_level", sentence = sprintf(texts[["advisory_overall_level"]], factors),
                detail = message))
        }
        list(kind = "other", sentence = texts[["advisory_other"]], detail = message)
    })
}

#' The standard results of a trained model (KIT-063), each a pure read of its saved product: (1) the variance each
#' factor explains per data set and in all per data set; (2) the factor values per biological sample with the group,
#' and each factor's relation to the group (`integrationFactorAssociation`), without p-values where D10's claim limit
#' applies; (3) every weight per factor and data set with the feature's display name, its scaled weight and its rank by
#' absolute weight, and the `top_features` first; (4) the correlation of the factors; (5) MOFA2's warnings in words;
#' and the claim limit's facts.
#' @param value The model's product (`integrationModelProduct`, v2).
#' @param object The imported inputs state, hydrated (display names, corrections).
#' @param parameters The Results parameters (absent ones take the table's defaults).
#' @param manager The owner of the store that holds the model (KIT-148: each data set's weights are read from it, one
#'   data set at a time).
#' @return A list of tables and facts (see the body).
#' @noRd
integrationResults <- function(value, object, parameters = list(), texts = integrationGuiTexts(), manager) {
    integrationAssertModelProduct(value)
    factors <- setdiff(names(value$factors), "biological_sample")
    p <- integrationResolveResultsParameters(parameters, factors)
    shown <- p$factors_shown
    views <- unique(value$variance_explained$view)
    # (1) Variance explained, per factor and data set, and all factors together per data set.
    variance <- value$variance_explained[value$variance_explained$factor %in% shown, c("factor", "view", "r2_percent")]
    rownames(variance) <- NULL
    total <- value$variance_total[c("view", "r2_percent")]
    # (2) Factor values with the design, and each factor's relation to the group.
    samples <- value$samples
    factor_values <- data.frame(biological_sample = value$factors$biological_sample,
        group = as.character(samples[[p$design_covariate]][match(value$factors$biological_sample, samples$biological_sample)]),
        value$factors[factors], check.names = FALSE, stringsAsFactors = FALSE)
    corrected <- integrationGroupCorrectedViews(object, views)
    association <- integrationFactorAssociation(factor_values, shown, tested = !any(corrected), adjust = p$p_adjust_method)
    # (3) Weights with display names, scaled and ranked per factor and data set.
    weights <- do.call(rbind, lapply(views, function(view) {
        w <- integrationModelViewTable(manager, value, "weights", view)
        features <- object@features[[view]]
        display <- if (!is.null(features) && "display_name" %in% names(features))
            as.character(features$display_name[match(w$feature_id, as.character(features$feature_id))]) else rep(NA_character_, nrow(w))
        display[is.na(display) | !nzchar(display)] <- w$feature_id[is.na(display) | !nzchar(display)]
        do.call(rbind, lapply(shown, function(name) {
            x <- w[[name]]
            largest <- max(abs(x), na.rm = TRUE)
            rank <- integer(length(x))
            rank[order(-abs(x), seq_along(x))] <- seq_along(x)
            data.frame(factor = name, view = view, rank = rank, feature_id = w$feature_id, display_name = display, weight = x,
                weight_scaled = if (is.finite(largest) && largest > 0) x / largest else 0, stringsAsFactors = FALSE)
        }))
    }))
    weights <- weights[order(match(weights$factor, shown), match(weights$view, views), weights$rank, method = "radix"), ,
        drop = FALSE]
    rownames(weights) <- NULL
    top <- weights[weights$rank <= p$top_features, , drop = FALSE]
    rownames(top) <- NULL
    # (4) The factors' correlation over the biological samples.
    correlation <- stats::cor(as.matrix(value$factors[factors]), use = "pairwise.complete.obs")
    # The claim limit's facts: each omic's correction in words, and the data sets corrected with the groups.
    corrections <- vapply(object@inputs, function(entry) integrationCorrectionText(entry$manifest, texts), character(1))
    list(parameters = p, factors = factors, shown = shown, views = views, variance = variance, variance_total = total,
        factor_values = factor_values, association = association, weights = weights, top_weights = top,
        correlation = correlation, advisories = integrationModelAdvisories(value, texts), notes = value$notes,
        shared_samples = value$shared_samples,
        claim = list(group_corrected = names(corrected)[corrected], corrections = corrections))
}

# ------------------------------------------------------------------------------------------------
# The figures (ggplot2, the package's theme_bw), each from `integrationResults()`'s tables.
# ------------------------------------------------------------------------------------------------

#' Variance explained per factor and data set, as a heatmap with the percentages written in.
#' @noRd
integrationPlotVariance <- function(results, texts = integrationGuiTexts()) {
    data <- results$variance
    data$factor <- factor(data$factor, levels = rev(results$shown))
    data$view <- factor(data$view, levels = results$views)
    ggplot2::ggplot(data, ggplot2::aes(x = .data$view, y = .data$factor, fill = .data$r2_percent)) +
        ggplot2::geom_tile(colour = "white") +
        ggplot2::geom_text(ggplot2::aes(label = sprintf("%.1f", .data$r2_percent)), size = 3.5) +
        ggplot2::scale_fill_gradient(low = "white", high = "steelblue", limits = c(0, NA), name = texts[["plot_variance"]]) +
        ggplot2::labs(x = texts[["plot_dataset"]], y = texts[["plot_factor"]]) +
        ggplot2::theme_bw() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1), panel.grid = ggplot2::element_blank())
}

#' Variance explained by all factors together, per data set.
#' @noRd
integrationPlotVarianceTotal <- function(results, texts = integrationGuiTexts()) {
    data <- results$variance_total
    data$view <- factor(data$view, levels = results$views)
    ggplot2::ggplot(data, ggplot2::aes(x = .data$view, y = .data$r2_percent)) +
        ggplot2::geom_col(fill = "steelblue") +
        ggplot2::geom_text(ggplot2::aes(label = sprintf("%.1f", .data$r2_percent)), vjust = -0.3, size = 3.5) +
        ggplot2::labs(x = texts[["plot_dataset"]], y = texts[["plot_variance"]]) +
        ggplot2::theme_bw() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
}

#' The factors' values per biological sample against the groups: one panel per factor shown, a box per group and the
#' samples as points (jittered with a fixed seed, so a saved figure is the same each time).
#' @noRd
integrationPlotFactorsByGroup <- function(results, texts = integrationGuiTexts()) {
    values <- results$factor_values
    data <- do.call(rbind, lapply(results$shown, function(name) data.frame(factor = name, group = values$group,
        value = values[[name]], stringsAsFactors = FALSE)))
    data <- data[!is.na(data$value), , drop = FALSE]
    data$group[is.na(data$group)] <- ""
    data$factor <- factor(data$factor, levels = results$shown)
    ggplot2::ggplot(data, ggplot2::aes(x = .data$group, y = .data$value, colour = .data$group)) +
        ggplot2::geom_boxplot(outlier.shape = NA) +
        ggplot2::geom_point(position = ggplot2::position_jitter(width = 0.15, height = 0, seed = 1L)) +
        ggplot2::facet_wrap(ggplot2::vars(.data$factor), scales = "free_y") +
        ggplot2::labs(x = texts[["plot_group"]], y = texts[["plot_factor_value"]]) +
        ggplot2::theme_bw() +
        ggplot2::theme(legend.position = "none")
}

#' The `top_features` largest weights of one factor in one data set, labelled with the features' display names (a
#' name shared by two features carries its identifier).
#' @param factor,view The factor and the data set.
#' @noRd
integrationPlotWeights <- function(results, factor, view, texts = integrationGuiTexts()) {
    data <- results$top_weights[results$top_weights$factor == factor & results$top_weights$view == view, , drop = FALSE]
    scaled <- identical(results$parameters$weight_scale, "scaled")
    data$shown <- if (scaled) data$weight_scaled else data$weight
    shared <- duplicated(data$display_name) | duplicated(data$display_name, fromLast = TRUE)
    data$label <- ifelse(shared, sprintf("%s (%s)", data$display_name, data$feature_id), data$display_name)
    data$label <- factor(data$label, levels = rev(data$label))
    ggplot2::ggplot(data, ggplot2::aes(x = .data$shown, y = .data$label, fill = .data$shown > 0)) +
        ggplot2::geom_col() +
        ggplot2::geom_vline(xintercept = 0, colour = "grey40") +
        ggplot2::scale_fill_manual(values = c(`TRUE` = "firebrick", `FALSE` = "steelblue"), guide = "none") +
        ggplot2::labs(x = texts[[if (scaled) "plot_weight_scaled" else "plot_weight"]], y = NULL,
            title = sprintf(texts[["plot_weights_title"]], factor, view)) +
        ggplot2::theme_bw()
}

#' The factors' correlation over the biological samples, as a heatmap with the coefficients written in.
#' @noRd
integrationPlotFactorCorrelation <- function(results, texts = integrationGuiTexts()) {
    correlation <- results$correlation[results$shown, results$shown, drop = FALSE]
    data <- data.frame(row = rep(rownames(correlation), times = ncol(correlation)),
        column = rep(colnames(correlation), each = nrow(correlation)), r = as.vector(correlation), stringsAsFactors = FALSE)
    data$row <- factor(data$row, levels = rev(results$shown))
    data$column <- factor(data$column, levels = results$shown)
    ggplot2::ggplot(data, ggplot2::aes(x = .data$column, y = .data$row, fill = .data$r)) +
        ggplot2::geom_tile(colour = "white") +
        ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f", .data$r)), size = 3.5) +
        ggplot2::scale_fill_gradient2(low = "steelblue", mid = "white", high = "firebrick", limits = c(-1, 1),
            name = texts[["plot_correlation"]]) +
        ggplot2::labs(x = NULL, y = NULL) +
        ggplot2::theme_bw() +
        ggplot2::theme(panel.grid = ggplot2::element_blank())
}

# ------------------------------------------------------------------------------------------------
# The tables as the tab shows them, the headers in words.
# ------------------------------------------------------------------------------------------------

#' The variance table: one row per factor, one column per data set, percentages to two decimals.
#' @noRd
integrationResultsVarianceTable <- function(results, texts = integrationGuiTexts()) {
    variance <- results$variance
    table <- data.frame(results$shown, stringsAsFactors = FALSE)
    names(table) <- texts[["col_factor"]]
    for (view in results$views) {
        table[[view]] <- sprintf("%.2f", variance$r2_percent[match(paste(results$shown, view), paste(variance$factor, variance$view))])
    }
    table
}

#' The relation of the factors to the groups, in words: a p-value not given under the claim limit, or not computable.
#' @noRd
integrationResultsAssociationTable <- function(results, texts = integrationGuiTexts()) {
    a <- results$association
    number <- function(x, format) ifelse(is.na(x), "", sprintf(format, x))
    p <- function(x) ifelse(a$status == "tested", ifelse(is.na(x), "", format.pval(x, digits = 3)),
        ifelse(a$status == "not_tested", texts[["cell_not_tested"]], texts[["cell_cannot_test"]]))
    table <- data.frame(a$factor, a$samples, a$groups, number(a$r2_percent, "%.1f"), number(a$f_statistic, "%.2f"),
        ifelse(is.na(a$df_group), "", sprintf("%d, %d", a$df_group, a$df_residual)), p(a$p_value), p(a$p_adjusted),
        stringsAsFactors = FALSE)
    names(table) <- unname(texts[c("col_factor", "col_samples_tested", "col_groups_n", "col_group_r2", "col_f", "col_df",
        "col_p", "col_p_adjusted")])
    table
}

#' The top features of one factor in one data set.
#' @noRd
integrationResultsWeightsTable <- function(results, factor, view, texts = integrationGuiTexts()) {
    top <- results$top_weights[results$top_weights$factor == factor & results$top_weights$view == view, , drop = FALSE]
    table <- data.frame(top$rank, top$display_name, top$feature_id, sprintf("%.4g", top$weight), sprintf("%.3f", top$weight_scaled),
        stringsAsFactors = FALSE)
    names(table) <- unname(texts[c("col_rank", "col_feature", "col_feature_id", "col_weight", "col_weight_scaled")])
    table
}

#' The factor correlation as a table, to two decimals.
#' @noRd
integrationResultsCorrelationTable <- function(results, texts = integrationGuiTexts()) {
    correlation <- results$correlation[results$shown, results$shown, drop = FALSE]
    table <- data.frame(rownames(correlation), matrix(sprintf("%.2f", correlation), nrow(correlation)), stringsAsFactors = FALSE)
    names(table) <- c(texts[["col_factor"]], colnames(correlation))
    table
}

#' The claim limit beside the factor-against-group figure (KIT-060 D10), in words: always that the factors are
#' descriptive; where a data set was corrected with the groups as replicate sets, that no p-value is given and why;
#' each omic's correction.
#' @return A character vector: the sentences, then one line per omic.
#' @noRd
integrationClaimLimitText <- function(results, texts = integrationGuiTexts()) {
    claim <- results$claim
    sentences <- texts[["claim_descriptive"]]
    if (length(claim$group_corrected)) {
        sentences <- c(sentences, sprintf(texts[["claim_group_corrected"]], paste(claim$group_corrected, collapse = ", ")))
    }
    c(sentences, sprintf(texts[["claim_correction"]], tools::toTitleCase(names(claim$corrections)), unname(claim$corrections)))
}

# ------------------------------------------------------------------------------------------------
# Saving the results to the project (the owner's action-output publication).
# ------------------------------------------------------------------------------------------------

#' Write the results' tables and figures into the staged folders (KIT-063): tables tab-separated in
#' `mofa_tables_dir`, figures in `mofa_plots_dir` in the table's format and size; every file checked
#' (`checkWorkflowFigureFiles`: a figure the device did not write stops the save by name).
#' @param results `integrationResults()`'s value.
#' @param staged The staged folders by role (`mofa_tables_dir`, `mofa_plots_dir`).
#' @return The display: the files written, by role.
#' @noRd
integrationWriteResults <- function(results, staged, texts = integrationGuiTexts()) {
    p <- results$parameters
    write <- function(table, name) {
        path <- file.path(staged$mofa_tables_dir, name)
        utils::write.table(table, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
        path
    }
    correlation <- data.frame(factor = rownames(results$correlation), results$correlation, check.names = FALSE,
        stringsAsFactors = FALSE, row.names = NULL)
    advisories <- data.frame(kind = vapply(results$advisories, `[[`, "", "kind"),
        sentence = vapply(results$advisories, `[[`, "", "sentence"),
        mofa2_text = gsub("[\t\r\n]+", " ", vapply(results$advisories, `[[`, "", "detail")), stringsAsFactors = FALSE)
    tables <- c(write(results$variance, "variance_explained.tsv"), write(results$variance_total, "variance_explained_total.tsv"),
        write(results$factor_values, "factor_values.tsv"), write(results$association, "factor_group_association.tsv"),
        write(results$weights, "weights.tsv"), write(results$top_weights, "top_weights.tsv"),
        write(correlation, "factor_correlation.tsv"), write(advisories, "model_warnings.tsv"))
    safe <- function(text) gsub("[^A-Za-z0-9_-]", "_", text)
    save <- function(plot, name) {
        path <- file.path(staged$mofa_plots_dir, paste0(safe(name), ".", p$figure_format))
        ggplot2::ggsave(path, plot, width = p$figure_width, height = p$figure_height, units = "in", dpi = p$figure_dpi)
        path
    }
    figures <- c(save(integrationPlotVariance(results, texts), "variance_explained"),
        save(integrationPlotVarianceTotal(results, texts), "variance_explained_total"),
        save(integrationPlotFactorsByGroup(results, texts), "factors_by_group"),
        save(integrationPlotFactorCorrelation(results, texts), "factor_correlation"))
    for (factor in results$shown) for (view in results$views) {
        figures <- c(figures, save(integrationPlotWeights(results, factor, view, texts), sprintf("weights_%s_%s", factor, view)))
    }
    checkWorkflowFigureFiles(c(tables, figures), "the integration results' files were not written")
    list(completed = TRUE, tables = basename(tables), figures = basename(figures))
}

#' Save a trained model's results to the project (KIT-063) through the owner's action-output publication: staged under
#' `<scripts>/.workflow-export`, drawn in a short folder (`withWorkflowShortStaging`; a Windows path over 259 characters
#' stops the save by name), journaled, then copied into `results/integration/mofa/{tables,plots}`. Its identity is the
#' model's product and the Results parameters: the same model and parameters saved again republish the same files.
#' @param manager The integration owner (the inputs state at its head).
#' @param model_reference The model product's reference.
#' @param results `integrationResults()`'s value for that model.
#' @param paths `project_dirs$integration` (`mofa_tables_dir`, `mofa_plots_dir`, `source_dir`, `results_dir`).
#' @return The display: the files written, by role, and the folders.
#' @noRd
integrationSaveResults <- function(manager, model_reference, results, paths, texts = integrationGuiTexts()) {
    destinations <- list(mofa_tables_dir = paths$mofa_tables_dir %||% file.path(paths$results_dir, "mofa", "tables"),
        mofa_plots_dir = paths$mofa_plots_dir %||% file.path(paths$results_dir, "mofa", "plots"))
    for (dir in destinations) dir.create(dir, recursive = TRUE, showWarnings = FALSE)
    identity <- digest::digest(list(schema = "multischolar.integration_results.v1", model = model_reference$entity_id,
        parameters = results$parameters), algo = "sha256")
    display <- withWorkflowStateScope(manager, function() {
        # The publication stages under the project's scripts folder, not inside the results folders.
        if (is.null(manager$.context$experiment_paths) && is.character(paths$source_dir)) {
            manager$.context$experiment_paths <- list(source_dir = paths$source_dir)
        }
        publishWorkflowActionOutputs(manager, getWorkflowStateReference(manager), INTEGRATION_RESULTS_LABEL, destinations,
            function(staging) withWorkflowShortStaging(staging, function(staged) integrationWriteResults(results, staged, texts),
                "the integration results could not be copied into the project's staging folder"),
            identity = identity)
    })
    # KIT-142: the saved results are kept with their model, for the analysis record.
    integrationRecordResultsPublication(manager, model_reference, results,
        list(tables = file.path(destinations$mofa_tables_dir, display$tables),
            figures = file.path(destinations$mofa_plots_dir, display$figures)), identity)
    c(display, list(tables_dir = destinations$mofa_tables_dir, plots_dir = destinations$mofa_plots_dir))
}
