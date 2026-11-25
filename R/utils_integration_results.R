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
# KIT-052). Every choice (how many features, which factors, the statistic, figure size and format) is a declared
# parameter with its default and the decision behind it in one table, `integrationResultsParameterSpec()`. Enrichment on the weights is KIT-141, the record and report KIT-142 (a saved
# publication is kept with its model for the record, `integrationRecordResultsPublication`).

#' The label of the results' publication on the owner's action-output journal.
#' @noRd
INTEGRATION_RESULTS_LABEL <- "integration-results"

#' The Results tab's parameters with their defaults, each with the decision that settled it (`settled_by`; KIT-240:
#' every default kept by the owner's decision VAL-221 of 2026-11-25, the values unchanged), in the form of
#' `integrationMofaParameterSpec()`.
#' @return Named list by parameter: `default`, `choices` (NULL for a number or a free value), `settled_by`.
#' @noRd
integrationResultsParameterSpec <- function() {
    kept <- "VAL-221 (decided 2026-11-25: the default kept)"
    list(
        # Which factors the figures and tables show: "all", or factor names (Factor1, ...).
        factors_shown = list(default = "all", choices = NULL, type = "factors", settled_by = kept),
        # How many features per factor and data set are listed and drawn (the starter workbook's plot_weights: 10).
        top_features = list(default = 10L, choices = NULL, type = "whole", settled_by = kept),
        # The weights' plot: scaled by the largest absolute weight of the factor in the data set (MOFA2's
        # plot_weights(scale = TRUE)), or as estimated.
        weight_scale = list(default = "scaled", choices = c("scaled", "raw"), settled_by = kept),
        # The design column the factors are related to (the integration design holds the group only).
        design_covariate = list(default = "group", choices = "group", settled_by = kept),
        # The factor-design statistic: one-way analysis of variance of each factor on the group, R^2 as the effect;
        # KIT-151 (plan item E): or the Kruskal-Wallis rank test, epsilon-squared as the effect. ANOVA stays the default.
        association_test = list(default = "anova", choices = c("anova", "kruskal_wallis"), settled_by = kept),
        # The p-values' adjustment across the factors tested.
        p_adjust_method = list(default = "BH", choices = c("BH", "holm", "bonferroni", "none"), settled_by = kept),
        figure_format = list(default = "png", choices = c("png", "pdf"), settled_by = kept),
        figure_width = list(default = 8, choices = NULL, type = "number", settled_by = kept),
        figure_height = list(default = 6, choices = NULL, type = "number", settled_by = kept),
        figure_dpi = list(default = 300L, choices = NULL, type = "whole", settled_by = kept),
        # KIT-152 (audit SCI-17; ledger row "integration: the values the top-feature scatter and heatmap show"): the
        # values handed to the model (the aligned table: runs averaged, the model's features, any per-feature scaling;
        # before MOFA2 centres and scales the data sets), the only choice.
        feature_values = list(default = "model_input", choices = "model_input", settled_by = kept),
        # KIT-152 (SCI-17): the heatmap's colours, display only: each feature centred on its mean, or also scaled to
        # standard deviation 1.
        heatmap_scaling = list(default = "centre", choices = c("centre", "centre_scale"), settled_by = kept),
        # KIT-152 (INT-9, audit SCI-15; ledger row "integration: UMAP and t-SNE maps of the factors"): the maps are off
        # unless chosen; uwot's default number of neighbours; a perplexity valid from 17 biological samples; one seed.
        umap = list(default = FALSE, choices = NULL, type = "flag", settled_by = kept),
        umap_n_neighbors = list(default = 15L, choices = NULL, type = "whole", settled_by = kept),
        tsne = list(default = FALSE, choices = NULL, type = "flag", settled_by = kept),
        tsne_perplexity = list(default = 5, choices = NULL, type = "number", settled_by = kept),
        map_seed = list(default = 20261119L, choices = NULL, type = "whole", settled_by = kept))
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
    # KIT-152: the maps' switches and settings (their validity for the number of samples is checked when a map is drawn).
    for (name in c("umap", "tsne")) {
        if (!(is.logical(out[[name]]) && length(out[[name]]) == 1L && !is.na(out[[name]]))) {
            invalid(sprintf("%s must be TRUE or FALSE", name))
        }
    }
    out$umap_n_neighbors <- positive("umap_n_neighbors", TRUE)
    out$tsne_perplexity <- positive("tsne_perplexity", FALSE)
    out$map_seed <- positive("map_seed", TRUE)
    out
}

#' Which data sets were corrected for unwanted variation with the groups as replicate sets (KIT-060 D10): a corrected
#' matrix whose correction names no technical-replicate structure (`replicate_structure` or
#' `ruv_replicate_structure`), except data sets the correction lists as not corrected. An unrecorded structure on a
#' corrected matrix counts as the groups (the cautious reading). KIT-153: the matrix the model used decides; with
#' `input_matrix = uncorrected` no data set was corrected, with `ruv_corrected` each data set the omic corrected was; a
#' model recorded before KIT-153 (no `input_matrix`) is read by its manifest's role, as before.
#' @param object The imported inputs state (`inputs` hydrated).
#' @param views The model's data sets.
#' @param input_matrix The model's `input_matrix`, or NULL for a model recorded before KIT-153.
#' @return Named logical by data set.
#' @noRd
integrationGroupCorrectedViews <- function(object, views, input_matrix = NULL) {
    out <- stats::setNames(rep(FALSE, length(views)), views)
    if (identical(input_matrix, "uncorrected")) return(out)
    for (entry in object@inputs) {
        manifest <- entry$manifest
        mine <- intersect(entry$views, views)
        used_corrected <- if (is.null(input_matrix)) identical(manifest$matrix_role, "corrected") else
            identical(input_matrix, "ruv_corrected")
        if (!length(mine) || !used_corrected || !is.list(manifest$correction)) next
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
#' percent) is the effect size, the F test's p-value is reported raw and adjusted across every factor of the model
#' (KIT-151, INT-20: the family is the model's factors, whichever are shown; `adjusted_over` says how many). Without a
#' test (`tested = FALSE`, D10) the p-values are NA; with fewer than two groups or no residual degree of freedom
#' nothing is computed. KIT-151 (plan item E): `test = "kruskal_wallis"` uses the Kruskal-Wallis rank test instead
#' (`stats::kruskal.test`, its chi-squared approximation with the tie correction), its H and degrees of freedom, and
#' epsilon-squared, H / (n - 1) (in percent), as the effect size; the claim limit applies alike.
#' @param values `data.frame(biological_sample, group, <factor columns>)`.
#' @param factors Every factor of the model: all are tested and the p-values adjusted over them.
#' @param tested Whether p-values are given (FALSE under D10's claim limit).
#' @param adjust The `p.adjust` method.
#' @param shown The factors whose rows are returned (default every factor), in the order of `factors`.
#' @param test `anova` or `kruskal_wallis`.
#' @return `data.frame(factor, samples, groups, test, r2_percent, f_statistic, df_group, df_residual,
#'   epsilon_squared_percent, h_statistic, p_value, p_adjusted, adjusted_over, status)`; the ANOVA's columns are NA under
#'   Kruskal-Wallis and its columns NA under ANOVA; `status` is `tested`, `not_tested` (D10) or `cannot_test`.
#' @noRd
integrationFactorAssociation <- function(values, factors, tested = TRUE, adjust = "BH", shown = factors, test = "anova") {
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
                df_group = NA_integer_, df_residual = NA_integer_, epsilon_squared_percent = NA_real_, h_statistic = NA_real_,
                p_value = NA_real_, status = "cannot_test", stringsAsFactors = FALSE))
        }
        if (identical(test, "kruskal_wallis")) {
            kw <- stats::kruskal.test(y, group)
            h <- unname(kw$statistic)
            return(data.frame(factor = name, samples = n, groups = k, r2_percent = NA_real_, f_statistic = NA_real_,
                df_group = as.integer(kw$parameter), df_residual = NA_integer_, epsilon_squared_percent = 100 * h / (n - 1),
                h_statistic = h, p_value = if (isTRUE(tested)) kw$p.value else NA_real_,
                status = if (isTRUE(tested)) "tested" else "not_tested", stringsAsFactors = FALSE))
        }
        table <- stats::anova(stats::lm(y ~ group))
        between <- table[["Sum Sq"]][[1L]]
        within <- table[["Sum Sq"]][[2L]]
        data.frame(factor = name, samples = n, groups = k, r2_percent = 100 * between / (between + within),
            f_statistic = table[["F value"]][[1L]], df_group = as.integer(table[["Df"]][[1L]]),
            df_residual = as.integer(table[["Df"]][[2L]]), epsilon_squared_percent = NA_real_, h_statistic = NA_real_,
            p_value = if (isTRUE(tested)) table[["Pr(>F)"]][[1L]] else NA_real_,
            status = if (isTRUE(tested)) "tested" else "not_tested", stringsAsFactors = FALSE)
    })
    out <- do.call(rbind, rows)
    out$p_adjusted <- NA_real_
    done <- out$status == "tested" & !is.na(out$p_value)
    if (any(done)) out$p_adjusted[done] <- stats::p.adjust(out$p_value[done], method = adjust)
    out$adjusted_over <- sum(done)
    # KIT-151 (INT-20): adjusted over every factor, then the factors shown.
    out <- out[out$factor %in% shown, , drop = FALSE]
    out$test <- rep(test, nrow(out))
    rownames(out) <- NULL
    out[c("factor", "samples", "groups", "test", "r2_percent", "f_statistic", "df_group", "df_residual",
        "epsilon_squared_percent", "h_statistic", "p_value", "p_adjusted", "adjusted_over", "status")]
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
#' and the claim limit's facts. KIT-152 adds (6) the top features' values in every biological sample
#' (`integrationTopFeatureValues`) with the scale they are on, for the scatter and the heatmap, and (7) the UMAP and
#' t-SNE maps of the factors when chosen (`integrationFactorMap`; a map refused for its settings is kept as its refusal,
#' so the other results still show).
#' @param value The model's product (`integrationModelProduct`, v2).
#' @param object The imported inputs state, hydrated (display names, corrections).
#' @param parameters The Results parameters (absent ones take the table's defaults).
#' @param manager The owner of the store that holds the model (KIT-148: each data set's weights are read from it, one
#'   data set at a time).
#' @param model_parameters The model's resolved parameters (for the per-feature scaling the values are on); NULL or
#'   empty takes those of the model accepted on the owner's head.
#' @return A list of tables and facts (see the body).
#' @noRd
integrationResults <- function(value, object, parameters = list(), texts = integrationGuiTexts(), manager,
    model_parameters = NULL) {
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
    # KIT-153: the claim limit follows the matrix the model used (its `input_matrix`).
    if (!length(model_parameters)) model_parameters <- manager$.context$snapshot$exports[[INTEGRATION_MODEL_OPERATION]]$parameters
    input_matrix <- model_parameters$input_matrix
    corrected <- integrationGroupCorrectedViews(object, views, input_matrix)
    # KIT-151 (INT-20): every factor tested and adjusted together; the factors shown kept.
    association <- integrationFactorAssociation(factor_values, factors, tested = !any(corrected), adjust = p$p_adjust_method,
        shown = shown, test = p$association_test)
    # (3) Weights with display names, scaled and ranked per factor and data set. KIT-152: each data set's feature order
    # (the weights' rows, which are the aligned values' rows) is kept to find the top features' rows.
    feature_rows <- list()
    weights <- do.call(rbind, lapply(views, function(view) {
        w <- integrationModelViewTable(manager, value, "weights", view)
        feature_rows[[view]] <<- w$feature_id
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
    # KIT-151 (INT-14): each factor's agreement across the training seeds (NULL for a model trained before KIT-151).
    stability <- value$seeds$stability
    if (!is.null(stability)) {
        stability <- stability[stability$factor %in% shown, , drop = FALSE]
        rownames(stability) <- NULL
    }
    # The claim limit's facts: each omic's correction in words, and the data sets corrected with the groups.
    corrections <- vapply(object@inputs, function(entry) integrationCorrectionText(entry$manifest, texts), character(1))
    # (6) KIT-152: the top features' values, and the scale each data set's values are on.
    feature_values <- integrationTopFeatureValues(manager, value, top, feature_rows, factor_values, factors)
    # (7) KIT-152: the maps, when chosen.
    maps <- list()
    for (method in c("umap", "tsne")[c(p$umap, p$tsne)]) {
        maps[[method]] <- tryCatch(integrationFactorMap(factor_values, factors, method, p, texts),
            integration_map_invalid = function(condition) list(refusal = conditionMessage(condition)))
    }
    list(parameters = p, factors = factors, shown = shown, views = views, variance = variance, variance_total = total,
        factor_values = factor_values, association = association, weights = weights, top_weights = top,
        correlation = correlation, advisories = integrationModelAdvisories(value, texts), notes = value$notes,
        shared_samples = value$shared_samples,
        claim = list(group_corrected = names(corrected)[corrected], corrections = corrections,
            uncorrected = identical(input_matrix, "uncorrected")),
        feature_values = feature_values, value_scales = integrationFeatureValueScales(object, views, model_parameters$scaling),
        maps = maps, seeds = value$seeds$seeds, stability = stability,
        # KIT-240 (VAL-221 (1)): the reported model's rank among the seeds by final ELBO (NULL for a model trained before).
        elbo_rank = value$seeds$elbo_rank)
}

# ------------------------------------------------------------------------------------------------
# KIT-152: the top features' values and the maps of the factors, from the stored tables only (no MOFA2 object).
# ------------------------------------------------------------------------------------------------

#' The scale of the values handed to the model, per data set (KIT-152; SCI-17): `feature_zscore` when the model scaled
#' each feature, else the scale the omic saved (`log2`, `log2_cpm`), `unknown` when it saved none.
#' @param object The imported inputs state.
#' @param views The model's data sets.
#' @param scaling The model's `scaling` parameter.
#' @return `data.frame(view, scale)`.
#' @noRd
integrationFeatureValueScales <- function(object, views, scaling) {
    saved <- character()
    for (entry in object@inputs) saved[entry$views] <- as.character(integrationGuiScalar(entry$manifest$scale))
    scale <- if (identical(scaling, "feature_zscore")) rep("feature_zscore", length(views)) else unname(saved[views])
    scale[is.na(scale) | !nzchar(scale)] <- "unknown"
    data.frame(view = views, scale = scale, stringsAsFactors = FALSE)
}

#' Rows of one data set's aligned values, read by position (KIT-152; audit ENG-13): the store reads at most 10,000 rows
#' and 64 columns at a time (`readWorkflowSelectedOrdinals`), so the rows are read in pages of 10,000 and the biological
#' samples in blocks of 64; each page's feature identifiers are read back and must be the ones asked for.
#' @param manager The owner of the store that holds the model.
#' @param value The model product (v2).
#' @param view The data set.
#' @param rows The rows' positions (1-based, in the model's feature order).
#' @param features The features expected at those rows.
#' @return A matrix: one row per position (named by feature), one column per biological sample.
#' @noRd
integrationModelAlignedRows <- function(manager, value, view, rows, features) {
    reference <- integrationModelTableReference(value, "aligned_values", view)
    withWorkflowStateScope(manager, function() {
        store <- manager$.context$store
        scope <- newWorkflowTableReadScope(store, reference, reader = "duckdb")
        on.exit(releaseWorkflowTableReadScope(scope), add = TRUE)
        samples <- setdiff(scope$metadata$columns, "feature_id")
        out <- matrix(NA_real_, length(rows), length(samples), dimnames = list(features, samples))
        for (page in split(seq_along(rows), ceiling(seq_along(rows) / 10000L))) {
            read <- readWorkflowSelectedOrdinals(store, reference, "feature_id", rows[page], scope = scope)$feature_id
            if (!identical(as.character(read), features[page])) {
                integrationRefuse(sprintf("the stored values of the data set %s do not match its weights", view),
                    "integration_model_outdated")
            }
            for (block in split(samples, ceiling(seq_along(samples) / 64L))) {
                part <- readWorkflowSelectedOrdinals(store, reference, block, rows[page], scope = scope)
                out[page, block] <- as.matrix(part[block])
            }
        }
        out
    })
}

#' The top features' values in every biological sample (KIT-152; SCI-17): for each factor shown and data set, the
#' `top_features` features ranked by absolute weight (`top_weights`), each with its value as handed to the model (the
#' stored aligned table, read by keyed rows: never a whole data set) and the factor's value, one row per feature and
#' sample; missing values stay missing.
#' @param manager,value The store's owner and the model product.
#' @param top The top weights (`integrationResults()`'s `top_weights`).
#' @param feature_rows Per data set, its feature identifiers in the model's order (the weights' rows).
#' @param factor_values `integrationResults()`'s `factor_values`.
#' @param factors The model's factors.
#' @return `data.frame(factor, view, rank, feature_id, display_name, weight, biological_sample, group, factor_value,
#'   value)`, in the order of `top` and then of the samples.
#' @noRd
integrationTopFeatureValues <- function(manager, value, top, feature_rows, factor_values, factors) {
    samples <- factor_values$biological_sample
    blocks <- lapply(unique(top$view), function(view) {
        features <- unique(top$feature_id[top$view == view])
        m <- integrationModelAlignedRows(manager, value, view, match(features, feature_rows[[view]]), features)
        rownames(m) <- paste(view, features, sep = "\r")
        m[, samples, drop = FALSE]
    })
    values <- do.call(rbind, blocks)
    i <- rep(seq_len(nrow(top)), each = length(samples))
    j <- rep(seq_along(samples), times = nrow(top))
    factor_matrix <- as.matrix(factor_values[factors])
    out <- data.frame(factor = top$factor[i], view = top$view[i], rank = top$rank[i], feature_id = top$feature_id[i],
        display_name = top$display_name[i], weight = top$weight[i], biological_sample = samples[j],
        group = factor_values$group[j], factor_value = factor_matrix[cbind(j, match(top$factor[i], factors))],
        value = values[cbind(match(paste(top$view[i], top$feature_id[i], sep = "\r"), rownames(values)), j)],
        stringsAsFactors = FALSE)
    rownames(out) <- NULL
    out
}

#' The heatmap's values for one factor and data set (KIT-152; SCI-17), display only: the top features' values, each
#' feature centred on its mean over the samples where it is observed and, with `centre_scale`, divided by its standard
#' deviation there (a feature with no spread is left centred); features in order of rank, samples in order of the
#' factor's value.
#' @param results `integrationResults()`'s value.
#' @param factor,view The factor and the data set.
#' @return A matrix: features (named by display name and rank) by biological samples.
#' @noRd
integrationHeatmapValues <- function(results, factor, view) {
    data <- results$feature_values[results$feature_values$factor == factor & results$feature_values$view == view, ,
        drop = FALSE]
    features <- unique(data$feature_id[order(data$rank)])
    order_samples <- results$factor_values$biological_sample[order(results$factor_values[[factor]],
        results$factor_values$biological_sample, method = "radix")]
    m <- matrix(NA_real_, length(features), length(order_samples), dimnames = list(features, order_samples))
    m[cbind(match(data$feature_id, features), match(data$biological_sample, order_samples))] <- data$value
    m <- m - rowMeans(m, na.rm = TRUE)
    if (identical(results$parameters$heatmap_scaling, "centre_scale")) {
        spread <- apply(m, 1L, stats::sd, na.rm = TRUE)
        spread[!is.finite(spread) | spread == 0] <- 1
        m <- m / spread
    }
    m
}

#' A map of the factors (KIT-152; INT-9, audit SCI-15): UMAP (`uwot::umap`) or t-SNE (`Rtsne::Rtsne`) of the biological
#' samples with every factor observed, on all the model's factors; single-threaded, seeded (`withr::with_seed`, R's
#' default generators), every setting recorded. Refused in plain words (`integration_map_invalid`) before anything is
#' computed: fewer than two factors; for UMAP a number of neighbours below 2 or above the samples minus 1 (uwot would
#' lower a larger number to the number of samples itself, silently changing the setting; at least one other sample is
#' kept outside each neighbourhood); for t-SNE a perplexity whose triple is not below the samples minus 1 (Rtsne's own
#' limit, `3 * perplexity < nrow(X) - 1`). t-SNE is exact (theta 0), started from the first two principal components
#' of the factors scaled to a standard deviation of 1e-4 (Kobak and Berens 2019), with Rtsne's usual early
#' exaggeration (250 iterations). The software missing is refused in words too.
#' @param factor_values `integrationResults()`'s `factor_values`.
#' @param factors The model's factors.
#' @param method `umap` or `tsne`.
#' @param p The resolved Results parameters.
#' @return `list(settings, coordinates = data.frame(biological_sample, group, map_1, map_2), notes)`: `notes` holds any
#'   warning the software gave.
#' @noRd
integrationFactorMap <- function(factor_values, factors, method = c("umap", "tsne"), p, texts = integrationGuiTexts()) {
    method <- match.arg(method)
    refuse <- function(sentence) integrationGuiRefuse(sentence, "integration_map_invalid")
    x <- as.matrix(factor_values[factors])
    complete <- stats::complete.cases(x)
    x <- x[complete, , drop = FALSE]
    n <- nrow(x)
    if (length(factors) < 2L) refuse(sprintf(texts[["map_refused_factors"]], length(factors)))
    package <- if (method == "umap") "uwot" else "Rtsne"
    if (method == "umap" && (p$umap_n_neighbors < 2L || p$umap_n_neighbors > n - 1L)) {
        refuse(sprintf(texts[["map_refused_neighbors"]], p$umap_n_neighbors, n, n - 1L))
    }
    if (method == "tsne" && !(p$tsne_perplexity > 0 && 3 * p$tsne_perplexity < n - 1)) {
        refuse(sprintf(texts[["map_refused_perplexity"]], format(p$tsne_perplexity), n, format(round((n - 1) / 3, 2))))
    }
    if (!optionalPackageAvailable(package)) {
        refuse(sprintf(texts[["map_unavailable"]], texts[[paste0("map_", method)]], package))
    }
    requireOptionalPackage(package, "the maps of the integration factors")
    common <- list(method = method, package = package, version = as.character(utils::packageVersion(package)),
        samples = n, factors = length(factors), seed = p$map_seed, threads = 1L,
        random_generator = withr::with_seed(p$map_seed, paste(RNGkind(), collapse = ", ")))
    notes <- character()
    keep <- function(warning) {
        notes <<- c(notes, conditionMessage(warning))
        invokeRestart("muffleWarning")
    }
    if (method == "umap") {
        settings <- c(common, list(n_neighbors = p$umap_n_neighbors, n_components = 2L, metric = "euclidean", min_dist = 0.01,
            spread = 1, n_epochs = 500L, init = "spectral", n_sgd_threads = 1L))
        xy <- withCallingHandlers(withr::with_seed(p$map_seed, uwot::umap(x, n_neighbors = settings$n_neighbors,
            n_components = 2L, metric = "euclidean", min_dist = settings$min_dist, spread = settings$spread,
            n_epochs = settings$n_epochs, init = settings$init, n_threads = 1L, n_sgd_threads = 1L, verbose = FALSE)),
            warning = keep)
    } else {
        pcs <- stats::prcomp(x, center = TRUE, scale. = FALSE)$x[, 1:2, drop = FALSE]
        start <- pcs / stats::sd(pcs[, 1L]) * 1e-4
        settings <- c(common, list(perplexity = p$tsne_perplexity, dims = 2L, theta = 0, max_iter = 1000L,
            initialisation = "first two principal components, standard deviation 1e-4", stop_lying_iter = 250L,
            mom_switch_iter = 250L, exaggeration_factor = 12, eta = 200, normalize = TRUE, pca = FALSE))
        xy <- withCallingHandlers(withr::with_seed(p$map_seed, Rtsne::Rtsne(x, dims = 2L, perplexity = settings$perplexity,
            theta = 0, check_duplicates = FALSE, pca = FALSE, max_iter = settings$max_iter, Y_init = start,
            stop_lying_iter = 250L, mom_switch_iter = 250L, exaggeration_factor = 12, eta = 200, normalize = TRUE,
            num_threads = 1L, verbose = FALSE)$Y), warning = keep)
    }
    list(settings = settings, notes = notes,
        coordinates = data.frame(biological_sample = factor_values$biological_sample[complete],
            group = factor_values$group[complete], map_1 = unname(xy[, 1L]), map_2 = unname(xy[, 2L]),
            stringsAsFactors = FALSE))
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

#' KIT-152: the factor pair plot. Each pair of the factors shown in its own panel (columns the first factor, rows the
#' second), each biological sample a point coloured by its group; with fewer than two factors shown, a sentence saying so.
#' @noRd
integrationPlotFactorPairs <- function(results, texts = integrationGuiTexts()) {
    shown <- results$shown
    if (length(shown) < 2L) {
        return(ggplot2::ggplot() + ggplot2::annotate("text", x = 0, y = 0, label = texts[["plot_pairs_one"]]) +
            ggplot2::theme_void())
    }
    values <- results$factor_values
    pairs <- utils::combn(shown, 2L)
    data <- do.call(rbind, lapply(seq_len(ncol(pairs)), function(k) data.frame(x_factor = pairs[1L, k],
        y_factor = pairs[2L, k], x = values[[pairs[1L, k]]], y = values[[pairs[2L, k]]], group = values$group,
        stringsAsFactors = FALSE)))
    data$x_factor <- factor(data$x_factor, levels = shown)
    data$y_factor <- factor(data$y_factor, levels = shown)
    ggplot2::ggplot(data, ggplot2::aes(x = .data$x, y = .data$y, colour = .data$group)) +
        ggplot2::geom_point(size = 2) +
        ggplot2::facet_grid(ggplot2::vars(.data$y_factor), ggplot2::vars(.data$x_factor), scales = "free") +
        ggplot2::labs(x = texts[["plot_factor_value"]], y = texts[["plot_factor_value"]], colour = texts[["plot_group"]]) +
        ggplot2::theme_bw()
}

#' KIT-152 (SCI-17): the top features' values against the factor, one panel per feature (in order of rank, named by
#' rank and display name), each biological sample a point coloured by its group; the axis says which values are shown.
#' @param factor,view The factor and the data set.
#' @noRd
integrationPlotTopFeatureScatter <- function(results, factor, view, texts = integrationGuiTexts()) {
    data <- results$feature_values[results$feature_values$factor == factor & results$feature_values$view == view &
        !is.na(results$feature_values$value), , drop = FALSE]
    labels <- unique(data[order(data$rank), c("rank", "feature_id", "display_name")])
    labels$label <- sprintf("%d. %s", labels$rank, labels$display_name)
    shared <- duplicated(labels$label) | duplicated(labels$label, fromLast = TRUE)
    labels$label[shared] <- sprintf("%s (%s)", labels$label[shared], labels$feature_id[shared])
    data$panel <- factor(labels$label[match(data$feature_id, labels$feature_id)], levels = labels$label)
    scale <- results$value_scales$scale[match(view, results$value_scales$view)]
    ggplot2::ggplot(data, ggplot2::aes(x = .data$factor_value, y = .data$value, colour = .data$group)) +
        ggplot2::geom_point(size = 1.8) +
        ggplot2::facet_wrap(ggplot2::vars(.data$panel), scales = "free_y") +
        ggplot2::labs(x = sprintf(texts[["plot_scatter_x"]], factor), y = texts[[paste0("plot_values_", scale)]],
            colour = texts[["plot_group"]], title = sprintf(texts[["plot_scatter_title"]], factor, view),
            caption = texts[["plot_scatter_caption"]]) +
        ggplot2::theme_bw()
}

#' KIT-152 (SCI-17): the top features' heatmap, display only (`integrationHeatmapValues`): features in order of rank,
#' biological samples in order of the factor's value, the legend saying how the values were centred (and scaled).
#' @param factor,view The factor and the data set.
#' @noRd
integrationPlotTopFeatureHeatmap <- function(results, factor, view, texts = integrationGuiTexts()) {
    m <- integrationHeatmapValues(results, factor, view)
    data <- results$feature_values[results$feature_values$factor == factor & results$feature_values$view == view, ,
        drop = FALSE]
    names <- data$display_name[match(rownames(m), data$feature_id)]
    shared <- duplicated(names) | duplicated(names, fromLast = TRUE)
    names[shared] <- sprintf("%s (%s)", names[shared], rownames(m)[shared])
    tiles <- data.frame(feature = rep(names, times = ncol(m)), sample = rep(colnames(m), each = nrow(m)), value = as.vector(m),
        stringsAsFactors = FALSE)
    tiles$feature <- factor(tiles$feature, levels = rev(names))
    tiles$sample <- factor(tiles$sample, levels = colnames(m))
    ggplot2::ggplot(tiles, ggplot2::aes(x = .data$sample, y = .data$feature, fill = .data$value)) +
        ggplot2::geom_tile() +
        ggplot2::scale_fill_gradient2(low = "steelblue", mid = "white", high = "firebrick", midpoint = 0, na.value = "grey85") +
        ggplot2::labs(x = sprintf(texts[["plot_heatmap_x"]], factor), y = NULL,
            fill = texts[[paste0("plot_heatmap_", results$parameters$heatmap_scaling)]],
            title = sprintf(texts[["plot_scatter_title"]], factor, view), caption = texts[["plot_heatmap_caption"]]) +
        ggplot2::theme_bw() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5),
            panel.grid = ggplot2::element_blank())
}

#' KIT-152 (SCI-15): the caution beside the maps, in words: always that they are a picture with the number of samples
#' on them; per map shown, what it depends on (UMAP: the seed and the neighbours; t-SNE, started from the principal
#' components: the perplexity, not the seed); where a data set was corrected with the groups as replicate sets, that
#' groups separate by construction.
#' @param results `integrationResults()`'s value.
#' @param samples The biological samples on the map.
#' @param methods The maps shown (`umap`, `tsne`).
#' @return A character vector of sentences.
#' @noRd
integrationMapCaution <- function(results, samples, methods, texts = integrationGuiTexts()) {
    corrected <- results$claim$group_corrected
    c(sprintf(texts[["maps_caution"]], as.integer(samples)),
        unname(texts[paste0("maps_caution_", intersect(c("umap", "tsne"), methods))]),
        if (length(corrected)) sprintf(texts[["maps_caution_group_corrected"]], paste(corrected, collapse = ", ")))
}

#' KIT-152: one map of the factors (`integrationFactorMap`), each biological sample a point coloured by its group, its
#' settings in the title and the caution under it.
#' @param method `umap` or `tsne`.
#' @noRd
integrationPlotFactorMap <- function(results, method, texts = integrationGuiTexts()) {
    map <- results$maps[[method]]
    name <- texts[[paste0("map_", method)]]
    s <- map$settings
    details <- if (method == "umap") sprintf(texts[["map_settings_umap"]], as.integer(s$n_neighbors), as.integer(s$seed)) else
        sprintf(texts[["map_settings_tsne"]], format(s$perplexity), as.integer(s$seed))
    caution <- paste(strwrap(paste(integrationMapCaution(results, s$samples, method, texts), collapse = " "), 110L), collapse = "\n")
    ggplot2::ggplot(map$coordinates, ggplot2::aes(x = .data$map_1, y = .data$map_2, colour = .data$group)) +
        ggplot2::geom_point(size = 2.5) +
        ggplot2::labs(x = sprintf(texts[["plot_map_axis"]], name, 1L), y = sprintf(texts[["plot_map_axis"]], name, 2L),
            colour = texts[["plot_group"]], title = sprintf(texts[["plot_map_title"]], name, as.integer(s$factors),
                as.integer(s$samples), details), caption = caution) +
        ggplot2::theme_bw() +
        ggplot2::theme(plot.caption = ggplot2::element_text(hjust = 0))
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

#' Each factor's agreement across the training seeds (KIT-151; INT-14), to three decimals; NULL when not measured.
#' @noRd
integrationResultsStabilityTable <- function(results, texts = integrationGuiTexts()) {
    stability <- results$stability
    if (is.null(stability) || !nrow(stability) || all(is.na(stability$agreement_lowest))) return(NULL)
    table <- data.frame(stability$factor, sprintf("%.3f", stability$agreement_lowest), sprintf("%.3f", stability$agreement_median),
        stringsAsFactors = FALSE)
    names(table) <- unname(texts[c("col_factor", "col_agreement_lowest", "col_agreement_median")])
    table
}

#' The relation of the factors to the groups, in words: a p-value not given under the claim limit, or not computable.
#' @noRd
integrationResultsAssociationTable <- function(results, texts = integrationGuiTexts()) {
    a <- results$association
    number <- function(x, format) ifelse(is.na(x), "", sprintf(format, x))
    p <- function(x) ifelse(a$status == "tested", ifelse(is.na(x), "", format.pval(x, digits = 3)),
        ifelse(a$status == "not_tested", texts[["cell_not_tested"]], texts[["cell_cannot_test"]]))
    # KIT-151 (plan item E): the Kruskal-Wallis test's own effect size, statistic and degrees of freedom.
    if (identical(results$parameters$association_test, "kruskal_wallis")) {
        table <- data.frame(a$factor, a$samples, a$groups, number(a$epsilon_squared_percent, "%.1f"),
            number(a$h_statistic, "%.2f"), ifelse(is.na(a$df_group), "", sprintf("%d", a$df_group)), p(a$p_value),
            p(a$p_adjusted), stringsAsFactors = FALSE)
        names(table) <- unname(texts[c("col_factor", "col_samples_tested", "col_groups_n", "col_epsilon2", "col_h", "col_df",
            "col_p", "col_p_adjusted")])
        return(table)
    }
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
#' KIT-153: where the model used the uncorrected matrices, that unwanted variation correlated with the groups was not
#' removed; each omic's correction.
#' @return A character vector: the sentences, then one line per omic.
#' @noRd
integrationClaimLimitText <- function(results, texts = integrationGuiTexts()) {
    claim <- results$claim
    sentences <- texts[["claim_descriptive"]]
    if (length(claim$group_corrected)) {
        sentences <- c(sentences, sprintf(texts[["claim_group_corrected"]], paste(claim$group_corrected, collapse = ", ")))
    }
    if (isTRUE(as.logical(integrationGuiScalar(claim$uncorrected)))) sentences <- c(sentences, texts[["claim_uncorrected"]])
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
        write(correlation, "factor_correlation.tsv"), write(advisories, "model_warnings.tsv"),
        # KIT-152: the top features' values (the scatter's and, before display centring, the heatmap's).
        write(results$feature_values, "top_feature_values.tsv"),
        # KIT-151 (INT-14): each factor's agreement across the training seeds.
        if (!is.null(results$stability)) write(results$stability, "factor_stability.tsv"))
    maps <- Filter(function(map) is.null(map$refusal), results$maps)
    for (method in names(maps)) tables <- c(tables, write(maps[[method]]$coordinates, sprintf("map_%s.tsv", method)))
    if (length(maps)) {
        settings <- do.call(rbind, lapply(names(maps), function(method) {
            s <- maps[[method]]$settings
            data.frame(method = method, setting = names(s), value = vapply(s, function(x) paste(format(x), collapse = ", "), ""),
                stringsAsFactors = FALSE, row.names = NULL)
        }))
        tables <- c(tables, write(settings, "map_settings.tsv"))
    }
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
    # KIT-152: the pair plot, the top features' scatter and heatmap per factor and data set, the maps chosen.
    figures <- c(figures, save(integrationPlotFactorPairs(results, texts), "factor_pairs"))
    for (factor in results$shown) for (view in results$views) {
        figures <- c(figures, save(integrationPlotTopFeatureScatter(results, factor, view, texts),
            sprintf("top_features_%s_%s", factor, view)),
            save(integrationPlotTopFeatureHeatmap(results, factor, view, texts), sprintf("top_features_heatmap_%s_%s", factor, view)))
    }
    for (method in names(maps)) figures <- c(figures, save(integrationPlotFactorMap(results, method, texts), paste0("map_", method)))
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
    # KIT-152: a map chosen but refused for its settings stops the save with its sentence; nothing is saved without it.
    for (map in results$maps) if (!is.null(map$refusal)) integrationGuiRefuse(map$refusal, "integration_map_invalid")
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
