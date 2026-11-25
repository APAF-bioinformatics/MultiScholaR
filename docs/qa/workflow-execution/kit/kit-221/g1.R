# KIT-221 G1 and G3: the digests of one-comparison calls, to be run on the base side (the tests-only commit, whose R/
# is the base's) and on the head side, and compared. Usage: Rscript g1.R <out dir> <tag>
# Writes <out dir>/g1-<tag>.tsv: case, field, sha256 of serialize(value), rows, columns.
# Cases: limma, limpa (fake dpcDE), ruv (factors in the model, three groups relabelled), kernel (the application's call
# shape, a one-row data.frame), alias (G3: the workbooks' form, a one-row slice of a tibble through deAnalysisWrapperFunction).
# The fixtures are those of tests/testthat/test-kit221-da-several-comparisons.R (copied). Synthetic data only.

args <- commandArgs(trailingOnly = TRUE)
out_dir <- args[[1L]]
tag <- args[[2L]]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
suppressPackageStartupMessages(library(MultiScholaR))
ns <- asNamespace("MultiScholaR")

quiet <- function(expr) {
    value <- NULL
    utils::capture.output(value <- suppressWarnings(suppressMessages(expr)))
    value
}
with_seed <- function(seed, code) withr::with_seed(seed, code, .rng_kind = "Mersenne-Twister",
    .rng_normal_kind = "Inversion", .rng_sample_kind = "Rejection")

make_matrix <- function(sizes = c(A = 4L, B = 4L, C = 4L), seed = 200001L, proteins = 60L) {
    groups <- rep(names(sizes), sizes)
    x <- with_seed(seed, {
        level <- stats::rnorm(proteins, 22, 2)
        effect <- as.numeric(seq_len(proteins) %in% 11:16)
        level + outer(effect, as.numeric(groups == "B")) +
            matrix(stats::rnorm(proteins * length(groups), sd = 0.3), proteins)
    })
    dimnames(x) <- list(paste0("P", seq_len(proteins)), paste0("S", seq_along(groups)))
    x
}
make_design <- function(sizes = c(A = 4L, B = 4L, C = 4L)) {
    runs <- paste0("S", seq_len(sum(sizes)))
    data.frame(Run = runs, group = rep(names(sizes), sizes), replicates = runs, row.names = runs, stringsAsFactors = FALSE)
}
make_object <- function(x, design) {
    rownames(design) <- NULL
    ProteinQuantitativeData(protein_quant_table = data.frame(Protein.Ids = rownames(x), x, row.names = NULL,
        check.names = FALSE, stringsAsFactors = FALSE), design_matrix = design, sample_id = "Run", group_id = "group",
        technical_replicate_id = "replicates", protein_id_column = "Protein.Ids",
        protein_id_table = data.frame(Protein.Ids = rownames(x), gene_names = paste0("GENE_", rownames(x)),
            stringsAsFactors = FALSE), args = list())
}
limma_object <- function() make_object(make_matrix(), make_design())
limpa_object <- function() {
    x <- make_matrix()
    elist <- methods::new("EList", list(E = x,
        genes = data.frame(NPrec = rep(2L, nrow(x)), PropObs = rep(1, nrow(x)), row.names = rownames(x)),
        other = list(standard.error = matrix(0.2, nrow(x), ncol(x), dimnames = dimnames(x)),
            n.observations = matrix(2L, nrow(x), ncol(x), dimnames = dimnames(x)))))
    object <- limma_object()
    object@args <- list(limpa_dpc_quant_results = list(quantified_elist = elist, standard_errors = elist$other$standard.error))
    object
}
ruv_state <- function(seed = 20261118L, proteins = 300L, per_group = 6L, missing = 0.25) {
    fixture <- with_seed(seed, {
        groups <- rep(c("control", "treated"), each = per_group)
        n_samples <- length(groups)
        mu <- stats::rnorm(proteins, 22, 2)
        s <- sqrt(4 * 0.25^2 / stats::rchisq(proteins, 4))
        effect <- numeric(proteins)
        where <- sample.int(proteins, round(0.2 * proteins))
        effect[where] <- sample(rep(c(1, -1), length.out = length(where)))
        scale1 <- function(v) { v <- v - mean(v); v / stats::sd(v) }
        w <- scale1(stats::lm.fit(stats::model.matrix(~ factor(groups)), stats::rnorm(n_samples))$residuals)
        tau <- sqrt(3 * 4 * 0.25^2 / stats::qchisq(0.5, 4))
        x <- matrix(mu, proteins, n_samples) + outer(effect, as.numeric(groups == "treated")) +
            outer(stats::rnorm(proteins, 0, tau), w) + matrix(stats::rnorm(proteins * n_samples), proteins, n_samples) * s
        lod <- stats::uniroot(function(l) mean(stats::plogis(l - x)) - missing, c(0, 60))$root
        x[stats::runif(length(x)) < stats::plogis(lod - x)] <- NA_real_
        dimnames(x) <- list(sprintf("PG_%05d", seq_len(proteins)), sprintf("run_%03d", seq_len(n_samples)))
        keep <- rowSums(is.finite(x[, groups == "control"])) >= per_group / 2 &
            rowSums(is.finite(x[, groups == "treated"])) >= per_group / 2
        design <- data.frame(Run = colnames(x), group = rep(c("A", "B", "C"), each = n_samples / 3L),
            replicates = seq_len(n_samples), tech_rep_group = paste0("u", seq_len(n_samples)), stringsAsFactors = FALSE)
        list(x = x[keep, , drop = FALSE], design = design)
    })
    quant <- as.data.frame(c(list(Protein.Group = rownames(fixture$x)), as.data.frame(fixture$x)), check.names = FALSE,
        stringsAsFactors = FALSE)
    normalised <- methods::new("ProteinQuantitativeData", protein_quant_table = quant, protein_id_column = "Protein.Group",
        protein_id_table = data.frame(Protein.Group = rownames(fixture$x), stringsAsFactors = FALSE),
        design_matrix = fixture$design, sample_id = "Run", group_id = "group", technical_replicate_id = "tech_rep_group",
        args = list(globalParameters = list(workflow_type = "LFQ")))
    factors <- ns$computeProtRuvModelFactors(normalised, ruv_k = 1L, ruv_k_rule = "user", ruv_grouping_variable = "group")
    corrected <- suppressMessages(ns$ruvIII_C_Varying(ns$protRuvGroupReplicateObject(normalised),
        ruv_grouping_variable = "group", ruv_number_k = factors$k, ctrl = factors$controls))
    ns$attachProtRuvModelFactors(corrected, normalised, factors)
}

sha <- function(x) digest::digest(serialize(x, NULL), algo = "sha256", serialize = FALSE)
is_plot_column <- function(column) is.list(column) && length(column) && all(vapply(column, inherits, logical(1L), "ggplot"))

# One row per returned data frame, per plot's data, per comparison table of contrasts_results, and the model record.
digests_of <- function(case, result) {
    rows <- list()
    add <- function(field, x) rows[[length(rows) + 1L]] <<- data.frame(case = case, field = field, sha256 = sha(x),
        rows = if (is.data.frame(x)) nrow(x) else NA_integer_, columns = if (is.data.frame(x)) ncol(x) else NA_integer_,
        stringsAsFactors = FALSE)
    for (field in names(result)) {
        x <- result[[field]]
        if (inherits(x, "ggplot")) {
            add(paste0(field, "$data"), x$data)
        } else if (is.data.frame(x)) {
            add(field, x[, !vapply(x, is_plot_column, logical(1L)), drop = FALSE])
        } else if (identical(field, "contrasts_results")) {
            for (name in names(x$results)) add(paste0(field, "$results$", name), x$results[[name]])
        } else if (field %in% c("da_model", "qvalue_warnings")) {
            add(field, x)
        } else if (is.list(x) && !isS4(x) && length(x) && all(vapply(x, is.data.frame, logical(1L)))) {
            for (name in names(x)) add(paste0(field, "$", name), x[[name]])
        } else if (is.list(x) && !isS4(x) && is.data.frame(x$table)) {
            add(paste0(field, "$table"), x$table)
            if (inherits(x$plot, "ggplot")) add(paste0(field, "$plot$data"), x$plot$data)
        }
    }
    do.call(rbind, rows)
}

plots <- list(plotRle = function(...) ggplot2::ggplot(), plotPca = function(...) ggplot2::ggplot(),
    plotNumOfValuesNoLog = function(...) ggplot2::ggplot())
one <- data.frame(contrasts = "B_vs_A=groupB-groupA", stringsAsFactors = FALSE)
da <- function(object, contrasts_tbl, row_id = "Protein.Ids", trend = TRUE, robust = TRUE, ...) {
    differentialAbundanceAnalysis(object, contrasts_tbl = contrasts_tbl, formula_string = "~ 0 + group",
        da_q_val_thresh = 0.05, treat_lfc_cutoff = 0, eBayes_trend = trend, eBayes_robust = robust, args_row_id = row_id,
        qvalue_column = "fdr_qvalue", raw_pvalue_column = "raw_pvalue", ...)
}
with_plots <- function(code) testthat::with_mocked_bindings(code, plotRle = plots$plotRle, plotPca = plots$plotPca,
    plotNumOfValuesNoLog = plots$plotNumOfValuesNoLog, .package = "MultiScholaR")

tables <- list()
tables$limma <- digests_of("limma", with_plots(quiet(da(limma_object(), one))))
tables$limpa <- digests_of("limpa", with_plots(testthat::with_mocked_bindings(
    quiet(da(limpa_object(), one, trend = FALSE, robust = FALSE)),
    dpcDE = function(y, design, plot = TRUE, ...) limma::lmFit(y$E, design), .package = "limpa")))
tables$ruv <- digests_of("ruv", with_plots(quiet(da(ruv_state(), one, row_id = "Protein.Group"))))
# The application's kernel call shape (R/utils_workflow_dia_da.R:535-543): a one-row data.frame, stringsAsFactors = FALSE,
# the stored arguments recorded on the object first.
kernel_object <- limma_object()
kernel_arguments <- list(contrasts_tbl = one, formula_string = "~ 0 + group", da_q_val_thresh = 0.05, treat_lfc_cutoff = 0,
    eBayes_trend = TRUE, eBayes_robust = TRUE, args_group_pattern = "(\\d+)", args_row_id = "Protein.Ids",
    qvalue_column = "fdr_qvalue", raw_pvalue_column = "raw_pvalue")
for (section in c("differentialAbundanceAnalysis", "differentialAbundanceAnalysisHelper")) {
    kernel_object@args[[section]] <- kernel_arguments
}
tables$kernel <- digests_of("kernel", with_plots(quiet(differentialAbundanceAnalysis(kernel_object, contrasts_tbl = one,
    formula_string = "~ 0 + group", da_q_val_thresh = 0.05, treat_lfc_cutoff = 0, qvalue_column = "fdr_qvalue",
    raw_pvalue_column = "raw_pvalue"))))
# G3: the workbooks' form, the alias on a one-row slice of a tibble (DIA_workflow_experienced.rmd:1107-1109 at 176f9d3ff).
tbl <- tibble::tibble(contrasts = c("B_vs_A=groupB-groupA", "C_vs_A=groupC-groupA"))
tables$alias <- digests_of("alias", with_plots(quiet(deAnalysisWrapperFunction(limma_object(),
    contrasts_tbl = dplyr::slice(tbl, 1L), formula_string = "~ 0 + group", de_q_val_thresh = 0.05, treat_lfc_cutoff = 0,
    eBayes_trend = TRUE, eBayes_robust = TRUE, args_row_id = "Protein.Ids"))))

result <- do.call(rbind, tables)
utils::write.table(result, file.path(out_dir, paste0("g1-", tag, ".tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("g1 %s: %d digests (%s)\n", tag, nrow(result), paste(names(tables), collapse = ", ")))
