# KIT-076: the old-vs-new DA comparison under the ticket's frozen rule. `~ 0 + group`: identical. Any other design:
# absolute difference <= 1e-12 on numeric columns with identical calls, row sets and order. Anything beyond it is a stop.
# Prints aggregates only. Usage: Rscript docs/qa/workflow-execution/kit/kit-076/da-compare.R <v1 dir> <v2 dir>
args <- commandArgs(TRUE); old_dir <- args[[1L]]; new_dir <- args[[2L]]
tolerance <- 1e-12
violations <- character()
violate <- function(where, what) violations <<- c(violations, paste(where, what))
# Largest absolute difference; a differing NA pattern is infinite.
difference <- function(a, b) {
    a <- as.numeric(a); b <- as.numeric(b)
    if (length(a) != length(b) || !identical(is.na(a), is.na(b)) || !identical(is.nan(a), is.nan(b))) return(Inf)
    keep <- !is.na(a)
    if (!any(keep)) return(0)
    finite <- is.finite(a[keep]) & is.finite(b[keep])
    if (!identical(a[keep][!finite], b[keep][!finite])) return(Inf)
    if (!any(finite)) return(0)
    max(abs(a[keep][finite] - b[keep][finite]))
}
# One table against the rule; returns its largest numeric difference.
compare_table <- function(old, new, rule, where) {
    old <- as.data.frame(old); new <- as.data.frame(new)
    if (!identical(names(old), names(new))) { violate(where, "columns differ"); return(Inf) }
    if (nrow(old) != nrow(new)) { violate(where, "row count differs"); return(Inf) }
    largest <- 0
    for (column in names(old)) {
        if (is.numeric(old[[column]]) && is.numeric(new[[column]])) {
            gap <- difference(old[[column]], new[[column]])
            largest <- max(largest, gap)
            if (identical(rule, "identical") && !identical(old[[column]], new[[column]])) violate(where, paste0("numeric column not identical: ", column))
            if (identical(rule, "tolerance") && gap > tolerance) violate(where, sprintf("numeric column beyond 1e-12: %s (%g)", column, gap))
        } else if (!identical(old[[column]], new[[column]])) {
            violate(where, paste0("column differs: ", column))
        }
    }
    largest
}
fit_matrices <- c("coefficients", "stdev.unscaled", "t", "p.value", "lods")
fit_vectors <- c("sigma", "df.residual", "df.prior", "s2.prior", "s2.post", "df.total", "Amean", "var.prior", "treat.lfc")
compare_fit <- function(old_fit, new_fit, contrast, rule, where) {
    if (!identical(rownames(old_fit$coefficients), rownames(new_fit$coefficients))) violate(where, "fit row order differs")
    if (!identical(colnames(new_fit$coefficients), contrast)) violate(where, "saved fit is not one contrast")
    largest <- 0
    for (field in c(fit_matrices, fit_vectors)) {
        a <- old_fit[[field]]; b <- new_fit[[field]]
        if (is.null(a) != is.null(b)) { violate(where, paste0("fit field present in one build only: ", field)); next }
        if (is.null(a)) next
        if (is.matrix(a)) a <- a[, contrast]
        # The per-coefficient prior is a named vector over the fit's coefficients.
        if (identical(field, "var.prior")) a <- a[[contrast]]
        gap <- difference(a, b)
        largest <- max(largest, gap)
        if (identical(rule, "identical") && gap != 0) violate(where, sprintf("fit field differs: %s (%g)", field, gap))
        if (identical(rule, "tolerance") && gap > tolerance) violate(where, sprintf("fit field beyond 1e-12: %s (%g)", field, gap))
    }
    largest
}
workbook_sheets <- function(path) {
    sheets <- readxl::excel_sheets(path)
    stats::setNames(lapply(sheets, function(sheet) as.data.frame(readxl::read_excel(path, sheet = sheet, guess_max = 1e5))), sheets)
}
compare_files <- function(old_root, new_root, rule, where) {
    old_files <- sort(list.files(old_root, recursive = TRUE)); new_files <- sort(list.files(new_root, recursive = TRUE))
    if (!identical(old_files, new_files)) violate(where, sprintf("file sets differ: only old [%s], only new [%s]",
        paste(setdiff(old_files, new_files), collapse = ", "), paste(setdiff(new_files, old_files), collapse = ", ")))
    largest <- 0
    counts <- c(tsv = 0L, xlsx = 0L, figure_same = 0L, figure_different = 0L)
    for (file in intersect(old_files, new_files)) {
        old_path <- file.path(old_root, file); new_path <- file.path(new_root, file)
        if (grepl("\\.tsv$", file)) {
            counts[["tsv"]] <- counts[["tsv"]] + 1L
            same <- identical(digest::digest(file = old_path, algo = "sha256"), digest::digest(file = new_path, algo = "sha256"))
            if (identical(rule, "identical") && !same) violate(file.path(where, file), "TSV bytes differ")
            if (identical(rule, "tolerance") && !same) largest <- max(largest, compare_table(
                utils::read.delim(old_path, check.names = FALSE), utils::read.delim(new_path, check.names = FALSE), "tolerance", file.path(where, file)))
        } else if (grepl("\\.xlsx$", file)) {
            counts[["xlsx"]] <- counts[["xlsx"]] + 1L
            old_sheets <- workbook_sheets(old_path); new_sheets <- workbook_sheets(new_path)
            if (!identical(names(old_sheets), names(new_sheets))) { violate(file.path(where, file), "sheets differ"); next }
            for (sheet in names(old_sheets)) largest <- max(largest, compare_table(old_sheets[[sheet]], new_sheets[[sheet]], rule, file.path(where, file, sheet)))
        } else {
            # Figures are reported, not judged: the ticket's rule is on tables, and PDFs embed a timestamp.
            same <- identical(digest::digest(file = old_path, algo = "sha256"), digest::digest(file = new_path, algo = "sha256"))
            slot <- if (same) "figure_same" else "figure_different"
            counts[[slot]] <- counts[[slot]] + 1L
        }
    }
    list(largest = largest, counts = counts)
}
stopifnot(requireNamespace("readxl", quietly = TRUE))
for (omic in c("lipid", "metab")) for (scenario in c("screen", "multi", "multi_treat", "nonorth", "e2e")) {
    where <- paste0(omic, "/", scenario)
    old <- readRDS(file.path(old_dir, omic, scenario, "capture.rds")); new <- readRDS(file.path(new_dir, omic, scenario, "capture.rds"))
    stopifnot(identical(old$rule, new$rule), identical(old$contrasts, new$contrasts), grepl("v1$", old$schema), grepl("v2$", new$schema))
    rule <- old$rule
    if (!identical(old$direct, new$direct)) violate(where, "the public kernel's list is not identical")
    long_name <- if (identical(omic, "lipid")) "da_lipids_long" else "da_metabolites_long"
    largest <- compare_table(old$results[[long_name]], new$results[[long_name]], rule, paste0(where, "/long"))
    if (!identical(old$results$significant_counts, new$results$significant_counts)) violate(where, "significant counts differ")
    if (!identical(old$results$qvalue_warnings, new$results$qvalue_warnings)) violate(where, "q-value warnings differ")
    if (!identical(names(old$results$per_assay_results), names(new$results$per_assay_results))) violate(where, "assays differ")
    if (length(new$results$failed_contrasts)) violate(where, "a contrast failed")
    fit_largest <- 0; tables <- 0L
    for (assay in names(old$results$per_assay_results)) {
        for (contrast in old$contrasts$contrasts) {
            inner <- paste(where, assay, contrast, sep = "/")
            a <- old$results$per_assay_results[[assay]]; b <- new$results$per_assay_results[[assay]]
            largest <- max(largest, compare_table(a[a$comparison == contrast, ], b[b$comparison == contrast, ], rule, inner))
            tables <- tables + 1L
            largest <- max(largest, compare_table(old$results$contrasts_results[[assay]]$results[[contrast]],
                new$results$contrasts_results[[assay]]$results[[contrast]], rule, paste0(inner, "/kernel")))
            if (!identical(rownames(old$results$contrasts_results[[assay]]$results[[contrast]]),
                rownames(new$results$contrasts_results[[assay]]$results[[contrast]]))) violate(inner, "row order differs")
            fit_largest <- max(fit_largest, compare_fit(old$results$contrasts_results[[assay]]$fit.eb, new$fits[[contrast]][[assay]], contrast, rule, paste0(inner, "/fit")))
        }
    }
    files <- compare_files(file.path(old_dir, omic, scenario), file.path(new_dir, omic, scenario), rule, where)
    cat(sprintf("%-18s %-22s rule %-9s contrasts %d, contrast-assay tables %3d: largest difference tables %.3g, fits %.3g, exports %.3g; files tsv %d xlsx %d, figures same %d different %d (reported)\n",
        where, old$formula, rule, nrow(old$contrasts), tables, largest, fit_largest, files$largest, files$counts[["tsv"]], files$counts[["xlsx"]], files$counts[["figure_same"]], files$counts[["figure_different"]]))
}
cat(sprintf("\nviolations: %d\n", length(violations)))
for (v in utils::head(violations, 60L)) cat(" ", v, "\n")
cat(if (length(violations)) "RESULT stop: beyond the frozen rule\n" else "RESULT within the frozen rule on every capture\n")
