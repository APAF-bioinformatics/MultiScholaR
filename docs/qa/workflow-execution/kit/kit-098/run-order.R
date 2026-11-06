# KIT-098 planning item (KIT-096 run order): does aligning runs by name give `identical` results?
# The GUI's design keeps the counts file's run order (KIT-096); the notebook reorders to design-file order
# (`count_matrix[, design_matrix$Run]`), and KIT-040's corpus shuffles the count columns on purpose. This runs KIT-093's
# kernels (the functions the declarations call, KIT-095) twice on the bound trx_screen corpus: with the design in
# design-file order (the notebook's) and with the same design rows in the counts file's order (the GUI's). Each product of
# the second run is aligned to design-file order by run name only, then compared with `identical()` to the KIT-040 full
# reference product. No tolerance. In-process, read-only; the kernels are the installed build's.
# Usage: R_LIBS=<kit library> Rscript docs/qa/workflow-execution/kit/kit-098/run-order.R <corpus dir> <full reference rds>
suppressMessages(library(MultiScholaR))
args <- commandArgs(TRUE)
dir <- args[[1L]]
reference <- readRDS(args[[2L]])$products
ns <- asNamespace("MultiScholaR")
k <- function(name) get(name, ns)
counts <- vroom::vroom(file.path(dir, "counts.tsv"), show_col_types = FALSE)
design_file <- vroom::vroom(file.path(dir, "design_matrix.tab"), show_col_types = FALSE)
contrast <- k("trxParseContrastFile")(file.path(dir, "contrast_strings.tab"))
runs_counts <- setdiff(names(counts), c("GeneID", "GeneName"))
runs_design <- design_file$Run
design_gui <- design_file[match(runs_counts, design_file$Run), ]
pipeline <- function(design) {
    object <- k("createTranscriptomicsAssayData")(list(counts = counts), design, c("RPMI", "Sera"))
    x <- k("trxAssayInputs")(object)
    f <- k("trxExpressionFilter")(x$counts, x$design$group)
    n <- k("trxTmmLogCpm")(f$filtered_counts, x$design$group)
    g <- k("trxRuvGrid")(f$filtered_counts, n$log_cpm_tmm, x$design)
    s <- k("trxRuvs")(f$filtered_counts, n$log_cpm_tmm, x$design, g$control_genes, g$best$best_k)
    d <- k("trxDaContrast")(s$ruvs_counts, s$pdata, contrast)
    r <- k("trxRuviii")(n$log_cpm_tmm, x$design, g$control_genes, g$best$best_k)
    list(keep_expr = f$keep_expr, gene_vars_positive = f$keep_var, count_matrix_filtered = f$filtered_counts,
        filter_effective_lib_size = f$effective_lib_size, dge_filtered_samples = n$norm_factors_filtered,
        log_cpm_tmm = n$log_cpm_tmm, log_cpm_raw = n$log_cpm_raw, ordered_features_idx = g$prelim_order,
        max_acceptable_k = g$max_acceptable_k, df_summary = g$grid_summary, best_k = g$best$best_k,
        best_control_genes = g$control_genes, w_matrix = s$W, ruvs_counts = s$ruvs_counts, log_cpm_corrected = s$log_cpm_ruvs,
        design_final = d$design_final, qlf_table = d$qlf$table, top_genes = d$top_table,
        ruviii_M = r$M, ruviii_newY = r$newY, ruviii_ctl_mask = r$ctl_mask)
}
# Align a product to design-file order by run name. Run-named columns, rows and vector names are permuted by name. A
# product indexed by run without names (one value or row per run in the design's order: the effective library sizes,
# W, the design matrix, M) is permuted by the design's own Run column, and its automatic row names ("1".."n") are kept
# automatic. Nothing else changes; a product indexed by gene only is returned unchanged.
perm <- match(runs_design, runs_counts)
align <- function(x) {
    n <- length(runs_design)
    by_name <- function(names_) !is.null(names_) && setequal(names_, runs_design) && length(names_) == n
    automatic <- function(names_) is.null(names_) || identical(as.character(names_), as.character(seq_len(n)))
    if (is.matrix(x) || is.data.frame(x)) {
        if (by_name(rownames(x))) x <- x[runs_design, , drop = FALSE] else if (nrow(x) == n && automatic(rownames(x))) {
            keep <- rownames(x)
            x <- x[perm, , drop = FALSE]
            if (is.data.frame(x)) attr(x, "row.names") <- .set_row_names(n) else rownames(x) <- keep
        }
        if (by_name(colnames(x))) x <- x[, runs_design, drop = FALSE]
        return(x)
    }
    if (is.atomic(x) && by_name(names(x))) return(x[runs_design])
    if (is.atomic(x) && is.null(names(x)) && length(x) == n) return(x[perm])
    x
}
gui_order <- pipeline(design_gui)
file_order <- pipeline(design_file)
cat("corpus", dir, "\ncounts-file order != design-file order:", !identical(runs_counts, runs_design), "\n")
cat(sprintf("%-26s %-12s %-14s %-14s %s\n", "product", "run_indexed", "file_order", "gui_aligned", "gui_unaligned"))
rows <- list()
for (name in names(gui_order)) {
    expected <- reference[[name]]
    gui <- gui_order[[name]]
    aligned <- align(gui)
    runs_named <- !identical(aligned, gui) || !identical(align(file_order[[name]]), file_order[[name]])
    row <- data.frame(product = name, runs_named = runs_named, file_order = identical(file_order[[name]], expected),
        gui_aligned = identical(aligned, expected), gui_unaligned = identical(gui, expected))
    if (!row$gui_aligned && is.numeric(as.matrix(aligned)) && identical(dim(as.matrix(aligned)), dim(as.matrix(expected))))
        row$max_abs_diff <- max(abs(as.matrix(aligned) - as.matrix(expected))) else row$max_abs_diff <- NA_real_
    rows[[name]] <- row
    cat(sprintf("%-26s %-12s %-14s %-14s %-14s %s\n", name, runs_named, row$file_order, row$gui_aligned, row$gui_unaligned,
        if (is.na(row$max_abs_diff)) "" else format(row$max_abs_diff)))
}
table <- do.call(rbind, rows)
cat("\nfile order identical to the reference:", sum(table$file_order), "of", nrow(table),
    "\ncounts-file order, aligned by name, identical to the reference:", sum(table$gui_aligned), "of", nrow(table), "\n")
cat("R", R.version.string, "| BLAS", extSoftVersion()[["BLAS"]], "| LAPACK", La_library(), "\n")
