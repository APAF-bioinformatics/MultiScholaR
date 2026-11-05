# KIT-040: render the transcriptomics notebook (dev/distributable/plans/transcriptomics_workflow.rmd) unchanged on a
# synthetic count corpus and freeze its products as the reference later tickets match.
# - Preflight, before any render: refuses NA, negative or non-whole counts (kit_trx_invalid_counts; integer-valued doubles
#   are accepted), a missing design or contrast file (kit_trx_missing_input, so the notebook's positional design fallback
#   never runs) and any package the notebook's p_load names that is absent (kit_trx_missing_package, so the notebook never
#   installs at render time).
# - The prepared directory (outside the repository) holds a .here sentinel, the corpus under the notebook's own file names
#   and a byte copy of the notebook whose sha256 is checked.
# - The render runs in a fresh `Rscript --vanilla` (R_LIBS carried from this process), set.seed(20261106) (only the unseeded
#   geom_text_repel draws), rmarkdown::render(envir = e) with knit_root_dir, output_dir and intermediates_dir all set to the
#   prepared directory. A knitr option hook writes chunk-markers.tsv and, with --fail-before <chunk>, raises
#   kit_trx_forced_failure at the start of that chunk. The .rmd is never edited.
# - Products are read from `e` (never from the notebook's TSVs), each bound by the sha256 of serialize(x, NULL, version = 3).
# Writes to <out dir>: chunk-markers.tsv, render.log, products.rds (full products, labels, binding), fingerprints.rds (the
# fingerprint store), output-hashes.tsv, product-map.tsv, status.json. Exit status 0 only when the render completed and
# every product is present.
# Usage: R_LIBS=<kit library> Rscript --vanilla tools/profiling/trx_notebook_reference.R <corpus dir> <work dir> <out dir> [--fail-before <chunk>]

TRX_NOTEBOOK_SHA256 <- "44d5f643946212c18cee10b5ebf77b1e7ebbe8d22ca37b40bafd08b5a0da652a"
TRX_RENDER_SEED <- 20261106L
TRX_CEILING_SECONDS <- 900L

trxAbort <- function(message, class) {
    stop(structure(class = c(class, "kit_trx_error", "error", "condition"), list(message = message, call = NULL)))
}

trxRepoRoot <- function() {
    file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
    script <- if (length(file)) file[[1L]] else if (!is.null(sys.frame(1)$ofile)) sys.frame(1)$ofile else "tools/profiling/trx_notebook_reference.R"
    normalizePath(file.path(dirname(normalizePath(script, mustWork = TRUE)), "..", ".."), mustWork = TRUE)
}

trxNotebookSource <- function(repo = trxRepoRoot()) file.path(repo, "dev/distributable/plans/transcriptomics_workflow.rmd")

# The notebook's own p_load list, plus pacman (loaded with library()), EDASeq (newSeqExpressionSet, attached by RUVSeq)
# and the render packages.
trxNotebookPackages <- function(notebook) {
    text <- paste(readLines(notebook, warn = FALSE), collapse = "\n")
    calls <- regmatches(text, gregexpr("p_load\\(([^)]*)\\)", text))[[1L]]
    named <- trimws(unlist(strsplit(sub("^p_load\\((.*)\\)$", "\\1", calls), ",")))
    unique(c("pacman", named[nzchar(named)], "EDASeq", "rmarkdown", "knitr"))
}

trxPreflight <- function(corpus, packages = trxNotebookPackages(trxNotebookSource()), lib = .libPaths()) {
    inputs <- file.path(corpus, c("counts.tsv", "design_matrix.tab", "contrast_strings.tab"))
    absent <- inputs[!file.exists(inputs)]
    if (length(absent)) trxAbort(sprintf("missing input %s in %s", paste(basename(absent), collapse = ", "), corpus), "kit_trx_missing_input")
    counts <- utils::read.delim(inputs[[1L]], colClasses = "character", na.strings = character(), check.names = FALSE, quote = "")
    cells <- as.matrix(counts[, -(1:2), drop = FALSE])
    x <- suppressWarnings(as.numeric(cells))
    bad <- which(!is.finite(x) | x < 0 | x != round(x))
    if (length(bad)) {
        at <- arrayInd(bad[[1L]], dim(cells))
        trxAbort(sprintf("count '%s' for gene %s in run %s is not a whole non-negative number (%d invalid cells)",
            cells[at], counts[[1L]][at[[1L]]], colnames(cells)[at[[2L]]], length(bad)), "kit_trx_invalid_counts")
    }
    missing <- packages[!vapply(packages, function(p) nzchar(system.file(package = p, lib.loc = lib)), logical(1L))]
    if (length(missing)) trxAbort(sprintf("package(s) the notebook loads are not installed: %s", paste(missing, collapse = ", ")), "kit_trx_missing_package")
    invisible(list(genes = nrow(cells), runs = ncol(cells), packages = packages))
}

trxSha256 <- function(path) digest::digest(file = path, algo = "sha256")
trxProductSha256 <- function(x) digest::digest(serialize(x, NULL, version = 3L), algo = "sha256", serialize = FALSE)

trxPrepare <- function(corpus, work, notebook, repo) {
    work <- normalizePath(work, mustWork = FALSE)
    if (startsWith(paste0(work, "/"), paste0(repo, "/"))) trxAbort(sprintf("the work directory %s is inside the repository", work), "kit_trx_work_dir")
    if (dir.exists(work) && length(list.files(work, all.files = TRUE, no.. = TRUE))) trxAbort(sprintf("the work directory %s is not empty", work), "kit_trx_work_dir")
    dir.create(file.path(work, "data", "transcriptomics"), recursive = TRUE)
    dir.create(file.path(work, "scripts", "transcriptomics"), recursive = TRUE)
    file.create(file.path(work, ".here"))
    copy <- function(from, to) stopifnot(file.copy(from, to, copy.date = FALSE), identical(trxSha256(from), trxSha256(to)))
    copy(file.path(corpus, "counts.tsv"), file.path(work, "data", "transcriptomics", "GSE152840_gene_counts.txt"))
    copy(file.path(corpus, "design_matrix.tab"), file.path(work, "scripts", "transcriptomics", "design_matrix.tab"))
    copy(file.path(corpus, "contrast_strings.tab"), file.path(work, "scripts", "transcriptomics", "contrast_strings.tab"))
    rmd <- file.path(work, "transcriptomics_workflow.rmd")
    stopifnot(file.copy(notebook, rmd, copy.date = FALSE))
    if (!identical(trxSha256(rmd), TRX_NOTEBOOK_SHA256)) trxAbort(sprintf("the notebook copy's sha256 is %s, not %s", trxSha256(rmd), TRX_NOTEBOOK_SHA256), "kit_trx_notebook_mismatch")
    work
}

# --- Products (read from the render environment) --------------------------------------------------------------------

# Each product: the notebook objects it needs, the expression evaluated in the render environment, the rmd lines and a label.
trxProduct <- function(needs, expr, source, label = "", derived = FALSE) list(needs = needs, expr = expr, source = source, label = label, derived = derived)

trxPcaVariance <- function(p) {
    scores <- as.matrix(p$data[grepl("^PC[0-9]+$", names(p$data))])
    v <- apply(scores, 2L, stats::var)
    v / sum(v)
}

trxRleWhiskers <- function(p) {
    d <- p$data
    q1 <- tapply(d$rle_value, d$Run, stats::quantile, 0.25, na.rm = TRUE)
    q3 <- tapply(d$rle_value, d$Run, stats::quantile, 0.75, na.rm = TRUE)
    data.frame(Run = names(q1), q1 = unname(q1), q3 = unname(q3), iqr = unname(q3 - q1),
        upper_whisker = unname(q3 + 1.5 * (q3 - q1)), lower_whisker = unname(q1 - 1.5 * (q3 - q1)))
}

trxRleBox <- function(p) {
    d <- ggplot2::layer_data(p, 1L)
    d[, intersect(c("x", "group", "PANEL", "fill", "ymin", "lower", "middle", "upper", "ymax", "notchlower", "notchupper", "relvarwidth", "width"), names(d))]
}

trxProductDefinitions <- function() {
    defs <- list(
        keep_expr = trxProduct("keep_expr", quote(keep_expr), "rmd:170", "edgeR::filterByExpr mask on the unfiltered, design-ordered counts"),
        gene_vars = trxProduct("gene_vars", quote(gene_vars), "rmd:174", "per-gene count variance after the expression filter"),
        gene_vars_positive = trxProduct("gene_vars", quote(gene_vars > 0), "rmd:175", "the zero-variance filter's mask"),
        count_matrix_filtered = trxProduct("count_matrix_filtered", quote(count_matrix_filtered), "rmd:171-175", "the filtered counts (design order)"),
        dge_raw_samples = trxProduct("dge_raw", quote(dge_raw$samples), "rmd:166-167", "pre-filter TMM factors (the ones filterByExpr uses)"),
        dge_filtered_samples = trxProduct("dge_filtered", quote(dge_filtered$samples), "rmd:178-179", "post-filter TMM factors"),
        log_cpm_tmm = trxProduct("log_cpm_tmm", quote(log_cpm_tmm), "rmd:180", "logCPM after TMM on the filtered counts"),
        filter_effective_lib_size = trxProduct("dge_raw", quote(dge_raw$samples$lib.size * dge_raw$samples$norm.factors), "rmd:170",
            "derived: the effective library sizes filterByExpr.DGEList passes (lib.size * norm.factors of dge_raw)", derived = TRUE),
        filter_min_sample_size = trxProduct("design_matrix", quote({
            # base:: throughout: the notebook attaches conflicted, which refuses an ambiguous unqualified name here.
            f <- base::formals(edgeR:::filterByExpr.default); n <- base::tabulate(base::as.factor(design_matrix$group)); m <- base::min(n[n > 0L])
            if (m > f$large.n) m <- f$large.n + (m - f$large.n) * f$min.prop
            m
        }), "rmd:170", "derived: filterByExpr.default's MinSampleSize by edgeR's own rule (large.n, min.prop defaults)", derived = TRUE),
        filter_median_lib_size = trxProduct("dge_raw", quote(stats::median(dge_raw$samples$lib.size * dge_raw$samples$norm.factors)), "rmd:170",
            "derived: filterByExpr.default's MedianLibSize of the effective library sizes", derived = TRUE),
        tt_prelim = trxProduct("tt_prelim", quote(tt_prelim), "rmd:236", "preliminary (non-robust) QL test table, sort.by none"),
        ordered_features_idx = trxProduct("ordered_features_idx", quote(ordered_features_idx), "rmd:237", "control ranking: order(PValue, decreasing = TRUE)"),
        max_acceptable_k = trxProduct("max_acceptable_k", quote(max_acceptable_k), "rmd:242", "the k cap from the sample count"),
        grid_members = trxProduct(c("percentage_range", "grid_results"), quote({
            ok <- stats::setNames(grid_results, vapply(grid_results, function(r) as.character(r$percentage), ""))
            do.call(rbind, lapply(percentage_range, function(pct) {
                r <- ok[[as.character(pct)]]
                if (is.null(r)) return(data.frame(percentage = pct, status = "failed", n_controls = NA_real_, best_k = NA_integer_,
                    separation_score = NA_real_, k_penalty = NA_real_, composite_score = NA_real_))
                data.frame(percentage = r$percentage, status = "complete", n_controls = r$n_controls, best_k = r$best_k,
                    separation_score = r$separation_score, k_penalty = r$k_penalty, composite_score = r$composite_score)
            }))
        }), "rmd:247-281", "one row per percentage, unrounded; a member dropped by compact() shows as failed"),
        grid_cancor_data = trxProduct("grid_results", quote(stats::setNames(lapply(grid_results, function(r) r$cancor_plot$data),
            vapply(grid_results, function(r) paste0("pct_", r$percentage), ""))), "rmd:253", "ruv_cancorplot data per complete member"),
        grid_control_genes = trxProduct("grid_results", quote(stats::setNames(lapply(grid_results, function(r) r$control_genes),
            vapply(grid_results, function(r) paste0("pct_", r$percentage), ""))), "rmd:279", "inferred control genes per complete member"),
        df_summary = trxProduct("df_summary", quote(df_summary), "rmd:283-292", "the notebook's 4-dp rounded grid summary (display)"),
        best_res = trxProduct("best_res", quote(best_res[c("percentage", "n_controls", "best_k", "separation_score", "k_penalty", "composite_score", "ctl_idx", "control_genes")]),
            "rmd:295-296", "the selected member's fields (cancor_plot excluded; its data is grid_cancor_data)"),
        best_k = trxProduct("best_k", quote(best_k), "rmd:297", "the selected k"),
        best_control_genes = trxProduct("best_control_genes", quote(best_control_genes), "rmd:298",
            "inferred: the least significant genes of the preliminary fit, not certified negative controls"),
        w_matrix = trxProduct("w_matrix", quote(w_matrix), "rmd:418", "RUVs W with Run names"),
        ruvs_counts = trxProduct("set_ruvs", quote(BiocGenerics::counts(set_ruvs)), "rmd:413", "counts(set_ruvs)"),
        ruvs_counts_identical = trxProduct(c("set_ruvs", "count_matrix_filtered"), quote(identical(BiocGenerics::counts(set_ruvs), count_matrix_filtered)),
            "rmd:407-413", "identical(counts(set_ruvs), count_matrix_filtered)"),
        log_cpm_corrected = trxProduct("log_cpm_corrected", quote(log_cpm_corrected), "rmd:422-426", "visualisation only: removeBatchEffect on log_cpm_tmm with W"),
        log_cpm_raw = trxProduct("log_cpm_raw", quote(log_cpm_raw), "rmd:490", "filtered counts, no TMM (the notebook's \"raw\" diagnostic)"),
        design_final = trxProduct("design_final", quote(design_final), "rmd:532-533", "the final model matrix (columns after make.names)"),
        design_formula_str = trxProduct("design_formula_str", quote(design_formula_str), "rmd:529", "the final model formula string"),
        contrast_entry = trxProduct("contrast_entry", quote(contrast_entry), "rmd:542", "the parsed contrast line"),
        contrasts_matrix = trxProduct("contrasts_matrix", quote(contrasts_matrix), "rmd:547-548", "makeContrasts result, named"),
        fit_final = trxProduct("fit_final", quote(fit_final), "rmd:537", "the robust glmQLFit (DGEGLM)"),
        qlf = trxProduct("lrt", quote(lrt), "rmd:554", "the glmQLFTest result (DGELRT; the notebook's lrt)"),
        qlf_table = trxProduct("lrt", quote(lrt$table), "rmd:554", "the glmQLFTest table in gene order"),
        top_genes = trxProduct("top_genes", quote(top_genes), "rmd:555", "topTags(n = Inf, sort.by = PValue), full, in its order"),
        feature_meta = trxProduct("feature_meta", quote(feature_meta), "rmd:144-148", "gene annotation (GeneName falls back to GeneID)"),
        top_genes_annotated = trxProduct("top_genes_annotated", quote(top_genes_annotated), "rmd:558-562", "the annotated DE table"),
        top_up = trxProduct("top_up", quote(top_up), "rmd:569-573", "top 10 up at FDR < 0.05"),
        top_down = trxProduct("top_down", quote(top_down), "rmd:575-579", "top 10 down at FDR < 0.05"),
        top_de_combined = trxProduct("top_de_combined", quote(top_de_combined), "rmd:581-582", "top up and down with comparison"),
        ruviii_M = trxProduct("M", quote(M), "rmd:704", "RUVIII M from group (notebook default; replicate assumption unconfirmed)"),
        ruviii_ctl_mask = trxProduct("ctl_mask", quote(ctl_mask), "rmd:707-708", "RUVIII control mask (inferred controls)"),
        ruviii_newY = trxProduct("res_ruviii", quote(res_ruviii$newY), "rmd:711-712", "RUVIII, M from group (notebook default; replicate assumption unconfirmed)"),
        glimma_html_exists = trxProduct("glimma_html_path", quote(file.exists(glimma_html_path)), "rmd:678-693", "the Glimma outcome: whether the HTML file exists")
    )
    stages <- list(raw = c("p_pca_raw", "p_rle_raw", "rmd:492,496"), tmm = c("p_pca_tmm", "p_rle_tmm", "rmd:493,497"),
        ruvs = c("p_pca_ruvs", "p_rle_ruvs", "rmd:494,498"), post_ruviii = c("p_pca_post_ruviii", "p_rle_post_ruviii", "rmd:731-732"))
    for (stage in names(stages)) {
        pca <- as.name(stages[[stage]][[1L]]); rle <- as.name(stages[[stage]][[2L]]); src <- stages[[stage]][[3L]]
        defs[[paste0("pca_", stage, "_scores")]] <- trxProduct(stages[[stage]][[1L]], bquote(.(pca)$data), src, "the PCA plot's data (scores of every PC, with the design)")
        defs[[paste0("pca_", stage, "_variance")]] <- trxProduct(stages[[stage]][[1L]], bquote(trxPcaVariance(.(pca))), src,
            "every PC's variance fraction from the PCA plot's scores, unrounded (var(PC) / sum(var))", derived = TRUE)
        defs[[paste0("pca_", stage, "_axis_labels")]] <- trxProduct(stages[[stage]][[1L]], bquote(unlist(ggplot2::get_labs(.(pca))[c("x", "y")])), src,
            "the notebook's axis labels (5-dp rounded fractions from summary.prcomp)")
        defs[[paste0("rle_", stage, "_data")]] <- trxProduct(stages[[stage]][[2L]], bquote(.(rle)$data), src, "the RLE plot's data: the gene-median-centred matrix, long, with the design")
        defs[[paste0("rle_", stage, "_box")]] <- trxProduct(stages[[stage]][[2L]], bquote(trxRleBox(.(rle))), src, "the RLE box statistics from ggplot2::layer_data()")
        defs[[paste0("rle_", stage, "_ylim")]] <- trxProduct(stages[[stage]][[2L]], bquote(.(rle)$coordinates$limits$y), src, "the RLE plot's y limits (coord_cartesian)")
        defs[[paste0("rle_", stage, "_whiskers")]] <- trxProduct(stages[[stage]][[2L]], bquote(trxRleWhiskers(.(rle))), src,
            "derived: the per-run whisker limits recomputed from the RLE plot's data by the notebook's rule", derived = TRUE)
    }
    defs
}

trxReadProducts <- function(e) {
    defs <- trxProductDefinitions()
    helpers <- new.env(parent = e)
    for (f in c("trxPcaVariance", "trxRleWhiskers", "trxRleBox")) assign(f, get(f), envir = helpers)
    products <- list(); absent <- character(); errors <- list()
    for (name in names(defs)) {
        def <- defs[[name]]
        if (!all(vapply(def$needs, exists, logical(1L), envir = e, inherits = FALSE))) { absent <- c(absent, name); next }
        value <- tryCatch(eval(def$expr, envir = helpers), error = function(err) err)
        if (inherits(value, "error")) { absent <- c(absent, name); errors[[name]] <- conditionMessage(value); next }
        products[name] <- list(value)
    }
    labels <- lapply(defs, function(d) d[c("source", "label", "derived")])
    list(products = products, absent = absent, errors = errors, labels = labels)
}

# --- The render child -----------------------------------------------------------------------------------------------

trxStamp <- function() format(Sys.time(), "%Y-%m-%dT%H:%M:%OS3%z")

trxRenderChild <- function(work, out, fail_before = NA_character_) {
    setwd(work)
    markers <- file.path(out, "chunk-markers.tsv")
    state <- new.env()
    state$rows <- data.frame(label = character(), start = character(), end = character(), stringsAsFactors = FALSE)
    flush <- function() utils::write.table(state$rows, markers, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
    flush()
    knitr::opts_chunk$set(kit_trx_marker = TRUE)
    knitr::opts_hooks$set(kit_trx_marker = function(options) {
        if (identical(options$label, fail_before)) trxAbort(sprintf("forced failure before chunk %s", fail_before), "kit_trx_forced_failure")
        state$rows[nrow(state$rows) + 1L, ] <- list(options$label, trxStamp(), NA_character_)
        flush()
        options
    })
    knitr::knit_hooks$set(kit_trx_marker = function(before, options, envir) {
        if (!before) {
            state$rows$end[state$rows$label == options$label] <- trxStamp()
            flush()
        }
        NULL
    })
    e <- new.env(parent = globalenv())
    first_error <- NULL
    set.seed(TRX_RENDER_SEED)
    started <- trxStamp()
    result <- tryCatch(withCallingHandlers(
        rmarkdown::render(file.path(work, "transcriptomics_workflow.rmd"), envir = e, knit_root_dir = work, output_dir = work,
            intermediates_dir = work, quiet = TRUE),
        error = function(err) if (is.null(first_error)) first_error <<- err), error = function(err) err)
    finished <- trxStamp()
    read <- trxReadProducts(e)
    rendered <- !inherits(result, "error")
    complete <- rendered && !length(read$absent)
    si <- utils::sessionInfo()
    saveRDS(list(schema = "kit.trx_reference.v1", complete = complete, rendered = rendered,
        error = if (!is.null(first_error)) list(class = class(first_error), message = conditionMessage(first_error)) else NULL,
        products = read$products, absent = read$absent, product_errors = read$errors, labels = read$labels,
        started = started, finished = finished, fail_before = fail_before,
        session = list(session_info = utils::capture.output(print(si)), blas = si$BLAS, lapack = si$LAPACK, la_version = La_version(),
            r = R.version.string, lib_paths = .libPaths(),
            packages = vapply(si$otherPkgs, function(p) p$Version, ""), loaded = vapply(si$loadedOnly, function(p) p$Version, ""))),
        file.path(out, "child.rds"), version = 3L)
    quit(save = "no", status = if (complete) 0L else 1L)
}

# --- The parent: preflight, prepare, render, bind -------------------------------------------------------------------

trxOutputHashes <- function(work) {
    files <- list.files(work, recursive = TRUE, all.files = TRUE)
    inputs <- c(".here", "data/transcriptomics/GSE152840_gene_counts.txt", "scripts/transcriptomics/design_matrix.tab",
        "scripts/transcriptomics/contrast_strings.tab", "transcriptomics_workflow.rmd")
    files <- sort(setdiff(files, inputs))
    data.frame(path = files, bytes = file.size(file.path(work, files)),
        sha256 = vapply(file.path(work, files), trxSha256, "", USE.NAMES = FALSE),
        class = ifelse(grepl("\\.(csv|tsv)$", files), "stable", "recorded_only"), stringsAsFactors = FALSE)
}

trxDescribe <- function(x) {
    d <- dim(x)
    data.frame(class = paste(class(x), collapse = "/"), dims = if (is.null(d)) as.character(length(x)) else paste(d, collapse = "x"),
        sha256 = trxProductSha256(x), stringsAsFactors = FALSE)
}

# KIT-093's kernel products (its kernel table) mapped to the notebook objects frozen here (consistency pass, 2026-11-06).
trxProductMap <- function() {
    row <- function(kernel, product, reference, object, lines, note = "") data.frame(kit093_kernel = kernel, kit093_product = product,
        reference_product = reference, notebook_object = object, rmd_lines = lines, note = note, stringsAsFactors = FALSE)
    rbind(
        row("trxExpressionFilter", "effective_lib_size", "filter_effective_lib_size", "dge_raw$samples$lib.size * norm.factors", "166-170", "derived by edgeR's rule"),
        row("trxExpressionFilter", "min_sample_size", "filter_min_sample_size", "(filterByExpr.default internal)", "170", "derived by edgeR's rule"),
        row("trxExpressionFilter", "median_lib_size", "filter_median_lib_size", "(filterByExpr.default internal)", "170", "derived by edgeR's rule"),
        row("trxExpressionFilter", "keep_expr", "keep_expr", "keep_expr", "170"),
        row("trxExpressionFilter", "keep_var", "gene_vars_positive", "gene_vars > 0", "174-175"),
        row("trxExpressionFilter", "filtered_counts", "count_matrix_filtered", "count_matrix_filtered", "171-175"),
        row("trxTmmLogCpm", "norm_factors_filtered", "dge_filtered_samples", "dge_filtered$samples", "178-179", "the norm.factors column"),
        row("trxTmmLogCpm", "log_cpm_tmm", "log_cpm_tmm", "log_cpm_tmm", "180"),
        row("trxTmmLogCpm", "log_cpm_raw", "log_cpm_raw", "log_cpm_raw", "490", "filtered counts, no TMM"),
        row("trxRuvGrid", "prelim_order", "ordered_features_idx", "ordered_features_idx", "237", "tt_prelim (236) is frozen as well"),
        row("trxRuvGrid", "max_acceptable_k", "max_acceptable_k", "max_acceptable_k", "242"),
        row("trxRuvGrid", "grid_members", "grid_members; grid_cancor_data; grid_control_genes", "grid_results", "247-281", "unrounded; complete rows of grid_members"),
        row("trxRuvGrid", "grid_failures", "grid_members", "grid_results (compact)", "281", "rows with status failed"),
        row("trxRuvGrid", "grid_summary", "df_summary", "df_summary", "283-292", "4-dp rounded"),
        row("trxRuvGrid", "best", "best_res", "best_res", "295-296"),
        row("trxRuvGrid", "control_genes", "best_control_genes", "best_control_genes", "298", "inferred"),
        row("trxRuvGrid", "control_source", "-", "(none)", "-", "a label, not a notebook object: the reference labels best_control_genes inferred"),
        row("trxRuvs", "W", "w_matrix", "w_matrix", "417-418"),
        row("trxRuvs", "ruvs_counts", "ruvs_counts", "counts(set_ruvs)", "413", "and ruvs_counts_identical"),
        row("trxRuvs", "log_cpm_ruvs", "log_cpm_corrected", "log_cpm_corrected", "421-426", "visualisation only"),
        row("trxRuvs", "control_source", "-", "(none)", "-", "a label, as above"),
        row("trxDaContrast", "design_final", "design_final", "design_final", "532-533"),
        row("trxDaContrast", "design_formula", "design_formula_str", "design_formula_str", "529"),
        row("trxDaContrast", "contrast_entry", "contrast_entry", "contrast_entry", "542"),
        row("trxDaContrast", "contrasts_matrix", "contrasts_matrix", "contrasts_matrix", "547-548"),
        row("trxDaContrast", "fit", "fit_final", "fit_final", "536-537"),
        row("trxDaContrast", "qlf", "qlf; qlf_table", "lrt", "554", "the notebook names it lrt"),
        row("trxDaContrast", "top_table", "top_genes", "top_genes", "555"),
        row("trxDaContrast", "method", "-", "(none)", "-", "a literal (edgeR_QL_robust), not a notebook object"),
        row("trxAnnotateDe", "feature_meta", "feature_meta", "feature_meta", "144-148"),
        row("trxAnnotateDe", "de_annotated", "top_genes_annotated", "top_genes_annotated", "558-562"),
        row("trxAnnotateDe", "top_de", "top_de_combined", "top_de_combined", "569-582", "top_up and top_down frozen as well"),
        row("trxRuviii", "newY", "ruviii_newY", "res_ruviii$newY", "711-712", "M from group (notebook default)"),
        row("trxRuviii", "M", "ruviii_M", "M", "704"),
        row("trxRuviii", "m_source", "-", "(none)", "-", "the notebook's M is from group; the reference label says so"),
        row("trxRuviii", "m_label", "-", "(none)", "-", "a label, as above"),
        row("trxRuviii", "control_genes", "best_control_genes; ruviii_ctl_mask", "best_control_genes, ctl_mask", "707-708"),
        row("trxRuviii", "control_source", "-", "(none)", "-", "a label, as above"),
        row("trxRuviii", "k", "best_k", "best_k", "297, 711"),
        row("trxPcaStats", "pca_scores", "pca_<stage>_scores", "p_pca_<stage>$data", "472-488, 492-494, 731", "stages raw, tmm, ruvs, post_ruviii"),
        row("trxPcaStats", "pca_variance", "pca_<stage>_variance; pca_<stage>_axis_labels", "summary(prcomp)$importance (labels)", "474", "unrounded from the plot's scores; the labels carry the notebook's rounded values"),
        row("trxRleStats", "rle_matrix", "rle_<stage>_data", "p_rle_<stage>$data", "441-447", "long form of the centred matrix, with the design"),
        row("trxRleStats", "rle_whiskers", "rle_<stage>_whiskers", "whisker_limits (local)", "449-458", "derived from the plot's data by the notebook's rule"),
        row("trxRleStats", "rle_ylim", "rle_<stage>_ylim", "ylim_vals (local)", "459-465", "the plot's coord_cartesian limits")
    )
}

trxRunReference <- function(corpus, work, out, fail_before = NA_character_) {
    repo <- trxRepoRoot()
    driver <- file.path(repo, "tools/profiling/trx_notebook_reference.R")
    notebook <- trxNotebookSource(repo)
    corpus <- normalizePath(corpus, mustWork = TRUE)
    dir.create(out, recursive = TRUE, showWarnings = FALSE); out <- normalizePath(out)
    if (!nzchar(Sys.getenv("R_LIBS"))) trxAbort("R_LIBS must name the kit library", "kit_trx_missing_package")
    trxPreflight(corpus, trxNotebookPackages(notebook))
    work <- trxPrepare(corpus, work, notebook, repo)
    rscript <- file.path(R.home("bin"), "Rscript")
    started <- Sys.time()
    status <- system2("timeout", c(TRX_CEILING_SECONDS, shQuote(rscript), "--vanilla", shQuote(driver), "--render-child", shQuote(work), shQuote(out),
        if (!is.na(fail_before)) c("--fail-before", shQuote(fail_before))), stdout = file.path(out, "render.log"), stderr = file.path(out, "render.log"))
    seconds <- as.numeric(difftime(Sys.time(), started, units = "secs"))
    child <- if (file.exists(file.path(out, "child.rds"))) readRDS(file.path(out, "child.rds")) else NULL
    receipt <- file.path(repo, "docs/qa/workflow-execution/kit/kit-040/library-receipt.json")
    binding <- list(corpus = corpus, corpus_manifest_sha256 = trxSha256(file.path(corpus, "manifest.json")),
        corpus_manifest = jsonlite::read_json(file.path(corpus, "manifest.json")),
        notebook_copy_sha256 = trxSha256(file.path(work, "transcriptomics_workflow.rmd")), driver_sha256 = trxSha256(driver),
        library = Sys.getenv("R_LIBS"), library_receipt_sha256 = if (file.exists(receipt)) trxSha256(receipt) else NA_character_,
        library_receipt = if (file.exists(receipt)) jsonlite::read_json(receipt) else NULL,
        render_seed = TRX_RENDER_SEED, ceiling_seconds = TRX_CEILING_SECONDS, work = work, fail_before = fail_before)
    outputs <- trxOutputHashes(work)
    utils::write.table(outputs, file.path(out, "output-hashes.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
    utils::write.table(trxProductMap(), file.path(out, "product-map.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
    complete <- identical(status, 0L) && isTRUE(child$complete)
    if (!is.null(child)) {
        products <- child$products
        prints <- if (length(products)) do.call(rbind, lapply(names(products), function(n) cbind(data.frame(product = n, stringsAsFactors = FALSE),
            trxDescribe(products[[n]]), source = child$labels[[n]]$source, label = child$labels[[n]]$label, derived = child$labels[[n]]$derived))) else NULL
        binding$session <- child$session
        binding$product_sha256 <- if (!is.null(prints)) stats::setNames(as.list(prints$sha256), prints$product) else list()
        reference <- list(schema = "kit.trx_reference.v1", complete = complete, rendered = child$rendered, error = child$error,
            absent = child$absent, products = products, labels = child$labels, binding = binding, outputs = outputs,
            started = child$started, finished = child$finished, seconds = seconds, exit_status = status)
        saveRDS(reference, file.path(out, "products.rds"), version = 3L)
        store <- list(schema = "kit.trx_fingerprint.v1", complete = complete, incomplete = !complete, absent = child$absent,
            corpus = binding$corpus_manifest$profile, samples = binding$corpus_manifest$samples, corpus_manifest_sha256 = binding$corpus_manifest_sha256,
            notebook_copy_sha256 = binding$notebook_copy_sha256, driver_sha256 = binding$driver_sha256,
            library_receipt_sha256 = binding$library_receipt_sha256, products_rds_sha256 = trxSha256(file.path(out, "products.rds")),
            fingerprints = prints, stable_outputs = outputs[outputs$class == "stable", ], recorded_at = child$finished, seconds = seconds)
        saveRDS(store, file.path(out, "fingerprints.rds"), version = 3L)
    }
    status_record <- list(complete = complete, exit_status = status, timed_out = identical(status, 124L), seconds = seconds,
        error = child$error, absent = child$absent, fail_before = fail_before, markers = nrow(utils::read.delim(file.path(out, "chunk-markers.tsv"))),
        corpus_manifest_sha256 = binding$corpus_manifest_sha256, notebook_copy_sha256 = binding$notebook_copy_sha256, driver_sha256 = binding$driver_sha256)
    jsonlite::write_json(status_record, file.path(out, "status.json"), auto_unbox = TRUE, pretty = TRUE, null = "null", digits = NA)
    status_record
}

if (sys.nframe() == 0L) {
    args <- commandArgs(TRUE)
    if (length(args) >= 1L && identical(args[[1L]], "--render-child")) {
        fail <- if (length(args) >= 5L && identical(args[[4L]], "--fail-before")) args[[5L]] else NA_character_
        trxRenderChild(args[[2L]], args[[3L]], fail)
    }
    if (length(args) != 3L && !(length(args) == 5L && identical(args[[4L]], "--fail-before"))) {
        stop("usage: R_LIBS=<kit library> Rscript --vanilla trx_notebook_reference.R <corpus dir> <work dir> <out dir> [--fail-before <chunk>]")
    }
    result <- trxRunReference(args[[1L]], args[[2L]], args[[3L]], if (length(args) == 5L) args[[5L]] else NA_character_)
    cat(sprintf("TRX_REFERENCE complete=%s exit=%s seconds=%.1f absent=%d error=%s\n", result$complete, result$exit_status, result$seconds,
        length(result$absent), if (is.null(result$error)) "none" else paste(result$error$class, collapse = "/")))
    quit(save = "no", status = if (isTRUE(result$complete)) 0L else 1L)
}
