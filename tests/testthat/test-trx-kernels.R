# KIT-093: the notebook stages as side-effect-free kernels, each product identical to its KIT-040 reference product
# (docs/qa/workflow-execution/kit/kit-040/reference-<corpus>.rds; criteria in docs/qa/workflow-execution/kit/kit-093/
# parity-criteria.md). A product is identical when the sha256 of serialize(x, NULL, version = 3) equals the frozen one.

TRX_STORE_SHA256 <- c(
    trx_edge_12 = "f460b376008219c667fdca32bd1184813ff1f5b8ba8c53d31f45c63b728cb12c",
    trx_edge_16 = "ac118c0ada3411ad5510efdd385e0ae08a4b6b3cecfb6c9a4ba436c75a7ecf59",
    trx_screen = "dd058d7bd344e606899f95dac67786ef8a60df473ffede1e71918bdca0d66840")
TRX_PER_GROUP <- c(trx_edge_12 = 3L, trx_edge_16 = 4L, trx_screen = 6L)

trxTestRoot <- function() normalizePath(testthat::test_path("..", ".."), mustWork = FALSE)

skipWithoutTrxReference <- function() {
    root <- trxTestRoot()
    testthat::skip_if_not(file.exists(file.path(root, "docs/qa/workflow-execution/kit/kit-040/reference-trx_screen.rds")) &&
        file.exists(file.path(root, "tools/profiling/kit_omics_corpus.R")), "the KIT-040 references are not shipped (R CMD check)")
}

trxSha <- function(x) digest::digest(serialize(x, NULL, version = 3L), algo = "sha256", serialize = FALSE)

# The store, checked against KIT-040's recorded hash, and the corpus regenerated and checked against the store's manifest hash.
trxTestReference <- function(corpus, env = parent.frame()) {
    path <- file.path(trxTestRoot(), "docs/qa/workflow-execution/kit/kit-040", sprintf("reference-%s.rds", corpus))
    testthat::expect_identical(digest::digest(file = path, algo = "sha256"), TRX_STORE_SHA256[[corpus]])
    store <- readRDS(path)
    h <- new.env(parent = globalenv())
    sys.source(file.path(trxTestRoot(), "tools/profiling/kit_omics_corpus.R"), h)
    dir <- withr::local_tempdir(.local_envir = env)
    h$writeKitTrxCorpus(dir, h$kitOmicsRecipe("trx_screen", TRX_PER_GROUP[[corpus]]))
    testthat::expect_identical(digest::digest(file = file.path(dir, "manifest.json"), algo = "sha256"), store$corpus_manifest_sha256)
    list(store = store, dir = dir, corpus = corpus)
}

# Read as the notebook reads (rmd:121, 136); the class does the rest.
trxTestInput <- function(dir) {
    list(counts = vroom::vroom(file.path(dir, "counts.tsv"), show_col_types = FALSE),
        design = vroom::vroom(file.path(dir, "design_matrix.tab"), show_col_types = FALSE),
        contrast = trxParseContrastFile(file.path(dir, "contrast_strings.tab")))
}

trxTestPipeline <- function(input, group_id = "group") {
    object <- createTranscriptomicsAssayData(list(counts = input$counts), input$design, c("RPMI", "Sera"), roles = list(group_id = group_id))
    x <- trxAssayInputs(object)
    f <- trxExpressionFilter(x$counts, x$design$group)
    n <- trxTmmLogCpm(f$filtered_counts, x$design$group)
    g <- trxRuvGrid(f$filtered_counts, n$log_cpm_tmm, x$design)
    s <- trxRuvs(f$filtered_counts, n$log_cpm_tmm, x$design, g$control_genes, g$best$best_k)
    d <- trxDaContrast(s$ruvs_counts, s$pdata, input$contrast)
    a <- trxAnnotateDe(d$top_table, object@count_data[[1L]], colnames(d$contrasts_matrix))
    r <- trxRuviii(n$log_cpm_tmm, x$design, g$control_genes, g$best$best_k)
    stages <- list(raw = n$log_cpm_raw, tmm = n$log_cpm_tmm, ruvs = s$log_cpm_ruvs, post_ruviii = t(r$newY))
    list(inputs = x, filter = f, norm = n, grid = g, ruvs = s, da = d, annot = a, ruviii = r,
        pca = lapply(stages, trxPcaStats, design = x$design), rle = lapply(stages, trxRleStats, design = x$design))
}

# The driver's grid_members shape (one row per percentage; a failed member NA but for its percentage and status).
trxTestGridMembers <- function(grid, percentage_range = seq(5, 30, by = 5)) {
    ok <- stats::setNames(grid$grid_members, vapply(grid$grid_members, function(r) as.character(r$percentage), ""))
    do.call(rbind, lapply(percentage_range, function(pct) {
        r <- ok[[as.character(pct)]]
        if (is.null(r)) return(data.frame(percentage = pct, status = "failed", n_controls = NA_real_, best_k = NA_integer_,
            separation_score = NA_real_, k_penalty = NA_real_, composite_score = NA_real_))
        data.frame(percentage = r$percentage, status = "complete", n_controls = r$n_controls, best_k = r$best_k,
            separation_score = r$separation_score, k_penalty = r$k_penalty, composite_score = r$composite_score)
    }))
}

trxBestFields <- c("percentage", "n_controls", "best_k", "separation_score", "k_penalty", "composite_score", "ctl_idx", "control_genes")

# Reference product name -> kernel product (kit-040/product-map.tsv).
trxKernelProducts <- function(p) {
    members <- p$grid$grid_members
    pct <- vapply(members, function(r) paste0("pct_", r$percentage), "")
    out <- list(
        filter_effective_lib_size = p$filter$effective_lib_size, filter_min_sample_size = p$filter$min_sample_size,
        filter_median_lib_size = p$filter$median_lib_size, keep_expr = p$filter$keep_expr, gene_vars_positive = p$filter$keep_var,
        count_matrix_filtered = p$filter$filtered_counts, dge_filtered_samples = p$norm$norm_factors_filtered,
        log_cpm_tmm = p$norm$log_cpm_tmm, log_cpm_raw = p$norm$log_cpm_raw,
        ordered_features_idx = p$grid$prelim_order, max_acceptable_k = p$grid$max_acceptable_k, grid_members = trxTestGridMembers(p$grid),
        grid_cancor_data = stats::setNames(lapply(members, `[[`, "cancor_data"), pct),
        grid_control_genes = stats::setNames(lapply(members, `[[`, "control_genes"), pct),
        df_summary = p$grid$grid_summary, best_res = p$grid$best[trxBestFields], best_control_genes = p$grid$control_genes,
        w_matrix = p$ruvs$W, ruvs_counts = p$ruvs$ruvs_counts, log_cpm_corrected = p$ruvs$log_cpm_ruvs,
        design_final = p$da$design_final, design_formula_str = p$da$design_formula, contrast_entry = p$da$contrast_entry,
        contrasts_matrix = p$da$contrasts_matrix, fit_final = p$da$fit, qlf = p$da$qlf, qlf_table = p$da$qlf$table, top_genes = p$da$top_table,
        feature_meta = p$annot$feature_meta, top_genes_annotated = p$annot$de_annotated, top_de_combined = p$annot$top_de,
        ruviii_newY = p$ruviii$newY, ruviii_M = p$ruviii$M, ruviii_ctl_mask = p$ruviii$ctl_mask, best_k = p$ruviii$k)
    for (stage in names(p$pca)) {
        out[[paste0("pca_", stage, "_scores")]] <- p$pca[[stage]]$pca_scores
        out[[paste0("pca_", stage, "_variance")]] <- p$pca[[stage]]$pca_variance
        out[[paste0("rle_", stage, "_data")]] <- p$rle[[stage]]$rle_matrix
        out[[paste0("rle_", stage, "_whiskers")]] <- p$rle[[stage]]$rle_whiskers
        out[[paste0("rle_", stage, "_ylim")]] <- p$rle[[stage]]$rle_ylim
    }
    out
}

# Reference products that are not kernel products: renderer products (test-trx-publication.R) and notebook intermediates.
TRX_NOT_KERNEL <- c("gene_vars", "dge_raw_samples", "tt_prelim", "ruvs_counts_identical", "top_up", "top_down", "glimma_html_exists")

expectTrxProducts <- function(products, reference, what = "") {
    prints <- reference$store$fingerprints
    for (name in names(products)) {
        expect_identical(trxSha(products[[name]]), prints$sha256[prints$product == name], info = paste(reference$corpus, what, name))
    }
}

test_that("KIT-093: every kernel product is identical to its KIT-040 reference product on the three corpora", {
    skipWithoutTrxReference()
    for (corpus in c("trx_screen", "trx_edge_12", "trx_edge_16")) {
        reference <- trxTestReference(corpus)
        p <- trxTestPipeline(trxTestInput(reference$dir))
        products <- trxKernelProducts(p)
        renderer <- grepl("_(layers|annotation|axis_labels|box)$|^(cca_|volcano_|bar_|glimma_anno|summary_text)", reference$store$fingerprints$product)
        expect_setequal(c(names(products), TRX_NOT_KERNEL), reference$store$fingerprints$product[!renderer])
        expect_identical(p$grid$max_acceptable_k, if (corpus == "trx_edge_12") 3L else 4L)
        expectTrxProducts(products, reference)
        expect_identical(p$ruvs$ruvs_counts, p$filter$filtered_counts, info = corpus)
    }
})

test_that("KIT-093: two consecutive runs on the screen corpus are identical; shuffled columns and another group_id change nothing", {
    skipWithoutTrxReference()
    reference <- trxTestReference("trx_screen")
    input <- trxTestInput(reference$dir)
    expect_identical(trxTestPipeline(input), trxTestPipeline(input))
    edge <- trxTestReference("trx_edge_12")
    input <- trxTestInput(edge$dir)
    runs <- setdiff(names(input$counts), c("GeneID", "GeneName"))
    shuffled <- input
    shuffled$counts <- input$counts[, c("GeneName", rev(runs), "GeneID")]
    expect_false(identical(names(shuffled$counts), names(input$counts)))
    expectTrxProducts(trxKernelProducts(trxTestPipeline(shuffled)), edge, "shuffled")
    renamed <- input
    names(renamed$design)[names(renamed$design) == "group"] <- "condition"
    expectTrxProducts(trxKernelProducts(trxTestPipeline(renamed, group_id = "condition")), edge, "group_id condition")
})

test_that("KIT-093: a constructed tie selects the first member; a forced member failure is recorded and leaves the winner", {
    member <- function(pct, score) list(percentage = pct, composite_score = score)
    expect_identical(trxSelectGridWinner(list(member(5, 0.5), member(10, 0.9), member(15, 0.9))), 2L)
    expect_identical(trxSelectGridWinner(list(member(5, 0.9), member(10, 0.9))), 1L)
    skipWithoutTrxReference()
    reference <- trxTestReference("trx_screen")
    x <- trxAssayInputs(createTranscriptomicsAssayData(list(counts = trxTestInput(reference$dir)$counts),
        trxTestInput(reference$dir)$design, c("RPMI", "Sera")))
    f <- trxExpressionFilter(x$counts, x$design$group)
    n <- trxTmmLogCpm(f$filtered_counts, x$design$group)
    failing <- max(10, floor(nrow(n$log_cpm_tmm) * 30 / 100))
    real <- trxCancor
    local_mocked_bindings(trxCancor = function(Y, X, ctl) {
        if (sum(ctl) == failing) stop(structure(class = c("kit093_forced_failure", "error", "condition"), list(message = "forced", call = NULL)))
        real(Y, X, ctl)
    })
    g <- trxRuvGrid(f$filtered_counts, n$log_cpm_tmm, x$design)
    expect_identical(g$grid_failures, data.frame(percentage = 30, condition_class = "kit093_forced_failure"))
    expect_identical(vapply(g$grid_members, `[[`, numeric(1), "percentage"), c(5, 10, 15, 20, 25))
    prints <- reference$store$fingerprints
    expect_identical(trxSha(g$best[trxBestFields]), prints$sha256[prints$product == "best_res"])
    expect_identical(trxSha(g$control_genes), prints$sha256[prints$product == "best_control_genes"])
    failed <- trxTestGridMembers(g)
    expect_identical(failed$status, c(rep("complete", 5L), "failed"))
    local_mocked_bindings(trxCancor = function(Y, X, ctl) stop("forced"))
    expect_error(trxRuvGrid(f$filtered_counts, n$log_cpm_tmm, x$design), class = "trx_grid_no_member")
    expect_error(trxRuvGrid(f$filtered_counts, n$log_cpm_tmm, x$design), class = "trx_error")
})

test_that("KIT-093: technical-replicate RUVIII equals the direct ruv::RUVIII call; labels and FDR (BH) are set", {
    skipWithoutTrxReference()
    reference <- trxTestReference("trx_screen")
    p <- trxTestPipeline(trxTestInput(reference$dir))
    tech <- utils::read.delim(file.path(reference$dir, "design.tsv"))
    design <- p$inputs$design
    design$tech_rep_group <- tech$tech_rep_group[match(design$Run, tech$sample)]
    k <- p$grid$best$best_k
    r <- trxRuviii(p$norm$log_cpm_tmm, design, p$grid$control_genes, k, m_source = "technical_replicate", technical_replicate_id = "tech_rep_group")
    direct <- ruv::RUVIII(Y = t(p$norm$log_cpm_tmm), M = ruv::replicate.matrix(design$tech_rep_group), ctl = p$ruviii$ctl_mask,
        k = k, return.info = TRUE)$newY
    expect_identical(r$newY, direct)
    expect_false(identical(r$newY, p$ruviii$newY))
    expect_identical(r$m_source, "technical_replicate")
    expect_identical(r$m_label, "M from technical replicates")
    expect_identical(p$ruviii$m_source, "group")
    expect_identical(p$ruviii$m_label, "M from treatment group (notebook default)")
    for (product in list(p$grid, p$ruvs, p$ruviii, r)) expect_identical(product$control_source, "inferred_least_significant")
    expect_true("FDR" %in% names(p$da$top_table))
    expect_false(any(grepl("q[._]?val", names(p$da$top_table), ignore.case = TRUE)))
    expect_identical(p$da$top_table$FDR, stats::p.adjust(p$da$top_table$PValue, method = "BH"))
    expect_identical(p$da$method, "edgeR_QL_robust")
    # The refusals: absent or NA column, a set spanning two groups (the biological `replicates` number), no set of two.
    refuse <- function(design, column) {
        expect_error(trxRuviii(p$norm$log_cpm_tmm, design, p$grid$control_genes, k, "technical_replicate", column), class = "trx_ruviii_replicates_invalid")
    }
    refuse(design, "no_such_column")
    refuse(design, NA_character_)
    with_na <- design; with_na$tech_rep_group[[3]] <- NA
    refuse(with_na, "tech_rep_group")
    refuse(design, "replicates")
    singles <- design; singles$tech_rep_group <- singles$Run
    refuse(singles, "tech_rep_group")
})

test_that("KIT-093: the preflight refusals are raised by class, all under trx_error", {
    skipWithoutTrxReference()
    group <- factor(rep(c("RPMI", "Sera"), each = 3L), levels = c("RPMI", "Sera"))
    flat <- matrix(0, nrow = 10L, ncol = 6L, dimnames = list(sprintf("g%02d", 1:10), sprintf("s%d", 1:6)))
    flat[1L, ] <- 1
    expect_error(trxExpressionFilter(flat, group), class = "trx_filter_empty")
    expect_error(trxExpressionFilter(flat, group), class = "trx_error")
    expect_error(trxExpressionFilter(flat * 0.5, group), class = "trx_counts_fractional")
    reference <- trxTestReference("trx_edge_12")
    p <- trxTestPipeline(trxTestInput(reference$dir))
    few <- p$filter$filtered_counts[1:9, ]
    expect_error(trxRuvGrid(few, p$norm$log_cpm_tmm[1:9, ], p$inputs$design), class = "trx_grid_too_few_features")
    confounded <- p$ruvs$pdata
    confounded$W_1 <- as.numeric(confounded$group == "Sera")
    expect_error(trxDaContrast(p$ruvs$ruvs_counts, confounded, p$da$contrast_entry), class = "trx_da_rank_deficient")
    expect_error(trxDaContrast(p$ruvs$ruvs_counts, p$ruvs$pdata, "Sera_vs_RPMI groupSera-groupRPMI"), class = "trx_contrast_unparseable")
    expect_error(trxDaContrast(p$ruvs$ruvs_counts, p$ruvs$pdata, NA_character_), class = "trx_contrast_unparseable")
    header_only <- withr::local_tempfile(lines = "contrasts")
    expect_error(trxParseContrastFile(header_only), class = "trx_contrast_unparseable")
    expect_error(trxDaContrast(p$ruvs$ruvs_counts, p$ruvs$pdata, "x=groupLPS-groupRPMI"), class = "trx_contrast_not_estimable")
    expect_error(trxDaContrast(p$ruvs$ruvs_counts, p$ruvs$pdata, "x=groupSera-groupSera"), class = "trx_contrast_not_estimable")
    expect_error(trxDaContrast(p$ruvs$ruvs_counts, p$ruvs$pdata, "x=groupSera-groupSera"), class = "trx_error")
})
