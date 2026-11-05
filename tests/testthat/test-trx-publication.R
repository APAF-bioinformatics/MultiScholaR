# KIT-093: the notebook's publication renderers, each figure compared with its KIT-040 reference by layer_data() of every
# layer and by its plotted data (the KIT-040 addendum of 2026-11-06), the tables and the summary text as the notebook writes
# them. Criteria: docs/qa/workflow-execution/kit/kit-093/parity-criteria.md. The reference helpers repeat those of
# test-trx-kernels.R (each test file runs on its own).

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


# The driver's RLE box columns of layer_data(p, 1).
trxTestRleBox <- function(p) {
    d <- ggplot2::layer_data(p, 1L)
    d[, intersect(c("x", "group", "PANEL", "fill", "ymin", "lower", "middle", "upper", "ymax", "notchlower", "notchupper", "relvarwidth", "width"), names(d))]
}

trxPanelProducts <- function(panels, stages, suffix = "") {
    out <- list()
    for (stage in stages) {
        pca <- panels[[paste0("pca_", stage)]]; rle <- panels[[paste0("rle_", stage)]]
        out[[paste0("pca_", stage, "_layers", suffix)]] <- pca$layers
        out[[paste0("pca_", stage, "_scores", suffix)]] <- pca$data
        out[[paste0("pca_", stage, "_axis_labels", suffix)]] <- unlist(ggplot2::get_labs(pca$plot)[c("x", "y")])
        out[[paste0("rle_", stage, "_layers", suffix)]] <- rle$layers
        out[[paste0("rle_", stage, "_data", suffix)]] <- rle$data
        out[[paste0("rle_", stage, "_box", suffix)]] <- trxTestRleBox(rle$plot)
    }
    out
}

# Reference product name -> renderer product (kit-040/product-map.tsv, renderer rows).
trxRendererProducts <- function(p) {
    contrast_name <- colnames(p$da$contrasts_matrix)
    cca <- trxCcaPanels(p$grid$best, p$grid$grid_summary, nrow(p$norm$log_cpm_tmm))
    norm <- trxNormStageGrid(p$norm$log_cpm_raw, p$norm$log_cpm_tmm, p$ruvs$log_cpm_ruvs, p$inputs$design, p$grid$best$best_k)
    ruviii <- trxRuviiiGrid(p$norm$log_cpm_tmm, p$ruviii$newY, p$inputs$design)
    volcano <- trxVolcano(p$annot$de_annotated, contrast_name)
    bar <- trxSignificanceBar(p$annot$de_annotated, contrast_name)
    tables <- trxDeTables(p$annot$de_annotated, p$annot$top_de)
    glimma <- trxGlimmaInputs(p$annot$de_annotated, p$da$qlf, contrast_name)
    c(list(cca_cancor_layers = cca$cancor$layers, cca_cancor_data = cca$cancor$data, cca_tuning_layers = cca$tuning$layers,
        cca_tuning_data = cca$tuning$data, cca_grid_layers = cca$grid$layers, cca_optimal_parameters = cca$workbook$Optimal_Parameters,
        df_summary = cca$workbook$Grid_Optimization, norm_grid_layers = norm$grid$layers, norm_grid_annotation = norm$grid$annotation,
        ruviii_grid_layers = ruviii$grid$layers, ruviii_grid_annotation = ruviii$grid$annotation,
        volcano_layers = volcano$layers, volcano_data = volcano$data, bar_layers = bar$layers, bar_data = bar$data,
        glimma_anno_aligned = glimma$anno, top_genes_annotated = tables$de_table, top_de_combined = tables$top_de),
        trxPanelProducts(norm$panels, c("raw", "tmm", "ruvs")), trxPanelProducts(ruviii$panels, "post_ruviii"),
        trxPanelProducts(ruviii$panels, "tmm", " (ruviii grid)"),
        list(summary_text = trxTestSummaryText(trxSummaryText(p$grid$best, nrow(p$norm$log_cpm_tmm), p$da$design_formula,
            p$da$contrast_entry, p$annot$de_annotated, p$grid$grid_summary))))
}

# Written as the notebook writes it (rmd:756-775), read back without the timestamp line.
trxTestSummaryText <- function(summary) {
    path <- withr::local_tempfile()
    writeLines(summary$lines, path)
    utils::write.table(summary$grid_table, file = path, append = TRUE, sep = "\t", row.names = FALSE, quote = FALSE)
    text <- readLines(path, warn = FALSE)
    text[!startsWith(text, "Timestamp: ")]
}

expectTrxRendererProducts <- function(products, reference) {
    prints <- reference$store$fingerprints
    for (name in names(products)) {
        expect_identical(trxSha(products[[name]]), prints$sha256[prints$product == sub(" .*$", "", name)], info = paste(reference$corpus, name))
    }
}

test_that("KIT-093: every renderer product is identical to its KIT-040 reference product on the three corpora", {
    skipWithoutTrxReference()
    for (corpus in c("trx_screen", "trx_edge_12", "trx_edge_16")) {
        reference <- trxTestReference(corpus)
        p <- trxTestPipeline(trxTestInput(reference$dir))
        products <- trxRendererProducts(p)
        added <- reference$store$addendum$added
        expect_true(all(setdiff(added, "glimma_anno") %in% names(products)), info = corpus)
        expectTrxRendererProducts(products, reference)
        # The DE tables written with the notebook's writers (rmd:566, 583) match the stable output hashes.
        tables <- trxDeTables(p$annot$de_annotated, p$annot$top_de)
        stable <- reference$store$stable_outputs
        csv <- withr::local_tempfile(fileext = ".csv"); tsv <- withr::local_tempfile(fileext = ".tsv")
        utils::write.csv(tables$de_table, csv, row.names = FALSE)
        vroom::vroom_write(tables$top_de, tsv)
        expect_identical(digest::digest(file = csv, algo = "sha256"), stable$sha256[basename(stable$path) == "de_results_Sera_vs_RPMI.csv"], info = corpus)
        expect_identical(digest::digest(file = tsv, algo = "sha256"), stable$sha256[basename(stable$path) == "top_de_genes.tsv"], info = corpus)
    }
})

test_that("KIT-093: the renderers return figures with their data, write no files, label FDR (BH) and keep zero-count bars", {
    skipWithoutTrxReference()
    reference <- trxTestReference("trx_edge_12")
    p <- trxTestPipeline(trxTestInput(reference$dir))
    before <- list.files(tempdir(), recursive = TRUE)
    volcano <- trxVolcano(p$annot$de_annotated, "Sera_vs_RPMI")
    expect_s3_class(volcano$plot, "ggplot")
    expect_identical(volcano$layers, lapply(seq_along(volcano$plot$layers), function(i) ggplot2::layer_data(volcano$plot, i)))
    expect_match(paste(deparse(ggplot2::get_labs(volcano$plot)$y), collapse = ""), "BH")
    none <- p$annot$de_annotated
    none$FDR <- 1
    bar <- trxSignificanceBar(none, "Sera_vs_RPMI")
    expect_identical(as.character(bar$data$significance), c("Upregulated", "Downregulated", "Not Significant"))
    expect_identical(bar$data$n, c(0L, 0L, nrow(none)))
    glimma <- trxGlimmaInputs(p$annot$de_annotated, p$da$qlf, "Sera_vs_RPMI")
    expect_identical(rownames(glimma$anno), rownames(p$da$qlf$table))
    expect_identical(glimma$display_columns, c("GeneID", "GeneName", "DisplaySymbol"))
    expect_identical(c(glimma$lfc, glimma$fdr), c(1, 0.05))
    summary <- trxSummaryText(p$grid$best, nrow(p$norm$log_cpm_tmm), p$da$design_formula, p$da$contrast_entry, p$annot$de_annotated, p$grid$grid_summary)
    expect_match(summary$lines[[2]], "^Timestamp: ")
    expect_identical(list.files(tempdir(), recursive = TRUE), before)
})
