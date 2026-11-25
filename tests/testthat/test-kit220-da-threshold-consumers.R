# KIT-220 (VAL-191): every published significant count, figure and gene list follows the DA threshold the analysis ran
# at. Proteomics: the first count table, the published NumSig table (with its `q_threshold` column), the report's count
# sentence and the interactive volcano's status (the user's q and fold-change cutoff). Transcriptomics: the published
# volcano, significance bar, top_de_genes.tsv, the summary line and the interactive status. E-7: the DA tab's summary
# states the run's threshold and says when the box differs from it.
# Every expected count is computed here from the returned long table at the threshold under test, never hard-coded, so
# each row holds on any fixture that has rows between the thresholds (ticket 10, stop rule 2; G0 counts them).

# ---- proteomics fixture: graded effects, so q values fall between 0.01, 0.05 and 0.10 -------------------------------

kit220ProtObject <- function(seed = 20261125L, proteins = 400L, per_group = 4L) {
    set.seed(seed)
    groups <- rep(c("control", "treated"), each = per_group)
    n <- length(groups)
    effect <- numeric(proteins)
    moved <- seq_len(round(0.4 * proteins))
    effect[moved] <- seq(0.3, 2.4, length.out = length(moved)) * rep(c(1, -1), length.out = length(moved))
    x <- matrix(stats::rnorm(proteins, 22, 2), proteins, n) + outer(effect, as.numeric(groups == "treated")) +
        matrix(stats::rnorm(proteins * n, 0, 0.5), proteins, n)
    ids <- sprintf("PG%04d", seq_len(proteins))
    samples <- sprintf("run_%02d", seq_len(n))
    dimnames(x) <- list(ids, samples)
    quant <- data.frame(Protein.Group = ids, as.data.frame(x), check.names = FALSE, stringsAsFactors = FALSE)
    rownames(quant) <- NULL # the DA kernel moves the identifier column to row names (instrument fix, run 45498530109)
    design <- data.frame(Run = samples, group = groups, replicates = seq_len(n),
        tech_rep_group = paste0("u", seq_len(n)), stringsAsFactors = FALSE)
    methods::new("ProteinQuantitativeData", protein_quant_table = quant, protein_id_column = "Protein.Group",
        protein_id_table = data.frame(Protein.Group = ids, stringsAsFactors = FALSE), design_matrix = design,
        sample_id = "Run", group_id = "group", technical_replicate_id = "tech_rep_group",
        args = list(globalParameters = list(workflow_type = "LFQ")))
}

kit220Contrast <- "treated_vs_control=grouptreated-groupcontrol"

kit220ProtDa <- function(object, qvalue, lfc = 0) {
    result <- NULL
    utils::capture.output(result <- suppressWarnings(suppressMessages(differentialAbundanceAnalysis(object,
        contrasts_tbl = data.frame(contrasts = kit220Contrast, stringsAsFactors = FALSE), formula_string = "~ 0 + group",
        da_q_val_thresh = qvalue, treat_lfc_cutoff = lfc, eBayes_trend = TRUE, eBayes_robust = TRUE,
        args_row_id = "Protein.Group", qvalue_column = "fdr_qvalue", raw_pvalue_column = "raw_pvalue"))))
    result
}

#' Up and down at threshold t by the published (sign) rule, from a long table.
kit220SignCounts <- function(long, t) {
    q <- long$fdr_qvalue
    c(up = sum(q < t & long$log2FC >= 0, na.rm = TRUE), down = sum(q < t & long$log2FC < 0, na.rm = TRUE))
}

#' Up and down of a count table (`status`, `counts`).
kit220TableCounts <- function(table) {
    table <- as.data.frame(table)
    c(up = sum(table$counts[table$status == "Significant and Up"]),
        down = sum(table$counts[table$status == "Significant and Down"]))
}

#' The rows of a long table between two thresholds (stop rule 2's count).
kit220Between <- function(long, low, high) sum(long$fdr_qvalue >= low & long$fdr_qvalue < high, na.rm = TRUE)

#' The report template's FormatDESummaryTable, taken from the template's source and evaluated on its own.
kit220FormatDESummaryTable <- function(template) {
    template_dir <- system.file("reports", "proteomics", package = "MultiScholaR")
    if (!nzchar(template_dir)) template_dir <- file.path(rprojroot::find_package_root_file(), "inst", "reports", "proteomics")
    lines <- readLines(file.path(template_dir, template), warn = FALSE)
    start <- grep("^FormatDESummaryTable <- function", lines)
    stopifnot(length(start) == 1L)
    exprs <- NULL
    for (end in seq(start, length(lines))) {
        exprs <- tryCatch(parse(text = lines[start:end]), error = function(e) NULL)
        if (!is.null(exprs) && length(exprs) == 1L) break
    }
    env <- new.env(parent = globalenv())
    env$`%>%` <- magrittr::`%>%`
    eval(exprs[[1L]], env)
    env$FormatDESummaryTable
}

#' The source of an R file of the package, parsed.
kit220PackageFile <- function(...) {
    root <- tryCatch(rprojroot::find_package_root_file(), error = function(e) NULL)
    if (is.null(root)) testthat::skip("the package source is not available")
    file.path(root, ...)
}

#' Every call to `fun` in a parsed file, as calls.
kit220CallsTo <- function(path, fun) {
    found <- list()
    walk <- function(e) {
        if (is.call(e)) {
            if (identical(e[[1L]], as.name(fun))) found[[length(found) + 1L]] <<- e
            for (i in seq_along(e)) if (!identical(e[[i]], quote(expr = )) && !is.pairlist(e[[i]])) walk(e[[i]])
        } else if (is.expression(e)) for (i in seq_along(e)) walk(e[[i]])
    }
    walk(parse(path, keep.source = FALSE))
    found
}

# ---- T-P1: the first count table follows the threshold ----------------------------------------------------------

test_that("KIT-220 T-P1: the first count table counts at the run's threshold (0.01 and 0.10), as the second table does", {
    object <- kit220ProtObject()
    for (t in c(0.01, 0.10)) {
        result <- kit220ProtDa(object, t)
        long <- result$da_proteins_long
        expected <- kit220SignCounts(long, t)
        expect_identical(kit220TableCounts(result$num_sig_da_molecules_first_go$table), expected,
            info = sprintf("t = %s; rows between t and 0.05: %d", t, kit220Between(long, min(t, 0.05), max(t, 0.05))))
        expect_identical(kit220TableCounts(result$num_sig_da_molecules), expected, info = paste("second table, t =", t))
    }
})

# ---- T-P2: TREAT at 0.5 ------------------------------------------------------------------------------------------

test_that("KIT-220 T-P2: with TREAT at 0.5 the first count table counts at the run's threshold; the tab counts by its rule", {
    object <- kit220ProtObject()
    t <- 0.01
    result <- kit220ProtDa(object, t, lfc = 0.5)
    long <- result$da_proteins_long
    first <- kit220TableCounts(result$num_sig_da_molecules_first_go$table)
    expect_identical(first, kit220SignCounts(long, t))
    # The tab (list branch, box set to the run's values): its own rule, log2FC beyond the cutoff.
    text <- diaDaSummaryStatsView(list(da_proteins_long = long), unique(long$comparison)[[1L]], t, 0.5)
    tab_up <- sum(long$fdr_qvalue < t & long$log2FC > 0.5, na.rm = TRUE)
    tab_down <- sum(long$fdr_qvalue < t & long$log2FC < -0.5, na.rm = TRUE)
    expect_match(text, sprintf("Up-regulated: %d", tab_up), fixed = TRUE)
    expect_match(text, sprintf("Down-regulated: %d", tab_down), fixed = TRUE)
    # The two rules agree when no row below t has |log2FC| <= 0.5 (TREAT; ticket 2.3).
    inside <- sum(long$fdr_qvalue < t & abs(long$log2FC) <= 0.5, na.rm = TRUE)
    if (inside == 0L) expect_identical(unname(first), c(tab_up, tab_down))
})

# ---- T-P3 and T-PAR (proteomics): the published table, its column and the report sentence --------------------------

kit220Publish <- function(object, result, t, root) {
    da_dir <- file.path(root, "da")
    graphs <- file.path(root, "graphs")
    object@args$outputDaResultsAllContrasts <- list(uniprot_tbl = data.frame(), da_output_dir = da_dir,
        publication_graphs_dir = graphs, file_prefix = "da_proteins", args_row_id = "Protein.Group",
        gene_names_column = "gene_names", uniprot_id_column = "Entry", da_q_val_thresh = t, fdr_column = "fdr_qvalue",
        log2fc_column = "log2FC")
    utils::capture.output(done <- suppressWarnings(suppressMessages(outputDaResultsAllContrasts(object,
        list(treated_vs_control = result), uniprot_tbl = data.frame(), da_output_dir = da_dir,
        publication_graphs_dir = graphs, file_prefix = "da_proteins", args_row_id = "Protein.Group",
        gene_names_column = "gene_names", uniprot_id_column = "Entry"))))
    file.path(graphs, "NumSigDaMolecules", "da_proteins_num_sig_da_molecules.tab")
}

test_that("KIT-220 T-P3: the published count table counts at the run's threshold and states it; the report sentence too", {
    object <- kit220ProtObject()
    t <- 0.01
    result <- kit220ProtDa(object, t)
    expected <- kit220SignCounts(result$da_proteins_long, t)
    path <- kit220Publish(object, result, t, withr::local_tempdir())
    expect_true(file.exists(path))
    table <- utils::read.delim(path, check.names = FALSE)
    expect_identical(kit220TableCounts(table), expected)
    expect_true("q_threshold" %in% names(table))
    if ("q_threshold" %in% names(table)) expect_equal(unique(table$q_threshold), t)
    for (template in c("DIANN_report.rmd", "LFQ_report.rmd", "TMT_report.rmd")) {
        sentence <- kit220FormatDESummaryTable(template)(path)
        expect_match(sentence, "(q-value < 0.01)", fixed = TRUE, info = template)
        expect_match(sentence, sprintf("there were %d upregulated and %d downregulated proteins", expected[["up"]],
            expected[["down"]]), fixed = TRUE, info = template)
    }
})

test_that("KIT-220 T-PAR (proteomics): at 0.01, the tab's counts (box = the run's threshold) equal the published table's", {
    object <- kit220ProtObject()
    t <- 0.01
    result <- kit220ProtDa(object, t)
    long <- result$da_proteins_long
    text <- diaDaSummaryStatsView(list(da_proteins_long = long), unique(long$comparison)[[1L]], t, 0)
    table <- utils::read.delim(kit220Publish(object, result, t, withr::local_tempdir()), check.names = FALSE)
    published <- kit220TableCounts(table)
    expect_match(text, sprintf("Significant (q < 0.010): %d", sum(published)), fixed = TRUE)
    expect_match(text, sprintf("Up-regulated: %d", published[["up"]]), fixed = TRUE)
    expect_match(text, sprintf("Down-regulated: %d", published[["down"]]), fixed = TRUE)
})

# ---- T-P5: the proteomics interactive volcano's status ------------------------------------------------------------

test_that("KIT-220 T-P5: the proteomics interactive volcano colours by the user's q and fold-change cutoff", {
    captured <- NULL
    testthat::local_mocked_bindings(glimmaXY = function(...) {
        captured <<- list(...)
        htmlwidgets::createWidget("glimmaXY", list(), package = "Glimma")
    }, .package = "Glimma")
    rows <- data.frame(Protein.Ids = paste0("P", 1:8), comparison = "T_vs_C",
        log2FC = c(0.2, -0.2, 0.8, -0.8, 1.5, -1.5, 0.8, 1.5), fdr_qvalue = c(rep(0.01, 6), 0.2, 0.2),
        raw_pvalue = c(rep(0.001, 6), 0.1, 0.1), stringsAsFactors = FALSE)
    for (lfc in c(0, 0.5)) {
        captured <- NULL
        suppressMessages(generateProtDAVolcanoPlotGlimma(list(da_proteins_long = rows), selected_contrast = "T_vs_C",
            args_row_id = "Protein.Ids", da_q_val_thresh = 0.05, lfc_thresh = lfc))
        expected <- as.integer(sign(rows$log2FC) * (rows$fdr_qvalue < 0.05 & abs(rows$log2FC) >= lfc))
        expect_identical(captured$status, expected, info = paste("lfc =", lfc))
    }
})

# ---- T-T1, T-T2 and T-PAR (transcriptomics) ---------------------------------------------------------------------

#' The trx_edge_12 owner after a DA at `t` (KIT-192 T6's path) and its DA publication, with every figure captured.
kit220TrxPublished <- function(t, env = parent.frame()) {
    reference <- trxDeclReference("trx_edge_12", env = env)
    fixture <- trxDeclarationsOwner(withr::local_tempdir(.local_envir = env), reference$dir)
    withr::defer(closeDIAStateAdapter(fixture$owner), envir = env)
    owner <- fixture$owner
    figures <- new.env(parent = emptyenv())
    real_save <- trxSaveFigure
    testthat::local_mocked_bindings(workflowExecutionBuildIdentity = function() strrep("f", 64),
        finishWorkflowExecutionJob = function(owner, status) {
            owner$job <- NULL
            list(status = status)
        },
        trxSaveFigure = function(plot, directory, name, width, height) {
            assign(name, plot, envir = figures)
            real_save(plot, directory, name, width, height)
        }, .package = "MultiScholaR", .env = env)
    suppressWarnings(suppressMessages(runTrxReferenceExpressionFilter(owner, list())))
    commitWorkflowOmicsStoreStep(owner, "trx_qc_finalize", "trx_qc_complete", list(), "QC processing complete - ready for normalization")
    suppressWarnings(suppressMessages(runTrxReferenceRuvs(owner, NULL, list(), design_parameters = trxRuvDesignParameters("none"))))
    suppressWarnings(suppressMessages(runTrxReferenceDaAnalysis(owner, getWorkflowStateReference(owner),
        list(da_q_val_thresh = t, treat_lfc_cutoff = 0), fixture$workflow$contrasts_tbl)))
    display <- suppressWarnings(suppressMessages(runWorkflowOmicsPublicationInProcess(transcriptomicsDeclarationBinding, owner,
        "trx_da_publication", list(), fixture$paths)))
    results <- withWorkflowStateScope(owner, function() {
        manifest <- resolveWorkflowStateValue(owner$.context, owner$.context$snapshot$da_results)
        r <- readWorkflowOmicsDaResults(transcriptomicsDeclarationBinding, owner, manifest)
        r$theObject <- readTrxWorkflowState(owner, manifest$source)
        r
    })
    list(fixture = fixture, display = display, figures = figures, results = results)
}

test_that("KIT-220 T-T1: the transcriptomics publication draws, counts and lists at the run's FDR threshold", {
    skipWithoutTrxDeclReference()
    t <- 0.01
    run <- kit220TrxPublished(t)
    expect_true(isTRUE(run$display$saved))
    volcanoes <- grep("^volcano_plot_", ls(run$figures), value = TRUE)
    bars <- grep("^summary_barplot_", ls(run$figures), value = TRUE)
    expect_gt(length(volcanoes), 0L)
    for (name in volcanoes) {
        plot <- get(name, envir = run$figures)
        hline <- which(vapply(plot$layers, function(layer) inherits(layer$geom, "GeomHline"), logical(1)))
        expect_equal(ggplot2::layer_data(plot, hline[[1L]])$yintercept, -log10(t), info = name)
        data <- plot$data
        expected <- ifelse(data$FDR < t & data$logFC > 0, "Upregulated", ifelse(data$FDR < t & data$logFC < 0,
            "Downregulated", "Not Significant"))
        expect_identical(as.character(data$significance), expected, info = name)
    }
    for (name in bars) {
        counts <- get(name, envir = run$figures)$data
        volcano <- get(sub("^summary_barplot_", "volcano_plot_", name), envir = run$figures)$data
        expect_identical(as.integer(counts$n[counts$significance == "Upregulated"]),
            sum(volcano$FDR < t & volcano$logFC > 0), info = name)
        expect_identical(as.integer(counts$n[counts$significance == "Downregulated"]),
            sum(volcano$FDR < t & volcano$logFC < 0), info = name)
    }
    top <- utils::read.delim(file.path(run$fixture$paths$da_output_dir, "top_de_genes.tsv"), check.names = FALSE)
    expect_true(all(top$FDR < t))
    long <- run$results$da_genes_long
    for (comparison in unique(long$comparison)) {
        rows <- long[long$comparison == comparison, , drop = FALSE]
        expect_identical(nrow(top[top$comparison == comparison, , drop = FALSE]),
            as.integer(min(10L, sum(rows$FDR < t & rows$logFC > 0)) + min(10L, sum(rows$FDR < t & rows$logFC < 0))),
            info = comparison)
    }
})

test_that("KIT-220 T-T1 (summary line): the RUV summary's count line states and uses the run's threshold", {
    best <- list(percentage = 5, n_controls = 10L, best_k = 1L, separation_score = 0.5, composite_score = 0.5)
    lines <- trxSummaryText(best, 100L, "~ 0 + group + W_1", "B_vs_A", data.frame(FDR = c(0.005, 0.02, 0.2)),
        data.frame(), fdr_thresh = 0.01)$lines
    expect_true("Differentially Expressed Genes (FDR < 0.01): 1" %in% lines)
    default <- trxSummaryText(best, 100L, "~ 0 + group + W_1", "B_vs_A", data.frame(FDR = c(0.005, 0.02, 0.2)),
        data.frame())$lines
    expect_true("Differentially Expressed Genes (FDR < 0.05): 2" %in% default)
})

test_that("KIT-220 T-T2 and T-PAR (transcriptomics): the interactive status is edgeR's at the run's FDR; tab = published", {
    skipWithoutTrxDeclReference()
    t <- 0.01
    run <- kit220TrxPublished(t)
    results <- run$results
    entries <- results$contrasts_results[[1L]]$results
    entry <- names(entries)[[1L]]
    qlf <- entries[[entry]]$qlf
    captured <- NULL
    testthat::local_mocked_bindings(glimmaVolcano = function(x, ...) {
        captured <<- list(...)
        "kit220-page"
    }, .package = "Glimma")
    assay <- names(results$contrasts_results)[[1L]]
    for (threshold in c(t, 0.05)) {
        captured <- NULL
        trxVolcanoWidget(results, entry, assay, da_q_val_thresh = threshold)
        expect_identical(captured$status, edgeR::decideTests.DGEExact(qlf, adjust.method = "BH", p.value = threshold),
            info = paste("threshold", threshold))
    }
    # At 0.05 the status is Glimma's documented default object.
    expect_identical(captured$status, edgeR::decideTests.DGEExact(qlf))
    # T-PAR: the tab's counts (the stored significant column, made at the run's threshold) equal the published bar's.
    long <- results$da_genes_long
    # The tab names each row by its comparison (completeOmicDaRun, D2).
    if (is.null(long$friendly_name)) long$friendly_name <- long$comparison
    for (name in grep("^summary_barplot_", ls(run$figures), value = TRUE)) {
        comparison <- sub("^summary_barplot_", "", name)
        rows <- long[long$comparison == comparison, , drop = FALSE]
        counts <- get(name, envir = run$figures)$data
        text <- omicDaSummaryText(long, comparison, "All", t, omicGuiConfig("transcriptomics")$da$summary)
        expect_match(text, sprintf("Up-regulated: %d", as.integer(counts$n[counts$significance == "Upregulated"])),
            fixed = TRUE, info = name)
        expect_match(text, sprintf("Down-regulated: %d", as.integer(counts$n[counts$significance == "Downregulated"])),
            fixed = TRUE, info = name)
    }
})

# ---- T-PAR (lipidomics, metabolomics): already at the run's threshold -------------------------------------------

kit220FeatureTable <- function(t, lfc = 0) {
    long <- data.frame(feature = paste0("F", 1:8), logFC = c(1.2, -0.9, 0.6, -0.4, 2, -2, 0.3, -0.3),
        raw_pvalue = c(0.001, 0.002, 0.004, 0.006, 0.0001, 0.0002, 0.2, 0.3),
        fdr_qvalue = c(0.005, 0.008, 0.02, 0.03, 0.001, 0.002, 0.4, 0.5),
        comparison = "groupB-groupA", friendly_name = "groupB-groupA", numerator = "groupB", denominator = "groupA",
        assay = "LCMS_Pos", intensity.S1.groupA = 10, intensity.S2.groupB = 12, stringsAsFactors = FALSE)
    long$significant <- ifelse(long$fdr_qvalue < t & abs(long$logFC) >= lfc, ifelse(long$logFC > 0, "Up", "Down"), "NS")
    long
}

test_that("KIT-220 T-PAR (lipidomics, metabolomics): at 0.01 the tab's counts equal the published table's", {
    t <- 0.01
    testthat::local_mocked_bindings(
        generateLipidDAVolcanoStatic = function(...) ggplot2::ggplot(data.frame(x = 1, y = 1), ggplot2::aes(x, y)) + ggplot2::geom_point(),
        generateLipidDAHeatmap = function(...) NULL,
        generateMetabDAVolcanoStatic = function(...) ggplot2::ggplot(data.frame(x = 1, y = 1), ggplot2::aes(x, y)) + ggplot2::geom_point(),
        generateMetabDAHeatmap = function(...) NULL, .package = "MultiScholaR")
    cases <- list(
        lipidomics = list(fun = outputLipidDaResultsAllContrasts, long = "da_lipids_long", id = "lipid_id",
            file = "lipids_num_sig_de_molecules.tab"),
        metabolomics = list(fun = outputMetabDaResultsAllContrasts, long = "da_metabolites_long", id = "metabolite_id",
            file = "metabolites_num_sig_de_molecules.tab"))
    for (omic in names(cases)) {
        case <- cases[[omic]]
        long <- kit220FeatureTable(t)
        names(long)[names(long) == "feature"] <- case$id
        root <- withr::local_tempdir()
        utils::capture.output(suppressWarnings(suppressMessages(case$fun(stats::setNames(list(long), case$long),
            da_output_dir = file.path(root, "da"), publication_graphs_dir = file.path(root, "graphs"),
            da_q_val_thresh = t, lfc_threshold = 0))))
        table <- utils::read.delim(file.path(root, "graphs", "NumSigDeMolecules", case$file), check.names = FALSE)
        text <- omicDaSummaryText(long, "groupB-groupA", "All", t, omicGuiConfig(omic)$da$summary)
        expect_match(text, sprintf("Up-regulated: %d", sum(table$up_regulated)), fixed = TRUE, info = omic)
        expect_match(text, sprintf("Down-regulated: %d", sum(table$down_regulated)), fixed = TRUE, info = omic)
        expect_equal(unique(table$q_threshold), t, info = omic)
    }
})

# ---- E-7: the tab states the run's threshold and says when the box differs ----------------------------------------

test_that("KIT-220 E-7: the feature omics' summary states the run's threshold and notes a box that differs from it", {
    long <- kit220FeatureTable(0.05)
    summary <- omicGuiConfig("lipidomics")$da$summary
    # The box at 0.01, the run at 0.05: the counts are the run's, so the label is the run's, and the note says so.
    text <- omicDaSummaryText(long, "groupB-groupA", "All", 0.01, summary, runQValThresh = 0.05)
    expect_match(text, sprintf("Significant (%s < 0.050): 6", summary$significance_label), fixed = TRUE)
    expect_match(text, "differs from the threshold the analysis ran at", fixed = TRUE)
    expect_match(text, "0.010", fixed = TRUE)
    # The box equal to the run: today's text, no note.
    same <- omicDaSummaryText(long, "groupB-groupA", "All", 0.05, summary, runQValThresh = 0.05)
    expect_identical(same, omicDaSummaryText(long, "groupB-groupA", "All", 0.05, summary))
    expect_false(grepl("differs", same, fixed = TRUE))
    # The module gives the run's threshold to the summary.
    body <- paste(deparse(body(registerOmicDaTableOutputs)), collapse = "\n")
    expect_match(body, "runQValThresh = results$da_q_val_thresh", fixed = TRUE)
})

test_that("KIT-220 E-7 (proteomics, SD-6 (a)): the summary notes a fold-change box that differs from the run's cutoff", {
    testthat::local_mocked_bindings(isDIAResultManifest = function(x) TRUE,
        readDIAResultPage = function(...) list(counts = list(total = 10L, significant = 3L, up = 2L, down = 1L)),
        .package = "MultiScholaR")
    manifest <- list(parameters = list(qvalue = 0.05, lfc = 0))
    differs <- diaDaSummaryStatsView(manifest, "B_vs_A", 0.05, 1)
    expect_match(differs, "differs from the cutoff the analysis ran at", fixed = TRUE)
    same <- diaDaSummaryStatsView(manifest, "B_vs_A", 0.05, 0)
    expect_false(grepl("differs", same, fixed = TRUE))
})

test_that("KIT-220 SD-6 (a): the interactive page is built and found with the run's fold-change cutoff", {
    captured <- NULL
    testthat::local_mocked_bindings(generateProtDAVolcanoPlotGlimma = function(...) {
        captured <<- list(...)
        NULL
    }, .package = "MultiScholaR")
    diaDaVolcanoGlimmaView(list(parameters = list(qvalue = 0.05, lfc = 0.5)), "B_vs_A", 0.05, kit220ProtObject(proteins = 20L),
        tempdir(), 1)
    expect_identical(captured$lfc_thresh, 0.5)
    body <- paste(deparse(body(diaInteractiveVolcanoIdentity)), collapse = "\n")
    expect_match(body, "lfc = ", fixed = TRUE)
})

# ---- G5: every app and workbook count call names its threshold ---------------------------------------------------

test_that("KIT-220 G5: every printCountDaGenesTable call of the DA helper and the workbook route names q_val_thresh", {
    for (file in c("func_prot_da_model_methods.R", "func_prot_da_compat.R")) {
        calls <- kit220CallsTo(kit220PackageFile("R", file), "printCountDaGenesTable")
        expect_gt(length(calls), 0L)
        for (call in calls) expect_true("q_val_thresh" %in% names(as.list(call)), info = paste(file, deparse(call)[[1L]]))
    }
    trx <- kit220PackageFile("R", "utils_workflow_trx_declarations.R")
    for (fun in c("trxVolcano", "trxSignificanceBar", "trxAnnotateDe", "trxSummaryText")) {
        calls <- kit220CallsTo(trx, fun)
        expect_gt(length(calls), 0L)
        # trxDaFormat's annotation (the tab's long table) uses no top-genes list and keeps the default.
        if (fun == "trxAnnotateDe") calls <- Filter(function(call) !identical(call[[3L]], as.name("annotation_info")), calls)
        for (call in calls) expect_true("fdr_thresh" %in% names(as.list(call)), info = paste(fun, deparse(call)[[1L]]))
    }
})

# ---- r1.3 (implementation audit I05): the top-genes list at the given threshold -----------------------------------

test_that("KIT-220 T-T1 (top genes): trxAnnotateDe lists the top genes at the given FDR threshold", {
    # trxAnnotateDe takes the gene ids from the top table's row names (rownames_to_column, R/func_trx_da.R:102).
    top <- data.frame(logFC = c(2, 1.5, 1, -2, -1.5, -1), FDR = c(0.001, 0.02, 0.04, 0.002, 0.03, 0.2),
        row.names = paste0("G", 1:6))
    counts <- data.frame(GeneID = paste0("G", 1:6), GeneName = paste0("N", 1:6), stringsAsFactors = FALSE)
    at_001 <- trxAnnotateDe(top, counts, "B_vs_A", fdr_thresh = 0.01)$top_de
    expect_identical(sort(at_001$GeneID), c("G1", "G4"))
    expect_identical(sort(trxAnnotateDe(top, counts, "B_vs_A")$top_de$GeneID), c("G1", "G2", "G3", "G4", "G5"))
})
