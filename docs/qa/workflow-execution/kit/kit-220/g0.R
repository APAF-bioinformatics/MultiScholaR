# KIT-220 gate G0 (ticket 9 A1, 10 G0; stop rule 2). Evidence code under docs/qa: no package code. Runs on a
# GitHub-hosted runner (run-tests `scripts`) from the package root after the commit under test is installed: on the
# tests-only commit (base code) for G0, and on the fix commit for the "after" values (A1, G2).
#
# Usage: Rscript docs/qa/workflow-execution/kit/kit-220/g0.R <out dir>
#
# (1) versions.tsv: R, edgeR, Glimma, limma, qvalue; whether edgeR has decideTestsDGE and decideTests.DGEExact (E-3),
#     and Glimma's glimmaVolcano.DGELRT formals (printed to the console).
# (2) g0.tsv, proteomics rows: the new test file's synthetic fixture (its kit220* functions, evaluated alone) at
#     t = 0.01 and 0.10 (lfc 0) and t = 0.01 with TREAT at 0.5: the rows between 0.01 and 0.05 and between 0.05 and
#     0.10, the long table's counts at t, the first and second count tables, the tab's counts, the published table's
#     counts and q_threshold, the DIA-NN report sentence.
# (3) g0.tsv, transcriptomics rows: trx_edge_12 (KIT-192 T6's path) at t = 0.01 and 0.10: per contrast the rows
#     between the thresholds, the tab's counts, the published volcano's dashed line and coloured counts, the bar's
#     counts, top_de_genes.tsv's rows and largest FDR.
# (4) w1.tsv: the KIT-112 W1 fixture (recipe "small") at q 0.01 with TREAT at 0.5, per contrast: rows with
#     0.01 <= q < 0.05 (section 6: whether dia_qvalue_only's leaves move under SD-1 (A)).
# (5) feature.tsv: lipidomics and metabolomics published count tables at 0.01 on the test file's feature table: md5
#     and counts (G2: base against head).
# Each section is wrapped: a failure writes its message to errors.tsv and the script continues. Exit status 0.

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1L)
out_dir <- args[[1L]]
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
writeTsv <- function(x, name) utils::write.table(x, file.path(out_dir, name), sep = "\t", quote = FALSE,
    row.names = FALSE, na = "NA")
Sys.setenv(http_proxy = "http://127.0.0.1:9", https_proxy = "http://127.0.0.1:9", HTTP_PROXY = "http://127.0.0.1:9",
    HTTPS_PROXY = "http://127.0.0.1:9", no_proxy = "", NO_PROXY = "", MULTISCHOLAR_TEST_MODE = "true",
    MULTISCHOLAR_TRX_MANIFEST_NOTE = "1")
suppressPackageStartupMessages(library(MultiScholaR))
library(testthat)
ns <- asNamespace("MultiScholaR")
h <- new.env(parent = ns)
# Every helper, as testthat loads them (instrument fix, run 45498530109: the trx owner needs omicsOwnerSession).
for (helper in sort(list.files(file.path("tests", "testthat"), pattern = "^helper-.*[.]R$", full.names = TRUE)))
    tryCatch(sys.source(helper, envir = h), error = function(e) message(sprintf("helper %s not loaded: %s", helper,
        conditionMessage(e))))
# The test file's kit220* functions only (its test_that blocks are not run).
for (e in parse(file.path("tests", "testthat", "test-kit220-da-threshold-consumers.R"), keep.source = FALSE)) {
    if (is.call(e) && identical(e[[1L]], as.name("<-")) && startsWith(as.character(e[[2L]]), "kit220")) eval(e, h)
}
version_of <- function(p) tryCatch(as.character(utils::packageVersion(p)), error = function(e) "not installed")
errors <- data.frame(section = character(), message = character(), stringsAsFactors = FALSE)
section <- function(name, code) tryCatch(code, error = function(e) {
    message(sprintf("SECTION %s FAILED: %s", name, conditionMessage(e)))
    errors[nrow(errors) + 1L, ] <<- c(name, substr(conditionMessage(e), 1L, 500L))
    NULL
})
rows <- list()
add <- function(...) rows[[length(rows) + 1L]] <<- data.frame(..., stringsAsFactors = FALSE)

# ---- (1) versions -----------------------------------------------------------------------------------------------
section("versions", {
    edgeR_ns <- asNamespace("edgeR")
    versions <- data.frame(item = c("R", "edgeR", "Glimma", "limma", "qvalue", "edgeR::decideTestsDGE exists",
        "edgeR exports decideTests.DGEExact", "edgeR exports decideTests.DGELRT"),
        value = c(R.version.string, vapply(c("edgeR", "Glimma", "limma", "qvalue"), version_of, ""),
            exists("decideTestsDGE", envir = edgeR_ns, inherits = FALSE),
            "decideTests.DGEExact" %in% getNamespaceExports("edgeR"),
            "decideTests.DGELRT" %in% getNamespaceExports("edgeR")), stringsAsFactors = FALSE)
    writeTsv(versions, "versions.tsv")
    print(versions)
    print(args(utils::getFromNamespace("glimmaVolcano.DGELRT", "Glimma")))
    print(args(edgeR::decideTests.DGEExact))
})

# ---- (2) proteomics ---------------------------------------------------------------------------------------------
section("proteomics", local({
    object <- h$kit220ProtObject()
    for (setting in list(c(0.01, 0), c(0.10, 0), c(0.01, 0.5))) {
        t <- setting[[1L]]
        lfc <- setting[[2L]]
        result <- h$kit220ProtDa(object, t, lfc)
        long <- result$da_proteins_long
        sign <- h$kit220SignCounts(long, t)
        first <- h$kit220TableCounts(result$num_sig_da_molecules_first_go$table)
        second <- h$kit220TableCounts(result$num_sig_da_molecules)
        tab <- ns$diaDaSummaryStatsView(list(da_proteins_long = long), unique(long$comparison)[[1L]], t, lfc)
        path <- h$kit220Publish(object, result, t, tempfile("kit220-g0-"))
        table <- utils::read.delim(path, check.names = FALSE)
        published <- h$kit220TableCounts(table)
        sentence <- h$kit220FormatDESummaryTable("DIANN_report.rmd")(path)
        add(omic = "proteomics", fixture = "kit220ProtObject", contrast = unique(long$comparison)[[1L]], t = t, lfc = lfc,
            rows = nrow(long), between_001_005 = h$kit220Between(long, 0.01, 0.05),
            between_005_010 = h$kit220Between(long, 0.05, 0.10),
            below_t_lfc_inside = sum(long$fdr_qvalue < t & abs(long$log2FC) <= lfc, na.rm = TRUE),
            long_up = sign[["up"]], long_down = sign[["down"]], first_up = first[["up"]], first_down = first[["down"]],
            second_up = second[["up"]], second_down = second[["down"]],
            tab = gsub("\n", " | ", tab), published_up = published[["up"]], published_down = published[["down"]],
            published_q_threshold = if ("q_threshold" %in% names(table)) paste(unique(table$q_threshold), collapse = ",") else "absent",
            published_md5 = unname(tools::md5sum(path)), extra = sentence)
    }
}))

# ---- (3) transcriptomics ----------------------------------------------------------------------------------------
section("transcriptomics", local({
    for (t in c(0.01, 0.10)) local({
        frame <- environment()
        run <- h$kit220TrxPublished(t, env = frame)
        long <- run$results$da_genes_long
        top_path <- file.path(run$fixture$paths$da_output_dir, "top_de_genes.tsv")
        top <- tryCatch(utils::read.delim(top_path, check.names = FALSE), error = function(e) data.frame(FDR = numeric(),
            comparison = character()))
        for (name in grep("^volcano_plot_", ls(run$figures), value = TRUE)) {
            comparison <- sub("^volcano_plot_", "", name)
            plot <- get(name, envir = run$figures)
            hline <- which(vapply(plot$layers, function(layer) inherits(layer$geom, "GeomHline"), logical(1)))
            line <- ggplot2::layer_data(plot, hline[[1L]])$yintercept
            classes <- table(factor(as.character(plot$data$significance), c("Upregulated", "Downregulated", "Not Significant")))
            bar <- get(sub("^volcano_plot_", "summary_barplot_", name), envir = run$figures)$data
            rows_c <- long[long$comparison == comparison, , drop = FALSE]
            top_c <- top[top$comparison == comparison, , drop = FALSE]
            add(omic = "transcriptomics", fixture = "trx_edge_12", contrast = comparison, t = t, lfc = 0,
                rows = nrow(rows_c), between_001_005 = sum(rows_c$FDR >= 0.01 & rows_c$FDR < 0.05),
                between_005_010 = sum(rows_c$FDR >= 0.05 & rows_c$FDR < 0.10), below_t_lfc_inside = NA,
                long_up = sum(rows_c$FDR < t & rows_c$logFC > 0), long_down = sum(rows_c$FDR < t & rows_c$logFC < 0),
                first_up = NA, first_down = NA, second_up = NA, second_down = NA,
                tab = sprintf("stored significant: Up %d, Down %d", sum(rows_c$significant == "Up"), sum(rows_c$significant == "Down")),
                published_up = as.integer(bar$n[bar$significance == "Upregulated"]),
                published_down = as.integer(bar$n[bar$significance == "Downregulated"]),
                published_q_threshold = sprintf("volcano line 10^-%.4f; coloured Up %d Down %d", line, classes[[1L]], classes[[2L]]),
                published_md5 = unname(tools::md5sum(top_path)),
                extra = sprintf("top_de_genes rows %d, largest FDR %s", nrow(top_c),
                    if (nrow(top_c)) format(max(top_c$FDR), digits = 4) else "NA"))
        }
    })
}))

# ---- (4) the KIT-112 W1 fixture ---------------------------------------------------------------------------------
section("w1", local({
    object <- h$makeWorkflowCheckpointProtein(h$readWorkflowCheckpointRecipe("small", repo_root = getwd()), repo_root = getwd())
    w1 <- list()
    for (contrast in c("treated.vs.control=grouptreated-groupcontrol", "control.vs.treated=groupcontrol-grouptreated")) {
        result <- NULL
        utils::capture.output(result <- suppressWarnings(suppressMessages(differentialAbundanceAnalysis(object,
            contrasts_tbl = data.frame(contrasts = contrast, stringsAsFactors = FALSE), formula_string = "~ 0 + group",
            da_q_val_thresh = 0.01, treat_lfc_cutoff = 0.5, eBayes_trend = TRUE, eBayes_robust = TRUE,
            args_row_id = object@protein_id_column, qvalue_column = "fdr_qvalue", raw_pvalue_column = "raw_pvalue"))))
        long <- result$da_proteins_long
        w1[[length(w1) + 1L]] <- data.frame(case = "dia_qvalue_only", contrast = contrast, rows = nrow(long),
            between_001_005 = h$kit220Between(long, 0.01, 0.05), below_001 = sum(long$fdr_qvalue < 0.01, na.rm = TRUE),
            stringsAsFactors = FALSE)
    }
    writeTsv(do.call(rbind, w1), "w1.tsv")
    print(do.call(rbind, w1))
}))

# ---- (5) lipidomics and metabolomics published tables (G2) --------------------------------------------------------
section("feature", local({
    frame <- environment()
    testthat::local_mocked_bindings(
        generateLipidDAVolcanoStatic = function(...) ggplot2::ggplot(data.frame(x = 1, y = 1), ggplot2::aes(x, y)) + ggplot2::geom_point(),
        generateLipidDAHeatmap = function(...) NULL,
        generateMetabDAVolcanoStatic = function(...) ggplot2::ggplot(data.frame(x = 1, y = 1), ggplot2::aes(x, y)) + ggplot2::geom_point(),
        generateMetabDAHeatmap = function(...) NULL, .package = "MultiScholaR", .env = frame)
    cases <- list(lipidomics = list(fun = outputLipidDaResultsAllContrasts, long = "da_lipids_long", id = "lipid_id",
            file = "lipids_num_sig_de_molecules.tab"),
        metabolomics = list(fun = outputMetabDaResultsAllContrasts, long = "da_metabolites_long", id = "metabolite_id",
            file = "metabolites_num_sig_de_molecules.tab"))
    feature <- list()
    for (omic in names(cases)) {
        case <- cases[[omic]]
        long <- h$kit220FeatureTable(0.01)
        names(long)[names(long) == "feature"] <- case$id
        root <- tempfile("kit220-g0-feature-")
        utils::capture.output(suppressWarnings(suppressMessages(case$fun(stats::setNames(list(long), case$long),
            da_output_dir = file.path(root, "da"), publication_graphs_dir = file.path(root, "graphs"),
            da_q_val_thresh = 0.01, lfc_threshold = 0))))
        path <- file.path(root, "graphs", "NumSigDeMolecules", case$file)
        table <- utils::read.delim(path, check.names = FALSE)
        feature[[omic]] <- data.frame(omic = omic, md5 = unname(tools::md5sum(path)), up = sum(table$up_regulated),
            down = sum(table$down_regulated), q_threshold = paste(unique(table$q_threshold), collapse = ","),
            stringsAsFactors = FALSE)
    }
    writeTsv(do.call(rbind, feature), "feature.tsv")
    print(do.call(rbind, feature))
}))

if (length(rows)) {
    g0 <- do.call(rbind, rows)
    writeTsv(g0, "g0.tsv")
    print(g0[, setdiff(names(g0), c("tab", "extra"))])
}
writeTsv(errors, "errors.tsv")
print(errors)
