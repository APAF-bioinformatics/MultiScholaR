# KIT-204 (VAL-129 (1), VAL-116 (1)): a proteomics DA row whose raw p-value is missing was not tested. The tab's summary
# counts it as "Not tested (no value in a compared group): n" (always shown, n may be 0) and the count tables give it the
# status "Not tested" (only when n > 0) instead of "Not significant". Calls to functions whose formals this ticket
# extends pass only the formals the function has, so the tests-only commit fails at a value, not "unused argument".
# Synthetic data only; builders are copied from test-kit200-prot-da-no-value-proteins.R and test-workflow-saved-da.R.

k204Statuses <- c("Not significant", "Significant and Up", "Significant and Down")
k204Significant <- c("Significant and Up", "Significant and Down")
k204NotTestedLine <- "Not tested (no value in a compared group): %d"

k204Quiet <- function(expr) {
    value <- NULL
    utils::capture.output(value <- suppressWarnings(suppressMessages(expr)))
    value
}

k204Call <- function(f, args) do.call(f, args[intersect(names(args), names(formals(f)))])

k204WithSeed <- function(seed, code) withr::with_seed(seed, code, .rng_kind = "Mersenne-Twister",
    .rng_normal_kind = "Inversion", .rng_sample_kind = "Rejection")

# kit200Matrix (KIT-200 section 4), copied: P1..P60, levels N(22, 2^2), residual sd 0.3, P11..P16 with a group-B effect.
k204Matrix <- function(sizes = c(A = 4L, B = 4L, C = 4L), seed = 200001L, proteins = 60L) {
    groups <- rep(names(sizes), sizes)
    x <- k204WithSeed(seed, {
        level <- stats::rnorm(proteins, 22, 2)
        effect <- as.numeric(seq_len(proteins) %in% 11:16)
        level + outer(effect, as.numeric(groups == "B")) +
            matrix(stats::rnorm(proteins * length(groups), sd = 0.3), proteins)
    })
    dimnames(x) <- list(paste0("P", seq_len(proteins)), paste0("S", seq_along(groups)))
    x
}

k204Design <- function(sizes = c(A = 4L, B = 4L, C = 4L)) {
    runs <- paste0("S", seq_len(sum(sizes)))
    data.frame(Run = runs, group = rep(names(sizes), sizes), replicates = runs, row.names = runs,
        stringsAsFactors = FALSE)
}

# The fixture: `P5` without a value in group C (untested in the C contrasts, tested in B-A); `empty` adds `P_empty`
# (no value in any of the 12 samples, KIT-200's untested protein); `complete` leaves every value.
k204Values <- function(complete = FALSE, empty = FALSE) {
    x <- k204Matrix()
    if (!complete) x["P5", paste0("S", 9:12)] <- NA_real_
    if (empty) x <- rbind(x, P_empty = rep(NA_real_, ncol(x)))
    x
}

k204Object <- function(x, design = k204Design()) {
    rownames(design) <- NULL
    ProteinQuantitativeData(protein_quant_table = data.frame(Protein.Ids = rownames(x), x, row.names = NULL,
        check.names = FALSE, stringsAsFactors = FALSE), design_matrix = design, sample_id = "Run", group_id = "group",
        technical_replicate_id = "replicates", protein_id_column = "Protein.Ids",
        protein_id_table = data.frame(Protein.Ids = rownames(x), gene_names = paste0("GENE_", rownames(x)),
            stringsAsFactors = FALSE), args = list())
}

# The QC figures the helper draws before the fit are replaced (they are not counted, stored or published).
k204MockFigures <- function(env = parent.frame()) {
    testthat::local_mocked_bindings(plotRle = function(...) ggplot2::ggplot(), plotPca = function(...) ggplot2::ggplot(),
        plotNumOfValuesNoLog = function(...) ggplot2::ggplot(), .package = "MultiScholaR", .env = env)
}

# One contrast per call, as the DA kernel calls the method (R/utils_workflow_dia_da.R:308, :316).
k204Da <- function(object, contrast, qvalue = 0.05) {
    k204Quiet(differentialAbundanceAnalysis(object,
        contrasts_tbl = data.frame(contrasts = contrast, stringsAsFactors = FALSE), formula_string = "~ 0 + group",
        da_q_val_thresh = qvalue, treat_lfc_cutoff = 0, eBayes_trend = TRUE, eBayes_robust = TRUE,
        args_row_id = "Protein.Ids", qvalue_column = "fdr_qvalue", raw_pvalue_column = "raw_pvalue"))
}

# A status's count in a count table (NA when the status is absent), as a double.
k204Count <- function(table, status) {
    at <- which(as.character(table$status) == status)
    if (!length(at)) NA_real_ else as.numeric(sum(table$counts[at]))
}

# The drawn statuses of a barplot: the rows of its data, every one drawn by its one layer.
k204Drawn <- function(plot) {
    expect_identical(nrow(ggplot2::layer_data(plot)), nrow(plot$data))
    sort(unique(as.character(plot$data$status)))
}

# --- The owner (savedDaW1Owner's shape, copied) and the tab's DA click in process ------------------------------------

k204Contrasts <- function(names = c("B.vs.A", "C.vs.A")) {
    formats <- c(B.vs.A = "groupB-groupA", C.vs.A = "groupC-groupA", C.vs.B = "groupC-groupB")
    data.frame(contrasts = names, full_format = paste0(names, "=", formats[names]), friendly_names = names,
        stringsAsFactors = FALSE)
}

k204Owner <- function(object, envir = parent.frame()) {
    root <- normalizePath(withr::local_tempdir(.local_envir = envir))
    paths <- list(source_dir = file.path(root, "source"), da_output_dir = file.path(root, "da"),
        publication_graphs_dir = file.path(root, "graphs"))
    for (path in paths) dir.create(path)
    manager <- newDIAStateAdapter(newWorkflowEntityStore(root, "kit204", "DIA"), WorkflowState$new(audit_enabled = FALSE))
    manager$.context$reference_states <- TRUE
    manager$setWorkflowType("DIA")
    manager$.context$experiment_paths <- paths
    manager$saveState("correlation_filtered", object, list(), "fixture")
    workflow <- list2env(list(state_manager = manager, tab_status = list()))
    withr::defer(closeDIAStateAdapter(workflow$state_manager), envir = envir)
    da <- new.env(parent = emptyenv())
    da$current_s4_object <- getDIAStateReference(manager)
    list(root = root, paths = paths, workflow = workflow, da = da, manager = manager)
}

k204Click <- function(fixture, contrasts = k204Contrasts()) {
    fixture$workflow$contrasts_tbl <- contrasts
    fixture$da$contrasts_available <- contrasts$contrasts
    settings <- list(formula_string = "~ 0 + group", da_q_val_thresh = 0.05, treat_lfc_cutoff = 0)
    k204Quiet(runDIAReferenceDaAnalysis(settings, fixture$da, fixture$workflow, fixture$paths, contrasts, NULL,
        publish = FALSE))
}

k204Summary <- function(fixture, manifest, contrast) {
    withDIAStateScope(fixture$manager, function() diaDaSummaryStatsView(manifest, contrast, 0.05, 0))
}

k204Long <- function(fixture, manifest, contrast) {
    withDIAStateScope(fixture$manager, function()
        as.data.frame(readDIAContrastOutput(manifest$results[[contrast]], "da_proteins_long")$da_proteins_long))
}

# The summary's lines with the "Not tested" line where the ticket places it (after "Down-regulated", before
# "Fold-change cutoff"); `n` NULL means the line must be absent.
k204ExpectSummaryLine <- function(text, n, label) {
    lines <- strsplit(text, "\n", fixed = TRUE)[[1L]]
    expect_true(startsWith(lines[[4L]], "Down-regulated: "), info = label)
    if (is.null(n)) {
        expect_false(any(startsWith(lines, "Not tested")), info = label)
        expect_true(startsWith(lines[[5L]], "Fold-change cutoff: "), info = label)
    } else {
        expect_identical(lines[5L], sprintf(k204NotTestedLine, as.integer(n)), info = label)
        expect_true(isTRUE(startsWith(lines[6L], "Fold-change cutoff: ")), info = label)
    }
}

# --- N1, N4 (a): countStatDaGenes ----------------------------------------------------------------------------------

k204N1Table <- function() data.frame(log2FC = c(1.5, 0.8, -1.2, 0.1, NA, NA),
    fdr_qvalue = c(0.001, 0.01, 0.02, 0.6, NA, NA), raw_pvalue = c(0.0001, 0.002, 0.004, 0.5, NA, NA))

k204CountStat <- function(table) k204Call(countStatDaGenes, list(data = table, lfc_thresh = 0, q_val_thresh = 0.05,
    log_fc_column = quote(log2FC), q_value_column = quote(fdr_qvalue), p_value_column = quote(raw_pvalue)))

test_that("KIT-204 N1: countStatDaGenes counts rows without a p-value as Not tested", {
    counts <- k204CountStat(k204N1Table())
    expect_identical(k204Count(counts, "Not significant"), 1)
    expect_identical(k204Count(counts, "Not tested"), 2)
    expect_identical(as.character(counts$status), c(k204Statuses, "Not tested"))
    expect_identical(counts$counts, c(1L, 2L, 1L, 2L))
})

test_that("KIT-204 N4 (a) (guard): without an untested row the table and its types are today's", {
    tested <- k204N1Table()[1:4, ]
    counts <- k204CountStat(tested)
    expect_identical(as.character(counts$status), k204Statuses)
    expect_identical(counts$counts, c(1L, 2L, 1L))
    no_down <- k204CountStat(tested[-3L, ])
    expect_identical(as.character(no_down$status), k204Statuses)
    expect_identical(no_down$counts, c(1, 2, 0))
})

# --- N2, N2b, N3, N4 (b): the S4 method's two count tables and barplots ------------------------------------------

# For one contrast's result: the expected "Not significant" (finite p, not q < 0.05; a missing q is not significant, as
# the package's last case_when branch has it) and "Not tested" (p missing).
k204Expected <- function(result) {
    table <- result$contrasts_results_table
    tested <- !is.na(table$raw_pvalue)
    significant <- !is.na(table$fdr_qvalue) & table$fdr_qvalue < 0.05
    list(not_significant = as.numeric(sum(tested & !significant)),
        not_tested = as.numeric(sum(!tested)))
}

k204ExpectCounts <- function(result, untested, label) {
    expected <- k204Expected(result)
    expect_identical(expected$not_tested, untested, info = paste(label, "precondition: rows without a p-value"))
    for (name in c("num_sig_da_molecules_first_go", "num_sig_da_molecules")) {
        table <- if (identical(name, "num_sig_da_molecules")) result[[name]] else result[[name]]$table
        info <- paste(label, name)
        expect_identical(k204Count(table, "Not significant"), expected$not_significant, info = info)
        expect_identical(k204Count(table, "Not tested"), if (untested > 0) untested else NA_real_, info = info)
    }
}

test_that("KIT-204 N2: the method's count tables count P5 as Not tested in C-A only", {
    k204MockFigures()
    object <- k204Object(k204Values())
    k204ExpectCounts(k204Da(object, "B.vs.A=groupB-groupA"), 0, "B-A")
    k204ExpectCounts(k204Da(object, "C.vs.A=groupC-groupA"), 1, "C-A")
})

test_that("KIT-204 N2b: KIT-200's protein without any value is counted as Not tested", {
    k204MockFigures()
    object <- k204Object(k204Values(empty = TRUE))
    k204ExpectCounts(k204Da(object, "B.vs.A=groupB-groupA"), 1, "B-A")
    k204ExpectCounts(k204Da(object, "C.vs.A=groupC-groupA"), 2, "C-A")
    # The summary of the tab (N5's route) on the same values.
    fixture <- k204Owner(k204Object(k204Values(empty = TRUE)))
    manifest <- k204Click(fixture)
    k204ExpectSummaryLine(k204Summary(fixture, manifest, "B.vs.A"), 1, "N2b summary B-A")
    k204ExpectSummaryLine(k204Summary(fixture, manifest, "C.vs.A"), 2, "N2b summary C-A")
})

test_that("KIT-204 N3 (i): C-B's every-status barplot draws Not tested; its significant-only plots do not", {
    k204MockFigures()
    result <- k204Da(k204Object(k204Values()), "C.vs.B=groupC-groupB")
    significant <- sum(result$num_sig_da_molecules$counts[result$num_sig_da_molecules$status %in% k204Significant])
    expect_gt(significant, 0) # precondition (an instrument defect if it fails): C-B has a significant protein
    expect_identical(k204Expected(result)$not_tested, 1) # precondition: P5 untested in C-B
    skip_if(significant == 0, "precondition failed: no significant protein in C-B")
    for (plot in list(result$num_sig_da_molecules_first_go$plot, result$num_sig_da_genes_barplot_only_significant)) {
        expect_true(all(k204Drawn(plot) %in% k204Significant))
    }
    expect_true("Not tested" %in% k204Drawn(result$num_sig_da_genes_barplot_with_not_significant))
})

test_that("KIT-204 N3 (ii) (guard): untested and nothing significant draws no barplot, as before", {
    k204MockFigures()
    result <- k204Da(k204Object(k204Values()), "C.vs.A=groupC-groupA", qvalue = 1e-12)
    statuses <- as.character(result$num_sig_da_molecules$status)
    expect_false(any(statuses %in% k204Significant)) # precondition: nothing significant at q < 1e-12
    expect_identical(k204Expected(result)$not_tested, 1) # precondition: P5 untested
    expect_length(grep("^num_sig_da_genes_barplot_", names(result), value = TRUE), 0L)
})

test_that("KIT-204 N4 (b) (guard): on complete values the method's two count tables are today's", {
    k204MockFigures()
    object <- k204Object(k204Values(complete = TRUE))
    for (contrast in c("B.vs.A=groupB-groupA", "C.vs.A=groupC-groupA")) {
        result <- k204Da(object, contrast)
        expect_identical(k204Expected(result)$not_tested, 0, info = contrast) # precondition
        nested <- list(result$contrasts_results_table)
        names(nested) <- names(result$contrasts_results$results)[1L]
        first_go <- k204Quiet(printCountDaGenesTable(list(nested), list("RUV applied")))
        expect_identical(result$num_sig_da_molecules_first_go$table, first_go$table, info = contrast)
        # Today's computation (R/func_prot_da_model_methods.R:597-606 at dd769671), copied.
        today <- k204Quiet(result$significant_rows |>
            dplyr::mutate(status = dplyr::case_when(
                fdr_qvalue >= 0.05 ~ "Not significant",
                log2FC >= 0 & fdr_qvalue < 0.05 ~ "Significant and Up",
                log2FC < 0 & fdr_qvalue < 0.05 ~ "Significant and Down",
                TRUE ~ "Not significant"
            )) |>
            dplyr::group_by(comparison, status) |>
            dplyr::summarise(counts = dplyr::n()) |>
            dplyr::ungroup())
        expect_identical(result$num_sig_da_molecules, today, info = contrast)
    }
})

# --- N5, N9, N3 (iii): the tab's summary and the publication on the owner -------------------------------------------

k204FormatDESummaryTable <- function() {
    lines <- readLines(system.file("reports", "proteomics", "DIANN_report.rmd", package = "MultiScholaR"), warn = FALSE,
        encoding = "UTF-8")
    defined <- grep("^FormatDESummaryTable <- function", lines)
    stopifnot(length(defined) == 1L)
    start <- max(grep("^```\\{r", lines[seq_len(defined)]))
    end <- defined + which(lines[-seq_len(defined)] == "```")[[1L]]
    chunk <- tempfile(fileext = ".Rmd")
    writeLines(lines[start:end], chunk, useBytes = TRUE)
    code <- tempfile(fileext = ".R")
    suppressMessages(knitr::purl(chunk, output = code, quiet = TRUE, documentation = 0))
    exprs <- parse(code, keep.source = FALSE, encoding = "UTF-8")
    at <- which(vapply(exprs, function(e) is.call(e) && identical(e[[1L]], as.name("<-")) &&
        identical(e[[2L]], as.name("FormatDESummaryTable")), logical(1)))
    env <- new.env(parent = asNamespace("MultiScholaR"))
    eval(exprs[[at]], env)
    env$FormatDESummaryTable
}

test_that("KIT-204 N5: the tab's summary always shows the Not tested line, from the store and from a list", {
    k204MockFigures()
    fixture <- k204Owner(k204Object(k204Values()))
    manifest <- k204Click(fixture)
    k204ExpectSummaryLine(k204Summary(fixture, manifest, "C.vs.A"), 1, "manifest C-A")
    k204ExpectSummaryLine(k204Summary(fixture, manifest, "B.vs.A"), 0, "manifest B-A")
    long <- k204Long(fixture, manifest, "C.vs.A")
    contrast <- unique(as.character(long$comparison))[[1L]]
    k204ExpectSummaryLine(diaDaSummaryStatsView(list(da_proteins_long = long), contrast, 0.05, 0), 1, "list C-A")
    without <- long[setdiff(names(long), "raw_pvalue")]
    k204ExpectSummaryLine(diaDaSummaryStatsView(list(da_proteins_long = without), contrast, 0.05, 0), NULL,
        "list without raw_pvalue")
})

test_that("KIT-204 N9 and N3 (iii): the published table has a Not tested row; the figure and report do not change", {
    k204MockFigures()
    fixture <- k204Owner(k204Object(k204Values()))
    # C-B carries the group-B effect, so the NumSig figure is drawn (N3 (iii)'s precondition).
    manifest <- k204Click(fixture, k204Contrasts(c("B.vs.A", "C.vs.A", "C.vs.B")))
    saved <- new.env(parent = emptyenv())
    saved$plots <- list()
    testthat::local_mocked_bindings(ggsave = function(filename, plot = ggplot2::last_plot(), ...) {
        saved$plots[[basename(filename)]] <- plot
        invisible(filename)
    }, .package = "ggplot2")
    out <- file.path(fixture$root, "publication")
    done <- k204Quiet(withDIAStateScope(fixture$manager, function() {
        object <- readDIAExecutionComponents(fixture$manager, manifest$source, fields = "design_matrix")
        # The arguments produceDIADaPublication records on the object (R/utils_workflow_dia_da.R:586-590).
        object@args$outputDaResultsAllContrasts <- list(uniprot_tbl = data.frame(),
            da_output_dir = file.path(out, "da"), publication_graphs_dir = file.path(out, "graphs"),
            file_prefix = "da_proteins", args_row_id = manifest$protein_id_column, gene_names_column = "gene_names",
            uniprot_id_column = "Entry", da_q_val_thresh = 0.05, fdr_column = "fdr_qvalue", log2fc_column = "log2FC")
        outputDaResultsAllContrasts(object, manifest$results, uniprot_tbl = data.frame(),
            da_output_dir = file.path(out, "da"), publication_graphs_dir = file.path(out, "graphs"),
            file_prefix = "da_proteins", args_row_id = manifest$protein_id_column, gene_names_column = "gene_names",
            uniprot_id_column = "Entry")
    }))
    expect_true(isTRUE(done))
    path <- file.path(out, "graphs", "NumSigDaMolecules", "da_proteins_num_sig_da_molecules.tab")
    expect_true(file.exists(path))
    table <- as.data.frame(readr::read_tsv(path, show_col_types = FALSE))
    at <- function(contrast) table[startsWith(as.character(table$comparison), contrast), , drop = FALSE]
    # N9: C-A's row (failing first); B-A has none.
    expect_identical(k204Count(at("C.vs.A"), "Not tested"), 1)
    expect_identical(k204Count(at("B.vs.A"), "Not tested"), NA_real_)
    # N3 (iii) (guard): the NumSig figure draws the two significant statuses only. Precondition: it was drawn (C-B has a
    # significant protein).
    figure <- saved$plots[["da_proteins_num_sig_da_molecules.png"]]
    expect_false(is.null(figure))
    if (!is.null(figure)) expect_true(all(k204Drawn(figure) %in% k204Significant))
    # N9 (guard): the report's sentence on this table equals its sentence on the table without the Not tested rows.
    format_summary <- k204FormatDESummaryTable()
    stripped <- tempfile(fileext = ".tab")
    readr::write_tsv(table[as.character(table$status) != "Not tested", , drop = FALSE], stripped)
    sentence <- format_summary(path)
    expect_false(startsWith(sentence, "Error"))
    expect_identical(sentence, format_summary(stripped))
})

# --- N6, N7: the complete fixture's summary and the display-selection identity ----------------------------------------

test_that("KIT-204 N6: on complete values the summary is today's five lines with the line at 0", {
    k204MockFigures()
    fixture <- k204Owner(k204Object(k204Values(complete = TRUE)))
    manifest <- k204Click(fixture, k204Contrasts("B.vs.A"))
    long <- k204Long(fixture, manifest, "B.vs.A")
    significant <- !is.na(long$fdr_qvalue) & long$fdr_qvalue < 0.05
    expected <- paste(sprintf("Total genes: %d", nrow(long)),
        sprintf("Significant (q < %.3f): %d", 0.05, sum(significant)),
        sprintf("Up-regulated: %d", sum(significant & long$log2FC > 0, na.rm = TRUE)),
        sprintf("Down-regulated: %d", sum(significant & long$log2FC < 0, na.rm = TRUE)),
        sprintf(k204NotTestedLine, 0L), sprintf("Fold-change cutoff: %.2f", 0), sep = "\n")
    expect_identical(k204Summary(fixture, manifest, "B.vs.A"), expected)
})

test_that("KIT-204 N7: a display selection saved under the v1 identity is not served", {
    k204MockFigures()
    fixture <- k204Owner(k204Object(k204Values()))
    manifest <- k204Click(fixture, k204Contrasts("C.vs.A"))
    manager <- fixture$manager
    counts <- withDIAStateScope(manager, function() {
        reference <- selectDIAContrastReference(manifest, "C.vs.A")$tables$da_proteins_long
        store <- manager$.context$store
        # A real selection at another cutoff lends its row index; the marker product carries total = -1.
        borrowed <- resolveDIAResultSelection(manager, reference, 0.5, 0)
        identity <- newSavedWorkflowScientificIdentity(reference$entity_id, "full-testing-population",
            digest::digest(list(0.05, 0), algo = "sha256"), "DA-display-selection-v1", manager$.context$software)
        seeded <- list(rows = borrowed$rows, counts = data.frame(total = -1, significant = 0, up = 0, down = 0))
        saveWorkflowConsumerEntity(store, newSavedWorkflowScientificProduct(seeded, "DA-display-selection"),
            basename(tempfile("kit204-seed-")), identity, "scientific_product")
        readDIAResultPage(manifest, "C.vs.A", limit = 1L)$counts
    })
    long <- k204Long(fixture, manifest, "C.vs.A")
    expect_identical(counts$total, as.double(nrow(long)))
    expect_identical(as.numeric(counts$not_tested), 1)
})

# --- N8: the lipid and metabolomics callers ------------------------------------------------------------------------

test_that("KIT-204 N8 (guard): printCountDaGenesTable without a p-value column keeps three statuses", {
    # R/func_lipid_s4_da_plot_methods.R:11 and R/func_metab_s4_da_results.R:80 pass no p-value column and
    # `formula_string = NA`; the table is given in the nested shape the shared functions name their comparisons from.
    table <- data.frame(logFC = c(1.5, 0.8, -1.2, 0.1, NA, NA), fdr_qvalue = c(0.001, 0.01, 0.02, 0.6, NA, NA),
        raw_pvalue = c(0.0001, 0.002, 0.004, 0.5, NA, NA))
    result <- k204Quiet(printCountDaGenesTable(list_of_da_tables = list(list(`groupB-groupA` = table)),
        list_of_descriptions = list("groupB-groupA"), formula_string = NA))
    expect_identical(sort(as.character(result$table$status)), sort(k204Statuses))
    expect_identical(k204Count(result$table, "Not significant"), 3)
    expect_identical(k204Count(result$table, "Significant and Up"), 2)
    expect_identical(k204Count(result$table, "Significant and Down"), 1)
    expect_true(all(k204Drawn(result$plot) %in% k204Significant))
})
