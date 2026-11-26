# KIT-236: the lipidomics and metabolomics reports (and part B, the four proteomics reports) print no finding the
# analysis did not compute, name the methods as they were run, and name a file or figure only when it is there (ticket
# section 5; the owner's decision OD-1, methods ledger VAL-245). Each lipid/metab case goes the real path: a
# study_parameters.txt in the writer's shape, a design_matrix.tab, optionally an analysis_record.json, the files the
# report looks for, then the template's own Main in an environment created as RenderReport creates it
# (new.env(parent = asNamespace("MultiScholaR")), the Project Setup chunk's paths, %||% and a d66_log stub). The
# proteomics cases run the template's AnalyzeSamples and its Format* functions. G1 and the W rows compare with the base
# templates (tests/testdata/kit236/base, verbatim at the branch point 29065062). The W rows record the sentences that
# wait for KIT-242 (the normalisation sentence D8 and the internal-standard sentence D19): unchanged on purpose.

k236Spec <- list(
    lipid = list(dir = "lipidomics", file = "lipidomics_report.rmd", omic = "lipidomics", prefix = "lipid", word = "lipids"),
    metab = list(dir = "metabolomics", file = "metabolomics_report.rmd", omic = "metabolomics", prefix = "metab",
        word = "metabolites"),
    LFQ = list(dir = "proteomics", file = "LFQ_report.rmd", omic = "proteomics", method = "FormatMethodDetailsLFQ"),
    TMT = list(dir = "proteomics", file = "TMT_report.rmd", omic = "proteomics", method = "FormatMethodDetailsTMT"),
    DIANN = list(dir = "proteomics", file = "DIANN_report.rmd", omic = "proteomics", method = "FormatMethodDetails"),
    DIANN_limpa = list(dir = "proteomics", file = "DIANN_limpa_report.rmd", omic = "proteomics",
        method = "FormatMethodDetails"))
k236Omics <- c("lipid", "metab")
k236ProtNames <- c("LFQ", "TMT", "DIANN", "DIANN_limpa")

#' A template file: the package's (head) or the stored base copy
k236Path <- function(name, side = "head") {
    spec <- k236Spec[[name]]
    if (identical(side, "head")) system.file("reports", spec$dir, spec$file, package = "MultiScholaR") else
        testthat::test_path("..", "testdata", "kit236", "base", spec$file)
}

#' The chunks of a template: label, header line and code lines
k236Chunks <- function(path) {
    lines <- readLines(path, warn = FALSE, encoding = "UTF-8")
    starts <- grep("^```\\{r", lines)
    lapply(starts, function(s) {
        e <- s + which(grepl("^```[[:space:]]*$", lines[-seq_len(s)]))[[1L]]
        list(label = trimws(sub("^```\\{r[[:space:],]*([^,}]*).*$", "\\1", lines[[s]])), header = lines[[s]],
            code = if (e > s + 1L) lines[(s + 1L):(e - 1L)] else character())
    })
}

#' The template's functions (every `name <- function` of every chunk) in one environment, with the report's paths
k236Env <- function(name, root, side = "head") {
    spec <- k236Spec[[name]]
    env <- new.env(parent = asNamespace("MultiScholaR"))
    for (chunk in k236Chunks(k236Path(name, side))) {
        for (e in parse(text = chunk$code, keep.source = FALSE)) {
            if (is.call(e) && identical(e[[1L]], as.name("<-")) && is.call(e[[3L]]) &&
                identical(e[[3L]][[1L]], as.name("function"))) eval(e, env)
        }
    }
    env$`%||%` <- function(a, b) if (!is.null(a)) a else b
    env$d66_log <- function(...) invisible(NULL)
    env$tibble <- tibble::tibble
    env$params <- list(omic_type = spec$omic, experiment_label = "kit236", workflow_name = "kit236")
    env$report_source_dir <- file.path(root, "scripts")
    env$report_results_summary_dir <- file.path(root, "summary")
    env$report_qc_figures_dir <- file.path(root, "summary", "QC_figures")
    env$report_publication_figures_dir <- file.path(root, "summary", "Publication_figures")
    env$report_publication_tables_dir <- file.path(root, "summary", "Publication_tables")
    for (d in c(env$report_source_dir, env$report_qc_figures_dir, env$report_publication_figures_dir,
        env$report_publication_tables_dir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
    env
}

k236Norm <- function(x) trimws(gsub("[[:space:]]+", " ", paste(as.character(x), collapse = " ")))
k236Has <- function(text, phrase) grepl(k236Norm(phrase), k236Norm(text), fixed = TRUE)
k236Touch <- function(path) {
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    file.create(path)
}

# --- Inputs: the writer's file shape (the GUI lanes' file, G0 run 45577606908), the design table, the record ---------

k236Lines <- function(method = "cyclicloess", offset = 0, log = TRUE, itsd = TRUE, itsd_block = character(),
    assays = c("LCMS_Pos", "LCMS_Neg"), de_extra = character(), ui_block = character(), ruv = "Not Applied") {
    c("Study Parameters", "================", "",
        "Basic Information:", "-----------------", "Workflow Name: kit236", "Description: KIT-236 test", "",
        "Version Information:", "--------------------", "Repository: MultiScholaR", "Branch: N/A", "Commit: abc1234", "",
        "Assay Information:", "-------------------", paste("* Number of Assays:", length(assays)),
        paste("* Assay Names:", paste(assays, collapse = ", ")), "",
        "ITSD Normalization:", "-------------------", paste("* Applied:", if (itsd) "Yes" else "No"),
        if (itsd) c("* Method Type: average_centered", "* Aggregation: median", "* ITSD Removed After Normalization: Yes"),
        itsd_block, "",
        "Log Transformation:", "------------------", paste("* Applied:", if (log) "Yes" else "No"),
        if (log) (if (offset == 0) "* Transform: log2" else paste("* Offset:", offset)), "",
        if (!is.null(method)) c("Between-Sample Normalization:", "-----------------------------", paste("* Method:", method), ""),
        "RUV-III Batch Correction:", "-------------------------", paste("* Status:", ruv), "",
        "Workflow Parameters:", "-------------------", "",
        "  Parameters from Final S4 Object:", "  --------------------------------",
        "[deAnalysisParameters]", "  formula_string = ~ 0 + group", de_extra, ui_block, "",
        "Contrasts:", "----------", "   KO_vs_WT=groupKO-groupWT")
}

k236Runs <- c("KO_1", "KO_2", "KO_3", "WT_1", "WT_2", "WT_3")
k236Design <- function(tech = k236Runs, runs = k236Runs, groups = sub("_[0-9]+$", "", runs)) {
    data.frame(Run = runs, group = groups, replicates = as.integer(sub("^.*_", "", runs)), tech_rep_group = tech,
        stringsAsFactors = FALSE)
}

k236Param <- function(name, value) list(name = name, source = "provenance", value = encodeWorkflowRecordValue(value),
    provenance = "default")

#' One record step: dimensions (features summed over the assays; columns = samples plus the id columns), outcomes
k236Step <- function(order, operation, rows, parameters = list(), columns = 8L, outcomes = list()) {
    list(order = order, operation = operation, dimensions = list(rows = rows, columns = columns),
        outcomes = c(list(rows = rows, columns = columns), outcomes),
        parameters = unname(Map(k236Param, names(parameters), parameters)))
}

#' A record in the application's shape (v2 lipid/metab; `prefix` NULL gives a proteomics record without steps)
k236Record <- function(prefix = "lipid", method = "cyclicloess", offset = 0, lfc = 0, treat_applied = FALSE,
    block_reason = "no_duplicates", contrasts = list(c("groupKO-groupWT", "KO_vs_WT")), tech = k236Runs, runs = k236Runs,
    fdr = NULL, pearson = TRUE, assays = c("LCMS_Pos", "LCMS_Neg"), steps = NULL, failures = NULL) {
    step <- function(order, operation, rows, parameters = list()) k236Step(order, operation, rows, parameters)
    if (is.null(steps) && !is.null(prefix)) {
        steps <- list(
            list(order = 1L, operation = "import", dimensions = list(rows = NULL, columns = NULL), parameters = list()),
            step(2L, "design_builder", 20L, list(technical_replicate_id = "tech_rep_group")),
            step(3L, paste0(prefix, "_intensity_filter"), 17L),
            step(4L, paste0(prefix, "_itsd_normalisation"), 17L),
            step(5L, paste0(prefix, "_log2"), 17L, list(offset = offset)),
            step(6L, paste0(prefix, "_between_sample_normalisation"), 17L, list(normalisation_method = method)))
        if (pearson) steps <- c(steps, list(step(7L, paste0(prefix, "_pearson_selection"), 17L,
            list(min_pearson_correlation_threshold = 0.5))))
    }
    model <- list(model = "limma", formula_used = "~ 0 + group", eBayes_trend = TRUE, eBayes_robust = TRUE)
    if (!is.null(treat_applied)) model$treat_applied <- treat_applied
    if (!is.null(block_reason)) model$block_reason <- block_reason
    results <- unlist(lapply(contrasts, function(x) lapply(if (is.null(prefix)) NA else assays, function(a) {
        entry <- list(contrast = x[[1L]], entity_id = "e1", rows = 9L)
        if (!is.na(a)) entry$assay <- a
        if (!is.null(fdr)) entry$fdr_method <- unname(fdr[[if (is.na(a)) 1L else a]])
        entry
    })), recursive = FALSE)
    list(schema = "multischolar.analysis_record.v1",
        software = list(package = "MultiScholaR", version = "0.9.0", build = "b1"),
        design = list(sample_id = "Run", group_id = "group", technical_replicate_id = "tech_rep_group",
            rows = lapply(seq_along(runs), function(i) list(Run = runs[[i]], group = sub("_[0-9]+$", "", runs[[i]]),
                tech_rep_group = tech[[i]]))),
        steps = steps %||% list(),
        analyses = list(differential_abundance = list(operation = if (is.null(prefix)) "DA-contrast-output" else
            paste0(prefix, "_da_contrast"),
            parameters = list(k236Param("formula", "~ 0 + group"), k236Param("qvalue", 0.05), k236Param("lfc", lfc)),
            contrasts = lapply(contrasts, function(x) list(contrasts = x[[1L]], friendly_names = x[[2L]],
                full_format = paste0(x[[2L]], "=", x[[1L]]))),
            model = model, results = results,
            failures = lapply(failures %||% list(), function(x) list(contrast = x[[1L]], assay = x[[2L]], error = x[[3L]])))))
}

#' One lipid/metab report through its Main: the report content, whitespace-normalised (`$norm`) and as printed (`$raw`)
k236Case <- function(name, lines = k236Lines(), design = k236Design(), record = NULL, contrast_strings = NULL,
    volcano = FALSE, numsig = FALSE, da_xlsx = FALSE, composite = TRUE, numsig_table = NULL, side = "head") {
    spec <- k236Spec[[name]]
    root <- withr::local_tempdir(pattern = "kit236-")
    env <- k236Env(name, root, side)
    writeLines(enc2utf8(lines), file.path(env$report_source_dir, "study_parameters.txt"), useBytes = TRUE)
    utils::write.table(design, file.path(env$report_source_dir, "design_matrix.tab"), sep = "\t", quote = FALSE,
        row.names = FALSE)
    if (!is.null(record)) jsonlite::write_json(record, file.path(env$report_source_dir, "analysis_record.json"),
        auto_unbox = TRUE, null = "null", digits = NA)
    if (!is.null(contrast_strings)) writeLines(contrast_strings, file.path(env$report_source_dir, "contrast_strings.tab"))
    if (volcano) k236Touch(file.path(env$report_publication_figures_dir, "Volcano_Plots", "a_volcano.png"))
    if (numsig) k236Touch(file.path(env$report_publication_figures_dir, "NumSigDeMolecules", "num_sig.png"))
    if (da_xlsx) k236Touch(file.path(env$report_publication_tables_dir, paste0("DA_results_", spec$omic, ".xlsx")))
    if (composite) k236Touch(file.path(env$report_qc_figures_dir, paste0(spec$omic, "_composite_QC_figure.png")))
    if (!is.null(numsig_table)) utils::write.table(numsig_table, file.path(env$report_publication_tables_dir,
        paste0("da_", spec$omic, "_num_sig_de_molecules.tab")), sep = "\t", quote = FALSE, row.names = FALSE)
    result <- suppressWarnings(env$Main())
    raw <- lapply(result$content, function(x) paste(as.character(x), collapse = ""))
    # The Results caption with RUV applied, from the same numbers (FormatResults called as Main calls it)
    args <- list(list(), ruv_applied = TRUE)
    if ("numbers" %in% names(formals(env$FormatResults))) args$numbers <- result$numbers
    list(raw = raw, norm = lapply(raw, k236Norm), ruv_results = k236Norm(do.call(env$FormatResults, args)))
}

#' One proteomics report's AnalyzeSamples and Format* functions
k236Prot <- function(name, design = k236Design(), record = NULL, lines = NULL, side = "head") {
    spec <- k236Spec[[name]]
    root <- withr::local_tempdir(pattern = "kit236-prot-")
    env <- k236Env(name, root, side)
    if (!is.null(record)) jsonlite::write_json(record, file.path(env$report_source_dir, "analysis_record.json"),
        auto_unbox = TRUE, null = "null", digits = NA)
    study <- list()
    if (!is.null(lines)) {
        path <- file.path(env$report_source_dir, "study_parameters.txt")
        writeLines(enc2utf8(lines), path, useBytes = TRUE)
        study <- env$ReadStudyParamsNew(path)
    }
    numbers <- suppressWarnings(env$AnalyzeSamples(design, study))$numbers
    if (!length(numbers)) stop("AnalyzeSamples returned no numbers for ", name)
    list(sample = k236Norm(env$FormatSampleDetails(numbers)), method = k236Norm(env[[spec$method]](numbers)),
        de = k236Norm(env$FormatDEDetails(numbers)))
}

#' The lines of a printed paragraph after the line equal to `heading`, up to the next blank line
k236Paragraph <- function(text, heading) {
    lines <- strsplit(text, "\n", fixed = TRUE)[[1L]]
    i <- which(trimws(lines) == heading)
    if (!length(i)) return(NA_character_)
    j <- i[[1L]] + 1L
    out <- character()
    while (j <= length(lines) && nzchar(trimws(lines[[j]]))) {
        out <- c(out, lines[[j]])
        j <- j + 1L
    }
    paste(out, collapse = "\n")
}

#' The printed paragraph whose first line starts with `start`
k236ParagraphFrom <- function(text, start) {
    lines <- strsplit(text, "\n", fixed = TRUE)[[1L]]
    i <- which(startsWith(lines, start))
    if (!length(i)) return(NA_character_)
    j <- i[[1L]]
    out <- character()
    while (j <= length(lines) && nzchar(trimws(lines[[j]]))) {
        out <- c(out, lines[[j]])
        j <- j + 1L
    }
    paste(out, collapse = "\n")
}

# --- Class S: stated findings ----------------------------------------------------------------------------------------

test_that("KIT-236 T1: Comments are a field for the analyst, not fixed findings", {
    for (name in k236Omics) {
        got <- k236Case(name, record = k236Record(k236Spec[[name]]$prefix))
        expect_false(k236Has(got$norm$comments, "All samples passed initial QC criteria"), info = name)
        expect_true(k236Has(got$norm$comments, "[To be completed by the analyst"), info = name)
    }
})

test_that("KIT-236 T1c: the computed counts print with a record only (features per step; samples kept by the filter)", {
    for (name in k236Omics) {
        with_record <- k236Case(name, record = k236Record(k236Spec[[name]]$prefix))
        expect_true(k236Has(with_record$norm$comments, paste("Features remaining after each step (total across the 2 assays):",
            "import and design 20; intensity filter 17; internal-standard normalisation 17; log2 transformation 17;",
            "between-sample normalisation 17; sample-correlation filter 17.")), info = name)
        expect_true(k236Has(with_record$norm$comments,
            "6 of 6 samples were kept by the sample-correlation filter (threshold 0.5)."), info = name)
        no_filter <- k236Case(name, record = k236Record(k236Spec[[name]]$prefix, pearson = FALSE))
        expect_false(k236Has(no_filter$norm$comments, "sample-correlation filter (threshold"), info = name)
        no_record <- k236Case(name)
        expect_false(k236Has(no_record$norm$comments, "Features remaining"), info = name)
        expect_false(k236Has(no_record$norm$comments, "samples were kept"), info = name)
    }
})

test_that("KIT-236 T2: Opinions are a field for the analyst", {
    lipid <- k236Case("lipid")
    expect_false(k236Has(lipid$norm$opinions, "data quality is sufficient for differential analysis"))
    expect_true(k236Has(lipid$norm$opinions, "[To be completed by the analyst"))
    metab <- k236Case("metab")
    expect_false(k236Has(metab$norm$opinions, "Technical reproducibility meets industry standards"))
    expect_true(k236Has(metab$norm$opinions, "[To be completed by the analyst"))
})

#' The Results caption with RUV applied (FormatResults with the case's numbers)
k236RuvResults <- function(name, ...) k236Case(name, ...)$ruv_results

test_that("KIT-236 T3: with RUV applied the caption describes the panels and draws no conclusion", {
    findings <- c("demonstrate the effectiveness of batch effect removal", "Initial log2 data showed notable batch effects",
        "batches merged effectively", "show progressive improvement in data quality",
        "high and consistent across all normalization methods", "Cyclic loess normalization (centre column")
    for (name in k236Omics) {
        text <- k236RuvResults(name)
        for (phrase in findings) expect_false(k236Has(text, phrase), info = paste(name, phrase))
        expect_true(k236Has(text, "this report draws no conclusion from them"), info = name)
        expect_true(k236Has(text, "after RUVIII-C correction"), info = name)
        expect_true(k236Has(text, "one row per assay and view"), info = name)
        expect_false(k236Has(text, "(panels a"), info = name)
    }
})

test_that("KIT-236 T3i: the left column names internal-standard normalisation only when a standard was found", {
    for (name in k236Omics) {
        none <- k236RuvResults(name)
        expect_true(k236Has(none, "this report draws no conclusion from them"), info = name)
        expect_false(k236Has(none, "before internal-standard normalisation"), info = name)
        found <- k236RuvResults(name, lines = k236Lines(itsd_block = c("* ITSD Features Per Assay:",
            "    - LCMS_Pos (2 features): IS_a, IS_b")))
        expect_true(k236Has(found, "on the intensities before internal-standard normalisation and log2 transformation (left)"),
            info = name)
    }
})

test_that("KIT-236 T4: with RUV not applied the caption describes the panels and says RUV was not applied", {
    for (name in k236Omics) {
        got <- k236Case(name)
        expect_false(k236Has(got$norm$results, "After cyclic loess normalization, technical variation was reduced"), info = name)
        expect_false(k236Has(got$norm$results, "Normalized data shows more consistent median values near zero"), info = name)
        expect_true(k236Has(got$norm$results, "RUV-III correction was not applied to this dataset."), info = name)
        expect_true(k236Has(got$norm$results, "after cyclic loess normalisation (normalizeCyclicLoess) (right)"), info = name)
    }
})

test_that("KIT-236 T5, T6n, T6r: the caption names the normalisation that ran, none, or that it is not recorded", {
    for (name in k236Omics) {
        prefix <- k236Spec[[name]]$prefix
        quantile <- k236Case(name, lines = k236Lines(method = "quantile"), record = k236Record(prefix, method = "quantile"))
        expect_false(k236Has(quantile$norm$results, "Cyclic loess normalization"), info = name)
        expect_true(k236Has(quantile$norm$results, "quantile normalisation (normalizeQuantiles)"), info = name)
        none <- k236Case(name, lines = k236Lines(method = "none"), record = k236Record(prefix, method = "none"))
        expect_true(k236Has(none$norm$results, "after the normalisation step, at which no normalisation was applied"), info = name)
        unrecorded <- k236Case(name, lines = k236Lines(method = NULL))
        expect_true(k236Has(unrecorded$norm$results, "after between-sample normalisation (the method is not recorded)"),
            info = name)
    }
})

test_that("KIT-236 T11, T11b, T11g: the sample details count the analysed design's technical replicates", {
    shared <- c("KO_1", "KO_1", "KO_3", "WT_1", "WT_2", "WT_3")
    for (name in k236Omics) {
        prefix <- k236Spec[[name]]$prefix
        record <- k236Case(name, record = k236Record(prefix, tech = shared))
        expect_false(k236Has(record$norm$sample_details, "with no technical replicates"), info = name)
        expect_true(k236Has(record$norm$sample_details, paste("The analysis included 6 measured runs from 5 biological",
            "samples; 1 biological sample was measured more than once (technical replicates), giving 1 additional runs.")),
            info = name)
        file <- k236Case(name, design = k236Design(tech = shared))
        expect_true(k236Has(file$norm$sample_details, "The design table lists 6 measured runs from 5 biological samples"),
            info = name)
        none <- k236Case(name, record = k236Record(prefix))
        expect_true(k236Has(none$norm$sample_details,
            "The analysis included 6 samples, each a separate biological sample (no technical replicates were declared in the design)."),
            info = name)
    }
})

# --- Class M: methods ------------------------------------------------------------------------------------------------

test_that("KIT-236 T7, T7b: zeros are described as they were handled", {
    for (name in k236Omics) {
        offset <- k236Case(name, lines = k236Lines(offset = 1))
        expect_false(k236Has(offset$norm$method_details, "were coded as NA"), info = name)
        expect_true(k236Has(offset$norm$method_details, "zeros were transformed with the offset"), info = name)
        plain <- k236Case(name)
        expect_false(k236Has(plain$norm$method_details, "were coded as NA"), info = name)
        expect_true(k236Has(plain$norm$method_details, "a value of zero was read as not detected"), info = name)
    }
})

test_that("KIT-236 T8: no integration step and no cross-assay check are described", {
    for (name in k236Omics) {
        got <- k236Case(name)
        expect_false(k236Has(got$norm$method_details, "before integration"), info = name)
        expect_false(k236Has(got$norm$method_details, "cross-assay consistency verified"), info = name)
        expect_true(k236Has(got$norm$method_details, paste("Each assay was processed and tested for differential abundance",
            "separately; the results of the 2 assays are reported together.")), info = name)
        # KIT-236 I09: today a sample the correlation filter removes in one assay stays in the others (KIT-255).
        expect_false(k236Has(got$norm$method_details, "removed from every assay"), info = name)
        expect_true(k236Has(got$norm$method_details, "The assays were required to contain the same samples."), info = name)
    }
})

test_that("KIT-236 T9, T9b: the comparisons are those of the contrast table", {
    three <- k236Design(runs = c("A_1", "A_2", "B_1", "B_2", "C_1", "C_2"), tech = c("A_1", "A_2", "B_1", "B_2", "C_1", "C_2"))
    for (name in k236Omics) {
        prefix <- k236Spec[[name]]$prefix
        record <- k236Record(prefix, contrasts = list(c("groupA-groupB", "A_vs_B")), runs = three$Run, tech = three$Run)
        got <- k236Case(name, design = three, record = record)
        expect_false(k236Has(got$norm$de_details, "comparing each pair of experimental groups"), info = name)
        expect_true(k236Has(got$norm$de_details, "the comparisons defined in the contrast table (1 comparison: A_vs_B)"),
            info = name)
        file <- k236Case(name, contrast_strings = "groupKO-groupWT")
        expect_true(k236Has(file$norm$de_details, "(1 comparison: KO_vs_WT)"), info = name)
    }
})

test_that("KIT-236 T10, T10b: the FDR sentence states the Benjamini-Hochberg fallback (rule, or per comparison)", {
    for (name in k236Omics) {
        prefix <- k236Spec[[name]]$prefix
        rule <- k236Case(name, record = k236Record(prefix))
        expect_true(k236Has(rule$norm$de_details,
            "False discovery rates were estimated as q-values (Storey, 2002); where they cannot be computed, Benjamini-Hochberg adjusted p-values are used."),
            info = name)
        # KIT-242 outcome 5's shape: fdr_method per contrast and assay.
        planted <- k236Case(name, record = k236Record(prefix, fdr = c(LCMS_Pos = "qvalue", LCMS_Neg = "BH")))
        expect_true(k236Has(planted$norm$de_details, paste("except for KO_vs_WT in LCMS_Neg, where q-values could not be",
            "computed and Benjamini-Hochberg adjusted p-values were used instead.")), info = name)
    }
})

test_that("KIT-236 T15: no SOP, Salesforce or test-condition placeholder prints as content; the version does", {
    for (name in k236Omics) {
        got <- k236Case(name, record = k236Record(k236Spec[[name]]$prefix))
        all_text <- k236Norm(unlist(got$norm))
        expect_false(k236Has(all_text, "SOP xxx"), info = name)
        expect_false(k236Has(all_text, "Salesforce"), info = name)
        expect_false(k236Has(all_text, "Comments on test specific conditions"), info = name)
        expect_true(k236Has(got$norm$method_details,
            "Data analyses were performed with MultiScholaR 0.9.0 (build b1); SOP reference: [to be completed by the analyst]."),
            info = name)
        expect_true(k236Has(got$norm$method_details, "[to be completed by the analyst]"), info = name)
    }
})

test_that("KIT-236 T17, T17b, T17g, T17n: the fold-change threshold and the significance rule as tested", {
    ui <- c("  User Interface Parameters:", "  -------------------------", "[Differential Expression UI Parameters]",
        "  q_value_threshold = 0.05", "  log_fold_change_cutoff = 0.5", "  treat_enabled = TRUE")
    for (name in k236Omics) {
        prefix <- k236Spec[[name]]$prefix
        treat <- k236Case(name, record = k236Record(prefix, lfc = 0.5, treat_applied = NULL))
        expect_false(k236Has(treat$norm$de_details, "test for any deviation from zero log-fold-change"), info = name)
        expect_true(k236Has(treat$norm$de_details, "log2 fold-change threshold of 0.5"), info = name)
        expect_true(k236Has(treat$norm$de_details, "absolute log2 fold change of at least 0.5"), info = name)
        # I10: the block the application writes ('Differential Expression UI Parameters'), no record.
        written <- k236Case(name, lines = k236Lines(ui_block = ui))
        expect_true(k236Has(written$norm$de_details, "log2 fold-change threshold of 0.5"), info = name)
        zero <- k236Case(name, record = k236Record(prefix))
        expect_true(k236Has(zero$norm$de_details, "for any difference from zero log2 fold change"), info = name)
        unrecorded <- k236Case(name)
        expect_true(k236Has(unrecorded$norm$de_details,
            "Whether a fold-change threshold (TREAT) was used is not recorded for this analysis."), info = name)
    }
})

test_that("KIT-236 T18, T18g: the technical-replicate block is stated when the record holds it", {
    for (name in k236Omics) {
        prefix <- k236Spec[[name]]$prefix
        for (reason in c("technical_replicates", "no_duplicates; technical_replicates")) {
            got <- k236Case(name, record = k236Record(prefix, block_reason = reason))
            expect_true(k236Has(got$norm$de_details, "modelled as a block through limma's duplicateCorrelation"),
                info = paste(name, reason))
        }
        plain <- k236Case(name, record = k236Record(prefix))
        expect_false(k236Has(plain$norm$de_details, "duplicateCorrelation"), info = name)
    }
})

# --- Class F: files and figures --------------------------------------------------------------------------------------

#' The chunk headers and bodies of a template, by label
k236Chunk <- function(name, label) Filter(function(ch) identical(ch$label, label), k236Chunks(k236Path(name)))[[1L]]

test_that("KIT-236 T12, T12g: the volcano figure is referred to only when it is there; no warning line", {
    for (name in k236Omics) {
        word <- k236Spec[[name]]$word
        absent <- k236Case(name)
        expect_false(k236Has(absent$norm$de_details, "are shown in Figure 3"), info = name)
        present <- k236Case(name, volcano = TRUE)
        expect_true(k236Has(present$norm$de_details, paste("Volcano plots of differentially abundant", word,
            "are shown in Figure 3")), info = name)
        for (label in c("grid_layout", "Volcano Plot Caption")) {
            chunk <- k236Chunk(name, label)
            expect_true(grepl("eval=", chunk$header, fixed = TRUE), info = paste(name, label))
            expect_false(any(grepl("Warning:", chunk$code, fixed = TRUE)), info = paste(name, label))
        }
    }
})

test_that("KIT-236 T13: the results table is named as the application writes it, only when it exists", {
    for (name in k236Omics) {
        spec <- k236Spec[[name]]
        present <- k236Case(name, da_xlsx = TRUE)
        expect_false(k236Has(present$norm$de_details, "DE_lipids_results.xlsx"), info = name)
        expect_false(k236Has(present$norm$de_details, "DE_metabolites_results.xlsx"), info = name)
        expect_true(k236Has(present$norm$de_details, paste0("DA_results_", spec$omic, ".xlsx")), info = name)
        absent <- k236Case(name)
        expect_false(k236Has(absent$norm$de_details, "results table"), info = name)
    }
})

test_that("KIT-236 T14: Figure 1 and Figure 2 print only when their files exist; no warning line", {
    for (name in k236Omics) {
        for (label in c("QC_metrics", "NumSigDE_Plot", "NumSigDE Caption")) {
            chunk <- k236Chunk(name, label)
            expect_true(grepl("eval=", chunk$header, fixed = TRUE), info = paste(name, label))
            expect_false(any(grepl("Warning:", chunk$code, fixed = TRUE)), info = paste(name, label))
        }
        got <- k236Case(name, composite = FALSE)
        expect_true(k236Has(got$norm$results, "The quality-control figure was not produced for this analysis."), info = name)
    }
})

test_that("KIT-236 T19, T19g: the count sentence reads the table Copy writes (da_<omic>_...), per comparison and assay", {
    for (name in k236Omics) {
        word <- k236Spec[[name]]$word
        table <- data.frame(total = c(9L, 8L), significant = c(3L, 0L), up_regulated = c(2L, 0L), down_regulated = c(1L, 0L),
            assay = c("LCMS_Pos", "LCMS_Neg"), mode = c("pos", "neg"), contrast = "groupKO-groupWT", q_threshold = 0.05)
        got <- k236Case(name, numsig_table = table, numsig = TRUE, record = k236Record(k236Spec[[name]]$prefix))
        expect_true(k236Has(got$norm$de_summary, paste0("The number of significant differentially abundant ", word,
            " is summarised in Figure 2 (q-value < 0.05).")), info = name)
        expect_false(k236Has(got$norm$de_summary, "groupKO-groupWT"), info = name)
        expect_true(k236Has(got$norm$de_summary, paste0("For the comparison KO_vs_WT in LCMS_Pos, there were 2 upregulated ",
            "and 1 downregulated ", word, ".")), info = name)
        none <- k236Case(name)
        expect_identical(none$norm$de_summary, "", info = name)
    }
})

# --- G1 and the waiting sentences (W) ---------------------------------------------------------------------------------

test_that("KIT-236 G1: the Methods paragraphs outside section 1 print the base text, with and without a record", {
    for (name in k236Omics) {
        for (record in list(NULL, k236Record(k236Spec[[name]]$prefix))) {
            info <- paste(name, if (is.null(record)) "file only" else "record")
            head <- k236Case(name, record = record)$raw$method_details
            base <- k236Case(name, record = record, side = "base")$raw$method_details
            for (heading in c("Batch Effect Correction", "Quality Control", "ITSD Normalization")) {
                expect_identical(k236Paragraph(head, heading), k236Paragraph(base, heading), info = paste(info, heading))
            }
            expect_identical(k236ParagraphFrom(head, "Between-sample normalization was performed"),
                k236ParagraphFrom(base, "Between-sample normalization was performed"), info = info)
            expect_true(k236Has(head, "The dataset contains 2 assays: LCMS_Pos, LCMS_Neg."), info = info)
            # The log2 sentence (KIT-172, the owner's wording pass SUM-7) is unchanged; only the zeros sentence after it went.
            head_log <- k236Norm(k236Paragraph(head, "Data Preprocessing and Normalization"))
            base_log <- k236Norm(k236Paragraph(base, "Data Preprocessing and Normalization"))
            zeros <- paste("Missing values (zeros or non-detects) were coded as NA to distinguish true absence from missing",
                "measurements.")
            expect_identical(k236Norm(sub(zeros, "", head_log, fixed = TRUE)), k236Norm(sub(zeros, "", base_log, fixed = TRUE)),
                info = info)
        }
    }
})

test_that("KIT-236 W1, W3: the normalisation sentence of a method applied, and the internal-standard sentence, wait for KIT-242", {
    for (name in k236Omics) {
        prefix <- k236Spec[[name]]$prefix
        quantile <- k236Case(name, lines = k236Lines(method = "quantile"), record = k236Record(prefix, method = "quantile"))
        expect_true(k236Has(quantile$norm$method_details,
            "Between-sample normalization was performed using the 'quantile' method from the 'limma' R package"), info = name)
        itsd <- k236Case(name)
        expect_true(k236Has(itsd$norm$method_details,
            "Internal standard (ITSD) normalization was applied using the average_centered method with median aggregation"),
            info = name)
    }
})

# --- Part B: the proteomics reports -----------------------------------------------------------------------------------

k236ProtContrasts <- list(c("groupKO-groupWT", "KO_vs_WT"))

test_that("KIT-236 B1: the proteomics sample details count the declared technical replicates", {
    shared <- c("KO_1", "KO_1", "KO_3", "WT_1", "WT_1", "WT_3")
    for (name in k236ProtNames) {
        got <- k236Prot(name, design = k236Design(tech = shared))
        expect_false(k236Has(got$sample, "with no technical replicates"), info = name)
        expect_true(k236Has(got$sample, paste("The design table lists 6 measured runs from 4 biological samples; 2 biological",
            "samples were measured more than once (technical replicates)")), info = name)
    }
})

test_that("KIT-236 B2 (I12): the proteomics 'none' normalisation says none was applied", {
    steps <- list(list(order = 1L, operation = "normalization-ruv", parameters = list(k236Param("norm_method", "none"))))
    for (name in k236ProtNames) {
        got <- k236Prot(name, record = k236Record(NULL, steps = steps, contrasts = k236ProtContrasts))
        expect_false(k236Has(got$method, "'none' method from the 'limma' R package"), info = name)
        expect_true(k236Has(got$method, "no between-sample normalisation was applied"), info = name)
    }
})

test_that("KIT-236 B3: the proteomics FDR sentence names a Benjamini-Hochberg comparison from the record", {
    for (name in k236ProtNames) {
        planted <- k236Prot(name, record = k236Record(NULL, contrasts = k236ProtContrasts, fdr = "BH"))
        expect_true(k236Has(planted$de, paste("except for KO_vs_WT, where q-values could not be computed and",
            "Benjamini-Hochberg adjusted p-values were used instead.")), info = name)
        rule <- k236Prot(name, record = k236Record(NULL, contrasts = k236ProtContrasts))
        expect_true(k236Has(rule$de, "where they cannot be computed, Benjamini-Hochberg adjusted p-values are used."),
            info = name)
    }
})

test_that("KIT-236 B4: the proteomics TREAT sentence as tested; DIA with limpa from the record's model only", {
    for (name in c("LFQ", "TMT", "DIANN")) {
        got <- k236Prot(name, record = k236Record(NULL, lfc = 0.5, treat_applied = NULL, contrasts = k236ProtContrasts))
        expect_false(k236Has(got$de, "test for any deviation from zero log-fold-change"), info = name)
        expect_true(k236Has(got$de, "log2 fold-change threshold of 0.5"), info = name)
        # The significance definition is KIT-220's and stays.
        expect_true(k236Has(got$de, "Significant differentially expressed proteins were defined as those with q-values less than"),
            info = name)
    }
    limpa <- k236Prot("DIANN_limpa", record = k236Record(NULL, lfc = 0.5, treat_applied = NULL, contrasts = k236ProtContrasts))
    expect_true(k236Has(limpa$de, "Whether a fold-change threshold (TREAT) was used is not recorded for this analysis."))
    limpa_zero <- k236Prot("DIANN_limpa", record = k236Record(NULL, lfc = 0.5, treat_applied = FALSE,
        contrasts = k236ProtContrasts))
    expect_true(k236Has(limpa_zero$de, "for any difference from zero log2 fold change"))
    # R1-N03: the non-limpa branch with a file cutoff of 0.5 and no record.
    lines <- c("Workflow Parameters:", "-------------------", "", "  Parameters from Final S4 Object:",
        "  --------------------------------", "[deAnalysisParameters]", "  formula_string = ~ 0 + group",
        "  treat_lfc_cutoff = 0.5", "  User Interface Parameters:", "  -------------------------",
        "[Differential Expression UI Parameters]", "  q_value_threshold = 0.05", "  log_fold_change_cutoff = 0.5",
        "  treat_enabled = TRUE", "")
    file_only <- k236Prot("DIANN_limpa", lines = lines)
    expect_false(k236Has(file_only$de, "TREAT method (McCarthy"))
    expect_true(k236Has(file_only$de, "Whether a fold-change threshold (TREAT) was used is not recorded for this analysis."))
})

test_that("KIT-236 B5: the proteomics technical-replicate block is stated when the record holds it", {
    for (name in c("LFQ", "TMT", "DIANN")) {
        got <- k236Prot(name, record = k236Record(NULL, block_reason = "technical_replicates", contrasts = k236ProtContrasts))
        expect_true(k236Has(got$de, "Technical replicates were modelled as a block through limma's duplicateCorrelation."),
            info = name)
        plain <- k236Prot(name, record = k236Record(NULL, contrasts = k236ProtContrasts))
        expect_false(k236Has(plain$de, "duplicateCorrelation"), info = name)
    }
})

test_that("KIT-236 B6: DIA with limpa no longer says technical replicates were averaged", {
    plain <- k236Prot("DIANN_limpa")
    expect_false(k236Has(plain$method, "Technical replicate averaging"))
    expect_false(k236Has(plain$method, "were averaged across replicates"))
    block <- k236Prot("DIANN_limpa", record = k236Record(NULL, block_reason = "technical_replicates",
        contrasts = k236ProtContrasts))
    expect_false(k236Has(block$method, "Technical replicate averaging"))
    expect_true(k236Has(block$method, "Technical replicates were modelled as a block through limma's duplicateCorrelation."))
})

# --- The correction pass (AUDIT-KIT-236-impl-20261125 I01, I02, I05 to I09, I12) -----------------------------------

test_that("KIT-236 T6, T6b (I12): the Methods say no normalisation was applied, or that the method is not recorded", {
    for (name in k236Omics) {
        prefix <- k236Spec[[name]]$prefix
        none <- k236Case(name, lines = k236Lines(method = "none"), record = k236Record(prefix, method = "none"))
        expect_false(k236Has(none$norm$method_details, "'none' method from"), info = name)
        expect_true(k236Has(none$norm$method_details, "No between-sample normalisation was applied."), info = name)
        unrecorded <- k236Case(name, lines = k236Lines(method = NULL))
        expect_false(k236Has(unrecorded$norm$method_details, "'cyclicloess' method from"), info = name)
        expect_true(k236Has(unrecorded$norm$method_details,
            "The between-sample normalisation method is not recorded for this analysis."), info = name)
    }
})

test_that("KIT-236 T1d (I02): the samples kept are counted from the filter's own step", {
    for (name in k236Omics) {
        prefix <- k236Spec[[name]]$prefix
        runs <- k236Runs[-6L]
        steps <- list(k236Step(2L, "design_builder", 20L, columns = 8L), k236Step(3L, paste0(prefix, "_intensity_filter"), 17L),
            k236Step(4L, paste0(prefix, "_between_sample_normalisation"), 17L, list(normalisation_method = "cyclicloess"),
                columns = 7L),
            k236Step(5L, paste0(prefix, "_pearson_selection"), 17L, list(min_pearson_correlation_threshold = 0.5), columns = 7L))
        earlier <- k236Case(name, record = k236Record(prefix, runs = runs, tech = runs, steps = steps))
        expect_true(k236Has(earlier$norm$comments, "5 of 5 samples were kept by the sample-correlation filter (threshold 0.5)."),
            info = name)
        steps[[3L]] <- k236Step(4L, paste0(prefix, "_between_sample_normalisation"), 17L, list(normalisation_method = "cyclicloess"))
        filtered <- k236Case(name, record = k236Record(prefix, runs = runs, tech = runs, steps = steps))
        expect_true(k236Has(filtered$norm$comments, "5 of 6 samples were kept by the sample-correlation filter (threshold 0.5)."),
            info = name)
    }
})

test_that("KIT-236 T1e (I05): the step list never prints an internal operation name", {
    for (name in k236Omics) {
        prefix <- k236Spec[[name]]$prefix
        steps <- list(k236Step(2L, "design_builder", 20L), k236Step(3L, paste0(prefix, "_intensity_filter"), 17L),
            k236Step(4L, paste0(prefix, "_ruv"), 17L), k236Step(5L, paste0(prefix, "_norm_complete"), 17L),
            k236Step(6L, paste0(prefix, "_reset"), 17L), k236Step(7L, "legacy", 17L))
        got <- k236Case(name, record = k236Record(prefix, steps = steps))
        expect_false(k236Has(got$norm$comments, "norm complete"), info = name)
        expect_false(k236Has(got$norm$comments, "legacy"), info = name)
        expect_false(k236Has(got$norm$comments, paste(prefix, "reset")), info = name)
        expect_true(k236Has(got$norm$comments, "intensity filter 17; RUV correction 17."), info = name)
    }
})

test_that("KIT-236 T20 (I01): under the declared RUV design the DA paragraph names the uncorrected matrix and the factors", {
    for (name in k236Omics) {
        prefix <- k236Spec[[name]]$prefix
        record <- k236Record(prefix)
        record$steps <- c(record$steps, list(k236Step(8L, paste0(prefix, "_ruv"), 17L, list(ruv_design = "ruv3_factors"),
            outcomes = list(ruv_design = "ruv3_factors", best_k = list(LCMS_Pos = 2L, LCMS_Neg = 3L)))))
        got <- k236Case(name, record = record)
        expect_false(k236Has(got$norm$de_details, "using the adjusted abundance matrix"), info = name)
        expect_true(k236Has(got$norm$de_details, paste("on the uncorrected abundance matrix, with each assay's RUV factors as",
            "covariates in the model (LCMS_Pos: 2 factors; LCMS_Neg: 3 factors)")), info = name)
        plain <- k236Case(name, record = k236Record(prefix))
        expect_true(k236Has(plain$norm$de_details, "using the adjusted abundance matrix"), info = name)
    }
})

test_that("KIT-236 T21, B7 (I06): at a fold-change threshold the p-values are named as TREAT's", {
    for (name in k236Omics) {
        prefix <- k236Spec[[name]]$prefix
        treat <- k236Case(name, record = k236Record(prefix, lfc = 0.5, treat_applied = TRUE))
        expect_false(k236Has(treat$norm$de_details, "p-values were calculated using the 'eBayes' function"), info = name)
        expect_true(k236Has(treat$norm$de_details, "p-values were calculated with limma's TREAT test"), info = name)
        zero <- k236Case(name, record = k236Record(prefix))
        expect_true(k236Has(zero$norm$de_details, "p-values were calculated using the 'eBayes' function"), info = name)
    }
    for (name in c("LFQ", "TMT", "DIANN")) {
        got <- k236Prot(name, record = k236Record(NULL, lfc = 0.5, treat_applied = TRUE, contrasts = k236ProtContrasts))
        expect_false(k236Has(got$de, "p-values were calculated using the 'eBayes' function"), info = name)
        expect_true(k236Has(got$de, "p-values were calculated with limma's TREAT test"), info = name)
    }
})

test_that("KIT-236 T19c (I07): at a fold-change threshold the count and volcano sentences state both criteria", {
    for (name in k236Omics) {
        word <- k236Spec[[name]]$word
        table <- data.frame(total = 9L, significant = 3L, up_regulated = 2L, down_regulated = 1L, assay = "LCMS_Pos",
            mode = "pos", contrast = "groupKO-groupWT", q_threshold = 0.05)
        got <- k236Case(name, numsig_table = table, numsig = TRUE, volcano = TRUE,
            record = k236Record(k236Spec[[name]]$prefix, lfc = 0.5, treat_applied = TRUE))
        both <- "(q-value < 0.05 and an absolute log2 fold change of at least 0.5)"
        expect_true(k236Has(got$norm$de_summary, paste0("is summarised in Figure 2 ", both, ".")), info = name)
        expect_true(k236Has(got$norm$de_details, paste("with a threshold of q-value < 0.05 and an absolute log2 fold change",
            "of at least 0.5.")), info = name)
    }
})

test_that("KIT-236 T22 (I09): an assay whose test failed is named, and the Methods do not say every assay was tested", {
    for (name in k236Omics) {
        prefix <- k236Spec[[name]]$prefix
        record <- k236Record(prefix, assays = "LCMS_Pos", failures = list(c("groupKO-groupWT", "LCMS_Neg", "too few samples")))
        got <- k236Case(name, record = record)
        expect_true(k236Has(got$norm$de_details,
            "The test could not be computed for KO_vs_WT in LCMS_Neg (too few samples)."), info = name)
        expect_true(k236Has(got$norm$method_details, paste("the results are reported together for the assays and comparisons",
            "where the test could be computed")), info = name)
        expect_false(k236Has(got$norm$method_details, "the results of the 2 assays are reported together"), info = name)
    }
})

test_that("KIT-255 (I01, I02): the comment counts each sample set missing once and says which comparisons were not tested", {
    for (name in k236Omics) {
        record <- k236Record(k236Spec[[name]]$prefix)
        record$steps[[7L]]$outcomes$correlation_excluded_per_assay <- list(
            LCMS_Neg = list(A_1 = list(best_correlation = -0.015, threshold = 0.5)),
            GCMS = list(A_1 = list(best_correlation = 0.1, threshold = 0.5)))
        record$steps[[7L]]$outcomes$correlation_not_tested_per_assay <- list(LCMS_Neg = list("B"))
        got <- k236Case(name, record = record)
        expect_true(k236Has(got$norm$comments, paste("6 of 6 samples were kept by the sample-correlation filter (threshold 0.5);",
            "1 sample was set missing in some assays only (A_1 in LCMS_Neg, GCMS).")), info = name)
        expect_true(k236Has(got$norm$comments, paste("Comparisons with group B in LCMS_Neg were not tested",
            "(no sample of the group has a value there).")), info = name)
    }
})
