# KIT-227: the report templates' RUV, grouping and consequent sentences state what was run, for every mode (ticket
# section 4, rows E1 to E10 and G1). Each proteomics, lipidomics and metabolomics case goes the real path: a
# study_parameters.txt, the template's own parser, its AnalyzeSamples (the analysis record overlay and, at the fix, the
# template's own read of the RUV step's parameters), then its FormatMethodDetails* and the `results` call of its Main.
# The transcriptomics cases evaluate the template's `ruv-text` chunk. G1 compares the default text with the base
# templates' chunks (tests/testdata/kit227/base, verbatim at cf8af245) byte for byte.

k227Templates <- list(
    LFQ = list(dir = "proteomics", file = "LFQ_report.rmd", parser = "New Hierarchical Parser",
        extract = "ExtractStudyParametersLFQ", method = "FormatMethodDetailsLFQ"),
    TMT = list(dir = "proteomics", file = "TMT_report.rmd", parser = "New Hierarchical Parser",
        extract = "ExtractStudyParametersTMT", method = "FormatMethodDetailsTMT"),
    DIANN = list(dir = "proteomics", file = "DIANN_report.rmd", parser = "New Hierarchical Parser",
        extract = "ExtractStudyParametersNew", method = "FormatMethodDetails"),
    DIANN_limpa = list(dir = "proteomics", file = "DIANN_limpa_report.rmd", parser = "New Hierarchical Parser",
        extract = "ExtractStudyParametersNew", method = "FormatMethodDetails"),
    lipid = list(dir = "lipidomics", file = "lipidomics_report.rmd", parser = "Hierarchical Parser",
        extract = "ExtractLipidStudyParameters", method = "FormatMethodDetails"),
    metab = list(dir = "metabolomics", file = "metabolomics_report.rmd", parser = "Hierarchical Parser",
        extract = "ExtractStudyParametersNew", method = "FormatMethodDetails"))

k227Proteomics <- c("LFQ", "TMT", "DIANN")
k227Omics <- c("lipid", "metab")

#' A template file: the package's (head) or the base copy of its chunks
k227Path <- function(dir, file, side) {
    if (identical(side, "head")) system.file("reports", dir, file, package = "MultiScholaR") else
        testthat::test_path("..", "testdata", "kit227", "base", file)
}

#' One chunk of a template, extracted with knitr::purl and parsed
k227Exprs <- function(path, label) {
    lines <- readLines(path, warn = FALSE, encoding = "UTF-8")
    start <- grep(paste0("^```\\{r ", label, "[,}]"), lines)
    stopifnot(length(start) == 1L)
    end <- start + which(lines[-seq_len(start)] == "```")[[1L]]
    chunk <- tempfile(fileext = ".Rmd")
    writeLines(lines[start:end], chunk, useBytes = TRUE)
    code <- tempfile(fileext = ".R")
    suppressMessages(knitr::purl(chunk, output = code, quiet = TRUE, documentation = 0))
    parse(code, keep.source = FALSE, encoding = "UTF-8")
}

#' The functions a template defines in its parser, method, result, sample-analysis and main chunks, in one environment
#' with the stubs of the chunks not evaluated (the import format, the D66 log, `tibble`)
k227Env <- function(name, side) {
    spec <- k227Templates[[name]]
    path <- k227Path(spec$dir, spec$file, side)
    env <- new.env(parent = globalenv())
    env$params <- list(workflow_name = "maxquant", omic_type = spec$dir, experiment_label = "kit227")
    env$d66_log <- function(...) invisible(NULL)
    env$tibble <- tibble::tibble
    for (label in c(spec$parser, "Method Section", "Result Section", "Sample Analysis Functions",
        "Main Execution Function")) {
        for (e in k227Exprs(path, label)) {
            if (is.call(e) && identical(e[[1L]], as.name("<-")) && is.call(e[[3L]]) &&
                identical(e[[3L]][[1L]], as.name("function"))) eval(e, env)
        }
    }
    env
}

#' The `results` element of Main's report content (the Results call as the template makes it)
k227ResultsCall <- function(env) {
    find <- function(expr) {
        if (!is.call(expr)) return(NULL)
        if (!is.null(names(expr)) && all(c("results", "method_details") %in% names(expr))) return(expr[["results"]])
        for (i in seq_along(expr)) if (is.call(expr[[i]])) {
            hit <- find(expr[[i]])
            if (!is.null(hit)) return(hit)
        }
        NULL
    }
    call <- find(body(env$Main))
    stopifnot(!is.null(call))
    call
}

#' An encoded step parameter, as the record stores it (`encodeWorkflowRecordValue`)
k227Param <- function(name, value) list(name = name, value = list(type = typeof(value), value = as.list(value),
    length = length(value)))

k227Step <- function(operation, parameters = list(), outcomes = list()) {
    list(operation = operation, parameters = unname(Map(k227Param, names(parameters), parameters)), outcomes = outcomes)
}

k227WriteRecord <- function(steps, path) {
    jsonlite::write_json(list(steps = steps), path, auto_unbox = TRUE, null = "null", digits = NA)
}

#' Run one case through a template: `list(method, results)`. `numbers` given: the method and Results calls on it
#' alone (AnalyzeSamples' error path).
k227Case <- function(name, file_lines = character(), record = NULL, side = "head", numbers = NULL) {
    spec <- k227Templates[[name]]
    env <- k227Env(name, side)
    dir <- withr::local_tempdir(pattern = "kit227-")
    env$report_source_dir <- dir
    if (is.null(numbers)) {
        params_path <- file.path(dir, "study_parameters.txt")
        # The lipid/metab method text needs the assay count (`num_assays > 1`).
        assays <- if (spec$dir %in% c("lipidomics", "metabolomics")) c("Assay Information:", "-----------------",
            "- Number of Assays: 1", "- Assay Names: LCMS_Pos", "")
        writeLines(c("Study Parameters", "================", "", "Basic Information:", "-----------------",
            "* Workflow Name: kit227", "", assays, file_lines), params_path, useBytes = TRUE)
        if (!is.null(record)) k227WriteRecord(record, file.path(dir, "analysis_record.json"))
        study <- env$ReadStudyParamsNew(params_path)
        numbers <- env$AnalyzeSamples(NULL, study)$numbers
        if (!length(numbers)) {
            # AnalyzeSamples returns list() on any error: surface the error outside its handler.
            extracted <- env[[spec$extract]](study)
            workflowRecordReportParameters(file.path(dir, "analysis_record.json"), extracted, study)
            stop("AnalyzeSamples returned no numbers for ", name)
        }
    }
    results_env <- new.env(parent = env)
    results_env$analysis_results <- list(numbers = numbers, table = data.frame())
    results_env$correlation_data <- list()
    list(method = paste(as.character(env[[spec$method]](numbers)), collapse = ""),
        results = paste(as.character(eval(k227ResultsCall(env), results_env)), collapse = ""))
}

k227Has <- function(text, pattern) grepl(pattern, text, fixed = TRUE)

#' KIT-236: the lines after the line equal to `heading` up to the next blank line followed by a heading line
k227Paragraph <- function(text, heading) {
    lines <- strsplit(text, "\n", fixed = TRUE)[[1L]]
    start <- which(lines == heading)
    if (!length(start)) return(NA_character_)
    end <- start[[1L]] + 1L
    while (end < length(lines) && !(lines[[end]] == "" && nzchar(lines[[end + 1L]]))) end <- end + 1L
    paste(lines[(start[[1L]] + 1L):(end - 1L)], collapse = "\n")
}

# --- File sections, in the writer's shape (R/func_general_filemgmt_study_parameter_support.R:230-317) ----------------

k227ProtSection <- function(percentage = "10.0%", k = "2", separation = "0.5000", composite = "0.4000", controls = "50",
    grouping = "group", metric = "max_difference", penalty = "0.5", adaptive = "TRUE", size = "12") {
    c("Automatic RUV Optimization Results:", "-----------------------------------",
        paste("* Best percentage:", percentage), paste("* Best k value:", k), paste("* Separation score:", separation),
        paste("* Composite score:", composite), paste("* Control genes:", controls),
        paste("* RUV grouping variable:", grouping), paste("* Separation metric:", metric),
        paste("* K penalty weight:", penalty), paste("* Adaptive penalty:", adaptive), paste("* Sample size:", size), "")
}
k227ProtAuto <- k227ProtSection()
# buildProtNormManualRuvResult: no scores, separation metric "manual", penalty NA, adaptive FALSE, no sample size.
k227ProtManual <- k227ProtSection(percentage = "5.0%", separation = "N/A", composite = "N/A", metric = "manual",
    penalty = "NA", adaptive = "FALSE", size = "N/A")
k227ProtSkip <- c("RUV-III Batch Correction:", "-------------------------", "* Status: Not Applied",
    "* Reason: User determined RUV was not appropriate due to dataset constraints", "")
k227ProtCounts <- c("Protein Filtering Summary:", "-------------------------", "* Proteins after QC filtering: 1000",
    "* Proteins after RUV filtering: 950", "")

# The lipid/metab parser reads '-' and '•' bullets only (not the writer's '*'), so the per-assay branch is reached from
# a '-' file; the '*' shape is the writer's, from which only the record carries the values.
k227OmicSection <- function(status = "Applied", grouping = "group", assays = "LCMS_Pos: k=2, controls=19") {
    c("RUV-III Batch Correction:", "-------------------------", paste("- Status:", status),
        if (!is.null(grouping)) paste("- Grouping Variable:", grouping),
        if (length(assays)) c("- Per-Assay Parameters:", paste0("    - ", assays)), "")
}
k227OmicWriterApplied <- c("RUV-III Batch Correction:", "-------------------------", "* Status: Applied",
    "* Grouping Variable: group", "* Per-Assay Parameters:", "    - LCMS_Pos: k=2, controls=19", "")
k227OmicWriterSkip <- c("RUV-III Batch Correction:", "-------------------------", "* Status: Not Applied", "")

k227ProtRecord <- function(mode = "automatic", percentage = 5, grouping = "group", applied = TRUE, k = 2,
    outcome_percentage = 10, count = 50) {
    parameters <- list(ruv_mode = mode, ruv_percentage = percentage, ruv_k = 3, grouping = grouping)
    list(k227Step("normalization-ruv", parameters[!vapply(parameters, is.null, logical(1))],
        list(ruv_applied = applied, k = k, negative_control_percentage = outcome_percentage,
            negative_control_count = count)))
}

k227OmicRecord <- function(name, mode = "automatic", grouping = "group", k = 2, percentage = 10, count = 19) {
    operation <- if (identical(name, "lipid")) "lipid_ruv" else "metab_ruv"
    list(k227Step(operation, list(ruv_mode = mode, ruv_grouping_variable = grouping, manual_k = 3, manual_percentage = 5),
        list(ruv_applied = TRUE, k = k, negative_control_percentage = percentage, negative_control_count = count)))
}

# --- E1, E1b: Manual -----------------------------------------------------------------------------------------------

test_that("KIT-227 E1: a manual RUV from the file is not described as an automatic optimisation (LFQ, TMT, DIA-NN)", {
    for (name in k227Proteomics) {
        text <- k227Case(name, k227ProtManual)$method
        expect_false(k227Has(text, "an automatic optimisation process was run"), info = name)
        expect_true(k227Has(text, "no automatic optimisation was run"), info = name)
        expect_true(k227Has(text, paste("For this study, no automatic optimisation was run: the analyst chose to use 5.0% of",
            "features as empirical negative controls (50 proteins) and to remove 2 unwanted components.")), info = name)
        expect_true(k227Has(text, "2 unwanted"), info = name)
    }
})

test_that("KIT-227 E1b: a manual RUV whose plot failed (record only) prints the analyst's percentage, never NA%", {
    record <- k227ProtRecord(mode = "manual", percentage = 5, outcome_percentage = NULL)
    for (name in k227Proteomics) {
        text <- k227Case(name, character(), record)$method
        expect_false(k227Has(text, "an automatic optimisation process was run"), info = name)
        expect_false(k227Has(text, "NA%"), info = name)
        expect_true(k227Has(text, "no automatic optimisation was run"), info = name)
        expect_true(k227Has(text, "5%"), info = name)
        expect_true(k227Has(text, "the analyst chose to use 5% of features as empirical negative controls (50 proteins)"),
            info = name)
    }
})

# --- E2, E6: Skip ----------------------------------------------------------------------------------------------------

test_that("KIT-227 E2: a skipped RUV is not attributed to dataset constraints (LFQ, TMT, DIA-NN, lipid, metab)", {
    for (name in k227Proteomics) {
        text <- k227Case(name, k227ProtSkip)$method
        expect_false(k227Has(text, "due to dataset constraints"), info = name)
        expect_true(k227Has(text, "was not applied"), info = name)
        expect_true(k227Has(text, paste("RUV-III batch correction was not applied. The analysis proceeded directly with",
            "normalised data for differential abundance analysis.")), info = name)
    }
    for (name in k227Omics) {
        for (lines in list(k227OmicSection(status = "Not Applied", grouping = NULL, assays = character()),
            k227OmicWriterSkip)) {
            text <- k227Case(name, lines)$method
            expect_false(k227Has(text, "due to dataset constraints"), info = name)
            expect_true(k227Has(text, "was not applied"), info = name)
            expect_true(k227Has(text, paste("RUV-III batch correction was not applied. The analysis proceeded directly with",
                "normalized data for differential abundance analysis.")), info = name)
        }
    }
})

test_that("KIT-227 E6: after a skipped RUV the count sentence does not follow RUV-III (LFQ, TMT, DIA-NN)", {
    for (name in k227Proteomics) {
        text <- k227Case(name, c(k227ProtSkip, k227ProtCounts))$method
        expect_false(k227Has(text, "Following RUV-III batch correction"), info = name)
        expect_true(k227Has(text, "was not applied"), info = name)
    }
})

# --- E3, E3b: nothing recorded ----------------------------------------------------------------------------------------

test_that("KIT-227 E3: with nothing recorded the report says so, in the method and the Results caption", {
    for (name in c(k227Proteomics, k227Omics)) {
        got <- k227Case(name, character())
        expect_false(k227Has(got$method, "not applied due to dataset constraints"), info = name)
        expect_true(k227Has(got$method, "Whether RUV-III batch correction was applied is not recorded for this analysis."),
            info = name)
        expect_true(k227Has(got$results, "is not recorded"), info = name)
        expect_true(k227Has(got$results, "Figure 1 shows four quality-control views at each stage of processing drawn for this analysis."),
            info = name)
    }
})

test_that("KIT-227 E3b: when AnalyzeSamples returned nothing the report says RUV is not recorded (LFQ, TMT, DIA-NN)", {
    for (name in k227Proteomics) {
        got <- k227Case(name, numbers = list())
        expect_false(k227Has(got$method, "auto-determined"), info = name)
        expect_true(k227Has(got$method, "is not recorded"), info = name)
        expect_true(k227Has(got$results, "is not recorded"), info = name)
    }
})

# --- E4, E7: DIA with limpa ---------------------------------------------------------------------------------------------

test_that("KIT-227 E4: DIA with limpa without a DPC engine: Skip says not applied, nothing recorded says so", {
    skip_case <- k227Case("DIANN_limpa", k227ProtSkip)
    expect_false(k227Has(skip_case$method, "an automatic optimisation process was run"))
    expect_false(k227Has(skip_case$method, "due to dataset constraints"))
    expect_true(k227Has(skip_case$method, "was not applied"))
    none <- k227Case("DIANN_limpa", character())
    expect_false(k227Has(none$method, "an automatic optimisation process was run"))
    expect_true(k227Has(none$method, "Whether RUV-III batch correction was applied is not recorded for this analysis."))
    expect_true(k227Has(none$results, "is not recorded"))
    expect_false(k227Has(none$results, "after RUVIII-C correction"))
})

test_that("KIT-227 E7: the limpa route's Results caption claims covariates only when they were the model's", {
    dpc <- k227Step("protein_dpc_quant", outcomes = list(engine = "fastlimpadpc", engine_version = "0.1.0"))
    correction <- k227Case("DIANN_limpa", character(), list(dpc, k227Step("protein_ruv_covariates",
        outcomes = list(k = 2, ruv_mode = "correction", percentage_as_neg_ctrl = 10, control_count = 100,
            replicate_structure = "group"))))
    expect_false(k227Has(correction$results, "carried the RUV-III factors as covariates"))
    expect_true(k227Has(correction$results,
        "after RUV-III correction (the corrected matrix entered the differential abundance model)"))
    skipped <- k227Case("DIANN_limpa", character(), list(dpc))
    expect_false(k227Has(skipped$results, "carried the RUV-III factors as covariates"))
    # KIT-266 (VAL-316 (4), declared edit): a limpa project with no covariates step has the two-column Figure 1, so its Results
    # paragraph names the two stages and their panel letters.
    expect_true(k227Has(skipped$results, paste0("Figure 1 shows four quality-control views at two stages of processing: ",
        "DPC-Quant's protein estimates (a complete log2 matrix; model-based where peptides were not detected) ",
        "(left column: panels a, c, e, g) and after the peptide count filter, with no between-sample normalisation applied, ",
        "before the sample-correlation filter (right column: b, d, f, h). RUV-III correction was not applied to this dataset.")))
    # The covariates mode keeps its caption.
    covariates <- k227Case("DIANN_limpa", character(), list(dpc, k227Step("protein_ruv_covariates",
        outcomes = list(k = 2, ruv_mode = "covariates", percentage_as_neg_ctrl = 10, control_count = 100,
            replicate_structure = "group"))))
    expect_true(k227Has(covariates$results, "carried the RUV-III factors as covariates"))
})

# --- E5, E9, E10: lipid/metab grouping and mode -------------------------------------------------------------------------

test_that("KIT-227 E5: lipid/metab with the grouping not recorded print 'not recorded', not 'NA'", {
    for (name in k227Omics) {
        text <- k227Case(name, k227OmicSection(grouping = NULL))$method
        expect_false(k227Has(text, "'NA' variable"), info = name)
        expect_true(k227Has(text, "not recorded"), info = name)
        expect_true(k227Has(text, paste("RUV correction was applied independently to each assay (the design column that",
            "defined the replicate samples is not recorded).")), info = name)
        # The writer's shape ('*' bullets, which the parser does not read) with a default record: the grouping is the
        # record's, not 'NA'.
        writer <- k227Case(name, k227OmicWriterApplied, k227OmicRecord(name))$method
        expect_false(k227Has(writer, "'NA' variable"), info = name)
        expect_true(k227Has(writer, "using the 'group' variable to preserve biological signal."), info = name)
    }
})

test_that("KIT-227 E9: lipid/metab grouped by a column other than group describe how it was used", {
    for (name in k227Omics) {
        text <- k227Case(name, k227OmicSection(grouping = "batch"))$method
        expect_false(k227Has(text, "'batch' variable to preserve biological signal"), info = name)
        expect_true(k227Has(text, "'batch'"), info = name)
        expect_true(k227Has(text, paste("RUV correction was applied independently to each assay. Samples with the same value",
            "of the design column 'batch' were treated as replicates, so the differences between those values were kept",
            "and other differences between samples could be removed as unwanted variation.")), info = name)
    }
})

test_that("KIT-227 E10: a manual lipid/metab RUV says the analyst set k and the controls", {
    for (name in k227Omics) {
        text <- k227Case(name, k227OmicSection(), k227OmicRecord(name, mode = "manual", percentage = 5))$method
        expect_true(k227Has(text, "set by the analyst"), info = name)
        expect_true(k227Has(text, paste("The number of unwanted components and the share of features used as controls were",
            "set by the analyst; no automatic optimisation was run.")), info = name)
    }
})

# --- E8: transcriptomics Manual ---------------------------------------------------------------------------------------

#' The transcriptomics `ruv-text` chunk's output, with the record read as the template reads it
k227TrxText <- function(outcomes, side = "head") {
    dir <- withr::local_tempdir(pattern = "kit227-")
    record_path <- file.path(dir, "analysis_record.json")
    k227WriteRecord(list(k227Step("trx_ruvs", outcomes = outcomes)), record_path)
    record <- jsonlite::read_json(record_path, simplifyVector = FALSE)
    env <- new.env(parent = globalenv())
    env$overlay <- workflowRecordReportParameters(record_path, list(), list())
    steps <- Filter(function(step) identical(step$operation, "trx_ruvs"), record$steps)
    env$ruv <- steps[[length(steps)]]$outcomes
    env$no_ruv <- identical(env$ruv$ruv_applied, FALSE)
    # The template's `value` (Project Setup chunk), verbatim.
    env$value <- function(x, digits = 4) if (is.null(x)) "not recorded" else if (is.numeric(x)) format(round(x, digits)) else as.character(x)
    exprs <- k227Exprs(k227Path("transcriptomics", "transcriptomics_report.rmd", side), "ruv-text")
    paste(utils::capture.output(for (e in exprs) eval(e, env)), collapse = "\n")
}

k227TrxOutcomes <- function(mode) list(ruv_mode = mode, percentage_min = 5, percentage_max = 30, percentage_step = 5,
    k_penalty_weight = 0.5, max_acceptable_k = 3, k = 2, negative_control_percentage = 10,
    negative_control_count = 1500, controls = "inferred", separation_score = 0.3, composite_score = 0.25,
    ruv_applied = TRUE)

test_that("KIT-227 E8: a manual transcriptomics RUVs is not described as a grid search", {
    text <- k227TrxText(k227TrxOutcomes("manual"))
    expect_false(k227Has(text, "choosing the maximum composite score"))
    expect_true(k227Has(text, "selection mode"))
    expect_true(k227Has(text, "least significant"))
    expect_true(k227Has(text, "inferred empirical negative control"))
    expect_true(k227Has(text, paste("The selection mode was manual: the analyst set k = 2 and 10 % of the genes, giving 1500",
        "inferred empirical negative controls; no grid of percentages was searched.")))
})

# --- G1: the default text is byte-identical to the base templates' -------------------------------------------------------

# KIT-272 (declared edit of G1, by substitution on the base side; the frozen base copies are untouched): the RUV-applied
# Results bullets and the Methods phrase about Figure 1's panels d to f changed by design (the owner's decisions VAL-316
# (7) and VAL-332). File-local copies of the constants of test-kit272-prot-box-panel-text.R, never sourced.
k272OldBullet <- "\u2022 Density plots (panels d-f): the distribution of protein intensities per sample."
k272NewBullet <- paste0("\u2022 Principal component scores by group (panels d-f): for each experimental group, a box plot of the ",
    "samples' scores on principal component 1 (upper plot) and principal component 2 (lower plot), the same scores that ",
    "panels a-c plot as points. The line in each box is the median, the box spans the first to the third quartile, the ",
    "whiskers reach the most extreme scores within 1.5 times the interquartile range of the box, and scores beyond them ",
    "are drawn as points. The notch spans the median plus and minus 1.58 times the interquartile range divided by the ",
    "square root of the number of samples in the group: as a rough guide, two groups whose notches do not overlap have ",
    "medians that differ at about the 5% level (McGill, Tukey and Larsen, 1978). With fewer than about ten samples in a ",
    "group the notch extends beyond the box and should not be read. Each stage is a separate analysis, so the scores of ",
    "different panels are not on a common scale and their signs are arbitrary.")
k272OldPca <- "\u2022 Principal component analysis (panels a-c): samples plotted on the first two principal components,\n  coloured by experimental group."
k272NewPca <- paste0("\u2022 Principal component analysis (panels a-c): samples plotted on the first two principal components of the ",
    "log2 protein values (centred, not scaled), computed at each stage from the tenth of the proteins with the largest ",
    "ratio of standard deviation to mean of their log2 values; coloured by experimental group.")
# KIT-280 (declared edit, ticket section 6): a report made without an analysis record no longer states the earlier
# filter sentence ("Intensities below the Nth percentile ... were treated as missing"): it states that the rule is not
# recorded. The base templates' rendered sentence (percentile 1; the line wrap differs by template) is substituted.
k280OldFilterSentence <- function(name) if (identical(name, "DIANN")) paste0("Intensities below the 1st percentile of all \n",
    "intensities in the table were treated as missing.") else paste0("Intensities below the \n",
    "1st percentile of all intensities in the table were treated as missing.")
k280NewFilterSentence <- paste("Whether intensities at or below the 1st percentile of all intensities counted as missing in",
    "this filter is not recorded for this analysis.")

k272OldMethod <- "density boxplots"
k272NewMethod <- "box plots of the principal component scores"

#' `text` with the one `old` replaced by `new`, after asserting that `old` is there exactly once
k272Substitute <- function(text, old, new) {
    expect_identical(lengths(regmatches(text, gregexpr(old, text, fixed = TRUE))), 1L, info = old)
    sub(old, new, text, fixed = TRUE)
}

test_that("KIT-227 G1: Automatic on the group column prints the base text byte for byte, without and with a record", {
    proteomics_cases <- list(
        file_only = list(lines = c(k227ProtAuto, k227ProtCounts), record = NULL),
        default_record = list(lines = c(k227ProtAuto, k227ProtCounts), record = k227ProtRecord()))
    for (name in c(k227Proteomics, "DIANN_limpa")) {
        for (case in names(proteomics_cases)) {
            spec <- proteomics_cases[[case]]
            head <- k227Case(name, spec$lines, spec$record)
            base <- k227Case(name, spec$lines, spec$record, side = "base")
            info <- paste(name, case)
            expect_true(k227Has(head$method, "an automatic optimisation process was run"), info = info)
            expect_true(k227Has(head$method, "provided the best separation of biological groups"), info = info)
            # KIT-236 (declared edit (a), extended from its G0c): DIA with limpa without a DPC engine no longer prints
            # the "Technical replicate averaging" paragraph; every other byte of the method text is the base's.
            base_method <- if (identical(name, "DIANN_limpa")) sub(paste0("(?s)\\*\\*Technical replicate averaging\\*\\*\n\n",
                ".*?\\(Oberg et al\\., 2008\\)\\.\n\n"), "", base$method, perl = TRUE) else base$method
            base_method <- k272Substitute(base_method, k272OldMethod, k272NewMethod)
            if (identical(case, "file_only") && name %in% k227Proteomics) base_method <- k272Substitute(base_method,
                k280OldFilterSentence(name), k280NewFilterSentence)
            expect_identical(head$method, base_method, info = info)
            base_results <- k272Substitute(base$results, k272OldBullet, k272NewBullet)
            base_results <- k272Substitute(base_results, k272OldPca, k272NewPca)
            expect_identical(head$results, base_results, info = info)
        }
    }
    for (name in k227Omics) {
        omic_cases <- list(
            file_only = list(lines = k227OmicSection(), record = NULL),
            default_record = list(lines = k227OmicSection(), record = k227OmicRecord(name)))
        for (case in names(omic_cases)) {
            spec <- omic_cases[[case]]
            head <- k227Case(name, spec$lines, spec$record)
            base <- k227Case(name, spec$lines, spec$record, side = "base")
            info <- paste(name, case)
            expect_true(k227Has(head$method, "RUVIII-C"), info = info)
            # KIT-236 (declared edit (a)): the lipid/metab Methods and Results outside the RUV paragraph changed by
            # design (KIT-236 section 1); the RUV paragraph is the base's byte for byte.
            expect_identical(k227Paragraph(head$method, "Batch Effect Correction"),
                k227Paragraph(base$method, "Batch Effect Correction"), info = info)
        }
    }
    trx_head <- k227TrxText(k227TrxOutcomes("automatic"))
    expect_true(k227Has(trx_head, "choosing the maximum composite score"))
    expect_identical(trx_head, k227TrxText(k227TrxOutcomes("automatic"), side = "base"))
})
